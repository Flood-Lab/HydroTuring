"""
    WflowSbmAdapter

HydroTuring adapter for wflow_sbm, the SBM land-hydrology model of Deltares' Wflow.jl, at
v1.0.4 (commit 82df720). The adapter builds a Wflow model for one lumped catchment, runs it
through Wflow's own time loop and reads what it needs from the model after every step.

The catchment
-------------
Wflow is a distributed model: a TOML configuration, a netCDF of static maps and a netCDF of
forcing. A lumped catchment is given to it as one representative cell, every flux and store
a depth over that cell, so the catchment's area enters nothing but the conversion of runoff
to discharge. The cell is a land cell, a river cell and the outlet (local drain direction 5,
a pit) at once. It sits in a 2 x 2 grid in metres (`cell_length_in_meter__flag`) with the
other three cells inactive, because Wflow reads the cell size from the spacing of the
coordinates and drops dimensions of length one.

Hillslope length has no lumped definition, so like every other parameter the catchment does
not name it comes from the model Deltares ships to test Wflow: the cell is the square with
the mean cell area of the Moselle schematisation as Wflow computes it (0.00833 degrees,
562 446 m2, a side of 749.96 m). Wflow gives a pit cell its diagonal as flow length and the
matching flow width, so the land and subsurface kinematic waves drain a hillslope of 1061 m.
The river runs the same length and is as wide as makes it cover the share of the cell that
the Moselle model's rivers cover of its catchment (length times width over its river cells,
over its active area: 1.183 %), 6.27 m. Neither changes with the catchment's area.

What the catchment names is used: the soil store's capacity, the canopy's, the snow
threshold and the degree-day factor. Wflow's soil water capacity is defined in its code as
soil thickness times (theta_s - theta_r), so the thickness is set to the catchment's soil
capacity divided by that porosity, and the soil column holds exactly `soil_capacity_mm`. The
maximum canopy storage is `canopy_capacity_mm`; the snowfall and melt thresholds are
`snow_threshold_degC`; the degree-day factor is `degree_day_factor_mm_per_C_day`. Nothing is
derived from the forcing.

Every other parameter has no lumped counterpart and comes from the Moselle schematisation
`staticmaps-moselle.nc` (wflow-artifacts v1.0.0): the median over its 50 063 active cells
(river parameters: its 5 809 river cells), with the canopy gap fraction as the median over
cells of Wflow's own exp(-Kext * LAI) at each cell's annual-mean LAI, held constant so that
nothing follows the calendar. Where the test model has no value, Wflow's documented default
applies.

Switched off, each because the probe's catchment has none of it: reservoirs and lakes,
glaciers, paddies, irrigation and all water demand, floodplains, lateral snow transport
(which at a pit would carry snow out of the map), the frozen-soil infiltration reduction,
open water outside the river (`WaterFrac` 0) and the leakage out of the bottom of the
saturated store (`MaxLeakage` 0, as in the test model).

The step
--------
Wflow takes the step as `timestepsecs` and converts its per-day parameters itself, so the
rows go in at the step the case gives, as depths per step, and fluxes come back as rates.
Wflow labels forcing at the end of its interval; a row stamped with the start of an
interval is written at start + step. The interception scheme is Wflow's choice by step:
Gash for steps of 23 hours or more, which carries no canopy store between steps, and the
modified Rutter model below that, which does.

What is reported
----------------
Fluxes, as rates in mm per day unless stated:

* `pr`       the forcing, echoed as the text it was given
* `evspsbl`  Wflow's total actual evapotranspiration (`actevap`: interception, soil
             evaporation, transpiration, open water)
* `mrro`     everything that leaves the cell as flow over the step, as a depth: the
             river's discharge at the outlet, plus the overland and lateral subsurface flow
             out of the outlet cell, which Wflow passes out of the map rather than into the
             river because the outlet has no cell downstream of it
* `dis`      `mrro` over the catchment's area, in m3/s
* `gwex`     minus the leakage out of the saturated store (zero: `MaxLeakage` is 0)

States, absolute, in mm:

* `mrso`     the whole SBM soil column: the unsaturated store over all layers plus the
             saturated store. SBM's saturated zone is the lower part of the same column
             below a pseudo water table, bounded by the same thickness and porosity, and it
             is the capacity of that column that `soil_capacity_mm` names
* `gw`       zero. The `sbm` model type has no store below the soil column (that is the
             `sbm_gwf` type); its lateral subsurface flow moves the saturated store above
* `snw`      dry snow plus the liquid water held in the pack
* `canopy`   the canopy storage (identically zero at a daily step, where Gash applies)
* `channel`  the water in the land and river kinematic waves: generated as runoff and not
             yet released at the outlet

The budget closes over these by construction of the mapping, and run.json records both the
adapter's per-step residual and Wflow's own mass-balance errors for each component. The
schematisation (static maps, forcing and TOML, well under a megabyte for a ten-year daily
case) is written to a temporary directory on /tmp and removed after the run; nothing else is
written outside /io/output. Wflow is deterministic; the request seed is recorded and
otherwise unused.
"""
module WflowSbmAdapter

using Dates: Dates, DateTime, Second, @dateformat_str
using JSON: JSON
using Logging: Logging, ConsoleLogger, with_logger
using NCDatasets: NCDataset, defDim, defVar
using PrecompileTools: @compile_workload, @setup_workload
using TOML: TOML
using Wflow: Wflow

const MODEL = Dict{String, Any}("name" => "wflow_sbm", "version" => "1.0.4-ht.2")
const WFLOW = Dict{String, Any}(
    "package" => "Wflow.jl", "version" => "1.0.4",
    "commit" => "82df72031511339d50fd9142fa159d0ec13e73c5", "model_type" => "sbm",
)
const COLUMNS = ("time", "pr", "evspsbl", "mrro", "dis", "gwex", "mrso", "snw", "canopy", "gw", "channel")
const STEP_SECONDS = Dict("PT1D" => 86400, "PT1H" => 3600, "PT15M" => 900, "PT5M" => 300, "PT1M" => 60)
const FILL = -9999.0
const TIME_UNITS = "seconds since 1900-01-01 00:00:00"

# Wflow's Moselle test schematisation (staticmaps-moselle.nc, wflow-artifacts v1.0.0).
# Parameters are medians over its active cells, or its river cells where marked; the file
# stores Float32 and six significant digits are kept. The geometry is computed as Wflow
# computes it for that grid (cell lengths from the 0.00833-degree spacing through its
# `lattometres`).
const MOSELLE = (
    cell_side_m = 749.964,                       # side of the square with the mean cell area, 562 446 m2
    river_area_share = 0.0118307,                # sum over river cells of length * width, over the active area
    theta_s = 0.435808,                          # thetaS, saturated water content
    theta_r = 0.165628,                          # thetaR, residual water content
    kv_0 = 256.983,                              # KsatVer, mm/day at the soil surface
    f = 0.00335494,                              # f, 1/mm, decline of KsatVer with depth
    c = (9.38966, 9.70355, 10.0999, 10.0956),   # c, Brooks-Corey exponent per layer
    soil_thickness = 2000.0,                     # SoilThickness, mm (only without a capacity, or as a sensitivity)
    rooting_depth = 387.0,                       # RootingDepth, mm
    khfrac = 100.0,                              # KsatHorFrac
    slope = 0.0741422,                           # Slope, land surface
    n_land = 0.4768,                             # N, Manning, land surface
    pathfrac = 0.0,                              # PathFrac, compacted area
    infiltcap_path = 5.0,                        # InfiltCapPath, mm/day
    rootdistpar = -500.0,                        # rootdistpar
    max_leakage = 0.0,                           # MaxLeakage, mm/day
    cf_soil = 0.038,                             # cf_soil
    tti = 2.0,                                   # TTI, snowfall interval, degC
    whc = 0.1,                                   # WHC, snow liquid holding capacity
    e_r = 0.11,                                  # EoverR, Gash
    gap_fraction = 0.305933,                     # median of exp(-Kext * annual-mean LAI)
    n_river = 0.03,                              # N_River (river cells)
    river_slope = 0.0017048,                     # RiverSlope (river cells)
    river_width = 30.0,                          # wflow_riverwidth, m (river cells); only the one-cell sensitivity uses it
    river_depth = 1.0,                           # RiverDepth, bankfull, m (river cells)
)
const WFLOW_DEFAULTS = (cmax = 1.0, cfmax = 3.75)  # documented defaults, used only as fallbacks

# Wflow standard name => the variable written to the static maps.
const STATIC_NAMES = Dict(
    "soil_water__saturated_volume_fraction" => "thetaS",
    "soil_water__residual_volume_fraction" => "thetaR",
    "soil_surface_water__vertical_saturated_hydraulic_conductivity" => "KsatVer",
    "soil_water__vertical_saturated_hydraulic_conductivity_scale_parameter" => "f",
    "soil_layer_water__brooks_corey_exponent" => "c",
    "soil__thickness" => "SoilThickness",
    "vegetation_root__depth" => "RootingDepth",
    "subsurface_water__horizontal_to_vertical_saturated_hydraulic_conductivity_ratio" => "KsatHorFrac",
    "land_surface__slope" => "Slope",
    "land_surface_water_flow__manning_n_parameter" => "N",
    "compacted_soil__area_fraction" => "PathFrac",
    "compacted_soil_surface_water__infiltration_capacity" => "InfiltCapPath",
    "soil_wet_root__sigmoid_function_shape_parameter" => "rootdistpar",
    "soil_water_saturated_zone_bottom__max_leakage_volume_flux" => "MaxLeakage",
    "soil_surface_water__infiltration_reduction_parameter" => "cf_soil",
    "atmosphere_air__snowfall_temperature_threshold" => "TT",
    "atmosphere_air__snowfall_temperature_interval" => "TTI",
    "snowpack__melting_temperature_threshold" => "TTM",
    "snowpack__degree_day_coefficient" => "Cfmax",
    "snowpack__liquid_water_holding_capacity" => "WHC",
    "vegetation_canopy_water__mean_evaporation_to_mean_precipitation_ratio" => "EoverR",
    "vegetation_canopy__gap_fraction" => "CanopyGapFraction",
    "vegetation_water__storage_capacity" => "Cmax",
    "land_water_covered__area_fraction" => "WaterFrac",
    "river__length" => "wflow_riverlength",
    "river__width" => "wflow_riverwidth",
    "river__slope" => "RiverSlope",
    "river_water_flow__manning_n_parameter" => "N_River",
    "river_bank_water__depth" => "RiverDepth",
)

# --- the record ----------------------------------------------------------------------------

struct Forcing
    time::Vector{String}
    pr_text::Vector{String}
    pr::Vector{Float64}
    tas::Vector{Float64}
    pet::Vector{Float64}
end

function read_forcing(path::AbstractString)
    lines = filter(line -> !isempty(strip(line)), readlines(path))
    header = [strip(h) for h in split(lines[1], ',')]
    col = Dict(h => i for (i, h) in enumerate(header))
    for name in ("time", "pr", "tas", "pet")
        haskey(col, name) || error("forcing.csv has no '$name' column")
    end
    rows = [split(line, ',') for line in lines[2:end]]
    number(name) = [parse(Float64, strip(r[col[name]])) for r in rows]
    return Forcing(
        [String(strip(r[col["time"]])) for r in rows],
        [String(strip(r[col["pr"]])) for r in rows],
        number("pr"), number("tas"), number("pet"),
    )
end

function parse_time(text::AbstractString)
    t = replace(strip(text), 'T' => ' ')
    length(t) == 10 && return DateTime(t, dateformat"yyyy-mm-dd")
    length(t) == 16 && return DateTime(t, dateformat"yyyy-mm-dd HH:MM")
    length(t) == 19 && return DateTime(t, dateformat"yyyy-mm-dd HH:MM:SS")
    error("unrecognised time stamp $(repr(text))")
end

# --- the catchment -------------------------------------------------------------------------

"River width that makes a river along the pit cell's diagonal cover the Moselle river-area share."
river_width_for_share(cell::Float64) = MOSELLE.river_area_share * cell^2 / (sqrt(2.0) * cell)

"The representative cell's side in metres, the soil capacity, and the static maps' values."
function catchment(static::AbstractDict)
    area_km2 = Float64(static["area_km2"])
    area_km2 > 0 || error("area_km2 must be positive, not $area_km2")
    cell = MOSELLE.cell_side_m
    theta_e = MOSELLE.theta_s - MOSELLE.theta_r
    capacity = haskey(static, "soil_capacity_mm") ? Float64(static["soil_capacity_mm"]) :
        MOSELLE.soil_thickness * theta_e
    threshold = Float64(get(static, "snow_threshold_degC", 0.0))
    maps = Dict{String, Any}(
        "wflow_subcatch" => 1.0,
        "wflow_ldd" => 5.0,
        "thetaS" => MOSELLE.theta_s,
        "thetaR" => MOSELLE.theta_r,
        "KsatVer" => MOSELLE.kv_0,
        "f" => MOSELLE.f,
        "c" => collect(MOSELLE.c),
        "SoilThickness" => capacity / theta_e,
        "RootingDepth" => MOSELLE.rooting_depth,
        "KsatHorFrac" => MOSELLE.khfrac,
        "Slope" => MOSELLE.slope,
        "N" => MOSELLE.n_land,
        "PathFrac" => MOSELLE.pathfrac,
        "InfiltCapPath" => MOSELLE.infiltcap_path,
        "rootdistpar" => MOSELLE.rootdistpar,
        "MaxLeakage" => MOSELLE.max_leakage,
        "cf_soil" => MOSELLE.cf_soil,
        "TT" => threshold,
        "TTI" => MOSELLE.tti,
        "TTM" => threshold,
        "Cfmax" => Float64(get(static, "degree_day_factor_mm_per_C_day", WFLOW_DEFAULTS.cfmax)),
        "WHC" => MOSELLE.whc,
        "EoverR" => MOSELLE.e_r,
        "CanopyGapFraction" => MOSELLE.gap_fraction,
        "Cmax" => Float64(get(static, "canopy_capacity_mm", WFLOW_DEFAULTS.cmax)),
        "WaterFrac" => 0.0,
        "wflow_riverlength" => sqrt(2.0) * cell,
        "wflow_riverwidth" => river_width_for_share(cell),
        "RiverSlope" => MOSELLE.river_slope,
        "N_River" => MOSELLE.n_river,
        "RiverDepth" => MOSELLE.river_depth,
    )
    return cell, capacity, maps
end

const GEOMETRY_OVERRIDES = ("one_cell_sqrt_area", "cell_side_m", "moselle_thickness", "rooting_depth_fraction")

"""
    apply_overrides!(maps, cell, static, capacity, overrides) -> cell

Sensitivity settings for the tables in README.md; empty on the contract path. Applied in
this order: `one_cell_sqrt_area` (the geometry of adapter version ht.1: one cell with the
catchment's own area, its river 30 m wide along the diagonal), `cell_side_m` (a
representative cell of another side, river share kept), `moselle_thickness` (the Moselle
model's 2000 mm soil thickness, with theta_s lowered so that the column still holds the
stated capacity), `rooting_depth_fraction` (rooting depth as that fraction of the soil
thickness), then any static-map value by its name.
"""
function apply_overrides!(maps::AbstractDict, cell::Float64, static::AbstractDict, capacity::Float64,
                          overrides::AbstractDict)
    if Float64(get(overrides, "one_cell_sqrt_area", 0.0)) != 0.0
        cell = sqrt(Float64(static["area_km2"]) * 1.0e6)
        maps["wflow_riverlength"] = sqrt(2.0) * cell
        maps["wflow_riverwidth"] = MOSELLE.river_width
    end
    if haskey(overrides, "cell_side_m")
        cell = Float64(overrides["cell_side_m"])
        maps["wflow_riverlength"] = sqrt(2.0) * cell
        maps["wflow_riverwidth"] = river_width_for_share(cell)
    end
    if Float64(get(overrides, "moselle_thickness", 0.0)) != 0.0
        maps["SoilThickness"] = MOSELLE.soil_thickness
        maps["thetaS"] = MOSELLE.theta_r + capacity / MOSELLE.soil_thickness
    end
    if haskey(overrides, "rooting_depth_fraction")
        maps["RootingDepth"] = Float64(overrides["rooting_depth_fraction"]) * maps["SoilThickness"]
    end
    for (name, value) in overrides
        name in GEOMETRY_OVERRIDES && continue
        maps[name] = value
    end
    return cell
end

coordinates(cell) = ([0.5 * cell, 1.5 * cell], [1.5 * cell, 0.5 * cell])

function on_cell(value::Real)
    grid = Array{Union{Missing, Float64}}(missing, 2, 2)
    grid[1, 1] = Float64(value)
    return grid
end

function write_static_maps(path::AbstractString, cell::Float64, maps::AbstractDict)
    x, y = coordinates(cell)
    NCDataset(path, "c") do ds
        defDim(ds, "x", 2)
        defDim(ds, "y", 2)
        defDim(ds, "layer", 4)
        defVar(ds, "x", x, ("x",))
        defVar(ds, "y", y, ("y",))
        defVar(ds, "layer", [1.0, 2.0, 3.0, 4.0], ("layer",))
        river = zeros(2, 2)
        river[1, 1] = 1.0
        defVar(ds, "wflow_river", river, ("x", "y"))
        for (name, value) in maps
            name == "c" && continue
            defVar(ds, name, on_cell(value), ("x", "y"); fillvalue = FILL)
        end
        c = Array{Union{Missing, Float64}}(missing, 2, 2, 4)
        c[1, 1, :] = maps["c"]
        defVar(ds, "c", c, ("x", "y", "layer"); fillvalue = FILL)
    end
    return path
end

function write_forcing(path::AbstractString, cell::Float64, stamps::Vector{DateTime},
                       series::Vector{Pair{String, Vector{Float64}}})
    x, y = coordinates(cell)
    n = length(stamps)
    NCDataset(path, "c") do ds
        defDim(ds, "x", 2)
        defDim(ds, "y", 2)
        defDim(ds, "time", n)
        defVar(ds, "x", x, ("x",))
        defVar(ds, "y", y, ("y",))
        defVar(ds, "time", stamps, ("time",);
            attrib = ["units" => TIME_UNITS, "calendar" => "standard"])
        for (name, values) in series
            grid = Array{Union{Missing, Float64}}(missing, 2, 2, n)
            grid[1, 1, :] = values
            defVar(ds, name, grid, ("x", "y", "time"); fillvalue = FILL)
        end
    end
    return path
end

function write_config(path::AbstractString, first_end::DateTime, last_end::DateTime, dt::Int)
    config = Dict{String, Any}(
        "dir_input" => ".",
        "dir_output" => ".",
        "time" => Dict{String, Any}(
            "calendar" => "standard",
            # Wflow's clock starts one step before the first right-labelled forcing time.
            "starttime" => first_end - Second(dt),
            "endtime" => last_end,
            "time_units" => TIME_UNITS,
            "timestepsecs" => dt,
        ),
        "logging" => Dict{String, Any}("loglevel" => "warn", "silent" => true),
        "model" => Dict{String, Any}(
            "type" => "sbm",
            "cold_start__flag" => true,
            "cell_length_in_meter__flag" => true,
            "snow__flag" => true,
            "snow_gravitational_transport__flag" => false,
            "glacier__flag" => false,
            "reservoir__flag" => false,
            "soil_infiltration_reduction__flag" => false,
            "soil_layer__thickness" => [100, 300, 800],
            "water_mass_balance__flag" => true,
        ),
        "input" => Dict{String, Any}(
            "path_forcing" => "forcing.nc",
            "path_static" => "staticmaps.nc",
            "basin__local_drain_direction" => "wflow_ldd",
            "river_location__mask" => "wflow_river",
            "subbasin_location__count" => "wflow_subcatch",
            "forcing" => Dict{String, Any}(
                "atmosphere_water__precipitation_volume_flux" => "precip",
                "land_surface_water__potential_evaporation_volume_flux" => "pet",
                "atmosphere_air__temperature" => "temp",
            ),
            "static" => Dict{String, Any}(STATIC_NAMES),
        ),
    )
    open(io -> TOML.print(io, config; sorted = true), path, "w")
    return path
end

# --- the model -----------------------------------------------------------------------------

"""
    simulate(forcing, static, timestep, workdir; overrides)

Run one case. `overrides` exists for the sensitivity tables in README.md and is empty on
the contract path (see `apply_overrides!`).
"""
function simulate(forcing::Forcing, static::AbstractDict, timestep::AbstractString,
                  workdir::AbstractString; overrides::AbstractDict = Dict{String, Any}())
    dt = STEP_SECONDS[timestep]
    dt_days = dt / 86400
    n = length(forcing.time)
    area_km2 = Float64(static["area_km2"])
    cell, capacity, maps = catchment(static)
    cell = apply_overrides!(maps, cell, static, capacity, overrides)

    stamps = [parse_time(t) + Second(dt) for t in forcing.time]
    write_static_maps(joinpath(workdir, "staticmaps.nc"), cell, maps)
    write_forcing(joinpath(workdir, "forcing.nc"), cell, stamps, [
        "precip" => forcing.pr .* dt_days,
        "temp" => forcing.tas,
        "pet" => forcing.pet .* dt_days,
    ])
    toml = write_config(joinpath(workdir, "wflow_sbm.toml"), stamps[1], stamps[end], dt)

    model = Wflow.Model(Wflow.Config(toml))
    Wflow.load_fixed_forcing!(model)

    soil = model.land.soil.variables
    canopy = model.land.interception.variables
    snow = model.land.snow.variables
    overland = model.routing.overland_flow.variables
    river = model.routing.river_flow.variables
    river_bc = model.routing.river_flow.boundary_conditions
    lateral = model.routing.subsurface_flow.variables
    balance = model.mass_balance
    land = model.domain.land.parameters
    area = land.area[1]
    mm(volume_m3) = volume_m3 / area * 1000.0

    columns = zeros(n, 9)  # evspsbl mrro dis gwex mrso snw canopy gw channel
    stored(i) = columns[i, 5] + columns[i, 6] + columns[i, 7] + columns[i, 8] + columns[i, 9]
    initial = soil.ustoredepth[1] + soil.satwaterdepth[1] + snow.snow_storage[1] + snow.snow_water[1] +
        canopy.canopy_storage[1] + mm(overland.storage[1] + river.storage[1])
    worst_residual, cumulative_residual = 0.0, 0.0
    wflow_errors = zeros(4)
    min_river_inflow = Inf

    for i in 1:n
        Wflow.run_timestep!(model; write_model_output = false)
        outflow = river.q_av[1] + overland.q_av[1] + lateral.ssf[1] / 86400.0  # m3/s out of the cell
        leakage = soil.actleakage[1]
        columns[i, 1] = soil.actevap[1] / dt_days
        columns[i, 2] = mm(outflow * 86400.0)
        columns[i, 3] = columns[i, 2] * area_km2 / 86.4
        columns[i, 4] = -leakage / dt_days
        columns[i, 5] = soil.ustoredepth[1] + soil.satwaterdepth[1]
        columns[i, 6] = snow.snow_storage[1] + snow.snow_water[1]
        columns[i, 7] = canopy.canopy_storage[1]
        columns[i, 8] = 0.0
        columns[i, 9] = mm(overland.storage[1] + river.storage[1])

        before = i == 1 ? initial : stored(i - 1)
        residual = forcing.pr[i] * dt_days - leakage - soil.actevap[1] - columns[i, 2] * dt_days -
            (stored(i) - before)
        worst_residual = max(worst_residual, abs(residual))
        cumulative_residual += residual
        wflow_errors[1] = max(wflow_errors[1], abs(balance.land_water_balance.error[1]))
        wflow_errors[2] = max(wflow_errors[2], abs(balance.routing.overland_water_balance.error[1]))
        wflow_errors[3] = max(wflow_errors[3], abs(balance.routing.river_water_balance.error[1]))
        wflow_errors[4] = max(wflow_errors[4], abs(balance.routing.subsurface_water_balance.error[1]))
        min_river_inflow = min(min_river_inflow, river_bc.inwater[1])
    end
    Wflow.close_files(model; delete_output = false)

    interception = model.land.interception isa Wflow.GashInterceptionModel ?
        "Gash (Wflow uses it for steps of 23 hours or more; no canopy store is carried between steps)" :
        "modified Rutter (Wflow uses it for steps under 23 hours; the canopy store is carried)"
    notes = Dict{String, Any}(
        "wflow" => WFLOW,
        "julia" => string(VERSION),
        "timestep" => timestep,
        "interception" => interception,
        "domain" => Dict{String, Any}(
            "cells" => "one representative cell in a 2 x 2 metre grid: land, river and outlet (ldd 5)",
            "cell_source" => "square with the mean cell area of Wflow's Moselle test model (0.00833 degrees, Wflow's lattometres)",
            "cell_side_m" => cell,
            "cell_area_m2" => area,
            "flow_length_m" => land.flow_length[1],
            "flow_width_m" => land.flow_width[1],
            "surface_flow_width_m" => land.surface_flow_width[1],
            "river_length_m" => maps["wflow_riverlength"],
            "river_width_m" => maps["wflow_riverwidth"],
            "river_fraction" => land.river_fraction[1],
            "catchment_area_km2" => area_km2,
            "catchment_area_enters" => "only dis = mrro * area_km2 / 86.4",
        ),
        "parameters" => Dict{String, Any}(
            "from_static_json" => Dict{String, Any}(
                "SoilThickness_mm" => maps["SoilThickness"],
                "soil_water_capacity_mm" => capacity,
                "Cmax_mm" => maps["Cmax"],
                "TT_TTM_degC" => maps["TT"],
                "Cfmax_mm_per_degC_day" => maps["Cfmax"],
            ),
            "from_moselle_test_model" => Dict{String, Any}(String(k) => v for (k, v) in pairs(MOSELLE)),
            "effective_rooting_depth_mm" => model.land.vegetation_parameters.rootingdepth[1],
            "static_attributes_unused" => sort([String(k) for k in keys(static) if !(k in (
                "area_km2", "soil_capacity_mm", "canopy_capacity_mm", "snow_threshold_degC",
                "degree_day_factor_mm_per_C_day"))]),
        ),
        "switched_off" => [
            "reservoirs and lakes", "glaciers", "water demand, irrigation and paddies",
            "floodplains", "lateral snow transport", "frozen-soil infiltration reduction",
            "open water outside the river (WaterFrac 0)", "leakage from the saturated store (MaxLeakage 0)",
        ],
        "reported" => Dict{String, Any}(
            "evspsbl" => "actevap: interception + soil evaporation + transpiration + open water",
            "mrro" => "river q_av at the outlet + overland q_av + lateral subsurface flow out of the outlet cell, as a depth over the cell",
            "dis" => "mrro over the catchment's area, m3/s",
            "gwex" => "minus the leakage from the saturated store (identically zero, MaxLeakage 0)",
            "mrso" => "unsaturated store (all layers) + saturated store: the SBM soil column",
            "gw" => "identically zero: the sbm model type has no store below the soil column",
            "snw" => "dry snow + liquid water in the pack",
            "canopy" => "canopy storage",
            "channel" => "land and river kinematic-wave storage, as a depth over the cell",
        ),
        "budget" => Dict{String, Any}(
            "adapter_max_step_residual_mm" => worst_residual,
            "adapter_cumulative_residual_mm" => cumulative_residual,
            "wflow_max_abs_error_land_mm" => wflow_errors[1],
            "wflow_max_abs_error_overland_m3s" => wflow_errors[2],
            "wflow_max_abs_error_river_m3s" => wflow_errors[3],
            "wflow_max_abs_error_subsurface_m3" => wflow_errors[4],
            "min_river_lateral_inflow_m3s" => min_river_inflow,
        ),
    )
    isempty(overrides) || (notes["overrides"] = Dict{String, Any}(String(k) => v for (k, v) in overrides))
    return columns, notes
end

# --- the contract --------------------------------------------------------------------------

function main(args::Vector{String})
    started = time()
    at = findfirst(==("--request"), args)
    (at === nothing || at == length(args)) && error("usage: ht_adapter.jl --request /io/request.json")
    request_path = abspath(args[at + 1])
    request = JSON.parsefile(request_path)
    io_dir = dirname(request_path)
    timestep = String(get(request, "timestep", "PT1D"))
    haskey(STEP_SECONDS, timestep) || error("unsupported timestep $(repr(timestep))")

    forcing = read_forcing(joinpath(io_dir, request["input"]["forcing"]))
    static = JSON.parsefile(joinpath(io_dir, request["input"]["static"]))
    if haskey(request, "n_steps") && Int(request["n_steps"]) != length(forcing.time)
        error("request asks for $(request["n_steps"]) steps; the forcing has $(length(forcing.time)) rows")
    end

    # The schematisation Wflow reads (static maps, forcing, TOML) goes to a scratch directory
    # on /tmp and is removed afterwards.
    workdir = mktempdir()
    columns, notes = try
        with_logger(ConsoleLogger(stderr, Logging.Warn)) do
            simulate(forcing, static, timestep, workdir)
        end
    finally
        rm(workdir; recursive = true, force = true)
    end

    table = joinpath(io_dir, request["output"]["table"])
    mkpath(dirname(table))
    open(table, "w") do io
        println(io, join(COLUMNS, ","))
        for i in eachindex(forcing.time)
            print(io, forcing.time[i], ",", forcing.pr_text[i])
            for j in 1:size(columns, 2)
                print(io, ",", columns[i, j])
            end
            println(io)
        end
    end

    notes["seed"] = get(request, "seed", nothing)
    notes["max_rss_mb"] = round(Sys.maxrss() / 2^20; digits = 1)
    run = Dict{String, Any}(
        "status" => "ok",
        "model" => MODEL,
        "n_steps" => length(forcing.time),
        "wall_seconds" => round(time() - started; digits = 2),
        "notes" => notes,
    )
    open(io -> write(io, JSON.json(run)), joinpath(io_dir, request["output"]["run"]), "w")
    return 0
end

# --- compiled when the image is built ------------------------------------------------------

"A small synthetic case with rain, a freeze and a thaw, for the precompile workload."
function synthetic_case(dir::AbstractString, timestep::AbstractString, n::Int)
    dt = STEP_SECONDS[timestep]
    format = timestep == "PT1D" ? dateformat"yyyy-mm-dd" : dateformat"yyyy-mm-dd HH:MM"
    mkpath(joinpath(dir, "input"))
    mkpath(joinpath(dir, "output"))
    open(joinpath(dir, "input", "forcing.csv"), "w") do io
        println(io, "time,pr,tas,pet")
        for k in 0:(n - 1)
            day = k * dt / 86400
            pr = k % 5 == 0 ? 24.0 : (k % 3 == 0 ? 2.5 : 0.0)
            tas = 4.0 - 9.0 * cos(2pi * day / 20)
            pet = max(0.0, 2.0 + 1.5 * sin(2pi * day / 20))
            println(io, Dates.format(DateTime(2000, 1, 1) + Second(dt * k), format), ",", pr, ",", tas, ",", pet)
        end
    end
    write(joinpath(dir, "input", "static.json"), JSON.json(Dict(
        "area_km2" => 250.0, "soil_capacity_mm" => 320.0, "canopy_capacity_mm" => 2.0,
        "degree_day_factor_mm_per_C_day" => 3.2, "baseflow_coefficient" => 0.006,
        "snow_threshold_degC" => 0.0, "latitude_deg" => 40.0,
    )))
    write(joinpath(dir, "request.json"), JSON.json(Dict(
        "case_id" => "warm-up", "seed" => 1, "timestep" => timestep, "n_steps" => n,
        "input" => Dict("forcing" => "input/forcing.csv", "static" => "input/static.json"),
        "output" => Dict("table" => "output/result.csv", "run" => "output/run.json"),
    )))
    return joinpath(dir, "request.json")
end

@setup_workload begin
    @compile_workload begin
        for (timestep, n) in (("PT1D", 60), ("PT1H", 72))
            dir = mktempdir()
            try
                main(["--request", synthetic_case(dir, timestep, n)])
            finally
                rm(dir; recursive = true, force = true)
            end
        end
    end
end

end # module
