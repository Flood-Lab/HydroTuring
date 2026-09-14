"""Add read-only, single-column thermal diagnostics to pinned Noah-MP source.

The two output files are written in the model process's working directory.
Use a serial, single-grid-cell HRLDAS case with soil_update_steps=1 and
OptSnowSoilTempTime=1. No native physical equations or parameter values change.
"""

from __future__ import annotations

import argparse
from pathlib import Path


THERMAL_COLUMNS = (
    "step,time_end_s,dt_s,grid_i,grid_j,opt_stc,snow_layers,"
    "g_top_w_m2,g_bottom_w_m2,t1_start_k,t1_thermal_end_k,t2_thermal_end_k,"
    "heat_capacity_area_j_m2_k,conductivity_1_w_m_k,z1_m,z2_m,"
    "penetrated_sw_w_m2,soil_moisture_1,soil_liquid_1,soil_ice_max_fraction,"
    "snow_swe_mm,snow_depth_m,tsoil_min_k"
)
COLUMN_COLUMNS = (
    "step,time_end_s,dt_s,grid_i,grid_j,soil_step_executed,"
    "t1_start_k,t1_final_k,soil_moisture_1_start,soil_moisture_1_end,"
    "soil_liquid_1_start,soil_liquid_1_end,soil_ice_max_fraction,"
    "snow_swe_mm,snow_layers,veg_fraction,tsoil_min_k,"
    "heat_latent_ground_w_m2,heat_ground_w_m2,heat_precip_advected_w_m2,"
    "evap_ground_mm_s,water_balance_error_mm,energy_balance_error_w_m2"
)


def replace_once(source: str, old: str, new: str) -> str:
    if source.count(old) != 1:
        raise ValueError(f"Expected one pinned-source insertion point: {old!r}")
    return source.replace(old, new, 1)


def fortran_header(columns: str, unit: str) -> str:
    chunks = [columns[i:i + 88] for i in range(0, len(columns), 88)]
    text = f"       write({unit},'(a)') &\n"
    text += " // &\n".join(f"          '{part}'" for part in chunks)
    return text + "\n"


def thermal_patch(source: str) -> str:
    source = replace_once(
        source, "! local variable\n",
        "! HYDROTURING_READONLY_DIAGNOSTICS: serial single-column validation.\n"
        "    integer, save :: ht_diag_unit, ht_diag_step = 0\n"
        "    real(kind=kind_noahmp), save :: ht_diag_time = 0.0_kind_noahmp\n"
        "    real(kind=kind_noahmp) :: ht_diag_t1_start, ht_diag_bottom\n\n"
        "! local variable\n",
    )
    source = replace_once(
        source, "    ! initialization\n",
        "    ht_diag_t1_start = TemperatureSoilSnow(1)\n\n"
        "    ! initialization\n",
    )
    diagnostic = """    ! HYDROTURING_READONLY_DIAGNOSTICS: direct Fourier flux, not a residual.
    if (ht_diag_step == 0) then
       open(newunit=ht_diag_unit, file='native_thermal_steps.csv', status='replace', action='write')
""" + fortran_header(THERMAL_COLUMNS, "ht_diag_unit") + """    endif
    ht_diag_step = ht_diag_step + 1
    ht_diag_time = ht_diag_time + TimeStep
    ht_diag_bottom = 2.0_kind_noahmp * noahmp%energy%state%ThermConductSoilSnow(1) * &
         (TemperatureSoilSnow(1) - TemperatureSoilSnow(2)) / &
         (-noahmp%config%domain%DepthSnowSoilLayer(2))
    write(ht_diag_unit,'(*(g0,:,","))') &
         ht_diag_step, ht_diag_time, TimeStep, &
         noahmp%config%domain%GridIndexI, noahmp%config%domain%GridIndexJ, &
         noahmp%config%nmlist%OptSnowSoilTempTime, NumSnowLayerNeg, &
         noahmp%energy%flux%HeatGroundTotMean, ht_diag_bottom, &
         ht_diag_t1_start, TemperatureSoilSnow(1), TemperatureSoilSnow(2), &
         (-noahmp%config%domain%DepthSnowSoilLayer(1)) * noahmp%energy%state%HeatCapacSoilSnow(1), &
         noahmp%energy%state%ThermConductSoilSnow(1), &
         noahmp%config%domain%DepthSnowSoilLayer(1), noahmp%config%domain%DepthSnowSoilLayer(2), &
         noahmp%energy%flux%RadSwPenetrateGrd(1), &
         noahmp%water%state%SoilMoisture(1), noahmp%water%state%SoilLiqWater(1), &
         maxval(noahmp%water%state%SoilMoisture - noahmp%water%state%SoilLiqWater), &
         noahmp%water%state%SnowWaterEquiv, noahmp%water%state%SnowDepth, &
         minval(TemperatureSoilSnow(1:NumSoilLayer))
    flush(ht_diag_unit)

"""
    return replace_once(
        source, "    ! deallocate local arrays to avoid memory leaks\n",
        diagnostic + "    ! deallocate local arrays to avoid memory leaks\n",
    )


def column_patch(source: str) -> str:
    declaration = "    type(noahmp_type), intent(inout) :: noahmp\n"
    source = replace_once(
        source, declaration,
        declaration + "\n"
        "    ! HYDROTURING_READONLY_DIAGNOSTICS: serial single-column validation.\n"
        "    integer, save :: ht_diag_unit, ht_diag_step = 0\n"
        "    real(kind=kind_noahmp), save :: ht_diag_time = 0.0_kind_noahmp\n"
        "    real(kind=kind_noahmp) :: ht_diag_t1_start, ht_diag_theta_start, ht_diag_liquid_start\n",
    )
    source = replace_once(
        source, "    call ProcessAtmosForcing(noahmp)\n",
        "    ht_diag_t1_start = noahmp%energy%state%TemperatureSoilSnow(1)\n"
        "    ht_diag_theta_start = noahmp%water%state%SoilMoisture(1)\n"
        "    ht_diag_liquid_start = noahmp%water%state%SoilLiqWater(1)\n\n"
        "    call ProcessAtmosForcing(noahmp)\n",
    )
    diagnostic = """
    ! HYDROTURING_READONLY_DIAGNOSTICS: end of every complete native column step.
    if (ht_diag_step == 0) then
       open(newunit=ht_diag_unit, file='native_column_steps.csv', status='replace', action='write')
""" + fortran_header(COLUMN_COLUMNS, "ht_diag_unit") + """    endif
    ht_diag_step = ht_diag_step + 1
    ht_diag_time = ht_diag_time + noahmp%config%domain%MainTimeStep
    write(ht_diag_unit,'(*(g0,:,","))') &
         ht_diag_step, ht_diag_time, noahmp%config%domain%MainTimeStep, &
         noahmp%config%domain%GridIndexI, noahmp%config%domain%GridIndexJ, &
         merge(1,0,noahmp%config%domain%FlagSoilProcess), &
         ht_diag_t1_start, noahmp%energy%state%TemperatureSoilSnow(1), &
         ht_diag_theta_start, noahmp%water%state%SoilMoisture(1), &
         ht_diag_liquid_start, noahmp%water%state%SoilLiqWater(1), &
         maxval(noahmp%water%state%SoilMoisture - noahmp%water%state%SoilLiqWater), &
         noahmp%water%state%SnowWaterEquiv, noahmp%config%domain%NumSnowLayerNeg, &
         noahmp%energy%state%VegFrac, &
         minval(noahmp%energy%state%TemperatureSoilSnow(1:noahmp%config%domain%NumSoilLayer)), &
         noahmp%energy%flux%HeatLatentGrd, noahmp%energy%flux%HeatGroundTot, &
         noahmp%energy%flux%HeatPrecipAdvSfc, noahmp%water%flux%EvapGroundNet, &
         noahmp%water%state%WaterBalanceError, noahmp%energy%state%EnergyBalanceError
    flush(ht_diag_unit)

"""
    return replace_once(
        source, "    call BalanceEnergyCheck(noahmp) \n",
        "    call BalanceEnergyCheck(noahmp) \n" + diagnostic,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path,
                        help="HRLDAS source directory containing noahmp/src")
    args = parser.parse_args()
    source_dir = args.source.resolve() / "noahmp" / "src"
    edits = []
    for name, patch in (
        ("SoilSnowTemperatureSolverMod.F90", thermal_patch),
        ("NoahmpMainMod.F90", column_patch),
    ):
        path = source_dir / name
        original = path.read_text(encoding="utf-8")
        if "HYDROTURING_READONLY_DIAGNOSTICS" in original:
            raise ValueError(f"Already instrumented: {path}")
        updated = patch(original)
        long_lines = [n for n, line in enumerate(updated.splitlines(), 1)
                      if len(line) > 132 and line not in original.splitlines()]
        if long_lines:
            raise ValueError(f"New Fortran lines exceed 132 columns in {name}: {long_lines}")
        edits.append((path, updated))
    for path, updated in edits:
        path.write_text(updated, encoding="utf-8", newline="\n")
        print(f"Instrumented {path}")


if __name__ == "__main__":
    main()
