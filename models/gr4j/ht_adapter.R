#!/usr/bin/env Rscript
# HydroTuring adapter for GR4J, from airGR 1.7.9.
#
# The physics is airGR's compiled kernel, called through RunModel_GR4J (RunModel_GR4H on an
# hourly row). This file is the units layer: rates in mm/day become depths per row on the way
# in and rates again on the way out, and the stores are read from what the kernel returned.
# README.md gives the reason for each choice below and how far it moves the results.

suppressPackageStartupMessages({
  library(airGR)
  library(jsonlite)
})

MODEL <- list(name = "gr4j", version = "1.7.9-4bb401c.2", airGR = as.character(packageVersion("airGR")))

# Length of one row in days, and the kernel that row runs. PT1D is GR4J and PT1H is GR4H.
STEPS <- list(PT1D = 1, PT1H = 1 / 24)

# The median of each column of Param_Sets_GR4J, the 27 generalist parameter sets airGR ships
# for use when a catchment cannot be calibrated (Andreassian et al. 2014, WRR 50(10),
# doi:10.1002/2013WR014761). They are per-day parameters and X4u is not yet X4.
GENERALIST <- list(X1 = 291.0, X2 = -1.1, X3 = 76.2, X4u = 2.15)

# airGR's help page for that data set: X4 = X4u / 5.995 * S^0.3, S in km2. S is held at the
# probes' 250 km2 rather than read from the case. Reading it would change the storm response
# with the area told, which is a statement about the adapter's scaling and not about GR4J,
# whose X4 is a fixed parameter of the catchment.
X4_REFERENCE_AREA_KM2 <- 250.0

parameters <- function(static) {
  x1_source <- "median of Param_Sets_GR4J"
  x1 <- GENERALIST$X1
  if (!is.null(static$soil_capacity_mm)) {
    x1 <- as.numeric(static$soil_capacity_mm)
    x1_source <- "static.json soil_capacity_mm"
  }
  list(
    X1 = x1, X1_source = x1_source,
    X2 = GENERALIST$X2,
    X3 = GENERALIST$X3,
    X4 = GENERALIST$X4u / 5.995 * X4_REFERENCE_AREA_KM2^0.3
  )
}

# The per-day parameters in the units of an hourly row. The routing store releases R^5 / (4 X3^4)
# per step when R << X3 and exchanges X2 (R/X3)^3.5 per step, so the same store drains at the same
# rate per unit time when X3 grows by 24^(1/4) and X2 shrinks by 24^(-1/8). X4 is a time, so it is
# counted in rows. X1 is a capacity and does not change. The percolation constant is the kernel's
# own: 9/4 in GR4J and 21/4 in GR4H, not ours to choose.
hourly_parameters <- function(p) {
  p$X2 <- p$X2 * 24^(-1 / 8)
  p$X3 <- p$X3 * 24^(1 / 4)
  p$X4 <- p$X4 * 24
  p
}

run_kernel <- function(forcing, static, timestep) {
  dt_days <- STEPS[[timestep]]
  daily <- parameters(static)
  p <- if (timestep == "PT1H") hourly_parameters(daily) else daily
  fun <- if (timestep == "PT1H") RunModel_GR4H else RunModel_GR4J

  # Rates to depths per row. The kernel's own unit is mm per its step.
  precip <- as.numeric(forcing$pr) * dt_days
  pet <- as.numeric(forcing$pet) * dt_days
  dates <- as.POSIXct(
    forcing$time, tz = "UTC",
    tryFormats = c("%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%dT%H:%M", "%Y-%m-%d")
  )
  inputs <- CreateInputsModel(FUN_MOD = fun, DatesR = dates, Precip = precip, PotEvap = pet, verbose = FALSE)
  # No warm-up period, so every row is reported, spinup included. The stores start at airGR's
  # default levels, 0.3 X1 and 0.5 X3 with both unit hydrographs empty, always read off the
  # per-day parameters. Taken from the hourly X3, which is larger by 24^(1/4), the hourly run
  # would start with 46 mm more water than the daily one and the resolution probe would
  # measure that instead of the step.
  prod0 <- 0.3 * daily$X1
  rout0 <- 0.5 * daily$X3
  initial <- CreateIniStates(FUN_MOD = fun, InputsModel = inputs, ProdStore = prod0, RoutStore = rout0, verbose = FALSE)
  options <- CreateRunOptions(
    FUN_MOD = fun, InputsModel = inputs, IndPeriod_Run = seq_along(precip),
    IndPeriod_WarmUp = 0L, IniStates = initial, Outputs_Sim = "all", warnings = FALSE, verbose = FALSE
  )
  out <- fun(InputsModel = inputs, RunOptions = options, Param = c(p$X1, p$X2, p$X3, p$X4))

  # The kernel carries the unit hydrographs' contents but returns them only for the last step,
  # so the store at every step is what entered them less what left: PR in, Q9 and Q1 out. That
  # is an account of the kernel's own fluxes, and it is held to the state the kernel returns.
  # The first ordinate of each hydrograph in StateEnd is what was already released this step and
  # is still sitting in the array until the next step overwrites it, so it is not counted.
  uh <- cumsum(out$PR - out$Q9 - out$Q1)
  uh_end <- sum(out$StateEnd$UH$UH1[-1]) + sum(out$StateEnd$UH$UH2[-1])
  uh_error <- abs(uh[length(uh)] - uh_end)
  if (uh_error > 1e-8) {
    stop(sprintf("unit-hydrograph account differs from the kernel's state by %.3g mm", uh_error))
  }

  n <- length(precip)
  residual <- sum(precip) - sum(out$AE) - sum(out$Qsim) + sum(out$AExch) -
    ((out$Prod[n] - prod0) + (out$Rout[n] - rout0) + uh[n])

  list(
    rows = data.frame(
      evspsbl = out$AE / dt_days,
      mrro = out$Qsim / dt_days,
      gwex = out$AExch / dt_days,
      mrso = out$Prod,
      gw = out$Rout,
      channel = uh
    ),
    meta = list(
      kernel = if (timestep == "PT1H") "RunModel_GR4H" else "RunModel_GR4J",
      parameters_per_day = daily,
      parameters_as_run = p[c("X1", "X2", "X3", "X4")],
      initial_state_mm = list(production_store = prod0, routing_store = rout0, unit_hydrographs = 0),
      checks = list(
        unit_hydrograph_vs_kernel_state_mm = uh_error,
        whole_record_budget_residual_mm = residual
      )
    )
  )
}

fmt <- function(x) sprintf("%.15g", x)

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  at <- match("--request", args)
  if (is.na(at) || at == length(args)) stop("usage: ht_adapter.R --request <request.json>")
  request_path <- normalizePath(args[at + 1])
  io <- dirname(request_path)
  request <- fromJSON(request_path, simplifyVector = TRUE)
  timestep <- request$timestep
  if (is.null(STEPS[[timestep]])) stop(sprintf("unsupported timestep %s", timestep))

  # Everything as text, so `pr` goes back exactly as it came.
  forcing <- read.csv(file.path(io, request$input$forcing), colClasses = "character", check.names = FALSE)
  static <- fromJSON(file.path(io, request$input$static), simplifyVector = TRUE)

  run <- run_kernel(forcing, static, timestep)
  rows <- run$rows
  area <- if (is.null(static$area_km2)) NA_real_ else as.numeric(static$area_km2)
  if (is.na(area)) stop("static.json has no area_km2, so `dis` cannot be reported")

  table <- data.frame(
    time = forcing$time,
    pr = forcing$pr,
    evspsbl = fmt(rows$evspsbl),
    mrro = fmt(rows$mrro),
    # mm/day over the catchment to m3/s.
    dis = fmt(rows$mrro * area / 86.4),
    gwex = fmt(rows$gwex),
    mrso = fmt(rows$mrso),
    snw = 0,
    canopy = 0,
    gw = fmt(rows$gw),
    channel = fmt(rows$channel)
  )
  write.csv(table, file.path(io, request$output$table), row.names = FALSE, quote = FALSE)

  meta <- c(
    list(status = "ok", model = MODEL, n_steps = nrow(table), timestep = timestep),
    run$meta,
    list(
      store_mapping = list(
        mrso = "production store (Prod)",
        gw = "routing store (Rout)",
        channel = "water inside the two unit hydrographs (cumulative PR - Q9 - Q1)"
      ),
      zero_stores = list(
        snw = "GR4J has no snow module: all precipitation is treated as liquid",
        canopy = "GR4J has no interception store: the P/E neutralisation holds no water between steps"
      ),
      not_reported = list(
        snm = "GR4J has no snow module, so nothing crosses a snow-module boundary",
        energy_fluxes = "GR4J has no energy balance: hfls, hfss, hfg and the radiation terms do not exist in it",
        stage = "GR4J has no channel geometry, so no water-surface elevation"
      )
    )
  )
  write(toJSON(meta, auto_unbox = TRUE, digits = NA, null = "null"), file.path(io, request$output$run))
  invisible(0L)
}

main()
