#!/usr/bin/env Rscript
# HydroTuring adapter for TOPMODEL (R package topmodel 0.7.5).
#
# The physics is the package's C core, unchanged, called through ht_driver.c, which is the same
# loop as c_topmodel() and also copies out the stores the package frees. This file is the units
# layer: mm/day to metres per row on the way in, metres per row to mm/day on the way out, and
# the stores in mm. README.md gives the reason for each choice and how far it moves the results.

suppressPackageStartupMessages(library(jsonlite))

MODEL <- list(name = "topmodel", version = "0.7.5-f60752c.2")

# Length of one row, in days and in the hours the core takes it in.
STEPS <- list(PT1D = list(days = 1, hours = 24), PT1H = list(days = 1 / 24, hours = 1))

# The datum of the saturated zone, in multiples of m. TOPMODEL's saturated zone is a deficit
# below saturation with no floor, and a store has to be reported as an amount. Transmissivity
# falls as exp(-S/m), so ten m down it is exp(-10) of its value at saturation and the store
# has effectively stopped draining. Only differences of it enter a budget.
GW_DATUM_IN_M <- 10

model_inputs <- function(static) {
  data("huagrahuma", package = "topmodel", envir = environment())
  p <- huagrahuma$parameters
  delay <- as.matrix(huagrahuma$delay)
  delay_source <- "huagrahuma$delay (the example catchment's own)"
  # The delay function is distance to the outlet against cumulative area. When the case gives
  # the channel, half the area lies within the centroid's distance and all of it within the
  # main channel's.
  if (!is.null(static$main_channel_length_km) && !is.null(static$centroid_channel_length_km)) {
    delay <- rbind(
      c(0, 0),
      c(1000 * as.numeric(static$centroid_channel_length_km), 0.5),
      c(1000 * as.numeric(static$main_channel_length_km), 1)
    )
    delay_source <- "main_channel_length_km and centroid_channel_length_km"
  }
  srmax_source <- "huagrahuma$parameters"
  if (!is.null(static$soil_capacity_mm)) {
    p[["Srmax"]] <- as.numeric(static$soil_capacity_mm) / 1000
    srmax_source <- "static.json soil_capacity_mm"
  }
  list(parameters = p, topidx = as.matrix(huagrahuma$topidx), delay = delay,
       delay_source = delay_source, srmax_source = srmax_source)
}

run_core <- function(forcing, static, timestep) {
  step <- STEPS[[timestep]]
  m <- model_inputs(static)
  p <- m$parameters
  p[["dt"]] <- step$hours

  # Rates in mm/day to metres per row, the unit the core's rain and ETp are in.
  rain <- as.numeric(forcing$pr) * step$days / 1000
  etp <- as.numeric(forcing$pet) * step$days / 1000
  n <- length(rain)

  dyn.load("/opt/topmodel/ht_topmodel.so")
  d <- .C("ht_topmodel", as.double(p), as.double(m$topidx), as.double(m$delay), as.double(rain), as.double(etp),
          as.integer(nrow(m$topidx)), as.integer(n), as.integer(nrow(m$delay)),
          Q = double(n), qo = double(n), qs = double(n), S = double(n), fex = double(n), Ea = double(n),
          Srz = double(n), Suz = double(n), ex = double(n), transit0 = double(1), weight = double(1),
          PACKAGE = "ht_topmodel")

  srmax <- p[["Srmax"]]
  datum <- GW_DATUM_IN_M * p[["m"]]
  # Water in the root zone is what its deficit leaves of the capacity, over the same class
  # weights the core uses, which do not quite sum to one on this table: see README.md.
  root <- d$weight * srmax - d$Srz
  gw <- d$Suz + (datum - d$S)
  # What the delay function holds: what was on its way at the start, plus what has entered it
  # (qo + qs, the core's total flow for the step), less what has reached the outlet.
  channel <- d$transit0 + cumsum(d$qo + d$qs) - cumsum(d$Q)

  # The whole-record budget over the reported stores, with the terms that explain it. The
  # initial stores are the core's own: Srz = Sr0 and Suz = 0 in every class, S_mean from qs0.
  root0 <- d$weight * (srmax - p[["Sr0"]])
  qss <- exp(p[["lnTe"]] + log(p[["dt"]]) - lambda_of(m$topidx))
  s0 <- -p[["m"]] * log(p[["qs0"]] * p[["dt"]] / qss)
  gw0 <- datum - s0
  storage_change <- (root[n] - root0) + (gw[n] - gw0) + (channel[n] - d$transit0)
  residual <- sum(rain) - sum(d$Ea) - sum(d$Q) - storage_change

  list(
    rows = data.frame(
      evspsbl = 1000 * d$Ea / step$days,
      mrro = 1000 * d$Q / step$days,
      mrso = 1000 * root,
      gw = 1000 * gw,
      channel = 1000 * channel
    ),
    meta = list(
      parameters = as.list(p),
      parameter_source = list(
        all = "huagrahuma$parameters, topmodel 0.7.5",
        Srmax = m$srmax_source,
        dt = "the case's step, in hours"
      ),
      topographic_index = "huagrahuma$topidx, 16 classes, topmodel 0.7.5",
      delay = list(source = m$delay_source, distance_m = m$delay[, 1], cumulative_area = m$delay[, 2]),
      gw_datum_mm = 1000 * datum,
      initial_state_mm = list(root_zone = 1000 * root0, unsaturated = 0, saturated_deficit = 1000 * s0,
                              delay_function = 1000 * d$transit0),
      class_weight_sum = d$weight,
      checks = list(
        # What the core removed as excess against what it reported as saturation-excess runoff.
        excess_removed_mm = 1000 * sum(d$ex),
        saturation_runoff_reported_mm = 1000 * (sum(d$qo) - sum(d$fex)),
        infiltration_excess_mm = 1000 * sum(d$fex),
        rain_mm = 1000 * sum(rain),
        whole_record_budget_residual_mm = 1000 * residual
      )
    )
  )
}

# The areal mean of the index, as misc.c's get_lambda takes it from the class table sorted from
# high to low. Needed only to state the core's initial deficit for the budget in run.json.
lambda_of <- function(topidx) {
  o <- order(topidx[, 1], decreasing = TRUE)
  atb <- topidx[o, 1]
  a <- topidx[o, 2]
  i <- seq_along(atb)[-1]
  sum(a[i] * (atb[i] + atb[i - 1]) / 2)
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
  if (is.null(static$area_km2)) stop("static.json has no area_km2, so `dis` cannot be reported")

  run <- run_core(forcing, static, timestep)
  rows <- run$rows
  table <- data.frame(
    time = forcing$time,
    pr = forcing$pr,
    evspsbl = fmt(rows$evspsbl),
    mrro = fmt(rows$mrro),
    # mm/day over the catchment to m3/s.
    dis = fmt(rows$mrro * as.numeric(static$area_km2) / 86.4),
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
        mrso = "root zone: class-weighted Srmax - Srz",
        gw = "unsaturated store Suz plus the saturated zone, datum - S_mean",
        channel = "water inside the delay function"
      ),
      zero_stores = list(
        snw = "TOPMODEL has no snow module: all precipitation is treated as liquid",
        canopy = "TOPMODEL has no interception store: rain goes directly to the root zone"
      ),
      not_reported = list(
        snm = "TOPMODEL has no snow module, so nothing crosses a snow-module boundary",
        gwex = "TOPMODEL has no exchange across the catchment boundary: its saturated zone drains only to qs",
        energy_fluxes = "TOPMODEL has no energy balance: hfls, hfss, hfg and the radiation terms do not exist in it",
        stage = "TOPMODEL's delay function has no channel geometry beyond distance, so no water-surface elevation"
      )
    )
  )
  write(toJSON(meta, auto_unbox = TRUE, digits = NA, null = "null"), file.path(io, request$output$run))
  invisible(0L)
}

main()
