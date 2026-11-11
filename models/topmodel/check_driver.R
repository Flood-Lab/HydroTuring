# Run at image build: the driver must return what the package returns, to the last bit we can
# ask of it, on the package's own example at its own 15-minute step, and at the two steps the
# adapter uses. A driver that drifted from topmodel() would no longer be the model's code.
library(topmodel)
dyn.load("/opt/topmodel/ht_topmodel.so")
data("huagrahuma", package = "topmodel")
h <- huagrahuma

check <- function(parameters, rain, etp) {
  ref <- topmodel(parameters, h$topidx, h$delay, rain, etp, verbose = TRUE)
  n <- length(rain)
  d <- .C("ht_topmodel", as.double(parameters), as.double(as.matrix(h$topidx)), as.double(as.matrix(h$delay)),
          as.double(rain), as.double(etp), as.integer(nrow(h$topidx)), as.integer(n), as.integer(nrow(h$delay)),
          Q = double(n), qo = double(n), qs = double(n), S = double(n), fex = double(n), Ea = double(n),
          Srz = double(n), Suz = double(n), ex = double(n), transit0 = double(1), weight = double(1),
          PACKAGE = "ht_topmodel")
  for (v in c("Q", "qo", "qs", "S", "fex", "Ea")) {
    stopifnot(max(abs(as.numeric(ref[[v]]) - d[[v]])) < 1e-12)
  }
  cat(sprintf("driver matches topmodel(): dt=%g h, %d rows\n", parameters[["dt"]], n))
}

check(h$parameters, h$rain, h$ETp)
for (dt in c(1, 24)) {
  p <- h$parameters
  p[["dt"]] <- dt
  k <- as.integer(dt / 0.25)
  groups <- ceiling(seq_along(h$rain) / k)
  check(p, as.numeric(tapply(h$rain, groups, sum)), as.numeric(tapply(h$ETp, groups, sum)))
}
