#!/usr/bin/env Rscript
# Evaluate the GSS model at the replication estimates.
#
# From the repository root, with the package installed:
#
#   Rscript validation/run_gss.R
#   Rscript validation/run_gss.R --threads 20
#   Rscript validation/run_gss.R --refit --threads 20
#
# The default call evaluates the shipped free coefficients with the
# package's Philox draws and does not start the optimiser. --refit
# continues from those values. Category counts are checked against the
# replication log before the likelihood is called.

args <- commandArgs(trailingOnly = TRUE)
threads <- 1L
refit <- FALSE
if ("--refit" %in% args) refit <- TRUE
thread_flag <- match("--threads", args)
if (!is.na(thread_flag)) threads <- as.integer(args[[thread_flag + 1L]])

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (file.exists("validation/gss_data.R")) {
  getwd()
} else if (length(file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."))
} else {
  getwd()
}
source(file.path(root, "validation/gss_data.R"))

data_path <- file.path(root, "validation/fortran/data/gss2.prn")
log_path <- file.path(root, "validation/fortran/src/gss_log.txt")
panel <- read_gss(data_path)
counts <- gss_category_counts(panel)
expected <- gss_frequencies()
for (name in names(expected)) {
  if (!identical(as.integer(counts[[name]]), expected[[name]])) {
    stop("Category counts for ", name, " do not match the replication log: ",
         paste(counts[[name]], collapse = " "), call. = FALSE)
  }
}
message(
  "People in the file: ", attr(panel, "n_people"),
  ". People in the likelihood: ", attr(panel, "n_included"), "."
)
if (!identical(attr(panel, "n_included"), 11252L)) {
  stop("Expected 11252 people in the likelihood.", call. = FALSE)
}

estimates <- read_gss_estimates(log_path)
if (length(estimates) != 164L) {
  stop("Expected 164 free coefficients, found ", length(estimates), ".", call. = FALSE)
}

suppressPackageStartupMessages(library(attitudedyn))
formula <- cbind(ideology, military, welfare, immigrants, party, environment) ~
  0 + `(Intercept)` + age + hs + college + adv + female + white + black +
  inc_mid + inc_high
control <- attitudedyn_control(
  maxit = if (refit) 200L else 0L,
  n_threads = threads,
  vcov = "none",
  use_gradient = FALSE,
  keep_scores = FALSE,
  seed = 1
)
fit <- attitudedyn(
  formula,
  data = panel,
  id = ~ id,
  wave = ~ wave,
  time_grid = gss_time_grid(),
  weights = weight,
  n_factors = 2L,
  draws = 10L,
  start = unname(estimates),
  control = control
)
cat(sprintf("free_coefficients %d\n", length(coef(fit))))
cat(sprintf("persons %d\n", fit$n_persons))
cat(sprintf("loglik %.16f\n", fit$loglik))
cat(sprintf("published_target %.2f\n", -111389.88))
if (refit) cat(sprintf("convergence %d\n", fit$convergence))
invisible(fit)
