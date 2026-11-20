#!/usr/bin/env Rscript

# Parameter recovery for a panel with the GSS wave dates.
#
# From the repository root, with the package installed:
#
#   Rscript validation/simulation_recovery.R
#   Rscript validation/simulation_recovery.R --threads 4 --replications 40
#
# Each replication draws a panel from known coefficients, fits it, and
# records the signed error and whether the 95 percent Wald interval covers
# the truth. The note written under notes/ summarises those replications.

args <- commandArgs(trailingOnly = TRUE)
number_of_threads <- 4L
number_of_replications <- 40L
thread_flag <- match("--threads", args)
if (!is.na(thread_flag)) number_of_threads <- as.integer(args[[thread_flag + 1L]])
replication_flag <- match("--replications", args)
if (!is.na(replication_flag)) number_of_replications <- as.integer(args[[replication_flag + 1L]])

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (file.exists("validation/simulation_recovery.R")) {
  getwd()
} else if (length(file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."))
} else {
  getwd()
}

suppressPackageStartupMessages(library(attitudedyn))

# The design is a short GSS-shaped panel: eight waves on the GSS dates,
# several ordinal attitudes, and one covariate besides the intercept.
# Forty replications at this size finish in a few minutes on four threads.
recovery_design <- function(number_of_replications, number_of_threads) {
  list(
    number_of_persons = 400L,
    number_of_waves = 8L,
    number_of_attitudes = 3L,
    number_of_factors = 1L,
    number_of_categories = c(5L, 3L, 3L),
    wave_dates = c(2, 3, 4, 5, 6, 7, 8, 10),
    number_of_draws = 8L,
    maximum_iterations = 100L,
    relative_tolerance = 1e-8,
    number_of_replications = as.integer(number_of_replications),
    number_of_threads = as.integer(number_of_threads),
    starting_seed = 20261120L
  )
}

# Map an unconstrained persistence coefficient onto (-1, 1).
natural_persistence <- function(unconstrained_coefficient) {
  2 * exp(unconstrained_coefficient) / (1 + exp(unconstrained_coefficient)) - 1
}

# Delta-method standard error. The derivative of the logistic map is
# 0.5 * (1 - rho^2).
persistence_standard_error <- function(natural_value, unconstrained_standard_error) {
  abs(0.5 * (1 - natural_value^2)) * unconstrained_standard_error
}

outcome_formula <- function(number_of_attitudes) {
  outcome_names <- paste0("y", seq_len(number_of_attitudes))
  outcome_call <- parse(text = paste0("cbind(", paste(outcome_names, collapse = ", "), ")"))[[1L]]
  stats::reformulate("x", response = outcome_call)
}

simulate_recovery_panel <- function(design, replication_index) {
  simulate_attitudes(
    n_persons = design$number_of_persons,
    n_waves = design$number_of_waves,
    n_attitudes = design$number_of_attitudes,
    n_factors = design$number_of_factors,
    n_categories = design$number_of_categories,
    seed = design$starting_seed + replication_index,
    wave_time = design$wave_dates,
    include_wave_effects = FALSE
  )
}

fit_recovery_replication <- function(panel, design) {
  attitudedyn(
    outcome_formula(design$number_of_attitudes),
    data = panel,
    id = ~ id,
    wave = ~ wave,
    time = ~ time,
    n_factors = design$number_of_factors,
    draws = design$number_of_draws,
    include_wave_effects = FALSE,
    control = attitudedyn_control(
      maxit = design$maximum_iterations,
      reltol = design$relative_tolerance,
      n_threads = design$number_of_threads,
      vcov = "opg",
      seed = 1L,
      trace = 0L
    )
  )
}

# One row per free coefficient. Coverage uses the unconstrained Wald interval.
replication_rows <- function(fit, truth, replication_index) {
  estimate <- coef(fit)
  standard_error <- sqrt(diag(vcov(fit)))
  truth_free <- truth[names(estimate)]
  critical <- stats::qnorm(0.975)
  lower <- estimate - critical * standard_error
  upper <- estimate + critical * standard_error
  persistence <- grepl("^persistence\\[", names(estimate))
  natural_truth <- natural_persistence(truth_free)
  natural_estimate <- natural_persistence(estimate)
  natural_standard_error <- persistence_standard_error(natural_estimate, standard_error)
  natural_lower <- natural_estimate - critical * natural_standard_error
  natural_upper <- natural_estimate + critical * natural_standard_error
  data.frame(
    replication = replication_index,
    parameter = names(estimate),
    truth = as.numeric(truth_free),
    estimate = as.numeric(estimate),
    standard_error = as.numeric(standard_error),
    covered = as.numeric(truth_free) >= lower & as.numeric(truth_free) <= upper,
    natural_truth = ifelse(persistence, natural_truth, NA_real_),
    natural_estimate = ifelse(persistence, natural_estimate, NA_real_),
    natural_covered = ifelse(
      persistence,
      natural_truth >= natural_lower & natural_truth <= natural_upper,
      NA
    ),
    converged = fit$convergence == 0L,
    stringsAsFactors = FALSE
  )
}

summarize_recovery <- function(detail) {
  parameters <- unique(detail$parameter)
  rows <- lapply(parameters, function(parameter_name) {
    block <- detail[detail$parameter == parameter_name, , drop = FALSE]
    error <- block$estimate - block$truth
    data.frame(
      parameter = parameter_name,
      truth = block$truth[[1L]],
      mean_estimate = mean(block$estimate),
      bias = mean(error),
      rmse = sqrt(mean(error^2)),
      coverage = mean(block$covered),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

summarize_natural_persistence <- function(detail) {
  block <- detail[grepl("^persistence\\[", detail$parameter), , drop = FALSE]
  parameters <- unique(block$parameter)
  rows <- lapply(parameters, function(parameter_name) {
    slice <- block[block$parameter == parameter_name, , drop = FALSE]
    error <- slice$natural_estimate - slice$natural_truth
    data.frame(
      parameter = parameter_name,
      truth = slice$natural_truth[[1L]],
      mean_estimate = mean(slice$natural_estimate),
      bias = mean(error),
      rmse = sqrt(mean(error^2)),
      coverage = mean(slice$natural_covered),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

format_recovery_table <- function(summary_frame) {
  lines <- c(
    "| Parameter | Truth | Mean estimate | Bias | RMSE | 95% coverage |",
    "| --- | ---: | ---: | ---: | ---: | ---: |"
  )
  formatted <- apply(summary_frame, 1L, function(row) {
    sprintf(
      "| %s | %.4f | %.4f | %.4f | %.4f | %.3f |",
      row[["parameter"]],
      as.numeric(row[["truth"]]),
      as.numeric(row[["mean_estimate"]]),
      as.numeric(row[["bias"]]),
      as.numeric(row[["rmse"]]),
      as.numeric(row[["coverage"]])
    )
  })
  paste(c(lines, formatted), collapse = "\n")
}

write_recovery_note <- function(path, design, detail, elapsed_seconds) {
  unconstrained <- summarize_recovery(detail)
  natural <- summarize_natural_persistence(detail)
  error <- detail$estimate - detail$truth
  natural_error <- detail$natural_estimate - detail$natural_truth
  natural_error <- natural_error[is.finite(natural_error)]
  note <- c(
    "# Parameter recovery",
    "",
    "Panels are drawn from `simulate_attitudes()` and refit with `attitudedyn()`.",
    "The wave dates are the GSS grid `c(2, 3, 4, 5, 6, 7, 8, 10)`.",
    "There are three ordinal attitudes, with 5, 3, and 3 categories, one",
    "standard normal covariate, an intercept, and one person factor.",
    "Wave effects are off, so each intercept is a free coefficient.",
    "",
    "The data-generating values are the ones `simulate_attitudes()` uses:",
    "diagonal persistence 0.45, off-diagonal persistence 0, slope 0.5,",
    "intercept 0.2, loading 0.5, and threshold increments `log(0.6)`.",
    "The slope, intercept, loading, and threshold rows below are on the",
    "unconstrained scale, which is the scale of `coef()` and of the Wald",
    "interval. Persistence is shown on both scales. The natural-scale",
    "interval uses the delta method, with standard error",
    "`|0.5 * (1 - rho^2)|` times the unconstrained standard error.",
    "",
    "## Design",
    "",
    "| Setting | Value |",
    "| --- | ---: |",
    sprintf("| Persons | %d |", design$number_of_persons),
    sprintf("| Waves | %d |", design$number_of_waves),
    sprintf("| Wave dates | %s |", paste(design$wave_dates, collapse = ", ")),
    sprintf("| Attitudes | %d |", design$number_of_attitudes),
    sprintf("| Categories | %s |", paste(design$number_of_categories, collapse = ", ")),
    sprintf("| Factors | %d |", design$number_of_factors),
    sprintf("| Draws (antithetic pairs) | %d |", design$number_of_draws),
    sprintf("| Maximum BFGS iterations | %d |", design$maximum_iterations),
    sprintf("| Replications | %d |", design$number_of_replications),
    sprintf("| Threads | %d |", design$number_of_threads),
    sprintf("| First simulation seed | %d |", design$starting_seed + 1L),
    "",
    "## Summary",
    "",
    sprintf("- Replications: %d", design$number_of_replications),
    sprintf(
      "- Converged: %d of %d",
      sum(tapply(detail$converged, detail$replication, function(flag) flag[[1L]])),
      design$number_of_replications
    ),
    sprintf("- Free coefficients: %d", length(unique(detail$parameter))),
    sprintf("- Mean signed bias, unconstrained: %.4f", mean(error)),
    sprintf("- RMSE, unconstrained: %.4f", sqrt(mean(error^2))),
    sprintf("- Coverage of unconstrained 95%% intervals: %.3f", mean(detail$covered)),
    sprintf("- Mean signed bias, natural persistence: %.4f", mean(natural_error)),
    sprintf("- RMSE, natural persistence: %.4f", sqrt(mean(natural_error^2))),
    sprintf(
      "- Coverage of natural-persistence 95%% intervals: %.3f",
      mean(detail$natural_covered, na.rm = TRUE)
    ),
    sprintf("- Elapsed seconds: %.1f", elapsed_seconds),
    "",
    "Bias is the mean of the estimate minus the data-generating value.",
    "RMSE is the square root of the mean squared error, pooled over",
    "coefficients and replications in the summary and computed within each",
    "coefficient in the tables. Coverage is the share of replications whose",
    "interval contains the data-generating value. Intervals are outer-product",
    "Wald intervals.",
    "",
    "## Unconstrained coefficients",
    "",
    format_recovery_table(unconstrained),
    "",
    "## Persistence on (-1, 1)",
    "",
    format_recovery_table(natural),
    "",
    "## Reproduce",
    "",
    "```bash",
    "Rscript validation/simulation_recovery.R --threads 4",
    "python python/validation/simulation_recovery.py --threads 4",
    "```",
    "",
    "The table above was written by `validation/simulation_recovery.R`.",
    "`python/validation/simulation_recovery.py` repeats the same design.",
    "A one-replication check of the bias, RMSE, and coverage calculations",
    "lives in `tests/testthat/test-recovery.R` and `python/tests/test_recovery.py`.",
    ""
  )
  writeLines(note, path)
}

design <- recovery_design(number_of_replications, number_of_threads)
message(
  "Recovery: ", design$number_of_replications, " replications, ",
  design$number_of_persons, " persons, ", design$number_of_threads, " threads."
)
started <- proc.time()
pieces <- vector("list", design$number_of_replications)
for (replication_index in seq_len(design$number_of_replications)) {
  panel <- simulate_recovery_panel(design, replication_index)
  fit <- fit_recovery_replication(panel, design)
  pieces[[replication_index]] <- replication_rows(
    fit, attr(panel, "theta_full"), replication_index
  )
  message(sprintf(
    "  %d/%d  converged=%d  loglik=%.2f",
    replication_index, design$number_of_replications,
    fit$convergence == 0L, fit$loglik
  ))
}
detail <- do.call(rbind, pieces)
elapsed_seconds <- (proc.time() - started)[["elapsed"]]
note_path <- file.path(root, "notes", "SIMULATION_RECOVERY.md")
write_recovery_note(note_path, design, detail, elapsed_seconds)
error <- detail$estimate - detail$truth
message(sprintf(
  "Bias %.4f  RMSE %.4f  coverage %.3f  elapsed %.1fs  wrote %s",
  mean(error), sqrt(mean(error^2)), mean(detail$covered), elapsed_seconds, note_path
))
