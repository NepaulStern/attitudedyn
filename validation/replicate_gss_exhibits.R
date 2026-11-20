#!/usr/bin/env Rscript
# GSS tables and the impulse-response figure, without a refit.
#
#   Rscript validation/replicate_gss_exhibits.R --threads 4
#   Rscript validation/replicate_gss_exhibits.R --fit validation/output/gss_refit.rds
#
# The default call evaluates the replication estimates (Philox, seed 1,
# 10 antithetic draws) and computes centred outer-product standard errors.
# It does not start the optimiser. --fit displays a saved refit instead
# of optimising it again.
#
# Tables and gss_impulse_response.png go to validation/output/. The figure
# is also copied to man/figures/ for the README. notes/GSS_EXHIBITS.md
# records the evaluation.

args <- commandArgs(trailingOnly = TRUE)
fit_flag <- match("--fit", args)
thread_flag <- match("--threads", args)
threads <- if (!is.na(thread_flag)) as.integer(args[[thread_flag + 1L]]) else 1L
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (file.exists("validation/gss_data.R")) {
  getwd()
} else if (length(file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."))
} else {
  getwd()
}
source(file.path(root, "validation/gss_data.R"))
suppressPackageStartupMessages(library(attitudedyn))

# Likelihood, scores, and post-estimation at a fixed coefficient vector.
fit_at_coefficients <- function(panel, coefficients, threads) {
  formula <- cbind(ideology, military, welfare, immigrants, party, environment) ~
    0 + `(Intercept)` + age + hs + college + adv + female + white + black +
    inc_mid + inc_high
  control <- attitudedyn_control(
    maxit = 0L,
    n_threads = threads,
    vcov = "opg_centered",
    use_gradient = TRUE,
    keep_scores = TRUE,
    seed = 1
  )
  attitudedyn(
    formula,
    data = panel,
    id = ~ id,
    wave = ~ wave,
    time_grid = gss_time_grid(),
    weights = weight,
    n_factors = 2L,
    draws = 10L,
    start = coefficients,
    control = control
  )
}

out_dir <- file.path(root, "validation", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
paths <- ensure_replication_files(root)
log_path <- paths$log
point <- "replication estimates (the X vector printed before iteration 8)"
started <- proc.time()
if (!is.na(fit_flag)) {
  fit_path <- args[[fit_flag + 1L]]
  if (!file.exists(fit_path)) {
    stop("No fitted object at ", fit_path, ".", call. = FALSE)
  }
  fit <- readRDS(fit_path)
  fit <- attitudedyn:::ensure_problem(fit)
  point <- paste0("saved fit in ", basename(fit_path))
} else {
  panel <- read_gss(paths$data)
  counts <- gss_category_counts(panel)
  expected <- gss_frequencies()
  for (name in names(expected)) {
    if (!identical(as.integer(counts[[name]]), expected[[name]])) {
      stop("Category counts for ", name, " do not match the replication log.", call. = FALSE)
    }
  }
  estimates <- read_gss_estimates(log_path)
  fit <- fit_at_coefficients(panel, unname(estimates), threads)
}
elapsed <- (proc.time() - started)[["elapsed"]]

comparison <- gss_comparison(fit, log_path)
coefficients <- comparison$coefficients
ratio <- abs(coefficients$estimate / coefficients$std.error)
coefficients$statistic <- coefficients$estimate / coefficients$std.error
coefficients$stars <- ifelse(
  !is.finite(ratio), "",
  ifelse(ratio > 2.576, "***",
    ifelse(ratio > 1.960, "**",
      ifelse(ratio > 1.645, "*", "")))
)
utils::write.csv(coefficients, file.path(out_dir, "gss_coefficients.csv"), row.names = FALSE)

persistence <- as.data.frame(comparison$persistence)
persistence <- cbind(attitude = rownames(comparison$persistence), persistence)
utils::write.csv(persistence, file.path(out_dir, "gss_persistence.csv"), row.names = FALSE)
published_persistence <- as.data.frame(comparison$published_persistence)
published_persistence <- cbind(
  attitude = rownames(comparison$published_persistence),
  published_persistence
)
utils::write.csv(
  published_persistence,
  file.path(out_dir, "gss_persistence_published.csv"),
  row.names = FALSE
)
utils::write.csv(comparison$variance, file.path(out_dir, "gss_variance.csv"), row.names = FALSE)

# Powers of the natural persistence matrix. The published path uses the
# replication rho, transformed onto (-1, 1), and the same recurrence.
horizon <- 8L
irf_from_matrix <- function(rho, names) {
  response <- attitudedyn:::attitudedyn_irf_cpp(
    as.numeric(t(rho)), nrow(rho), as.integer(horizon)
  )
  dimension <- nrow(rho)
  array <- array(response$response, dim = c(dimension, dimension, horizon + 1L))
  array <- aperm(array, c(3L, 2L, 1L))
  dimnames(array) <- list(
    horizon = paste0("h", seq.int(0, horizon)),
    to = names,
    from = names
  )
  list(response = array, spectral_radius = response$spectral_radius)
}
attitude_names <- rownames(comparison$persistence)
fitted_irf <- impulse_response(fit, horizon = horizon)
published_irf <- irf_from_matrix(comparison$published_persistence, attitude_names)
irf_rows <- vector("list", length(attitude_names)^2)
cursor <- 0L
for (receiver in attitude_names) {
  for (source in attitude_names) {
    cursor <- cursor + 1L
    irf_rows[[cursor]] <- data.frame(
      to = receiver,
      from = source,
      horizon = seq.int(0, horizon),
      estimate = as.numeric(fitted_irf$response[, receiver, source]),
      published = as.numeric(published_irf$response[, receiver, source]),
      stringsAsFactors = FALSE
    )
  }
}
irf_table <- do.call(rbind, irf_rows)
irf_table$difference <- irf_table$estimate - irf_table$published
utils::write.csv(irf_table, file.path(out_dir, "gss_impulse_response.csv"), row.names = FALSE)

png_path <- file.path(out_dir, "gss_impulse_response.png")
largest_gap <- max(abs(irf_table$difference), na.rm = TRUE)
show_replication <- is.finite(largest_gap) && largest_gap > 1e-3
grDevices::png(png_path, width = 1600, height = 1360, res = 140)
dimension <- length(attitude_names)
graphics::par(mfrow = c(dimension, dimension), mar = c(1.8, 2.2, 0.3, 0.3), oma = c(3.2, 3.4, 2.4, 0.6))
steps <- seq.int(0, horizon)
for (row in seq_along(attitude_names)) {
  for (column in seq_along(attitude_names)) {
    receiver <- attitude_names[[row]]
    source <- attitude_names[[column]]
    estimate <- fitted_irf$response[, receiver, source]
    reference <- published_irf$response[, receiver, source]
    graphics::plot(
      steps, estimate, type = "l", col = "#1f4e79", lwd = 2,
      ylim = range(c(estimate, if (show_replication) reference else estimate), na.rm = TRUE),
      xlab = "", ylab = "", axes = FALSE
    )
    graphics::axis(1, at = c(0, 4, 8), cex.axis = 0.7, padj = -1)
    if (column == 1L) graphics::axis(2, cex.axis = 0.7, las = 1)
    graphics::abline(h = 0, col = "#b0b0b0")
    if (show_replication) {
      graphics::lines(steps, reference, col = "#b85c38", lwd = 1.5, lty = 2)
    }
    if (row == dimension) graphics::mtext(source, side = 1, line = 1.5, cex = 0.65)
    if (column == 1L) graphics::mtext(receiver, side = 2, line = 2.6, cex = 0.65)
  }
}
graphics::mtext("GSS impulse responses", outer = TRUE, cex = 1.05)
graphics::mtext("Shock to", side = 1, outer = TRUE, line = 1.6, cex = 0.85)
graphics::mtext("Response of", side = 2, outer = TRUE, line = 1.8, cex = 0.85)
grDevices::dev.off()
# The README embeds this copy. validation/output/ is not part of the package.
figure_dir <- file.path(root, "man", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
file.copy(png_path, file.path(figure_dir, "gss_impulse_response.png"), overwrite = TRUE)

summary_text <- paste(utils::capture.output(summary(fit)), collapse = "\n")
writeLines(summary_text, file.path(out_dir, "gss_summary.txt"))

abs_coef <- abs(comparison$coefficients$difference)
abs_se <- abs(comparison$coefficients$se_difference)
abs_rho <- abs(comparison$persistence - comparison$published_persistence)
abs_var <- abs(comparison$variance$difference)
abs_irf <- abs(irf_table$difference)
finite_median <- function(x) stats::median(x[is.finite(x)])
finite_max <- function(x) max(x[is.finite(x)])
note <- c(
  "# GSS exhibits",
  "",
  paste0("Evaluation at the ", point, "."),
  "Draws are Philox, seed 1, with 10 antithetic pairs. The persistence matrix",
  "is on (-1, 1), and the impulse responses are its powers.",
  "",
  sprintf(
    "Log-likelihood %.10f. People in the likelihood: %d. Free coefficients: %d.",
    comparison$loglik, fit$n_persons, length(coef(fit))
  ),
  sprintf("Spectral radius %.6f.", fitted_irf$spectral_radius),
  sprintf("Threads: %d. Elapsed: %.1f seconds.", threads, elapsed),
  "",
  "Regenerate with `Rscript validation/replicate_gss_exhibits.R` or",
  "`python python/validation/replicate_gss_exhibits.py`.",
  "The General Social Survey file comes from the CC0 Harvard Dataverse deposit",
  "doi:10.7910/DVN/DUCHYS and is downloaded into `validation/fortran/vendor/`.",
  ""
)
writeLines(note, file.path(root, "notes", "GSS_EXHIBITS.md"))

summary_lines <- c(
  sprintf("loglik %.16f", comparison$loglik),
  sprintf("coefficient_abs_diff_max %.6g", finite_max(abs_coef)),
  sprintf("coefficient_abs_diff_median %.6g", finite_median(abs_coef)),
  sprintf("se_abs_diff_max %.6g", finite_max(abs_se)),
  sprintf("se_abs_diff_median %.6g", finite_median(abs_se)),
  sprintf("persistence_abs_diff_max %.6g", finite_max(abs_rho)),
  sprintf("persistence_abs_diff_median %.6g", finite_median(abs_rho)),
  sprintf("variance_abs_diff_max %.6g", finite_max(abs_var)),
  sprintf("variance_abs_diff_median %.6g", finite_median(abs_var)),
  sprintf("impulse_abs_diff_max %.6g", finite_max(abs_irf)),
  sprintf("impulse_abs_diff_median %.6g", finite_median(abs_irf)),
  sprintf("spectral_radius %.6f", fitted_irf$spectral_radius),
  sprintf("spectral_radius_published %.6f", published_irf$spectral_radius)
)
writeLines(summary_lines, file.path(out_dir, "gss_exhibit_comparison.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n", sep = "")
cat(sprintf("elapsed_seconds %.3f\n", elapsed))
cat(sprintf("threads %d\n", threads))
cat("Wrote exhibits to ", out_dir, "\n", sep = "")
