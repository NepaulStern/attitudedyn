#!/usr/bin/env Rscript
# Exhibits for the GSS application, from a fitted attitudedyn object.
#
#   Rscript validation/run_gss.R --threads 20 --refit
#   Rscript validation/replicate_gss_exhibits.R
#   Rscript validation/replicate_gss_exhibits.R --fit validation/output/gss_refit.rds
#
# Writes coefficient, persistence, variance, and impulse-response tables,
# plus an impulse-response plot, under validation/output/. That directory
# is gitignored. Each table is compared with the replication log.

args <- commandArgs(trailingOnly = TRUE)
fit_flag <- match("--fit", args)
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (file.exists("validation/gss_data.R")) {
  getwd()
} else if (length(file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."))
} else {
  getwd()
}
source(file.path(root, "validation/gss_data.R"))

fit_path <- if (!is.na(fit_flag)) {
  args[[fit_flag + 1L]]
} else {
  file.path(root, "validation", "output", "gss_refit.rds")
}
if (!file.exists(fit_path)) {
  stop(
    "No fitted object at ", fit_path, ". Run validation/run_gss.R --refit first.",
    call. = FALSE
  )
}

suppressPackageStartupMessages(library(attitudedyn))
fit <- readRDS(fit_path)
fit <- attitudedyn:::ensure_problem(fit)
log_path <- file.path(root, "validation/fortran/src/gss_log.txt")
comparison <- gss_comparison(fit, log_path)

out_dir <- file.path(root, "validation", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

coefficients <- comparison$coefficients
coefficients$stars <- ifelse(
  !is.finite(coefficients$std.error), "",
  ifelse(abs(coefficients$estimate / coefficients$std.error) > 2.576, "***",
    ifelse(abs(coefficients$estimate / coefficients$std.error) > 1.960, "**",
      ifelse(abs(coefficients$estimate / coefficients$std.error) > 1.645, "*", "")))
)
utils::write.csv(coefficients, file.path(out_dir, "gss_coefficients.csv"), row.names = FALSE)

persistence <- as.data.frame(comparison$persistence)
persistence <- cbind(attitude = rownames(comparison$persistence), persistence)
utils::write.csv(persistence, file.path(out_dir, "gss_persistence.csv"), row.names = FALSE)
published_persistence <- as.data.frame(comparison$published_persistence)
published_persistence <- cbind(attitude = rownames(comparison$published_persistence), published_persistence)
utils::write.csv(
  published_persistence,
  file.path(out_dir, "gss_persistence_published.csv"),
  row.names = FALSE
)
utils::write.csv(comparison$variance, file.path(out_dir, "gss_variance.csv"), row.names = FALSE)

horizon <- 8L
irf_array <- function(rho, names) {
  response <- attitudedyn:::attitudedyn_irf_cpp(
    as.numeric(t(rho)), nrow(rho), as.integer(horizon)
  )
  J <- nrow(rho)
  array <- array(response$response, dim = c(J, J, horizon + 1L))
  array <- aperm(array, c(3L, 2L, 1L))
  dimnames(array) <- list(
    horizon = paste0("h", seq.int(0, horizon)),
    to = names,
    from = names
  )
  list(response = array, spectral_radius = response$spectral_radius)
}
names <- rownames(comparison$persistence)
fitted_irf <- impulse_response(fit, horizon = horizon)
published_irf <- irf_array(comparison$published_persistence, names)
irf_rows <- list()
for (to in names) {
  for (from in names) {
    irf_rows[[length(irf_rows) + 1L]] <- data.frame(
      to = to,
      from = from,
      horizon = seq.int(0, horizon),
      estimate = as.numeric(fitted_irf$response[, to, from]),
      published = as.numeric(published_irf$response[, to, from]),
      stringsAsFactors = FALSE
    )
  }
}
irf_table <- do.call(rbind, irf_rows)
irf_table$difference <- irf_table$estimate - irf_table$published
utils::write.csv(irf_table, file.path(out_dir, "gss_impulse_response.csv"), row.names = FALSE)

png_path <- file.path(out_dir, "gss_impulse_response.png")
grDevices::png(png_path, width = 1600, height = 1400, res = 140)
J <- length(names)
graphics::par(mfrow = c(J, J), mar = c(1.6, 1.6, 1.2, 0.4), oma = c(1.5, 1.5, 2.2, 0.4))
h <- seq.int(0, horizon)
for (to in names) {
  for (from in names) {
    y <- range(fitted_irf$response[, to, from], published_irf$response[, to, from])
    graphics::plot(
      h, fitted_irf$response[, to, from], type = "l", col = "#1f4e79", lwd = 2,
      ylim = y, xlab = "", ylab = "", main = paste0(from, " -> ", to), cex.main = 0.7
    )
    graphics::lines(h, published_irf$response[, to, from], col = "#b85c38", lwd = 1.5, lty = 2)
    graphics::abline(h = 0, col = "grey70")
  }
}
graphics::mtext("Impulse response of the persistence matrix", outer = TRUE, cex = 1)
graphics::mtext("solid: refit    dashed: replication log", outer = TRUE, side = 1, cex = 0.8)
grDevices::dev.off()

abs_coef <- abs(comparison$coefficients$difference)
abs_rho <- abs(comparison$persistence - comparison$published_persistence)
abs_var <- abs(comparison$variance$difference)
abs_irf <- abs(irf_table$difference)
summary_lines <- c(
  sprintf("fit %s", fit_path),
  sprintf("loglik %.16f", comparison$loglik),
  sprintf("convergence %d", comparison$convergence),
  sprintf("coefficient_abs_diff_max %.6g", max(abs_coef)),
  sprintf("coefficient_abs_diff_median %.6g", stats::median(abs_coef)),
  sprintf("persistence_abs_diff_max %.6g", max(abs_rho)),
  sprintf("persistence_abs_diff_median %.6g", stats::median(abs_rho)),
  sprintf("variance_abs_diff_max %.6g", max(abs_var, na.rm = TRUE)),
  sprintf("variance_abs_diff_median %.6g", stats::median(abs_var, na.rm = TRUE)),
  sprintf("impulse_abs_diff_max %.6g", max(abs_irf, na.rm = TRUE)),
  sprintf("impulse_abs_diff_median %.6g", stats::median(abs_irf, na.rm = TRUE)),
  sprintf("spectral_radius_refit %.6f", fitted_irf$spectral_radius),
  sprintf("spectral_radius_published %.6f", published_irf$spectral_radius)
)
writeLines(summary_lines, file.path(out_dir, "gss_exhibit_comparison.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n", sep = "")
cat("Wrote exhibits to ", out_dir, "\n", sep = "")
