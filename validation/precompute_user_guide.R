#!/usr/bin/env Rscript
# Write inst/extdata/guide_* for the user-guide vignette.
# Run from a checkout, after the package is installed:
#   Rscript validation/precompute_user_guide.R
suppressPackageStartupMessages(library(attitudedyn))
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (file.exists("validation/gss_data.R")) {
  getwd()
} else if (length(file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."))
} else {
  getwd()
}
source(file.path(root, "validation/gss_data.R"))
out <- file.path(root, "inst/extdata")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
paths <- ensure_replication_files(root)

panel <- read_gss(paths$data)
estimates <- read_gss_estimates(paths$log)
counts <- gss_category_counts(panel)
stopifnot(identical(lapply(counts, as.integer), gss_frequencies()))

included <- panel$weight > 0
outcomes <- c("ideology", "military", "welfare", "immigrants", "party", "environment")
n_missing <- sapply(outcomes, function(nm) sum(is.na(panel[[nm]][included])))
waves_per_person <- table(panel$id[included])
gap_rows <- included & panel$`(Intercept)` == 0
# A person with an interior gap and a person with a missing item.
gap_id <- panel$id[which(gap_rows)[1]]
miss_id <- NULL
for (nm in outcomes) {
  hit <- which(included & is.na(panel[[nm]]))
  if (length(hit)) {
    miss_id <- panel$id[hit[1]]
    break
  }
}
zero_id <- panel$id[which(!included)[1]]
keep_ids <- unique(c(panel$id[1], gap_id, miss_id, zero_id))
excerpt <- panel[panel$id %in% keep_ids, ]
rownames(excerpt) <- NULL

layout_info <- list(
  n_rows_file = 21648L,
  n_people_file = attr(panel, "n_people"),
  n_included = attr(panel, "n_included"),
  n_rows_panel = nrow(panel),
  n_rows_included = sum(included),
  n_gap_rows = sum(gap_rows),
  n_missing = n_missing,
  wave_table = as.integer(table(panel$wave)),
  waves_per_person = as.integer(table(waves_per_person)),
  names_waves_per_person = as.integer(names(table(waves_per_person))),
  time_grid = gss_time_grid(),
  category_counts = counts
)
saveRDS(excerpt, file.path(out, "guide_excerpt.rds"), compress = "gzip")
saveRDS(layout_info, file.path(out, "guide_layout.rds"), compress = "gzip")
cat("layout people", layout_info$n_people_file, "included", layout_info$n_included,
    "panel rows", layout_info$n_rows_panel, "gap rows", layout_info$n_gap_rows, "\n")
cat("missing", paste(names(n_missing), n_missing, sep = "=", collapse = " "), "\n")
cat("waves per person", paste(layout_info$names_waves_per_person, layout_info$waves_per_person, sep = ":", collapse = " "), "\n")
cat("excerpt ids", paste(keep_ids, collapse = ","), "rows", nrow(excerpt), "\n")
print(excerpt)

started <- proc.time()
fit <- attitudedyn(
  cbind(ideology, military, welfare, immigrants, party, environment) ~
    0 + `(Intercept)` + age + hs + college + adv + female + white + black +
    inc_mid + inc_high,
  data = panel, id = ~ id, wave = ~ wave, time_grid = gss_time_grid(),
  weights = weight, n_factors = 2, draws = 10, start = unname(estimates),
  control = attitudedyn_control(maxit = 0, n_threads = 4, vcov = "opg_centered", seed = 1)
)
elapsed <- proc.time() - started
cat("FIT seconds", elapsed[["elapsed"]], "loglik", sprintf("%.16f", fit$loglik), "\n")

writeLines(capture.output(summary(fit)), file.path(out, "guide_summary.txt"))
writeLines(capture.output(glance(fit)), file.path(out, "guide_glance.txt"))
writeLines(capture.output(utils::head(tidy(fit), 8)), file.path(out, "guide_tidy_head.txt"))
writeLines(capture.output(persistence_matrix(fit)), file.path(out, "guide_persistence.txt"))
writeLines(capture.output(variance_decomposition(fit)), file.path(out, "guide_variance.txt"))
writeLines(capture.output(cutpoints(fit)), file.path(out, "guide_cutpoints.txt"))
vc <- vcov(fit)[1:3, 1:3]
writeLines(capture.output(print(vc, digits = 4)), file.path(out, "guide_vcov.txt"))

irf <- impulse_response(fit, horizon = 8)
saveRDS(list(response = irf$response, spectral_radius = irf$spectral_radius),
        file.path(out, "guide_irf.rds"), compress = "gzip")
cat("spectral", sprintf("%.16f", irf$spectral_radius), "\n")

fixed <- fit$layout[fit$layout$is_free == 0L, c("name", "block")]
cat("fixed count", nrow(fixed), "\n")
print(fixed)

cat("PREDICT start\n")
t0 <- proc.time()
pred_class <- predict(fit, type = "class", n_draws = 2, seed = 2)
pred_prob <- predict(fit, type = "prob", n_draws = 2, seed = 2)
cat("PREDICT seconds", (proc.time() - t0)[["elapsed"]], "class rows", nrow(pred_class), "\n")
writeLines(capture.output(utils::head(pred_class, 8)), file.path(out, "guide_predict_class.txt"))
writeLines(capture.output(utils::head(pred_prob$ideology, 6)), file.path(out, "guide_predict_prob.txt"))

cat("SIMULATE start\n")
t0 <- proc.time()
sim <- simulate(fit, nsim = 1, seed = 2)
cat("SIMULATE seconds", (proc.time() - t0)[["elapsed"]], "rows", nrow(sim), "\n")
writeLines(capture.output(utils::head(sim, 6)), file.path(out, "guide_simulate.txt"))

path <- tempfile(fileext = ".rds")
t0 <- proc.time()
saveRDS(fit, path)
fit2 <- readRDS(path)
same_loglik <- identical(fit2$loglik, fit$loglik)
same_pers <- isTRUE(all.equal(persistence_matrix(fit2), persistence_matrix(fit)))
cat("RELOAD seconds", (proc.time() - t0)[["elapsed"]], "bytes", file.info(path)$size,
    "same", same_loglik, same_pers, "\n")
reload_lines <- c(
  capture.output(print(fit2)),
  paste0("persistence identical after reload: ", same_pers)
)
writeLines(reload_lines, file.path(out, "guide_reload.txt"))

cat("DONE\n")
