simulation_frame <- function(panel, attitude_names, covariate_names) {
  n_attitudes <- panel$n_attitudes
  n_covariates <- panel$n_covariates
  n_rows <- length(panel$id)
  outcomes <- matrix(panel$outcome, n_rows, n_attitudes, byrow = TRUE)
  colnames(outcomes) <- attitude_names
  covariates <- matrix(panel$covariates, n_rows, n_covariates, byrow = TRUE)
  # Drop the intercept column from the returned data. The formula adds it back.
  keep <- covariate_names != "(Intercept)"
  covariates <- covariates[, keep, drop = FALSE]
  colnames(covariates) <- covariate_names[keep]
  data.frame(id = panel$id, wave = panel$wave, time = panel$time, covariates, outcomes,
             check.names = FALSE)
}

#' Simulate a panel from the structural model
#'
#' Draws covariates and ordinal outcomes from the latent-attitude recursion
#' in the replication likelihood. The returned data frame is ready for
#' [attitudedyn()]. Category 1 is the lowest category.
#'
#' @param n_persons Number of people.
#' @param n_waves Number of equally spaced waves, unless `wave_time` is set.
#' @param n_attitudes Number of ordinal outcomes, named `y1`, `y2`, ...
#' @param n_factors Number of person factors.
#' @param n_categories Categories per attitude. A single integer is recycled.
#' @param seed Seed for the counter-based generator.
#' @param wave_time Optional wave dates. The default is `1, 2, ..., n_waves`.
#' @param include_wave_effects If `TRUE`, the data-generating process includes
#'   wave effects. The default is `FALSE`, which matches the formula
#'   `cbind(y1, y2) ~ x` used in the examples.
#' @return A data frame with columns `id`, `wave`, `time`, `x`, and the outcomes.
#' @export
#'
#' @examples
#' simulate_attitudes(n_persons = 5, n_waves = 3, seed = 1)
simulate_attitudes <- function(n_persons = 100, n_waves = 4, n_attitudes = 2,
                               n_factors = 1, n_categories = 3,
                               seed = 1, wave_time = NULL,
                               include_wave_effects = FALSE) {
  if (length(n_categories) == 1L) n_categories <- rep(as.integer(n_categories), n_attitudes)
  if (is.null(wave_time)) wave_time <- as.numeric(seq_len(n_waves))
  attitude_names <- paste0("y", seq_len(n_attitudes))
  covariate_names <- c("(Intercept)", "x")
  layout <- attitudedyn_layout_cpp(
    n_waves = as.integer(n_waves), n_attitudes = as.integer(n_attitudes),
    n_covariates = 2L, n_factors = as.integer(n_factors),
    n_categories = as.integer(n_categories), wave_time = as.numeric(wave_time),
    include_wave_effects = isTRUE(include_wave_effects),
    free_contemporaneous = FALSE, fix_last_loading = n_factors >= 2,
    fix_intercepts = FALSE, intercept_columns = 0L,
    empty_wave_effect = integer(0),
    attitude_names = attitude_names, covariate_names = covariate_names,
    wave_names = as.character(seq_len(n_waves))
  )
  layout_df <- data.frame(name = layout$name, block = layout$block,
                          attitude = layout$attitude, index2 = layout$index2,
                          is_free = layout$is_free, stringsAsFactors = FALSE)
  theta <- default_theta(layout_df)
  slope <- layout_df$block == "beta" & layout_df$index2 == 2L
  intercept <- layout_df$block == "beta" & layout_df$index2 == 1L
  theta[slope] <- 0.5
  theta[intercept] <- 0.2
  diagonal <- layout_df$block == "persistence" & layout_df$attitude == layout_df$index2
  theta[diagonal] <- log((1 + 0.45) / (1 - 0.45))
  panel <- attitudedyn_simulate_cpp(
    n_waves = as.integer(n_waves), n_attitudes = as.integer(n_attitudes),
    n_covariates = 2L, n_factors = as.integer(n_factors),
    n_categories = as.integer(n_categories), wave_time = as.numeric(wave_time),
    include_wave_effects = isTRUE(include_wave_effects),
    free_contemporaneous = FALSE, fix_last_loading = n_factors >= 2,
    fix_intercepts = FALSE, intercept_columns = 0L, empty_wave_effect = integer(0),
    attitude_names = attitude_names, covariate_names = covariate_names,
    wave_names = as.character(seq_len(n_waves)),
    theta_full = theta, n_persons = as.integer(n_persons),
    seed = as.numeric(seed)
  )
  frame <- simulation_frame(panel, attitude_names, covariate_names)
  attr(frame, "theta_full") <- setNames(theta, layout$name)
  frame
}
