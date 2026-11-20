#' Persistence matrix on the natural scale
#'
#' Row `to` and column `from` is the coefficient that multiplies attitude
#' `from` in the recursion for attitude `to`, before it is raised to the gap
#' between waves. Values lie in (-1, 1).
#'
#' @param object A fitted [attitudedyn()] model.
#' @return A square numeric matrix.
#' @export
#'
#' @examples
#' panel <- simulate_attitudes(n_persons = 25, n_waves = 3, n_factors = 1, seed = 1)
#' fit <- attitudedyn(
#'   cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave, time = ~ time,
#'   n_factors = 1, draws = 3,
#'   control = attitudedyn_control(maxit = 0, n_threads = 1, vcov = "none")
#' )
#' persistence_matrix(fit)
persistence_matrix <- function(object) {
  object <- ensure_problem(object)
  natural <- attitudedyn_natural_cpp(object$cache$problem, object$theta_full)
  J <- object$inputs$n_attitudes
  matrix(natural$persistence, J, J, byrow = TRUE,
         dimnames = list(object$attitude_names, object$attitude_names))
}

#' Ordered cutpoints
#'
#' @param object A fitted [attitudedyn()] model.
#' @return A list of numeric vectors, one per attitude. Each vector contains
#'   the `K + 1` cutpoints, with the normalised values -100, 0 and 100 at the
#'   ends and between category 1 and category 2.
#' @examples
#' panel <- simulate_attitudes(n_persons = 20, n_waves = 3, n_factors = 1, seed = 1)
#' fit <- attitudedyn(
#'   cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave, time = ~ time,
#'   n_factors = 1, draws = 2,
#'   control = attitudedyn_control(maxit = 0, n_threads = 1, vcov = "none")
#' )
#' cutpoints(fit)
#' @export
cutpoints <- function(object) {
  object <- ensure_problem(object)
  natural <- attitudedyn_natural_cpp(object$cache$problem, object$theta_full)
  names(natural$cutpoints) <- object$attitude_names
  natural$cutpoints
}

#' Impulse responses of the persistence matrix
#'
#' The response at horizon 0 is the identity. Later horizons are powers of the
#' natural persistence matrix, which is the recursion over a gap of one time
#' unit. The path is that of a one-time innovation to the latent attitude.
#'
#' @param object A fitted model, or use [persistence_matrix()] yourself.
#' @param horizon Number of steps after the impact.
#' @return A list with the array `response` (horizon + 1 by attitude by
#'   attitude, indexed as receiving attitude then source attitude) and the
#'   spectral radius.
#' @examples
#' panel <- simulate_attitudes(n_persons = 20, n_waves = 3, n_factors = 1, seed = 1)
#' fit <- attitudedyn(
#'   cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave, time = ~ time,
#'   n_factors = 1, draws = 2,
#'   control = attitudedyn_control(maxit = 0, n_threads = 1, vcov = "none")
#' )
#' impulse_response(fit, horizon = 2)
#' @export
impulse_response <- function(object, horizon = 8L) {
  matrix <- persistence_matrix(object)
  response <- attitudedyn_irf_cpp(as.numeric(t(matrix)), nrow(matrix), as.integer(horizon))
  # t(matrix) because as.numeric is column-major and the core wants row-major.
  # persistence_matrix is already row = to, column = from, and as.numeric of
  # that matrix is column-major. The C++ routine wants row-major, which is
  # as.numeric(t(matrix)).
  J <- nrow(matrix)
  # The core stores each horizon row-major with the source attitude fastest
  # and the receiving attitude next, so the first R index is the source.
  array <- array(response$response, dim = c(J, J, horizon + 1L))
  array <- aperm(array, c(3L, 2L, 1L))
  dimnames(array) <- list(
    horizon = paste0("h", seq.int(0, horizon)),
    to = object$attitude_names,
    from = object$attitude_names
  )
  list(response = array, spectral_radius = response$spectral_radius)
}

#' Long-run influence of one attitude's innovation on another
#'
#' Under a stationary recursion `state = R %*% state_lag + innovation`, the
#' share is the fraction of the long-run variance of the receiving attitude
#' that comes from innovations in the source attitude.
#'
#' @param object A fitted model.
#' @return A matrix of shares, with a `variance` attribute holding the
#'   long-run variances. The matrix is empty when `R` is not stationary.
#' @examples
#' panel <- simulate_attitudes(n_persons = 20, n_waves = 3, n_factors = 1, seed = 1)
#' fit <- attitudedyn(
#'   cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave, time = ~ time,
#'   n_factors = 1, draws = 2,
#'   control = attitudedyn_control(maxit = 0, n_threads = 1, vcov = "none")
#' )
#' influence_shares(fit)
#' @export
influence_shares <- function(object) {
  matrix <- persistence_matrix(object)
  result <- attitudedyn_influence_cpp(as.numeric(t(matrix)), nrow(matrix))
  if (!isTRUE(result$stationary)) {
    cli::cli_warn("The persistence matrix is not stationary, so long-run shares are not defined.")
    return(structure(matrix(numeric(0), 0, 0), stationary = FALSE))
  }
  J <- nrow(matrix)
  shares <- matrix(result$share, J, J, byrow = TRUE,
                   dimnames = list(object$attitude_names, object$attitude_names))
  attr(shares, "variance") <- setNames(result$long_run_variance, object$attitude_names)
  attr(shares, "stationary") <- TRUE
  shares
}

#' Variance decomposition of the latent attitude
#'
#' Reports the four pieces printed by the replication program. `covariate`
#' is the variance of the wave effect plus the covariate index. `wave_effect`
#' is the variance of the wave effect on its own. `factor` is the sum of
#' squared loadings. `shock` averages, across included people, the
#' finite-horizon contribution built from the gap-scaled persistence
#' coefficients. Both variances use every wave after a person's first
#' observed wave through that person's last wave, including a wave with no
#' interview. The divisor is the number of person-waves.
#'
#' @param object A fitted model.
#' @return A data frame with one row per attitude. `total` is the sum of the
#'   four pieces. A `note` attribute describes the columns.
#' @examples
#' panel <- simulate_attitudes(n_persons = 20, n_waves = 3, n_factors = 1, seed = 1)
#' fit <- attitudedyn(
#'   cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave, time = ~ time,
#'   n_factors = 1, draws = 2,
#'   control = attitudedyn_control(maxit = 0, n_threads = 1, vcov = "none")
#' )
#' variance_decomposition(fit)
#' @export
variance_decomposition <- function(object) {
  object <- ensure_problem(object)
  component <- attitudedyn_variance_cpp(
    object$cache$problem, object$theta_full
  )
  out <- data.frame(
    attitude = object$attitude_names,
    covariate = component$covariate,
    wave_effect = component$wave_effect,
    factor = component$factor,
    shock = component$shock,
    total = component$total,
    stringsAsFactors = FALSE
  )
  attr(out, "note") <- component$note
  out
}

#' Wald test of linear restrictions
#'
#' @param object A fitted model.
#' @param hypothesis A named numeric vector. Names are unconstrained
#'   coefficient names and values are the values under the null.
#' @param type Covariance type passed to [vcov()].
#' @return A list with the chi-squared statistic, degrees of freedom and
#'   p-value.
#' @examples
#' panel <- simulate_attitudes(n_persons = 20, n_waves = 3, n_factors = 1, seed = 1)
#' fit <- attitudedyn(
#'   cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave, time = ~ time,
#'   n_factors = 1, draws = 2,
#'   control = attitudedyn_control(maxit = 0, n_threads = 1, vcov = "opg")
#' )
#' wald_test(fit, c("beta[y1,x]" = 0))
#' @export
wald_test <- function(object, hypothesis, type = NULL) {
  if (is.null(names(hypothesis)) || any(names(hypothesis) == "")) {
    cli::cli_abort("{.arg hypothesis} must be a named numeric vector.")
  }
  unknown <- setdiff(names(hypothesis), names(object$coefficients))
  if (length(unknown)) cli::cli_abort("Unknown coefficient{?s}: {.val {unknown}}.")
  estimates <- coef(object)[names(hypothesis)]
  variance <- vcov(object, type = type)[names(hypothesis), names(hypothesis), drop = FALSE]
  gap <- estimates - as.numeric(hypothesis)
  statistic <- as.numeric(t(gap) %*% solve(variance, gap))
  df <- length(hypothesis)
  list(statistic = statistic, df = df, p.value = stats::pchisq(statistic, df, lower.tail = FALSE))
}

#' Score test of a nested restriction
#'
#' The model in `object` must already impose the null: the coefficients being
#' tested are fixed at the null value, and the remaining free coefficients are
#' at the restricted maximum. The statistic is the Lagrange-multiplier form
#' `s2' (G^{-1})_{22} s2` using the outer product of the person scores.
#'
#' @param object A fit of the restricted model.
#' @param coefficients Names of coefficients that are fixed in `object` and
#'   are being tested.
#' @return A list with the chi-squared statistic, degrees of freedom and
#'   p-value.
#' @examples
#' panel <- simulate_attitudes(n_persons = 20, n_waves = 3, n_factors = 1, seed = 1)
#' restricted <- attitudedyn(
#'   cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave, time = ~ time,
#'   n_factors = 1, draws = 2, fix = c("beta[y1,x]" = 0),
#'   control = attitudedyn_control(maxit = 0, n_threads = 1, vcov = "none")
#' )
#' suppressWarnings(lm_test(restricted, "beta[y1,x]"))
#' @export
lm_test <- function(object, coefficients) {
  layout <- object$layout
  rows <- match(coefficients, layout$name)
  if (anyNA(rows)) {
    cli::cli_abort("Unknown coefficient{?s}: {.val {coefficients[is.na(rows)]}}.")
  }
  if (any(layout$is_free[rows] == 1L)) {
    cli::cli_abort("Each tested coefficient must be fixed in {.arg object}. Refit with {.arg fix}.")
  }
  object <- ensure_problem(object)
  free_rows <- which(layout$is_free == 1L)
  columns <- c(free_rows, rows)
  evaluation <- attitudedyn_eval_cpp(
    object$cache$problem, object$theta_full, FALSE,
    numerics_index(object$control$numerics), object$control$fortran_ghk,
    object$control$n_threads, TRUE, TRUE, FALSE, NULL, NULL
  )
  # The full-gradient evaluation differentiates every coefficient, in layout
  # order, only when grad_columns is every index. attitudedyn_eval_cpp does
  # that when theta_is_free is FALSE and want_gradient is TRUE.
  n_grad <- length(object$theta_full)
  scores <- matrix(evaluation$person_score, ncol = n_grad, byrow = TRUE)
  keep <- object$inputs$weight > 0
  scores <- scores[keep, columns, drop = FALSE]
  information <- crossprod(scores)
  total <- colSums(scores)
  n_free <- length(free_rows)
  tested <- seq_len(length(rows)) + n_free
  inverse <- solve(information)
  block <- inverse[tested, tested, drop = FALSE]
  score <- total[tested]
  statistic <- as.numeric(t(score) %*% block %*% score)
  if (n_free > 0L) {
    free_score <- total[seq_len(n_free)]
    free_scale <- sqrt(diag(information)[seq_len(n_free)])
    free_scale[free_scale == 0] <- 1
    if (max(abs(free_score) / free_scale) > 0.05) {
      cli::cli_warn("Free-coefficient scores are not close to zero. The statistic is the score test only at the restricted maximum.")
    }
  }
  df <- length(rows)
  list(statistic = statistic, df = df, p.value = stats::pchisq(statistic, df, lower.tail = FALSE))
}
