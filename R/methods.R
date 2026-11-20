#' @export
print.attitudedyn <- function(x, ...) {
  cat("Dynamic ordinal attitude model\n")
  cat("Log-likelihood: ", format(x$loglik, digits = 6), "\n", sep = "")
  cat("People: ", x$n_persons, "   Attitudes: ", x$inputs$n_attitudes,
      "   Waves: ", x$inputs$n_waves, "\n", sep = "")
  cat("Draws: ", x$draws, if (x$control$antithetic) " (plus antithetic twins)" else "", "\n", sep = "")
  cat("Free coefficients: ", length(x$coefficients), "\n", sep = "")
  if (x$convergence == 0L) cat("Converged.\n") else cat("Optimiser code: ", x$convergence, "\n", sep = "")
  invisible(x)
}

#' @export
coef.attitudedyn <- function(object, scale = c("unconstrained", "natural"), ...) {
  scale <- match.arg(scale)
  if (scale == "unconstrained") return(object$coefficients)
  natural_free_coefficients(object)
}

natural_free_coefficients <- function(object) {
  object <- ensure_problem(object)
  natural <- attitudedyn_natural_cpp(object$cache$problem, object$theta_full)
  J <- object$inputs$n_attitudes
  persistence <- matrix(natural$persistence, J, J, byrow = TRUE)
  out <- object$coefficients
  layout <- object$layout
  free <- layout[layout$is_free == 1L, , drop = FALSE]
  for (i in seq_len(nrow(free))) {
    if (free$block[i] == "persistence") {
      out[i] <- persistence[free$attitude[i], free$index2[i]]
    } else if (free$block[i] == "threshold") {
      cuts <- natural$cutpoints[[free$attitude[i]]]
      out[i] <- cuts[free$index2[i]]
    }
  }
  out
}

#' @export
nobs.attitudedyn <- function(object, ...) object$n_persons

#' @export
logLik.attitudedyn <- function(object, ...) {
  value <- object$loglik
  attr(value, "df") <- length(object$coefficients)
  attr(value, "nobs") <- nobs(object)
  class(value) <- "logLik"
  value
}

person_scores <- function(object) {
  object <- ensure_problem(object)
  scores <- object$scores
  if (is.null(scores)) {
    evaluation <- attitudedyn_eval_cpp(
      object$cache$problem, object$coefficients, TRUE,
      numerics_index(object$control$numerics), object$control$fortran_ghk,
      object$control$n_threads, TRUE, TRUE, FALSE, NULL, NULL
    )
    scores <- matrix(evaluation$person_score, ncol = length(object$coefficients), byrow = TRUE)
    colnames(scores) <- names(object$coefficients)
  }
  keep <- object$inputs$weight > 0
  scores[keep, , drop = FALSE]
}

outer_product_of_scores_covariance <- function(person_scores, center = FALSE) {
  # Invert sum (score)(score)' , or the same sum after subtracting column means.
  if (center) person_scores <- sweep(person_scores, 2L, colMeans(person_scores), check.margin = FALSE)
  information <- crossprod(person_scores)
  if (qr(information)$rank < ncol(information)) {
    cli::cli_warn("The outer product of the scores is singular. A ridge was added.")
    information <- information + diag(1e-8 * mean(diag(information)), ncol(information))
  }
  solve(information)
}

numerical_hessian_of_analytic_gradient <- function(object) {
  # Central differences of the analytic gradient. The step is
  # 1e-5 * max(1, |coefficient|). The returned matrix is symmetrised.
  object <- ensure_problem(object)
  coefficients <- object$coefficients
  n_coefficients <- length(coefficients)
  hessian <- matrix(0, n_coefficients, n_coefficients)
  for (column in seq_len(n_coefficients)) {
    step <- 1e-5 * max(1, abs(coefficients[column]))
    coefficient_up <- coefficient_down <- coefficients
    coefficient_up[column] <- coefficient_up[column] + step
    coefficient_down[column] <- coefficient_down[column] - step
    gradient_up <- attitudedyn_eval_cpp(
      object$cache$problem, coefficient_up, TRUE,
      numerics_index(object$control$numerics), object$control$fortran_ghk,
      object$control$n_threads, TRUE, FALSE, FALSE, NULL, NULL
    )$gradient
    gradient_down <- attitudedyn_eval_cpp(
      object$cache$problem, coefficient_down, TRUE,
      numerics_index(object$control$numerics), object$control$fortran_ghk,
      object$control$n_threads, TRUE, FALSE, FALSE, NULL, NULL
    )$gradient
    hessian[, column] <- (gradient_up - gradient_down) / (2 * step)
  }
  0.5 * (hessian + t(hessian))
}

#' @export
vcov.attitudedyn <- function(object, type = NULL, ...) {
  if (is.null(type)) type <- object$control$vcov
  type <- match.arg(type, c("opg", "opg_centered", "hessian", "sandwich", "none"))
  if (type == "none") {
    cli::cli_abort("This fit was asked not to compute a covariance. Call {.fn vcov} with {.code type}.")
  }
  cache <- object$cache
  if (!is.null(cache$vcov) && !is.null(cache$vcov[[type]])) return(cache$vcov[[type]])
  person_scores_matrix <- person_scores(object)
  inverse_negative_hessian <- NULL
  if (type %in% c("hessian", "sandwich")) {
    if (object$control$numerics != "stable" || object$contemporaneous != "identity") {
      cli::cli_abort("The numerical Hessian needs analytic gradients, hence stable numerics and an identity correlation.")
    }
    # H below is the Hessian of the log-likelihood. The covariance uses -H.
    hessian_of_log_likelihood <- numerical_hessian_of_analytic_gradient(object)
    inverse_negative_hessian <- tryCatch(solve(-hessian_of_log_likelihood), error = function(e) {
      cli::cli_abort("The Hessian could not be inverted.")
    })
  }
  estimate <- switch(
    type,
    opg = outer_product_of_scores_covariance(person_scores_matrix, center = FALSE),
    opg_centered = outer_product_of_scores_covariance(person_scores_matrix, center = TRUE),
    hessian = inverse_negative_hessian,
    sandwich = {
      # (-H)^{-1} (S'S) (-H)^{-1}, with S the uncentred person scores.
      outer_product <- crossprod(person_scores_matrix)
      inverse_negative_hessian %*% outer_product %*% inverse_negative_hessian
    }
  )
  colnames(estimate) <- rownames(estimate) <- names(object$coefficients)
  if (is.null(cache$vcov)) cache$vcov <- list()
  cache$vcov[[type]] <- estimate
  estimate
}

#' @export
confint.attitudedyn <- function(object, parm = NULL, level = 0.95, type = NULL, ...) {
  estimates <- coef(object)
  if (is.null(parm)) parm <- names(estimates)
  if (is.numeric(parm)) parm <- names(estimates)[parm]
  variance <- vcov(object, type = type)
  se <- sqrt(diag(variance))
  z <- stats::qnorm(1 - (1 - level) / 2)
  out <- cbind(estimates - z * se, estimates + z * se)
  colnames(out) <- paste(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2), trim = TRUE), "%")
  out[parm, , drop = FALSE]
}

#' Summarise a fitted attitude model
#'
#' @param object A fitted model.
#' @param vcov_type Covariance passed to [vcov()]. `"opg_centered"` centres
#'   the person scores and inverts their outer product. `"hessian"` inverts
#'   the negative numerical Hessian of the analytic gradient. `"sandwich"`
#'   is that inverse, times the uncentred outer product of the scores, times
#'   the inverse. The default is the covariance stored with the fit.
#' @param type The same choice as `vcov_type`.
#' @param ... Unused.
#' @return An object of class `"summary.attitudedyn"`.
#' @export
summary.attitudedyn <- function(object, vcov_type = NULL, type = NULL, ...) {
  if (is.null(vcov_type)) vcov_type <- type
  if (is.null(vcov_type)) vcov_type <- object$control$vcov
  vcov_type <- match.arg(vcov_type, c("opg", "opg_centered", "hessian", "sandwich", "none"))
  estimates <- coef(object)
  variance <- tryCatch(vcov(object, type = vcov_type), error = function(e) NULL)
  if (is.null(variance)) {
    se <- rep(NA_real_, length(estimates))
  } else {
    se <- sqrt(diag(variance))
  }
  statistic <- estimates / se
  p_value <- 2 * stats::pnorm(-abs(statistic))
  table <- cbind(Estimate = estimates, `Std. Error` = se, `z value` = statistic, `Pr(>|z|)` = p_value)
  natural <- coef(object, scale = "natural")
  object$summary <- list(
    coefficients = table,
    natural = natural,
    type = vcov_type
  )
  class(object) <- c("summary.attitudedyn", class(object))
  object
}

#' @export
print.summary.attitudedyn <- function(x, ...) {
  print.attitudedyn(x)
  cat("\nCoefficients (unconstrained scale):\n")
  stats::printCoefmat(x$summary$coefficients, signif.stars = TRUE, digits = 4)
  cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  cat("The 10% cutoff is 1.645.\n")
  cat("\nPersistence on (-1, 1), receiving attitude by source attitude:\n")
  print(round(persistence_matrix(x), 3))
  cat("\nCovariance type: ", x$summary$type, "\n", sep = "")
  invisible(x)
}

#' @export
formula.attitudedyn <- function(x, ...) x$formula

#' Predicted Likert probabilities for each person-wave
#'
#' Averages structural paths from the fitted model. Each path uses the
#' person's covariates, factors, and lagged latent attitude. The category
#' probability on a path is the normal interval probability around that
#' latent value, and the reported probability is the average across paths.
#' The same routine supplies `type = "class"` and `type = "latent"`.
#'
#' @param object A fitted [attitudedyn()] model.
#' @param type `"prob"` returns category probabilities. `"class"` returns
#'   the category with the largest probability. `"latent"` returns the
#'   average latent attitude before the innovation is added.
#' @param n_draws Number of structural paths per person. The default is
#'   `max(20, object$draws)`. These paths are not the likelihood draws.
#' @param seed Generator seed. The default is the fit seed plus one, so a
#'   later call with the same seed repeats the probabilities.
#' @param ... Unused. Included so the method matches [stats::predict()].
#' @return For `"prob"`, a named list of data frames, one per attitude.
#'   Each row is a person-wave from that person's first wave through the
#'   last, including a skipped wave inside the window, and the columns
#'   `p1` through `pK` are the probabilities of categories 1 through K.
#'   For `"class"`, a data frame with one integer column per attitude.
#'   For `"latent"`, a data frame with one numeric column per attitude.
#' @examples
#' panel <- simulate_attitudes(n_persons = 20, n_waves = 3, n_factors = 1, seed = 1)
#' fit <- attitudedyn(
#'   cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave, time = ~ time,
#'   n_factors = 1, draws = 2,
#'   control = attitudedyn_control(maxit = 0, n_threads = 1, vcov = "none")
#' )
#' head(predict(fit, type = "prob", n_draws = 2)$y1)
#' head(predict(fit, type = "class", n_draws = 2))
#' @export
predict.attitudedyn <- function(object, type = c("prob", "class", "latent"),
                                n_draws = NULL, seed = NULL, ...) {
  type <- match.arg(type)
  object <- ensure_problem(object)
  if (is.null(n_draws)) n_draws <- max(20L, object$draws)
  if (is.null(seed)) seed <- object$control$seed + 1
  prediction <- attitudedyn_predict_cpp(
    object$cache$problem, object$theta_full,
    as.integer(n_draws), as.numeric(seed), object$control$n_threads
  )
  inputs <- object$inputs
  J <- inputs$n_attitudes
  T <- inputs$n_waves
  Kmax <- prediction$max_categories
  # The C++ array includes every wave for every person, including people and
  # waves outside the window. Rebuild the long index in that same order.
  id_out <- integer(0)
  wave_out <- integer(0)
  keep <- logical(0)
  for (i in seq_len(inputs$n_persons)) {
    for (t in seq_len(T)) {
      inside <- t - 1L >= inputs$first_wave[i] && t - 1L <= inputs$last_wave[i]
      id_out <- c(id_out, inputs$person_ids[i])
      wave_out <- c(wave_out, inputs$wave_values[t])
      keep <- c(keep, inside && inputs$weight[i] > 0)
    }
  }
  if (type == "latent") {
    latent <- matrix(prediction$mean_latent, ncol = J, byrow = TRUE)
    colnames(latent) <- inputs$attitude_names
    out <- data.frame(id = id_out, wave = wave_out, latent, check.names = FALSE)
    return(out[keep, , drop = FALSE])
  }
  if (type == "class") {
    classes <- matrix(NA_integer_, length(id_out), J)
    colnames(classes) <- inputs$attitude_names
    for (j in seq_len(J)) {
      K <- inputs$n_categories[j]
      for (row in which(keep)) {
        person <- ((row - 1L) %/% T) + 1L
        wave <- ((row - 1L) %% T) + 1L
        base <- (((person - 1L) * T + (wave - 1L)) * J + (j - 1L)) * Kmax
        probs <- prediction$probability[base + seq_len(K)]
        classes[row, j] <- which.max(probs)
      }
    }
    out <- data.frame(id = id_out, wave = wave_out, classes, check.names = FALSE)
    return(out[keep, , drop = FALSE])
  }
  probabilities <- vector("list", J)
  names(probabilities) <- inputs$attitude_names
  for (j in seq_len(J)) {
    K <- inputs$n_categories[j]
    columns <- matrix(NA_real_, length(id_out), K)
    colnames(columns) <- paste0("p", seq_len(K))
    for (row in which(keep)) {
      person <- ((row - 1L) %/% T) + 1L
      wave <- ((row - 1L) %% T) + 1L
      base <- (((person - 1L) * T + (wave - 1L)) * J + (j - 1L)) * Kmax
      columns[row, ] <- prediction$probability[base + seq_len(K)]
    }
    probabilities[[j]] <- data.frame(id = id_out, wave = wave_out, columns, check.names = FALSE)[keep, , drop = FALSE]
  }
  probabilities
}

#' @export
fitted.attitudedyn <- function(object, ...) {
  probabilities <- predict(object, type = "prob", ...)
  observed <- vector("list", length(probabilities))
  names(observed) <- names(probabilities)
  inputs <- object$inputs
  T <- inputs$n_waves
  for (j in seq_along(probabilities)) {
    frame <- probabilities[[j]]
    category <- integer(nrow(frame))
    for (row in seq_len(nrow(frame))) {
      person <- match(frame$id[row], inputs$person_ids)
      wave <- match(frame$wave[row], inputs$wave_values)
      base <- ((person - 1L) * inputs$n_attitudes + (j - 1L)) * T
      category[row] <- inputs$category[base + wave]
    }
    value <- rep(NA_real_, nrow(frame))
    known <- category >= 1L
    value[known] <- frame[cbind(which(known), category[known] + 2L)]
    observed[[j]] <- data.frame(id = frame$id, wave = frame$wave, probability = value)
  }
  observed
}

#' Draw panels from a fitted attitude model
#'
#' Draws covariates and ordinal outcomes from the structural model at the
#' fitted coefficients. Intercept columns are filled with ones and the other
#' covariates are standard normal. This is a new panel of the same size as
#' the estimation sample, not a resample of its rows.
#'
#' @param object A fitted [attitudedyn()] model.
#' @param nsim Number of panels.
#' @param seed Seed of the first panel. Later panels use the following
#'   integers. The default is the fit seed plus 17. The seed does not change
#'   R's random-number generator.
#' @param ... Unused. Included so the method matches [stats::simulate()].
#' @return A data frame when `nsim` is 1, and a list of data frames
#'   otherwise. Columns are `id`, `wave`, `time`, the covariates other than
#'   the intercept, and the outcomes. Category 1 is the lowest category.
#'   The frame can be passed back to [attitudedyn()].
#' @examples
#' panel <- simulate_attitudes(n_persons = 20, n_waves = 3, n_factors = 1, seed = 1)
#' fit <- attitudedyn(
#'   cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave, time = ~ time,
#'   n_factors = 1, draws = 2,
#'   control = attitudedyn_control(maxit = 0, n_threads = 1, vcov = "none")
#' )
#' head(simulate(fit, nsim = 1, seed = 2))
#' @export
simulate.attitudedyn <- function(object, nsim = 1, seed = NULL, ...) {
  if (is.null(seed)) seed <- object$control$seed + 17
  simulations <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    panel <- attitudedyn_simulate_cpp(
      n_waves = object$inputs$n_waves,
      n_attitudes = object$inputs$n_attitudes,
      n_covariates = object$inputs$n_covariates,
      n_factors = object$n_factors,
      n_categories = object$n_categories,
      wave_time = object$time_grid,
      include_wave_effects = object$inputs$include_wave_effects,
      free_contemporaneous = identical(object$contemporaneous, "free"),
      fix_last_loading = object$control$fix_last_loading,
      fix_intercepts = object$control$fix_intercept && object$inputs$include_wave_effects,
      intercept_columns = object$inputs$intercept_columns,
      empty_wave_effect = object$inputs$empty_wave_effect,
      attitude_names = object$attitude_names,
      covariate_names = object$covariate_names,
      wave_names = object$inputs$wave_names,
      theta_full = object$theta_full,
      n_persons = object$n_persons,
      seed = as.numeric(seed) + s - 1
    )
    simulations[[s]] <- simulation_frame(panel, object$attitude_names, object$covariate_names)
  }
  if (nsim == 1L) simulations[[1L]] else simulations
}

#' Tidy a fitted attitude model
#'
#' `tidy()` returns one row per unconstrained coefficient. `glance()` returns
#' one row of fit statistics. Both follow the column names used by broom.
#'
#' @param x A fitted model.
#' @param conf.int Include Wald confidence limits.
#' @param conf.level Confidence level.
#' @param scale `"unconstrained"` matches [coef()]. `"natural"` transforms
#'   persistence coefficients onto (-1, 1) and reports delta-method standard
#'   errors. Threshold rows on the natural scale are the cutpoints; their
#'   standard errors stay on the unconstrained scale.
#' @param ... Unused.
#' @return `tidy()` returns a data frame with one row per free coefficient.
#'   `glance()` returns a one-row data frame of fit statistics.
#' @examples
#' panel <- simulate_attitudes(n_persons = 20, n_waves = 3, n_factors = 1, seed = 1)
#' fit <- attitudedyn(
#'   cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave, time = ~ time,
#'   n_factors = 1, draws = 2,
#'   control = attitudedyn_control(maxit = 0, n_threads = 1, vcov = "opg")
#' )
#' tidy(fit)
#' glance(fit)
#' @exportS3Method broom::tidy
tidy.attitudedyn <- function(x, conf.int = FALSE, conf.level = 0.95,
                             scale = c("unconstrained", "natural"), ...) {
  scale <- match.arg(scale)
  estimates <- coef(x, scale = scale)
  variance <- tryCatch(vcov(x), error = function(e) NULL)
  if (scale == "natural" && !is.null(variance)) {
    se <- natural_standard_errors(x, variance)
  } else if (!is.null(variance)) {
    se <- sqrt(diag(variance))
  } else {
    se <- rep(NA_real_, length(estimates))
  }
  out <- data.frame(
    term = names(estimates),
    estimate = as.numeric(estimates),
    std.error = as.numeric(se),
    statistic = as.numeric(estimates / se),
    stringsAsFactors = FALSE
  )
  out$p.value <- 2 * stats::pnorm(-abs(out$statistic))
  if (conf.int) {
    z <- stats::qnorm(1 - (1 - conf.level) / 2)
    out$conf.low <- out$estimate - z * out$std.error
    out$conf.high <- out$estimate + z * out$std.error
  }
  out
}

natural_standard_errors <- function(object, variance) {
  se <- sqrt(diag(variance))
  natural <- attitudedyn_natural_cpp(object$cache$problem, object$theta_full)
  J <- object$inputs$n_attitudes
  persistence <- matrix(natural$persistence, J, J, byrow = TRUE)
  free <- object$layout[object$layout$is_free == 1L, , drop = FALSE]
  for (i in seq_len(nrow(free))) {
    if (free$block[i] == "persistence") {
      rho <- persistence[free$attitude[i], free$index2[i]]
      se[i] <- abs(0.5 * (1 - rho^2)) * se[i]
    }
  }
  se
}

#' @rdname tidy.attitudedyn
#' @exportS3Method broom::glance
glance.attitudedyn <- function(x, ...) {
  data.frame(
    logLik = as.numeric(logLik(x)),
    AIC = AIC(x),
    BIC = BIC(x),
    nobs = nobs(x),
    n_free = length(x$coefficients),
    draws = x$draws,
    converged = x$convergence == 0L,
    stringsAsFactors = FALSE
  )
}
