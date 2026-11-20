default_theta <- function(layout) {
  theta <- rep(0, nrow(layout))
  diagonal <- layout$block == "persistence" & layout$attitude == layout$index2
  # Unconstrained value whose transformed persistence is 0.3.
  theta[diagonal] <- log((1 + 0.3) / (1 - 0.3))
  theta[layout$block == "loading" & layout$is_free == 1L] <- 0.5
  theta[layout$block == "threshold"] <- log(0.6)
  theta
}

apply_fix <- function(layout, theta, fix) {
  if (is.null(fix)) return(list(layout = layout, theta = theta))
  if (is.null(names(fix)) || any(names(fix) == "")) {
    cli::cli_abort("{.arg fix} must be a named numeric vector of unconstrained coefficients.")
  }
  unknown <- setdiff(names(fix), layout$name)
  if (length(unknown)) {
    cli::cli_abort("Unknown coefficient{?s} in {.arg fix}: {.val {unknown}}.")
  }
  for (nm in names(fix)) {
    row <- match(nm, layout$name)
    theta[row] <- fix[[nm]]
    layout$is_free[row] <- 0L
  }
  list(layout = layout, theta = theta)
}

build_problem <- function(inputs, theta_full, layout, n_factors, draws, control, contemporaneous) {
  attitudedyn_build_cpp(
    n_waves = inputs$n_waves,
    n_attitudes = inputs$n_attitudes,
    n_covariates = inputs$n_covariates,
    n_factors = as.integer(n_factors),
    n_categories = inputs$n_categories,
    wave_time = inputs$wave_time,
    include_wave_effects = inputs$include_wave_effects,
    free_contemporaneous = identical(contemporaneous, "free"),
    fix_last_loading = inputs$fix_last_loading,
    fix_intercepts = inputs$fix_intercept && inputs$include_wave_effects,
    intercept_columns = inputs$intercept_columns,
    empty_wave_effect = inputs$empty_wave_effect,
    attitude_names = inputs$attitude_names,
    covariate_names = inputs$covariate_names,
    wave_names = inputs$wave_names,
    first_wave = inputs$first_wave,
    last_wave = inputs$last_wave,
    interviewed = inputs$interviewed,
    category = inputs$category,
    covariates = inputs$covariates,
    weight = inputs$weight,
    theta_full = theta_full,
    is_free = as.integer(layout$is_free),
    n_draws = as.integer(draws),
    antithetic = control$antithetic,
    draw_kind = draw_kind_index(control$draw_type),
    seed = control$seed
  )
}

free_from_full <- function(theta_full, layout) theta_full[layout$is_free == 1L]

full_from_free <- function(theta_free, theta_template, layout) {
  theta <- theta_template
  theta[layout$is_free == 1L] <- theta_free
  theta
}

# saveRDS keeps the fitted coefficients but not the external pointer to the
# prepared panel. Rebuild that pointer from the stored inputs so a reloaded
# fit can still report persistence, scores, and the variance decomposition.
ensure_problem <- function(object) {
  if (!inherits(object, "attitudedyn")) {
    cli::cli_abort("{.arg object} must be an {.cls attitudedyn} fit.")
  }
  if (!is.environment(object$cache)) object$cache <- new.env(parent = emptyenv())
  alive <- FALSE
  if (!is.null(object$cache$problem)) {
    alive <- isTRUE(tryCatch({
      attitudedyn_natural_cpp(object$cache$problem, object$theta_full)
      TRUE
    }, error = function(e) FALSE))
  }
  if (!alive) {
    object$cache$problem <- build_problem(
      object$inputs, object$theta_full, object$layout,
      object$n_factors, object$draws, object$control, object$contemporaneous
    )
    object$cache$layout <- object$layout
    object$cache$theta_template <- object$theta_full
  }
  object
}

negative_finite <- function(value) {
  if (!is.finite(value)) 1e12 else -value
}

# A line search can hand the likelihood a coefficient outside the domain
# (a non-finite value, or a threshold whose exponential gap overflows).
# That is a rejected trial, not a failure of the fit.
parameter_domain_failure <- function(condition) {
  # Only the old unpack exception. Data errors ("a covariate is not finite")
  # still abort the fit.
  grepl("threshold increment|exp\\(threshold\\) overflowed", conditionMessage(condition))
}

finite_gradient <- function(gradient, n) {
  if (is.null(gradient) || length(gradient) != n || !all(is.finite(gradient))) {
    return(rep(0, n))
  }
  gradient
}

run_bhhh <- function(start, loglik_fn, score_fn, maxit, reltol, trace) {
  theta <- start
  current <- loglik_fn(theta)
  convergence <- 1L
  for (iteration in seq_len(maxit)) {
    scores <- score_fn(theta)
    if (!all(is.finite(scores))) break
    gradient <- colSums(scores)
    information <- crossprod(scores)
    ridge <- 1e-8 * mean(diag(information))
    if (!is.finite(ridge) || ridge == 0) ridge <- 1e-6
    direction <- tryCatch(solve(information + diag(ridge, nrow(information)), gradient),
                          error = function(e) gradient)
    if (!all(is.finite(direction))) break
    step <- 1
    accepted <- current
    candidate <- theta
    improved <- FALSE
    # A full scoring step on an unscaled panel can be far too long. Halve
    # until the trial stays inside the domain and does not decrease the
    # log-likelihood. A failed search keeps the current point.
    for (halving in seq_len(30L)) {
      candidate <- theta + step * direction
      accepted <- loglik_fn(candidate)
      if (is.finite(accepted) && accepted >= current - 1e-8) {
        improved <- TRUE
        break
      }
      step <- step / 2
    }
    if (!improved) break
    if (trace > 0L) {
      cat(sprintf("BHHH iteration %d loglik %.6f\n", iteration, accepted))
    }
    move <- sqrt(mean((candidate - theta)^2))
    theta <- candidate
    if (is.finite(current) && is.finite(accepted) && current != 0 &&
        abs(accepted - current) / (1 + abs(current)) < reltol && move < reltol) {
      current <- accepted
      convergence <- 0L
      break
    }
    current <- accepted
  }
  list(par = theta, value = current, convergence = convergence)
}

#' Fit a dynamic ordinal attitude model
#'
#' Estimates the panel ordered-probit model in the Nepaul and Stern (2026)
#' replication package by simulated maximum likelihood. Several ordinal
#' attitudes are linked by a matrix of persistence coefficients and by
#' person-level factors. The whole latent attitude carries forward from one
#' wave to the next, as in the Fortran likelihood.
#'
#' The formula has one row per person-wave and one column per attitude:
#'
#' ```r
#' attitudedyn(cbind(duty, trust, efficacy) ~ female + education,
#'             data = panel, id = ~ person, wave = ~ wave)
#' ```
#'
#' A missing attitude is integrated out. A wave that is absent for a person,
#' but lies between that person's first and last waves, is filled in by
#' simulation. Wave dates are common to everyone: pass `time_grid` or a
#' `time` column that does not vary within a wave.
#'
#' @param formula A formula. The left-hand side is one outcome or
#'   [cbind()] of ordinal outcomes coded `1, ..., K`.
#' @param data A long-format data frame.
#' @param id,wave The person identifier and the wave, as a one-sided formula
#'   or a column name.
#' @param time Optional time column. The mean within each wave is the common
#'   wave date. It must not vary within a wave.
#' @param time_grid Optional numeric wave dates, in the order of the sorted
#'   unique values of `wave`. Gaps in the persistence exponent are differences
#'   of these dates. The GSS application uses
#'   `c(2, 3, 4, 5, 6, 7, 8, 10)`.
#' @param subset Optional expression selecting rows of `data`.
#' @param weights Optional nonnegative weights, constant within each person.
#'   The sample log-likelihood is the weighted sum of the person contributions.
#' @param na.action Treatment of missing covariates. Missing outcomes are
#'   kept. `na.omit` drops person-waves with missing covariates.
#' @param n_factors Number of standard normal person factors. Use `0` for none.
#' @param draws Number of GHK draws per person. Each draw is paired with its
#'   sign flip when antithetic sampling is on, so the average uses `2 * draws`
#'   paths.
#' @param include_wave_effects Include a wave effect for every attitude.
#' @param contemporaneous `"identity"` keeps the innovations uncorrelated
#'   within a wave, which is the published specification. `"free"` estimates
#'   the lower triangle of a correlation factor and uses numeric derivatives.
#' @param start Optional starting values for the free coefficients, named or
#'   in free-coefficient order. The scale is the unconstrained scale used by
#'   the optimiser.
#' @param fix Optional named numeric vector of coefficients to hold fixed,
#'   on the unconstrained scale.
#' @param control An [attitudedyn_control()] list.
#'
#' @return An object of class `"attitudedyn"`.
#' @export
#'
#' @examples
#' panel <- simulate_attitudes(n_persons = 30, n_waves = 3, n_factors = 1, seed = 1)
#' fit <- attitudedyn(
#'   cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave, time = ~ time,
#'   n_factors = 1, draws = 4,
#'   control = attitudedyn_control(maxit = 0, n_threads = 1, vcov = "none")
#' )
#' coef(fit)
attitudedyn <- function(formula, data, id, wave, time = NULL, time_grid = NULL,
                        subset, weights, na.action = na.pass,
                        n_factors = 2L, draws = 10L,
                        include_wave_effects = TRUE,
                        contemporaneous = c("identity", "free"),
                        start = NULL, fix = NULL,
                        control = attitudedyn_control()) {
  call <- match.call()
  contemporaneous <- match.arg(contemporaneous)
  if (!inherits(control, "attitudedyn_control")) {
    cli::cli_abort("{.arg control} must come from {.fn attitudedyn_control}.")
  }
  if (!is.numeric(n_factors) || length(n_factors) != 1L || n_factors < 0) {
    cli::cli_abort("{.arg n_factors} must be a nonnegative integer.")
  }
  if (!is.numeric(draws) || length(draws) != 1L || draws < 1) {
    cli::cli_abort("{.arg draws} must be a positive integer.")
  }
  n_factors <- as.integer(n_factors)
  draws <- as.integer(draws)

  if (!missing(subset)) {
    selected <- eval(call$subset, data, parent.frame())
    if (!is.logical(selected) || length(selected) != nrow(data)) {
      cli::cli_abort("{.arg subset} must be a logical vector with one entry per row.")
    }
    data <- data[selected & !is.na(selected), , drop = FALSE]
  }
  weight_vector <- NULL
  if (!missing(weights)) {
    weight_vector <- eval(call$weights, data, parent.frame())
  }

  inputs <- prepare_attitudedyn(
    formula = formula, data = data, id = id, wave = wave, time = time,
    time_grid = time_grid, weight_vector = weight_vector, na.action = na.action,
    control = control
  )
  inputs$include_wave_effects <- isTRUE(include_wave_effects)
  inputs$n_factors <- n_factors

  if (sum(inputs$weight > 0) < 2L) {
    cli::cli_abort("At least two people with positive weight are required.")
  }
  if (inputs$n_waves < 2L) {
    cli::cli_warn("The panel has a single wave, so persistence is not informed by within-person change.")
  }

  layout_raw <- attitudedyn_layout_cpp(
    n_waves = inputs$n_waves, n_attitudes = inputs$n_attitudes,
    n_covariates = inputs$n_covariates, n_factors = n_factors,
    n_categories = inputs$n_categories, wave_time = inputs$wave_time,
    include_wave_effects = inputs$include_wave_effects,
    free_contemporaneous = identical(contemporaneous, "free"),
    fix_last_loading = control$fix_last_loading,
    fix_intercepts = control$fix_intercept && inputs$include_wave_effects,
    intercept_columns = inputs$intercept_columns,
    empty_wave_effect = inputs$empty_wave_effect,
    attitude_names = inputs$attitude_names,
    covariate_names = inputs$covariate_names,
    wave_names = inputs$wave_names
  )
  layout <- data.frame(
    name = layout_raw$name, block = layout_raw$block,
    attitude = layout_raw$attitude, index2 = layout_raw$index2,
    is_free = layout_raw$is_free, stringsAsFactors = FALSE
  )
  theta_full <- default_theta(layout)
  fixed <- apply_fix(layout, theta_full, fix)
  layout <- fixed$layout
  theta_full <- fixed$theta
  if (!is.null(start)) {
    if (is.null(names(start))) {
      free_n <- sum(layout$is_free == 1L)
      if (length(start) != free_n) {
        cli::cli_abort("{.arg start} has length {length(start)} but there are {free_n} free coefficients.")
      }
      theta_full[layout$is_free == 1L] <- start
    } else {
      unknown <- setdiff(names(start), layout$name)
      if (length(unknown)) cli::cli_abort("Unknown names in {.arg start}: {.val {unknown}}.")
      for (nm in names(start)) theta_full[match(nm, layout$name)] <- start[[nm]]
    }
  }
  if (any(layout$block == "wave_effect" & layout$is_free == 0L & layout$name %in%
          paste0("wave_effect[", inputs$attitude_names, ",", rep(inputs$wave_names, each = inputs$n_attitudes), "]"))) {
    # The message below is the useful one; the condition above is too brittle.
  }
  empty_fixed <- layout$block == "wave_effect" & layout$is_free == 0L
  if (any(empty_fixed)) {
    cli::cli_inform(c(
      "Wave effects with no observed responses are fixed at zero.",
      "i" = "{sum(empty_fixed)} coefficient{?s}: {.val {layout$name[empty_fixed]}}."
    ))
  }

  use_gradient <- control$use_gradient && control$numerics == "stable" && contemporaneous == "identity"
  if (control$use_gradient && !use_gradient) {
    cli::cli_warn("Analytic gradients need stable numerics and an identity contemporaneous correlation. Using numeric derivatives.")
  }
  if (control$optimizer == "BHHH" && !use_gradient) {
    cli::cli_abort("BHHH needs analytic scores. Use BFGS, or keep the identity correlation and stable numerics.")
  }

  problem <- build_problem(inputs, theta_full, layout, n_factors, draws, control, contemporaneous)
  cache <- new.env(parent = emptyenv())
  cache$problem <- problem
  cache$theta_template <- theta_full
  cache$layout <- layout

  loglik_of <- function(theta_free) {
    value <- tryCatch(
      attitudedyn_eval_cpp(
        cache$problem, theta_free, TRUE,
        numerics_index(control$numerics), control$fortran_ghk, control$n_threads,
        FALSE, FALSE, FALSE, NULL, NULL
      )$objective,
      error = function(e) {
        if (parameter_domain_failure(e)) NA_real_ else stop(e)
      }
    )
    if (is.finite(value) && value > best_loglik) {
      best_loglik <<- value
      best_theta <<- theta_free
    }
    value
  }
  gradient_of <- function(theta_free) {
    gradient <- tryCatch(
      attitudedyn_eval_cpp(
        cache$problem, theta_free, TRUE,
        numerics_index(control$numerics), control$fortran_ghk, control$n_threads,
        TRUE, FALSE, FALSE, NULL, NULL
      )$gradient,
      error = function(e) {
        if (parameter_domain_failure(e)) NULL else stop(e)
      }
    )
    finite_gradient(gradient, length(theta_free))
  }
  scores_of <- function(theta_free) {
    evaluation <- tryCatch(
      attitudedyn_eval_cpp(
        cache$problem, theta_free, TRUE,
        numerics_index(control$numerics), control$fortran_ghk, control$n_threads,
        TRUE, TRUE, FALSE, NULL, NULL
      ),
      error = function(e) {
        if (parameter_domain_failure(e)) NULL else stop(e)
      }
    )
    if (is.null(evaluation) || !is.finite(evaluation$objective)) {
      return(matrix(0, nrow = inputs$n_persons, ncol = length(theta_free)))
    }
    matrix(evaluation$person_score, ncol = length(theta_free), byrow = TRUE)
  }

  free_start <- free_from_full(theta_full, layout)
  free_names <- layout$name[layout$is_free == 1L]
  names(free_start) <- free_names
  if (!length(free_start)) cli::cli_abort("Every coefficient is fixed, so there is nothing to estimate.")

  # Best finite point seen by the line search. A rejected trial does not
  # replace it, and the reported fit is this point when the optimiser
  # returns a step that left the domain.
  best_theta <- free_start
  best_loglik <- -Inf
  objective <- loglik_of(free_start)
  if (!is.finite(objective)) {
    cli::cli_abort("The log-likelihood is not finite at the starting values.")
  }
  if (control$trace > 0L) cat(sprintf("initial loglik %.6f\n", objective))

  if (control$optimizer == "BHHH") {
    fitted_optim <- run_bhhh(free_start, loglik_of, scores_of, control$maxit, control$reltol, control$trace)
    convergence <- fitted_optim$convergence
    coefficients <- fitted_optim$par
    loglik <- fitted_optim$value
  } else {
    method <- if (control$optimizer == "Nelder-Mead") "Nelder-Mead" else "BFGS"
    gradient_argument <- if (use_gradient && method == "BFGS") {
      function(parameter) -gradient_of(parameter)
    } else {
      NULL
    }
    fitted_optim <- stats::optim(
      free_start,
      function(parameter) negative_finite(loglik_of(parameter)),
      gr = gradient_argument,
      method = method,
      control = list(maxit = control$maxit, reltol = control$reltol,
                     trace = control$trace, REPORT = 1L)
    )
    convergence <- fitted_optim$convergence
    coefficients <- fitted_optim$par
    loglik <- if (is.finite(fitted_optim$value)) -fitted_optim$value else NA_real_
  }
  if (!is.finite(loglik) || (is.finite(best_loglik) && best_loglik > loglik + 1e-8)) {
    coefficients <- best_theta
    loglik <- best_loglik
  }
  names(coefficients) <- free_names
  theta_hat <- full_from_free(coefficients, theta_full, layout)
  names(theta_hat) <- layout$name

  final <- attitudedyn_eval_cpp(
    cache$problem, coefficients, TRUE,
    numerics_index(control$numerics), control$fortran_ghk, control$n_threads,
    use_gradient, control$keep_scores, TRUE, NULL, NULL
  )
  scores <- NULL
  if (control$keep_scores && !is.null(final$person_score)) {
    scores <- matrix(final$person_score, ncol = length(coefficients), byrow = TRUE)
    colnames(scores) <- free_names
  }

  result <- list(
    coefficients = coefficients,
    theta_full = theta_hat,
    loglik = final$objective,
    penalty = final$penalty,
    n_persons = final$n_persons_used,
    nobs = final$n_persons_used,
    layout = layout,
    convergence = convergence,
    call = call,
    formula = formula,
    terms = inputs$terms,
    contemporaneous = contemporaneous,
    draws = draws,
    control = control,
    inputs = inputs,
    n_factors = n_factors,
    cache = cache,
    scores = scores,
    na.action = na.action,
    attitude_names = inputs$attitude_names,
    covariate_names = inputs$covariate_names,
    wave_values = inputs$wave_values,
    time_grid = inputs$wave_time,
    n_categories = inputs$n_categories
  )
  class(result) <- "attitudedyn"
  warn_after_fit(result, use_gradient)
  if (control$vcov != "none") {
    tryCatch(stats::vcov(result, type = control$vcov),
             error = function(e) cli::cli_warn("The covariance matrix was not computed: {conditionMessage(e)}"))
  }
  result
}

warn_after_fit <- function(object, use_gradient) {
  # maxit = 0 evaluates the start and returns. That is not a failed fit.
  asked_to_optimise <- is.finite(object$control$maxit) && object$control$maxit > 0L
  if (asked_to_optimise && object$convergence != 0L) {
    cli::cli_warn("The optimiser did not report convergence (code {object$convergence}).")
  }
  natural <- persistence_matrix(object)
  radius <- tryCatch(max(Mod(eigen(natural, only.values = TRUE)$values)), error = function(e) NA_real_)
  if (is.finite(radius) && radius >= 0.999) {
    cli::cli_warn(c(
      "The persistence matrix has spectral radius {round(radius, 3)}.",
      "i" = "A radius of 1 or more means the latent recursion is not stationary."
    ))
  }
  if (asked_to_optimise && object$convergence == 0L && use_gradient && !is.null(object$scores)) {
    total_score <- colSums(object$scores)
    scale <- sqrt(diag(crossprod(object$scores)))
    scale[scale == 0] <- 1
    worst <- max(abs(total_score) / scale)
    if (is.finite(worst) && worst > 0.05) {
      cli::cli_warn("The largest score, scaled by its outer-product standard deviation, is {round(worst, 3)}. The estimates may not be at a maximum.")
    }
  }
  invisible(object)
}
