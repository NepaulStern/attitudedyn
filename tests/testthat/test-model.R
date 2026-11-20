panel_small <- function() {
  simulate_attitudes(n_persons = 40, n_waves = 3, n_attitudes = 2,
                     n_factors = 1, n_categories = 3, seed = 4)
}

test_that("bad inputs are rejected with a useful error", {
  panel <- panel_small()
  expect_error(
    attitudedyn(cbind(y1, y2) ~ x, data = panel, id = ~ missing, wave = ~ wave),
    "missing"
  )
  shifted <- panel
  shifted$y1 <- shifted$y1 - 1L
  expect_error(
    attitudedyn(cbind(y1, y2) ~ x, data = shifted, id = ~ id, wave = ~ wave),
    "1, 2"
  )
  expect_error(
    attitudedyn(cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave, draws = 0),
    "draws"
  )
  broken <- panel
  broken$x[1] <- Inf
  expect_error(
    attitudedyn(cbind(y1, y2) ~ x, data = broken, id = ~ id, wave = ~ wave),
    "covariate"
  )
})

test_that("the log-likelihood does not depend on the thread count", {
  panel <- panel_small()
  quiet <- attitudedyn_control(maxit = 0L, n_threads = 1L, vcov = "none", trace = 0L)
  one <- attitudedyn(cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave,
                     time = ~ time, n_factors = 1L, draws = 4L, control = quiet)
  two <- attitudedyn(cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave,
                     time = ~ time, n_factors = 1L, draws = 4L,
                     control = attitudedyn_control(maxit = 0L, n_threads = 2L, vcov = "none"))
  expect_identical(one$loglik, two$loglik)
  expect_identical(coef(one), coef(two))
})

test_that("the analytic gradient matches a central difference", {
  panel <- panel_small()
  fit <- attitudedyn(cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave,
                     time = ~ time, n_factors = 1L, draws = 4L,
                     control = attitudedyn_control(maxit = 0L, n_threads = 1L, vcov = "none"))
  base <- attitudedyn:::attitudedyn_eval_cpp(
    fit$cache$problem, fit$coefficients, TRUE, 0L, FALSE, 1L, TRUE, FALSE, FALSE, NULL, NULL
  )
  step <- 1e-5
  up <- down <- fit$coefficients
  up[1] <- up[1] + step
  down[1] <- down[1] - step
  plus <- attitudedyn:::attitudedyn_eval_cpp(fit$cache$problem, up, TRUE, 0L, FALSE, 1L, FALSE, FALSE, FALSE, NULL, NULL)
  minus <- attitudedyn:::attitudedyn_eval_cpp(fit$cache$problem, down, TRUE, 0L, FALSE, 1L, FALSE, FALSE, FALSE, NULL, NULL)
  numeric <- (plus$objective - minus$objective) / (2 * step)
  expect_equal(base$gradient[1], numeric, tolerance = 1e-4)
})

test_that("an overflowing threshold is a non-finite objective, not an uncaught error", {
  panel <- panel_small()
  quiet <- attitudedyn_control(maxit = 0L, n_threads = 1L, vcov = "none")
  fit <- attitudedyn(
    cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave, time = ~ time,
    n_factors = 1L, draws = 3L, control = quiet
  )
  start <- coef(fit)
  thresh <- grep("^threshold\\[", names(start))
  expect_gte(length(thresh), 1L)

  extreme <- start
  extreme[thresh[[1L]]] <- 80
  evaluation <- attitudedyn:::attitudedyn_eval_cpp(
    fit$cache$problem, extreme, TRUE, 0L, FALSE, 1L, TRUE, FALSE, FALSE, NULL, NULL
  )
  expect_false(is.finite(evaluation$objective))
  expect_length(evaluation$gradient, length(extreme))
  expect_true(all(is.finite(evaluation$gradient)))

  nonfinite <- start
  nonfinite[thresh[[1L]]] <- Inf
  inf_eval <- attitudedyn:::attitudedyn_eval_cpp(
    fit$cache$problem, nonfinite, TRUE, 0L, FALSE, 1L, TRUE, FALSE, FALSE, NULL, NULL
  )
  expect_false(is.finite(inf_eval$objective))
  expect_true(all(is.finite(inf_eval$gradient)))

  expect_error(
    attitudedyn(
      cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave, time = ~ time,
      n_factors = 1L, draws = 3L, start = extreme, control = quiet
    ),
    "starting values"
  )
})

overflow_start <- function() {
  # At this start the threshold scores are positive and larger than the
  # distance to the exp() guard, so the first BFGS / L-BFGS-B / nlminb trial
  # steps past it. The panel is still small enough for a unit test.
  panel <- simulate_attitudes(
    n_persons = 120, n_waves = 3, n_attitudes = 2,
    n_factors = 1, n_categories = 3, seed = 4
  )
  quiet <- attitudedyn_control(maxit = 0L, n_threads = 1L, vcov = "none")
  base <- attitudedyn(
    cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave, time = ~ time,
    n_factors = 1L, draws = 2L, control = quiet
  )
  start <- coef(base)
  start[grep("^threshold\\[", names(start))] <- -8
  list(panel = panel, start = start, problem = base$cache$problem)
}

test_that("BFGS, L-BFGS-B and nlminb reject an overflowing threshold trial", {
  prepared <- overflow_start()
  objective <- function(theta) {
    value <- attitudedyn:::attitudedyn_eval_cpp(
      prepared$problem, theta, TRUE, 0L, FALSE, 1L, FALSE, FALSE, FALSE, NULL, NULL
    )$objective
    if (!is.finite(value)) 1e12 else -value
  }
  bfgs <- optim(prepared$start, objective, method = "BFGS", control = list(maxit = 4L))
  lbfgs <- optim(prepared$start, objective, method = "L-BFGS-B", control = list(maxit = 4L))
  bounded <- nlminb(prepared$start, objective, control = list(iter.max = 4L))
  expect_lt(bfgs$value, 1e11)
  expect_lt(lbfgs$value, 1e11)
  expect_lt(bounded$objective, 1e11)
  expect_true(all(is.finite(bfgs$par)))
  expect_true(all(bfgs$par[grep("^threshold\\[", names(prepared$start))] <= 50))

  fit <- suppressWarnings(attitudedyn(
    cbind(y1, y2) ~ x, data = prepared$panel, id = ~ id, wave = ~ wave, time = ~ time,
    n_factors = 1L, draws = 2L, start = prepared$start,
    control = attitudedyn_control(maxit = 4L, n_threads = 1L, vcov = "none")
  ))
  expect_true(is.finite(fit$loglik))
  expect_gt(fit$loglik, -1e6)
  expect_true(all(is.finite(coef(fit))))
  expect_true(all(coef(fit)[grep("^threshold\\[", names(coef(fit)))] <= 50))

  numeric_fit <- suppressWarnings(attitudedyn(
    cbind(y1, y2) ~ x, data = prepared$panel, id = ~ id, wave = ~ wave, time = ~ time,
    n_factors = 1L, draws = 2L, start = prepared$start,
    control = attitudedyn_control(
      maxit = 1L, n_threads = 1L, vcov = "none", use_gradient = FALSE
    )
  ))
  expect_true(is.finite(numeric_fit$loglik))
  expect_gt(numeric_fit$loglik, -1e6)
})

test_that("a short fit returns the methods the package promises", {
  panel <- panel_small()
  fit <- suppressWarnings(attitudedyn(
    cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave, time = ~ time,
    n_factors = 1L, draws = 4L,
    control = attitudedyn_control(maxit = 2L, n_threads = 1L, vcov = "opg")
  ))
  expect_s3_class(fit, "attitudedyn")
  expect_true(is.finite(logLik(fit)))
  expect_equal(nobs(fit), length(unique(panel$id)))
  expect_true(all(is.finite(coef(fit))))
  expect_equal(dim(vcov(fit)), c(length(coef(fit)), length(coef(fit))))
  s <- summary(fit)
  expect_s3_class(s, "summary.attitudedyn")
  expect_s3_class(tidy(fit), "data.frame")
  expect_s3_class(glance(fit), "data.frame")
  rho <- persistence_matrix(fit)
  expect_equal(dim(rho), c(2L, 2L))
  expect_true(all(abs(rho) < 1))
  irf <- impulse_response(fit, horizon = 2)
  expect_equal(dim(irf$response), c(3L, 2L, 2L))
  expect_equal(irf$response[1, , ], diag(2), tolerance = 1e-8, ignore_attr = TRUE)
  parts <- variance_decomposition(fit)
  expect_equal(nrow(parts), 2L)
  expect_equal(parts$total, parts$covariate + parts$wave_effect + parts$factor + parts$shock)
  expect_true(all(parts$shock >= 1))
  wald <- wald_test(fit, c("beta[y1,x]" = 0))
  expect_true(is.finite(wald$p.value))
  predicted <- predict(fit, type = "class", n_draws = 2)
  expect_true(all(predicted$y1 %in% 1:3))
  path <- tempfile(fileext = ".rds")
  saveRDS(fit, path)
  loaded <- readRDS(path)
  expect_equal(persistence_matrix(loaded), persistence_matrix(fit))
  expect_equal(variance_decomposition(loaded)$total, variance_decomposition(fit)$total)
})

test_that("centred outer product is the default, and hessian and sandwich are covariances", {
  expect_identical(attitudedyn_control()$vcov, "opg_centered")
  panel <- simulate_attitudes(
    n_persons = 80, n_waves = 3, n_attitudes = 2, n_factors = 1, n_categories = 3, seed = 4
  )
  truth <- attr(panel, "theta_full")
  fit <- attitudedyn(
    cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave, time = ~ time,
    n_factors = 1L, draws = 4L, start = truth,
    control = attitudedyn_control(maxit = 40L, n_threads = 1L, vcov = "none")
  )
  expect_identical(fit$convergence, 0L)
  centred <- vcov(fit, type = "opg_centered")
  hessian <- vcov(fit, type = "hessian")
  sandwich <- vcov(fit, type = "sandwich")
  expect_equal(dim(hessian), dim(centred))
  expect_equal(dim(sandwich), dim(centred))
  expect_equal(hessian, t(hessian))
  expect_equal(sandwich, t(sandwich))
  expect_true(all(diag(hessian) > 0))
  expect_true(all(diag(sandwich) > 0))

  kept_scores <- fit$scores[fit$inputs$weight > 0, , drop = FALSE]
  centred_scores <- sweep(kept_scores, 2L, colMeans(kept_scores), check.margin = FALSE)
  expect_equal(centred, solve(crossprod(centred_scores)))
  inverse_negative_hessian <- solve(-attitudedyn:::numerical_hessian_of_analytic_gradient(fit))
  outer_product <- crossprod(kept_scores)
  expect_equal(unname(hessian), unname(inverse_negative_hessian))
  expect_equal(
    unname(sandwich),
    unname(inverse_negative_hessian %*% outer_product %*% inverse_negative_hessian)
  )

  shown <- summary(fit, vcov_type = "hessian")
  expect_identical(shown$summary$type, "hessian")
  expect_equal(unname(shown$summary$coefficients[, "Std. Error"]), unname(sqrt(diag(hessian))))
  expect_identical(summary(fit, type = "sandwich")$summary$type, "sandwich")
  expect_identical(summary(fit)$summary$type, "none")
})
