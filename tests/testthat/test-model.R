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
  wald <- wald_test(fit, c("beta[y1,x]" = 0))
  expect_true(is.finite(wald$p.value))
  predicted <- predict(fit, type = "class", n_draws = 2)
  expect_true(all(predicted$y1 %in% 1:3))
})
