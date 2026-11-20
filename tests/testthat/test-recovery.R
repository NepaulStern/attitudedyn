# A tiny copy of validation/simulation_recovery.R.
# One short fit is enough to check bias, RMSE, and coverage of a 95% interval.

normal_critical_value <- function(level = 0.95) {
  stats::qnorm(1 - (1 - level) / 2)
}

recovery_record <- function(fit, truth) {
  estimate <- coef(fit)
  truth_free <- truth[names(estimate)]
  standard_error <- sqrt(diag(vcov(fit)))
  critical <- normal_critical_value()
  lower <- estimate - critical * standard_error
  upper <- estimate + critical * standard_error
  error <- estimate - truth_free
  data.frame(
    parameter = names(estimate),
    truth = as.numeric(truth_free),
    estimate = as.numeric(estimate),
    bias = as.numeric(error),
    rmse = abs(as.numeric(error)),
    covered = truth_free >= lower & truth_free <= upper,
    stringsAsFactors = FALSE
  )
}

test_that("a tiny panel reports finite bias, RMSE, and coverage", {
  panel <- simulate_attitudes(
    n_persons = 36, n_waves = 3, n_attitudes = 2,
    n_factors = 1, n_categories = 3, seed = 15
  )
  truth <- attr(panel, "theta_full")
  fit <- suppressWarnings(attitudedyn(
    cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave, time = ~ time,
    n_factors = 1, draws = 3, include_wave_effects = FALSE,
    control = attitudedyn_control(maxit = 4, n_threads = 1, vcov = "opg", reltol = 1e-6)
  ))
  record <- recovery_record(fit, truth)
  expect_true(all(is.finite(record$bias)))
  expect_true(all(record$rmse >= 0))
  expect_true(all(record$covered %in% c(TRUE, FALSE)))
  expect_gte(mean(record$covered), 0)
  expect_lte(mean(record$covered), 1)
  expect_equal(sqrt(mean(record$bias^2)), {
    error <- record$estimate - record$truth
    sqrt(mean(error^2))
  })
})

test_that("bias at the starting values equals the start minus the truth", {
  panel <- simulate_attitudes(
    n_persons = 20, n_waves = 3, n_attitudes = 2,
    n_factors = 1, n_categories = 3, seed = 6
  )
  truth <- attr(panel, "theta_full")
  fit <- attitudedyn(
    cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave, time = ~ time,
    n_factors = 1, draws = 2, include_wave_effects = FALSE,
    control = attitudedyn_control(maxit = 0, n_threads = 1, vcov = "opg")
  )
  # maxit = 0 stays at the default start. The slope on x starts at 0 and the
  # data-generating slope is 0.5, so the signed error is -0.5.
  record <- recovery_record(fit, truth)
  slope <- record[record$parameter == "beta[y1,x]", ]
  expect_equal(slope$truth, 0.5)
  expect_equal(slope$estimate, 0)
  expect_equal(slope$bias, -0.5)
  expect_equal(slope$rmse, 0.5)
})
