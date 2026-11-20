# Category probabilities and draws from a fitted model.

fitted_example <- function() {
  panel <- simulate_attitudes(
    n_persons = 25, n_waves = 3, n_attitudes = 2,
    n_factors = 1, n_categories = c(4L, 3L), seed = 9
  )
  fit <- attitudedyn(
    cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave, time = ~ time,
    n_factors = 1, draws = 3, include_wave_effects = FALSE,
    control = attitudedyn_control(maxit = 0, n_threads = 1, vcov = "none")
  )
  list(panel = panel, fit = fit)
}

test_that("category probabilities sum to one and the class is the mode", {
  prepared <- fitted_example()
  probabilities <- predict(prepared$fit, type = "prob", n_draws = 4, seed = 3)
  classes <- predict(prepared$fit, type = "class", n_draws = 4, seed = 3)
  expect_named(probabilities, c("y1", "y2"))
  expect_equal(nrow(probabilities$y1), nrow(prepared$panel))
  expect_equal(probabilities$y1$id, classes$id)
  expect_equal(probabilities$y1$wave, classes$wave)

  # Four categories on y1 and three on y2. Unused slots must not leak in.
  expect_named(probabilities$y1, c("id", "wave", "p1", "p2", "p3", "p4"))
  expect_named(probabilities$y2, c("id", "wave", "p1", "p2", "p3"))
  y1_probability <- as.matrix(probabilities$y1[, paste0("p", 1:4)])
  y2_probability <- as.matrix(probabilities$y2[, paste0("p", 1:3)])
  expect_equal(unname(rowSums(y1_probability)), rep(1, nrow(y1_probability)), tolerance = 1e-8)
  expect_equal(unname(rowSums(y2_probability)), rep(1, nrow(y2_probability)), tolerance = 1e-8)
  expect_true(all(y1_probability >= 0))
  expect_equal(classes$y1, max.col(y1_probability, ties.method = "first"))
  expect_equal(classes$y2, max.col(y2_probability, ties.method = "first"))
  expect_true(all(classes$y1 %in% 1:4))
  expect_true(all(classes$y2 %in% 1:3))
})

test_that("the same seed repeats the probabilities", {
  prepared <- fitted_example()
  first <- predict(prepared$fit, type = "prob", n_draws = 3, seed = 8)
  second <- predict(prepared$fit, type = "prob", n_draws = 3, seed = 8)
  expect_equal(first, second)
})

test_that("simulate draws a panel that can be fitted again", {
  prepared <- fitted_example()
  one <- simulate(prepared$fit, nsim = 1, seed = 4)
  two <- simulate(prepared$fit, nsim = 2, seed = 4)
  expect_s3_class(one, "data.frame")
  expect_type(two, "list")
  expect_length(two, 2)
  expect_equal(nrow(one), nrow(prepared$panel))
  expect_equal(one, two[[1]])
  expect_true(all(one$y1 %in% 1:4))
  expect_true(all(one$y2 %in% 1:3))
  expect_false("(Intercept)" %in% names(one))
  refit <- attitudedyn(
    cbind(y1, y2) ~ x, data = one, id = ~ id, wave = ~ wave, time = ~ time,
    n_factors = 1, draws = 2, include_wave_effects = FALSE,
    control = attitudedyn_control(maxit = 0, n_threads = 1, vcov = "none")
  )
  expect_true(is.finite(refit$loglik))
})
