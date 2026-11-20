#' Control parameters for [attitudedyn()]
#'
#' Settings for the random draws, the search, and the standard errors when
#' fitting The Nepaul-Stern Model of Attitude Stability and Constraint.
#'
#' @param antithetic Pair each draw with its sign-flipped twin, as in the
#'   original program. The reported draw count is the number of pairs.
#' @param draw_type `"philox"` is a counter-based pseudorandom generator keyed
#'   by person and draw. `"halton"` and `"sobol"` are Owen-scrambled
#'   quasi-random sequences, also keyed by person so that people do not share
#'   one sequence.
#' @param seed Non-negative seed. It does not change R's random-number state.
#' @param n_threads Number of threads used to loop over people. `1` is the
#'   default so examples and tests stay inside the limit expected by CRAN.
#'   Results do not depend on this value.
#' @param optimizer `"BFGS"` uses the analytic gradient. `"BHHH"` is a
#'   scoring step with the outer product of the person scores. `"Nelder-Mead"`
#'   is derivative-free.
#' @param maxit Maximum number of optimiser iterations.
#' @param reltol Relative tolerance on the log-likelihood, passed to
#'   [stats::optim()] for BFGS and Nelder-Mead.
#' @param use_gradient Use the analytic gradient. Turn this off only to
#'   compare with a numeric gradient. A free contemporaneous correlation has
#'   no analytic gradient; the fit then uses numeric derivatives and says so.
#' @param vcov Which covariance to compute and store with the fit.
#'   `"opg_centered"`, the default, centres the person scores and inverts
#'   their outer product, which is what the Fortran program did.
#'   `"opg"` is the same inverse without centring. `"hessian"` inverts the
#'   negative numerical Hessian of the analytic gradient. `"sandwich"` is
#'   that inverse, times the uncentred outer product of the scores, times
#'   the inverse. `"none"` skips it.
#' @param numerics `"stable"`, the default, uses a double-precision log of a
#'   normal interval probability and supplies analytic gradients. `"shared"`
#'   reproduces the original program's normal-probability routine for
#'   comparisons.
#' @param fortran_ghk Use the original within-wave conditioning when the
#'   contemporaneous correlation is not the identity. With an identity
#'   correlation the two conditionings agree.
#' @param fix_last_loading Fix the loading of the last attitude on the last
#'   factor at zero when there are two or more factors. This is the rotation
#'   normalisation used in the published estimates. It is not applied when
#'   there is only one factor, because that factor's scale is already fixed
#'   by the unit variance of the person effect.
#' @param fix_intercept Fix each attitude's intercept at zero when wave
#'   effects are included. The intercept and a full set of wave effects are
#'   not separately identified once a cutpoint is normalised at zero.
#' @param keep_scores Store the person score matrix. It is needed for the
#'   outer-product covariance and is about `n_persons * n_free` doubles.
#' @param trace If positive, print the log-likelihood at each iteration.
#'
#' @return A named list of class `"attitudedyn_control"`.
#' @export
#'
#' @examples
#' attitudedyn_control(n_threads = 1L, maxit = 20L)
attitudedyn_control <- function(antithetic = TRUE,
                                draw_type = c("philox", "halton", "sobol"),
                                seed = 1,
                                n_threads = 1L,
                                optimizer = c("BFGS", "BHHH", "Nelder-Mead"),
                                maxit = 200L,
                                reltol = 1e-8,
                                use_gradient = TRUE,
                                vcov = c("opg_centered", "opg", "hessian", "sandwich", "none"),
                                numerics = c("stable", "shared"),
                                fortran_ghk = FALSE,
                                fix_last_loading = TRUE,
                                fix_intercept = TRUE,
                                keep_scores = TRUE,
                                trace = 0L) {
  draw_type <- match.arg(draw_type)
  optimizer <- match.arg(optimizer)
  vcov <- match.arg(vcov)
  numerics <- match.arg(numerics)
  if (!is.finite(seed) || seed < 0) {
    cli::cli_abort("{.arg seed} must be a finite, non-negative number.")
  }
  if (!is.numeric(n_threads) || length(n_threads) != 1L || n_threads < 1) {
    cli::cli_abort("{.arg n_threads} must be a positive integer.")
  }
  if (!is.numeric(maxit) || length(maxit) != 1L || maxit < 0) {
    cli::cli_abort("{.arg maxit} must be a non-negative integer.")
  }
  if (!is.numeric(reltol) || length(reltol) != 1L || reltol <= 0) {
    cli::cli_abort("{.arg reltol} must be a positive number.")
  }
  structure(
    list(
      antithetic = isTRUE(antithetic),
      draw_type = draw_type,
      seed = as.numeric(seed),
      n_threads = as.integer(n_threads),
      optimizer = optimizer,
      maxit = as.integer(maxit),
      reltol = as.numeric(reltol),
      use_gradient = isTRUE(use_gradient),
      vcov = vcov,
      numerics = numerics,
      fortran_ghk = isTRUE(fortran_ghk),
      fix_last_loading = isTRUE(fix_last_loading),
      fix_intercept = isTRUE(fix_intercept),
      keep_scores = isTRUE(keep_scores),
      trace = as.integer(trace)
    ),
    class = "attitudedyn_control"
  )
}

draw_kind_index <- function(draw_type) {
  switch(draw_type, philox = 0L, halton = 1L, sobol = 2L, 0L)
}

numerics_index <- function(numerics) {
  if (identical(numerics, "shared")) 1L else 0L
}
