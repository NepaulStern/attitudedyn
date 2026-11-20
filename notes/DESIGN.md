# Design of attitudedyn

The package reimplements the authors' replication program in C++ and R. The likelihood is the GHK simulated likelihood in that program. One evaluation is the sum of the person contributions, minus one penalty term. With the published GSS specification the contemporaneous correlation is the identity, and the penalty is zero.

## The recursion

For each person, each draw, and the sign flip of that draw, the latent attitudes are updated from the person's first wave to the last. On a wave the contribution of the lag is the persistence matrix times the previous latent attitude. The covariate index, the wave effect, and the factor are part of that latent attitude, so they carry forward. A missing item adds an unrestricted standard normal shock and contributes nothing to the log-likelihood. An observed item contributes the log probability of its interval and then draws the shock from the truncated normal that lands in the interval.

The gap between wave `t` and wave `t-1` is `time_grid[t] - time_grid[t-1]`, using one date per wave, the same date for every person. The persistence entry that multiplies a lag is `sign(rho) * |rho|^gap`. For the GSS grid `c(2, 3, 4, 5, 6, 7, 8, 10)` the gaps are 1, 1, 1, 1, 1, 1, 2. The age imputation in the GSS reader uses the wave index, `2 * (k - j)`, which is a separate step from this exponent.

`rho` stored by `coef()` is the unconstrained coefficient. The value inside (-1, 1) is `2 * exp(c) / (1 + exp(c)) - 1`. `persistence_matrix()` reports that value. The row is the attitude that receives the lag. The column is the attitude that is lagged.

## Categories, factors, and what is held fixed

Outcomes are integers `1, ..., K`. The cut between category 1 and category 2 is 0. The remaining cuts are built from log-increments, so they stay ordered. The top and bottom of the real line are fixed at 100 and -100.

The published specification has two standard normal person factors, uncorrelated with each other and with the shocks. Shocks are uncorrelated within a wave. The loading of the last attitude on the last factor is fixed at 0 when there are two or more factors. When wave effects are included, each attitude's intercept is fixed at 0. A wave effect with no observed responses is fixed at 0.

## Draws and the gradient

The default generator is Philox, a counter-based generator. The counter is the person, the draw, and a stream number. Factor `f` uses stream `f`. The shock of attitude `a` at wave `t` uses stream `n_factors + a * n_waves + t`. Antithetic sampling flips the sign of the normal, which doubles the number of paths. Halton and Sobol are available; both are Owen-scrambled and keyed by the person, and the Sobol sequence skips the origin. The generator is not a bit-for-bit copy of the generator in the replication program. Two runs with the same seed, the same draws, and the same data produce the same log-likelihood, including across thread counts.

The default probability of an interval is the double-precision log of `Phi(b) - Phi(a)`. That setting, `numerics = "stable"`, supplies the analytic gradient. `numerics = "shared"` reproduces the original program's normal-probability routine for comparisons. Analytic gradients are available for `numerics = "stable"` only.

`contemporaneous = "free"` estimates a within-wave correlation factor. The analytic gradient is not provided for that case, and the fit says so.

## Reporting

`coef()` is the unconstrained vector, in the same order the replication program prints its free coefficients: wave effects, slopes, unconstrained persistence, loadings, the correlation factor, and threshold increments. Standard errors default to the centred outer product of the person scores (`vcov(fit, type = "opg_centered")`). `opg` skips the centring. `hessian` inverts the negative numerical Hessian of the analytic gradient. `sandwich` is that inverse, times the uncentred outer product of the scores, times the inverse. `summary(fit, vcov_type = ...)` selects the same covariance. Significance stars use the normal critical value 1.645 for 10 percent.

`variance_decomposition()` follows the replication program. The covariate piece is the variance of the wave effect plus the covariate index, the wave-effect piece is the variance of the wave effect, the factor piece is the sum of squared loadings, and the shock piece averages the finite-horizon persistence contribution over included people. Both variances skip each person's first observed wave and then keep every later wave in the window. `impulse_response()` reports powers of the persistence matrix.

## What is deliberately outside the package

The optimiser, the line search, and the linear-algebra routines that the replication program links are not part of the package. The comparison driver under `validation/fortran/` downloads them from the CC0 Harvard Dataverse deposit doi:10.7910/DVN/DUCHYS into a git-ignored folder and links them so that it can execute the original likelihood. The R package calls `stats::optim` for BFGS and Nelder-Mead, and a short scoring loop for BHHH.
