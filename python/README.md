# attitudedyn for Python

Python interface to the dynamic ordered-probit model in Nepaul and Stern
(2026, *Political Analysis*). The likelihood is the GHK simulated likelihood
from the authors' replication package. The numerical work is the C++17
library in `../src/core/`; this directory is a wrapper around that library and
does not reimplement it. There is no R dependency.

## Getting started

[`docs/user_guide.md`](docs/user_guide.md) is the walkthrough. It reads the
General Social Survey panel, fits the application, and reads the summary,
persistence, impulse responses, variance decomposition, predictions, and
simulated panels. The R guide is `vignettes/user-guide.Rmd`.

## Install

A C++17 compiler and CMake are required. From a checkout of the repository:

```bash
pip install ./python
```

`pip` compiles `src/core/` with `scikit-build-core` and pybind11. The same
sources build with Apple clang (OpenMP is not used) and with MSVC. Person
contributions are summed with the core's `std::thread` pool. The `threads`
argument changes the speed and not the value.

```python
import attitudedyn

panel = attitudedyn.simulate_attitudes(n_persons=50, n_waves=3, n_factors=1, seed=1)
result = attitudedyn.fit(
    panel,
    formula="cbind(y1, y2) ~ x",
    id_col="id",
    wave_col="wave",
    time_col="time",
    n_factors=1,
    draws=8,
    maxit=10,
    threads=1,
    vcov="opg",
)
print(result.summary())
result.to_dataframe().head()
result.persistence_matrix()
result.variance_decomposition()
result.impulse_response(horizon=4).to_dataframe()
```

`summary()` is the text table. `to_dataframe()` is one row per unconstrained
coefficient, with the standard error, z statistic, and p-value.
`persistence_matrix()` is on (-1, 1): the row receives the lag and the
column is lagged. `variance_decomposition()` reports the four pieces from
the replication program. `covariate` is the variance of the wave effect
plus the covariate index, `wave_effect` is the variance of the wave effect,
`factor` is the sum of squared loadings, and `shock` is the average
finite-horizon persistence contribution. `impulse_response()` returns
powers of the persistence matrix. Horizon 0 is the identity.

The General Social Survey walkthrough, with the output of that call, is
[`docs/user_guide.md`](docs/user_guide.md). Run it from a checkout so
`validation/fortran/data/gss2.prn` resolves. `maxit=0` evaluates `start`
and does not iterate.

Column lists work the same way when there is no formula. The formula helper
accepts `cbind(...)` and bare column names, including `` 0 + `(Intercept)` ``
to use a column of ones already in the data.

`result.params`, `result.bse`, and `result.pvalues` are named Series.
`result.vcov()` defaults to the centred outer product of the person scores
(`opg_centered`), which is the covariance printed with the replication
estimates. `opg` skips the centring. `hessian` inverts the negative
numerical Hessian of the analytic gradient. `sandwich` is that inverse,
times the uncentred outer product of the scores, times the inverse.
`summary(vcov_type=...)` selects the same covariance.

BFGS uses the same Nash iteration as `stats::optim(method = "BFGS")` in the
R package. L-BFGS-B uses SciPy. Both use the analytic gradient. A trial step
that is not finite, or whose threshold increment would overflow, is reported
to the optimiser as a large finite penalty with a zero gradient, so the line
search rejects it.

The GSS check is `python/validation/run_gss.py`. Build notes are in
`notes/PYTHON_GUIDE.md`.
