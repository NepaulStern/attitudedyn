# attitudedyn for Python

Python interface to the dynamic ordered-probit model in Nepaul and Stern
(2026, *Political Analysis*). The likelihood is the GHK simulated likelihood
from the authors' replication package. The numerical work is the C++17
library in `../core/`; this directory is a wrapper around that library and
does not reimplement it. There is no R dependency.

## Install

A C++17 compiler and CMake are required. From a checkout of the repository:

```bash
pip install ./python
```

`pip` compiles `core/` with `scikit-build-core` and pybind11. The same
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

The General Social Survey call has the same shape and does not iterate
when `maxit` is 0. Run it from a checkout so the data path resolves.

```python
from attitudedyn.gss import COVARIATES, OUTCOMES, TIME_GRID, read_gss, read_gss_estimates

panel = read_gss("validation/fortran/data/gss2.prn")
estimates = read_gss_estimates("validation/fortran/src/gss_log.txt")
gss = attitudedyn.fit(
    panel,
    outcomes=OUTCOMES,
    covariates=COVARIATES,
    id_col="id",
    wave_col="wave",
    time_grid=TIME_GRID,
    weights="weight",
    n_factors=2,
    draws=10,
    add_intercept=False,
    start=estimates.to_numpy(),
    threads=1,
    maxit=0,
    vcov="opg_centered",
    seed=1,
)
print(gss.summary())
gss.to_dataframe().head()
gss.persistence_matrix()
gss.variance_decomposition()
gss.impulse_response().to_dataframe()
```

Column lists work the same way when there is no formula. The formula helper
accepts `cbind(...)` and bare column names, including `` 0 + `(Intercept)` ``
to use a column of ones already in the data.

`result.params`, `result.bse`, and `result.pvalues` are named Series.
`result.vcov()` defaults to the outer product of the person scores (`opg`).
`opg_centered` centres those scores first, which is the covariance printed
with the replication estimates. `hessian` and `sandwich` are the other
covariance options.

BFGS uses the same Nash iteration as `stats::optim(method = "BFGS")` in the
R package. L-BFGS-B uses SciPy. Both use the analytic gradient. A trial step
that is not finite, or whose threshold increment would overflow, is reported
to the optimiser as a large finite penalty with a zero gradient, so the line
search rejects it.

The GSS check is `python/validation/run_gss.py`. Build notes are in
`notes/PYTHON_GUIDE.md`.
