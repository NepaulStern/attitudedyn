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
)
print(result.summary())
result.persistence_matrix()
result.variance_decomposition()
result.impulse_response().to_dataframe()
```

Column lists work the same way: pass `outcomes=` and `covariates=` instead
of `formula`. The formula helper accepts `cbind(...)` and bare column names,
including `` 0 + `(Intercept)` `` to use a column of ones already in the data.

`result.params`, `result.bse`, and `result.pvalues` are named Series.
`result.vcov()` defaults to the outer product of the person scores (`opg`).
`opg_centered`, `hessian`, and `sandwich` are the other covariance options.
BFGS uses the same Nash iteration as `stats::optim(method = "BFGS")` in the
R package. L-BFGS-B uses SciPy. Both use the analytic gradient. A trial step
that is not finite, or whose threshold increment would overflow, is reported
to the optimiser as a large finite penalty with a zero gradient, so the line
search rejects it.

The GSS check is `python/validation/run_gss.py`. Build notes are in
`notes/PYTHON_GUIDE.md`.
