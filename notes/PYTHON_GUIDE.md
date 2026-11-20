# The Python package, in plain language

`attitudedyn` can be called from Python. The numbers are computed by the
same C++ library the R package uses. Python prepares the panel, asks that
library for the log-likelihood and the analytic gradient, and hands the
gradient to `scipy.optimize`. Nothing in `python/` rewrites the recursion,
the draws, or the normal interval probability.

## What you install

From a checkout of this repository:

```bash
pip install ./python
```

That command needs Python 3.10 or newer, a C++17 compiler, and CMake. On
macOS the compiler is Apple clang. OpenMP is not required and is not used.
On Windows the compiler is MSVC. `pip` builds a wheel with
`scikit-build-core` and pybind11, and the wheel links the sources in
`src/core/` directly. R is not installed and not imported.

The Python dependencies are NumPy, SciPy, and pandas. Tests additionally
need pytest:

```bash
python -m pip install pytest
python -m pytest python/tests
```

## What each file is

- `python/pyproject.toml` names the package and the build backend.
- `python/CMakeLists.txt` compiles `src/core/src/*.cpp` together with
  `python/src/bindings.cpp`. The binding copies arrays into the types
  declared in `src/core/include/attitudedyn/api.hpp` and copies results back.
- `python/src/attitudedyn/` is the code a user imports: `fit`, the result
  object, simulation, and the GSS reader.
- `python/tests/` is the pytest suite.
- `python/validation/run_gss.py` evaluates the GSS application.

`src/core/` itself is unchanged by the wrapper. The R package and
`R CMD check` do not see `python/`, because `.Rbuildignore` excludes it.

## How a fit runs

`attitudedyn.fit` reads a long data frame, one row per person-wave. Outcomes
are integers `1, ..., K`. Missing outcomes stay missing and are integrated
out. People are stored in first-appearance order, which is the order of the
Philox (or Halton, or Sobol) streams. A wave that is missing inside a
person's first and last waves is a gap: the latent attitude is simulated
through it.

The default generator is Philox, seed 1, 10 antithetic draws. `threads`
splits the person loop. Each thread writes its own person contributions, and
the main thread adds those contributions in person order, so the
log-likelihood and the gradient do not depend on the thread count.

BFGS minimises the negative log-likelihood with the same Nash iteration
as `stats::optim(method = "BFGS")`, so the same start, analytic gradient,
and `reltol` take the same steps as the R package. L-BFGS-B uses SciPy.
The gradient is the analytic gradient from the core. If a trial coefficient
is not finite, or a threshold increment is large enough that `exp` of it
would overflow, the core returns a non-finite objective and a zero gradient.
The wrapper turns that objective into `1e12` so the line search can reject
the step and continue. The reported fit is the best finite point seen during
the search.

Standard errors default to the centred outer product of the person scores
(`vcov="opg_centered"`), which is the covariance printed with the
replication estimates. `opg` skips the centring. `hessian` inverts the
negative numerical Hessian of the analytic gradient. `sandwich` is that
inverse, times the uncentred outer product of the scores, times the
inverse. `summary(vcov_type=...)` selects the same covariance.

`AttitudeDynResult` exposes `params`, `bse`, `pvalues`, `conf_int`,
`summary`, `loglik`, `aic`, `bic`, `nobs`, `vcov`, `persistence_matrix`,
`variance_decomposition`, `impulse_response`, `predict`, `simulate`, and
`to_dataframe`. Persistence reported by `persistence_matrix` is on (-1, 1).
`params` stays on the unconstrained scale the optimiser uses.

## GSS

The replication sample is `validation/fortran/data/gss2.prn`. The shipped
free coefficients, standard errors, and variance decomposition are in
`validation/fortran/src/gss_log.txt`.

```bash
python python/validation/run_gss.py --threads 4
python python/validation/replicate_gss_exhibits.py --threads 4
python python/validation/run_gss.py --refit --threads 4
```

The first command evaluates the shipped coefficients. With Philox, seed 1,
and 10 antithetic draws the log-likelihood is `-111570.97022260929` on the
libm that produced that figure. Another math library can move the sum by a
few units in the last place; the test allows an absolute gap under `1e-9`.
`replicate_gss_exhibits.py` is that same evaluation plus the centred
outer-product standard errors, and it writes the coefficient table, the
persistence matrix, the variance decomposition, the impulse responses, and
`notes/GSS_EXHIBITS.md`. `--refit` starts BFGS from those coefficients and
runs until the relative tolerance is met. Pass `--maxit` only to impose a
cap.

The published Fortran objective at the same coefficients is not that
number. The package generator is Philox. The check above is the Philox
figure.
