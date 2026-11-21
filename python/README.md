# attitudedyn for Python

Python interface for The Nepaul-Stern Model of Attitude Stability and Constraint. The model is for several Likert-scale attitudes measured on the same people across survey waves. A person's latent attitude carries forward from one wave to the next, and the recorded answer is the Likert category that latent attitude falls into. Estimation uses a simulated likelihood.

The same model is available from R. The two give the same estimates. The project site is <https://nepaulstern.github.io/>.

## Getting started

[`docs/user_guide.md`](https://github.com/gabrielgggg/attitudedyn/blob/master/python/docs/user_guide.md) is the walkthrough. It reads the General Social Survey panel, fits the application, and reads the summary, the persistence matrix, the impulse responses, the variance decomposition, the predicted Likert categories, and a new panel drawn from the fit. The R guide is [`vignettes/user-guide.Rmd`](https://github.com/gabrielgggg/attitudedyn/blob/master/vignettes/user-guide.Rmd). A PDF of that guide is [`docs/attitudedyn_user_guide.pdf`](https://github.com/gabrielgggg/attitudedyn/blob/master/docs/attitudedyn_user_guide.pdf).

## Install

Installing from a checkout needs a C++17 compiler and CMake. The installer compiles the shared numerical code. You do not write any C++. There is no R dependency.

From a checkout of the repository:

```bash
pip install ./python
```

Person contributions are summed after the cores finish, in person order. The `threads` argument changes the speed and not the estimates.

## A first example

These lines simulate a small panel and fit The Nepaul-Stern Model for a few steps. `y1` and `y2` are Likert items. `x` is a covariate. `id` is the person, `wave` is the interview, and `time` is the date of the wave.

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

`summary()` is the text table. `to_dataframe()` is one row per free coefficient, with the standard error, z statistic, and p-value. `persistence_matrix()` is on (-1, 1): the row receives the carry-over from the earlier attitude in the column. `variance_decomposition()` reports four pieces of each latent attitude. `covariate` is the part from the covariates and the wave shifts together, `wave_effect` is the wave shift on its own, `factor` is the person-specific factor, and `shock` is the part from the random shocks carried forward by persistence. `impulse_response()` traces a one-time shock. Horizon 0 is the moment of the shock.

The General Social Survey walkthrough, with the output of that call, is [`docs/user_guide.md`](https://github.com/gabrielgggg/attitudedyn/blob/master/python/docs/user_guide.md). The survey file comes from the replication deposit on Harvard Dataverse (doi:10.7910/DVN/DUCHYS), dedicated to the public domain under CC0 1.0. From a checkout, `python3 validation/fortran/fetch_vendor.py` downloads it into `validation/fortran/vendor/data/gss2.prn`. `maxit=0` evaluates `start` and does not search.

You can pass column lists instead of a formula. The formula helper accepts `cbind(...)` and bare column names, including `` 0 + `(Intercept)` `` when the data already contain a column of ones.

`result.params`, `result.bse`, and `result.pvalues` are named Series. `result.vcov()` defaults to `opg_centered`, the standard errors stored with the replication estimates. That calculation uses the variation in each person's contribution to the slope of the likelihood, after subtracting the average. `opg` skips the subtraction. `hessian` uses the curvature of the likelihood. `sandwich` combines the curvature with the person-level variation. `summary(vcov_type=...)` selects the same calculation.

The default search is BFGS, the same steps `stats::optim(method = "BFGS")` uses in the R package. L-BFGS-B uses SciPy. Both use the slope of the likelihood. A trial step that is not finite, or that would make a threshold increment overflow, is handed back to the search as a large finite penalty with a zero slope, so the search rejects it.

The GSS check is `python/validation/run_gss.py`. Notes on the source layout are in [`notes/PYTHON_GUIDE.md`](https://github.com/gabrielgggg/attitudedyn/blob/master/notes/PYTHON_GUIDE.md).

## Acknowledgments

The C++/R/Python port, the documentation, and the deployment infrastructure were developed with Grok Bot and Cursor Cloud Agent.
