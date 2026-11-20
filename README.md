# attitudedyn

`attitudedyn` estimates The Nepaul-Stern Model of Attitude Stability and Constraint. The model is for several Likert-scale attitudes measured on the same people across survey waves. A person's latent attitude carries forward from one wave to the next, and the recorded answer is the Likert category that latent attitude falls into. Estimation uses a simulated likelihood.

The same estimates are available from R and from Python.

The project site is <https://nepaulstern.github.io/>.

## Install

Installing from source needs a C++17 compiler. That is the usual requirement for a package with compiled code. You do not write any C++.

From the source directory of a checkout:

```r
install.packages(".", repos = NULL, type = "source")
```

```bash
pip install ./python
```

## A first example

These lines simulate a small panel and fit The Nepaul-Stern Model for a few steps, so you can see the shape of the output. `y1` and `y2` are Likert items. `x` is a covariate. `id` is the person and `wave` is the interview.

```r
library(attitudedyn)
panel <- simulate_attitudes(n_persons = 50, n_waves = 3, n_factors = 1, seed = 1)
fit <- attitudedyn(
  cbind(y1, y2) ~ x,
  data = panel,
  id = ~ id,
  wave = ~ wave,
  time = ~ time,
  n_factors = 1,
  draws = 8,
  control = attitudedyn_control(maxit = 10, n_threads = 1)
)
summary(fit)
```

`simulate_attitudes()` builds the practice data. `attitudedyn()` estimates the model. `time = ~ time` uses the wave dates in that panel. `summary(fit)` prints the log-likelihood, the coefficients, their standard errors, and the persistence matrix. The persistence matrix says how strongly each earlier attitude pulls the later one.

The same call in Python:

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
```

`time_col="time"` tells Python which column holds the date of each wave. In the R call above, the simulated panel carries that date and the function uses it when you pass `time = ~ time`. The short R example omits it; equally spaced waves are then assumed. Pass `time = ~ time` when your file has its own wave dates.

## Where to read next

The user guide is the place to start. It reads the General Social Survey panel, fits the application, and shows how to read the estimates, the persistence matrix, the impulse responses, the variance decomposition, the predicted Likert categories, and a new panel drawn from the fit.

* R: `vignette("attitudedyn user guide", package = "attitudedyn")`. Source: [`vignettes/user-guide.Rmd`](vignettes/user-guide.Rmd). PDF: [`docs/attitudedyn_user_guide.pdf`](docs/attitudedyn_user_guide.pdf).
* Python: [`python/docs/user_guide.md`](python/docs/user_guide.md).

Notes on the source layout, for someone changing the package, are in `notes/R_ECOSYSTEM_GUIDE.md` and `notes/PYTHON_GUIDE.md`.

## The General Social Survey application

The guide is the walkthrough. The General Social Survey file comes from the replication deposit on Harvard Dataverse (doi:10.7910/DVN/DUCHYS), dedicated to the public domain under CC0 1.0. `python3 validation/fortran/fetch_vendor.py` downloads that file, the replication log, and the Fortran comparison sources into `validation/fortran/vendor/`. The file has 11,282 people and 11,252 of them enter the likelihood, with 164 free coefficients. Evaluating the replication estimates, with the default of 10 random draws and seed 1, gives log-likelihood -111570.9702226093. `summary()` prints that value rounded.

On 4 processor cores that evaluation takes a few seconds. Searching from those coefficients, with the iteration cap lifted, took 202.6 seconds in R and 276.3 seconds in Python. Both searches stopped successfully at log-likelihood -111545.0088160454.

```bash
Rscript validation/run_gss.R --refit --threads 4
python python/validation/run_gss.py --refit --threads 4
```

`impulse_response()` on the replication estimates reports a spectral radius of 0.636678. A spectral radius below 1 means a one-time shock to a latent attitude dies out. In the figure, each row is the attitude that responds and each column is the attitude that was shocked. The left edge of each small plot is the moment of the shock.

![GSS impulse responses](man/figures/gss_impulse_response.png)

`validation/replicate_gss_exhibits.R` and `python/validation/replicate_gss_exhibits.py` write the coefficient table, persistence matrix, variance decomposition, and impulse responses to `validation/output/`. `notes/GSS_EXHIBITS.md` records that evaluation. How the comparison driver is built is in `notes/VALIDATION.md`. The likelihood and the reporting functions are described in `notes/DESIGN.md`.

## Acknowledgments

The C++/R/Python port, the documentation, and the deployment infrastructure were developed with Grok Bot and Cursor Cloud Agent.
