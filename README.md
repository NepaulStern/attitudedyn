# attitudedyn

R package for the dynamic ordered-probit model of ordinal panel attitudes in Nepaul and Stern (2026, *Political Analysis*). The likelihood is a GHK simulated likelihood. The numerical work lives in a dependency-free C++ library under `src/core/`. The R package and the Python package under `python/` are wrappers around it. The latent attitude persists from wave to wave.

## Getting started

The user guide is the place to start. It installs the package, reads the General Social Survey panel, fits the application, and reads `summary()`, the persistence matrix, the impulse responses, the variance decomposition, `predict()`, and `simulate()`.

* R: `vignette("attitudedyn user guide", package = "attitudedyn")`. Source: [`vignettes/user-guide.Rmd`](vignettes/user-guide.Rmd). PDF: [`docs/attitudedyn_user_guide.pdf`](docs/attitudedyn_user_guide.pdf).
* Python: [`python/docs/user_guide.md`](python/docs/user_guide.md).

A C++17 compiler is required.

```r
# from the source directory of a checkout
install.packages(".", repos = NULL, type = "source")
```

```bash
pip install ./python
```

The smoke test below uses a simulated panel. The guide uses the GSS sample.

```r
library(attitudedyn)
panel <- simulate_attitudes(n_persons = 50, n_waves = 3, n_factors = 1, seed = 1)
fit <- attitudedyn(cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave,
                   n_factors = 1, draws = 8,
                   control = attitudedyn_control(maxit = 10, n_threads = 1))
summary(fit)
```

Build notes for someone new to R packages are in `notes/R_ECOSYSTEM_GUIDE.md`. Python build notes are in `notes/PYTHON_GUIDE.md`.

## GSS application

The guide in the section above is the walkthrough. The replication file is `validation/fortran/data/gss2.prn`: 11,282 people in the file and 11,252 in the likelihood, with 164 free coefficients. Evaluating the replication estimates (Philox, seed 1, 10 antithetic draws) gives log-likelihood -111570.9702226093. `summary()` prints that value rounded.

With 4 threads, that evaluation takes a few seconds. Dropping the iteration cap and searching from those coefficients took 202.6 seconds in R and 276.3 seconds in Python. Both stopped with convergence code 0 at log-likelihood -111545.0088160454.

```bash
Rscript validation/run_gss.R --refit --threads 4
python python/validation/run_gss.py --refit --threads 4
```

`impulse_response()` on the replication estimates reports a spectral radius of 0.636678. Rows of the figure receive the shock named in the column. Horizon 0 is the impact.

![GSS impulse responses](man/figures/gss_impulse_response.png)

`validation/replicate_gss_exhibits.R` and `python/validation/replicate_gss_exhibits.py` write the coefficient table, persistence matrix, variance decomposition, and impulse responses to `validation/output/`. `notes/GSS_EXHIBITS.md` records the median and maximum absolute differences from the replication log.
