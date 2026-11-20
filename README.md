# attitudedyn

R package for the dynamic ordered-probit model of ordinal panel attitudes in Nepaul and Stern (2026, *Political Analysis*). The likelihood is a GHK simulated likelihood. The numerical work lives in a dependency-free C++ library under `core/`; the R package is a wrapper around it.

`attitudedyn` reimplements the authors' replication package in C++ and R. The same C++ core is wrapped for Python under `python/`; see `notes/PYTHON_GUIDE.md`. The latent attitude persists from wave to wave.

```r
library(attitudedyn)
panel <- simulate_attitudes(n_persons = 50, n_waves = 3, n_factors = 1, seed = 1)
fit <- attitudedyn(cbind(y1, y2) ~ x, data = panel, id = ~ id, wave = ~ wave,
                   n_factors = 1, draws = 8,
                   control = attitudedyn_control(maxit = 10, n_threads = 1))
summary(fit)
```

Build notes for someone new to R packages are in `notes/R_ECOSYSTEM_GUIDE.md`. The GSS replication script is `validation/run_gss.R`.
