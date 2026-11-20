# Names in the replication program, and the names in this package

The replication program is fixed-form Fortran. The package keeps the same quantities and uses names that read as R. Indices in the Fortran names are 1-based. Indices in the C++ core are 0-based.

| Replication program | In the package | What it is |
| --- | --- | --- |
| `likely12` | the objective returned by the core | Sum of the person log-likelihoods, minus one penalty |
| `likobs12` | one person's contribution | Log mean of the simulated paths |
| `preplik12` | the latent update inside a path | Covariate index, wave effect, factor, and the lagged latent attitude |
| `probmake1` | `numerics = "shared"` | Interval probability with the original tail blend |
| `ghk1` | `fortran_ghk = TRUE` | Within-wave conditioning when the shocks are correlated |
| `fixef{a}{t}` | `wave_effect[attitude, wave]` | Wave effect for one attitude |
| `beta` | `beta[attitude, covariate]` | Slope, on the latent scale |
| `rho{row}{col}` | `persistence[receiver, source]` | Unconstrained persistence; the natural value is in (-1, 1) |
| `alamb{f}{a}` | `loading[factor, attitude]` | Factor loading |
| `tau` | `threshold[attitude, cut]` | Log-increment of a cutpoint |
| `nseg` | `n_categories` | Number of categories of each attitude |
| `nr` | `draws` | Draws per person, before the antithetic sign flip |
| `nfac` | `n_factors` | Person factors |
| `ntime`, `time` | `time_grid` | Common wave dates. GSS uses the first 8 entries |
| `nobs` in the parameter file | rows read from `gss2.prn` | 21648 person-waves |
| `nobs` after the file is collapsed | people in the file | 11282 |
| `incluv` | `weights` | 1 if the person enters the likelihood, 0 otherwise |
| `istartv`, `iendv` | first and last wave of a person | The window the recursion runs over |
| `iobsrvv` | whether a wave was in the file | A gap inside the window still carries the latent state |
| `iyv` | the outcome columns | Party identification is stored plus one |
| `xv` | the covariate columns | Constant, age, schooling, female, race, income |
| `icoef == 0` | `layout$is_free == 1` | The coefficient is estimated |
| `erv` | factor draws | Person, factor, draw |
| `urv` | shock draws | Person, draw, attitude, wave |
| `idata = 3` | `validation/run_gss.R` | The GSS application |
| `model = 2` | the recursion in the core | The latent attitude carries forward |

The GSS attitudes, in order, are ideology, military spending, welfare spending, immigrants, party identification, and the environment. The covariates, in order, are the constant, age, high school, college, advanced degree, female, white, black, income in the middle band, and high income.
