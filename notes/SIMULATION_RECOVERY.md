# Simulation recovery

80 replications at each of 200, 400, and 1600 people. Every replication converged (optimiser code 0). No covariance was singular. The run used 4 threads and took 354.4 seconds.

```bash
python validation/simulation_recovery.py --replications 80 --draws 8 --threads 4 --maxit 50 --n-persons 200 400 1600
```

## Design

`simulate_attitudes` with 4 waves, 2 attitudes, 1 factor, 3 categories, and no wave effects. The data-generating coefficients, on the unconstrained scale, are:

| Coefficient | Value |
| --- | ---: |
| `beta[y1,(Intercept)]`, `beta[y2,(Intercept)]` | 0.2 |
| `beta[y1,x]`, `beta[y2,x]` | 0.5 |
| `persistence[y1,y1]`, `persistence[y2,y2]` | 0.9694005571881035 |
| `persistence[y1,y2]`, `persistence[y2,y1]` | 0 |
| `loading[y1,f1]`, `loading[y2,f1]` | 0.5 |
| `threshold[y1,3]`, `threshold[y2,3]` | -0.5108256237659907 |

The diagonal persistence value is `log((1 + 0.45) / (1 - 0.45))`. The threshold value is `log(0.6)`.

Estimation is BFGS with the analytic gradient, relative tolerance `1e-8`, at most 50 iterations, started at those coefficients. Draws are Philox, seed 1, 8 antithetic pairs. Replication `r` uses data seed `1000 + r` for `r = 0, ..., 79`.

A 95 percent Wald interval is the estimate plus or minus `1.959963984540054` standard errors. Coverage is the number of free coefficients, across the 80 replications, whose interval contains the data-generating value, divided by 960 (`80 * 12` coefficients).

`opg_centered` centres the person scores and inverts their outer product. `hessian` inverts the negative numerical Hessian of the analytic gradient. The difference step is `1e-5 * max(1, abs(coefficient))`, and the matrix is symmetrised. `sandwich` is that inverse, times the uncentred outer product of the person scores, times the inverse.

## Coverage

| People | Seconds | opg_centered | hessian | sandwich |
| ---: | ---: | ---: | ---: | ---: |
| 200 | 34.5 | 877/960 = 0.913542 | 886/960 = 0.922917 | 907/960 = 0.944792 |
| 400 | 83.0 | 850/960 = 0.885417 | 869/960 = 0.905208 | 880/960 = 0.916667 |
| 1600 | 236.9 | 767/960 = 0.798958 | 793/960 = 0.826042 | 818/960 = 0.852083 |

Covered coefficients out of 80 replications:

| Coefficient | n=200 centred | n=200 hessian | n=200 sandwich | n=400 centred | n=400 hessian | n=400 sandwich | n=1600 centred | n=1600 hessian | n=1600 sandwich |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `beta[y1,(Intercept)]` | 76 | 76 | 76 | 76 | 76 | 76 | 71 | 73 | 76 |
| `beta[y1,x]` | 77 | 77 | 79 | 77 | 76 | 78 | 71 | 72 | 73 |
| `beta[y2,(Intercept)]` | 73 | 75 | 77 | 76 | 77 | 78 | 75 | 75 | 75 |
| `beta[y2,x]` | 78 | 78 | 76 | 74 | 75 | 75 | 72 | 72 | 72 |
| `persistence[y1,y1]` | 72 | 72 | 75 | 72 | 75 | 75 | 66 | 71 | 74 |
| `persistence[y1,y2]` | 73 | 75 | 75 | 66 | 69 | 71 | 56 | 56 | 60 |
| `persistence[y2,y1]` | 67 | 67 | 72 | 64 | 67 | 68 | 49 | 53 | 55 |
| `persistence[y2,y2]` | 74 | 75 | 76 | 73 | 74 | 74 | 73 | 74 | 77 |
| `loading[y1,f1]` | 67 | 69 | 74 | 62 | 66 | 69 | 41 | 47 | 49 |
| `loading[y2,f1]` | 70 | 70 | 73 | 61 | 66 | 68 | 46 | 52 | 61 |
| `threshold[y1,3]` | 75 | 76 | 78 | 76 | 75 | 75 | 72 | 72 | 71 |
| `threshold[y2,3]` | 75 | 76 | 76 | 73 | 73 | 73 | 75 | 76 | 75 |
