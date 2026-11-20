# Parameter recovery

Panels are drawn from `simulate_attitudes()` and refit with `attitudedyn()`.
The wave dates are the GSS grid `c(2, 3, 4, 5, 6, 7, 8, 10)`.
There are three ordinal attitudes, with 5, 3, and 3 categories, one
standard normal covariate, an intercept, and one person factor.
Wave effects are off, so each intercept is a free coefficient.

The data-generating values are the ones `simulate_attitudes()` uses:
diagonal persistence 0.45, off-diagonal persistence 0, slope 0.5,
intercept 0.2, loading 0.5, and threshold increments `log(0.6)`.
The slope, intercept, loading, and threshold rows below are on the
unconstrained scale, which is the scale of `coef()` and of the Wald
interval. Persistence is shown on both scales. The natural-scale
interval uses the delta method, with standard error
`|0.5 * (1 - rho^2)|` times the unconstrained standard error.

## Design

| Setting | Value |
| --- | ---: |
| Persons | 400 |
| Waves | 8 |
| Wave dates | 2, 3, 4, 5, 6, 7, 8, 10 |
| Attitudes | 3 |
| Categories | 5, 3, 3 |
| Factors | 1 |
| Draws (antithetic pairs) | 8 |
| Maximum BFGS iterations | 100 |
| Replications | 40 |
| Threads | 4 |
| First simulation seed | 20261121 |

## Summary

- Replications: 40
- Converged: 40 of 40
- Free coefficients: 23
- Mean signed bias, unconstrained: -0.0015
- RMSE, unconstrained: 0.0505
- Coverage of unconstrained 95% intervals: 0.897
- Mean signed bias, natural persistence: 0.0072
- RMSE, natural persistence: 0.0298
- Coverage of natural-persistence 95% intervals: 0.892
- Elapsed seconds: 87.4

Bias is the mean of the estimate minus the data-generating value.
RMSE is the square root of the mean squared error, pooled over
coefficients and replications in the summary and computed within each
coefficient in the tables. Coverage is the share of replications whose
interval contains the data-generating value. Intervals are outer-product
Wald intervals.

## Unconstrained coefficients

| Parameter | Truth | Mean estimate | Bias | RMSE | 95% coverage |
| --- | ---: | ---: | ---: | ---: | ---: |
| beta[y1,(Intercept)] | 0.2000 | 0.1942 | -0.0058 | 0.0431 | 0.825 |
| beta[y1,x] | 0.5000 | 0.4921 | -0.0079 | 0.0259 | 0.900 |
| beta[y2,(Intercept)] | 0.2000 | 0.1930 | -0.0070 | 0.0393 | 0.800 |
| beta[y2,x] | 0.5000 | 0.4862 | -0.0138 | 0.0272 | 0.925 |
| beta[y3,(Intercept)] | 0.2000 | 0.1909 | -0.0091 | 0.0426 | 0.800 |
| beta[y3,x] | 0.5000 | 0.4931 | -0.0069 | 0.0189 | 1.000 |
| persistence[y1,y1] | 0.9694 | 0.9559 | -0.0135 | 0.0580 | 0.925 |
| persistence[y1,y2] | 0.0000 | 0.0250 | 0.0250 | 0.0544 | 0.900 |
| persistence[y1,y3] | 0.0000 | 0.0452 | 0.0452 | 0.0694 | 0.775 |
| persistence[y2,y1] | 0.0000 | 0.0429 | 0.0429 | 0.0659 | 0.850 |
| persistence[y2,y2] | 0.9694 | 0.9442 | -0.0252 | 0.0647 | 1.000 |
| persistence[y2,y3] | 0.0000 | 0.0122 | 0.0122 | 0.0555 | 0.950 |
| persistence[y3,y1] | 0.0000 | 0.0367 | 0.0367 | 0.0639 | 0.825 |
| persistence[y3,y2] | 0.0000 | 0.0361 | 0.0361 | 0.0687 | 0.900 |
| persistence[y3,y3] | 0.9694 | 0.9257 | -0.0437 | 0.0692 | 0.900 |
| loading[y1,f1] | 0.5000 | 0.4904 | -0.0096 | 0.0479 | 0.850 |
| loading[y2,f1] | 0.5000 | 0.4963 | -0.0037 | 0.0432 | 0.900 |
| loading[y3,f1] | 0.5000 | 0.4852 | -0.0148 | 0.0430 | 0.925 |
| threshold[y1,3] | -0.5108 | -0.5232 | -0.0124 | 0.0358 | 0.950 |
| threshold[y1,4] | -0.5108 | -0.5220 | -0.0112 | 0.0477 | 0.925 |
| threshold[y1,5] | -0.5108 | -0.5226 | -0.0118 | 0.0480 | 0.950 |
| threshold[y2,3] | -0.5108 | -0.5275 | -0.0167 | 0.0393 | 0.950 |
| threshold[y3,3] | -0.5108 | -0.5306 | -0.0197 | 0.0448 | 0.900 |

## Persistence on (-1, 1)

| Parameter | Truth | Mean estimate | Bias | RMSE | 95% coverage |
| --- | ---: | ---: | ---: | ---: | ---: |
| persistence[y1,y1] | 0.4500 | 0.4443 | -0.0057 | 0.0234 | 0.925 |
| persistence[y1,y2] | 0.0000 | 0.0125 | 0.0125 | 0.0272 | 0.875 |
| persistence[y1,y3] | 0.0000 | 0.0226 | 0.0226 | 0.0347 | 0.775 |
| persistence[y2,y1] | 0.0000 | 0.0215 | 0.0215 | 0.0329 | 0.850 |
| persistence[y2,y2] | 0.4500 | 0.4396 | -0.0104 | 0.0261 | 1.000 |
| persistence[y2,y3] | 0.0000 | 0.0061 | 0.0061 | 0.0277 | 0.950 |
| persistence[y3,y1] | 0.0000 | 0.0183 | 0.0183 | 0.0319 | 0.825 |
| persistence[y3,y2] | 0.0000 | 0.0180 | 0.0180 | 0.0343 | 0.900 |
| persistence[y3,y3] | 0.4500 | 0.4322 | -0.0178 | 0.0282 | 0.925 |

## Reproduce

```bash
Rscript validation/simulation_recovery.R --threads 4
python python/validation/simulation_recovery.py --threads 4
```

The table above was written by `validation/simulation_recovery.R`.
`python/validation/simulation_recovery.py` repeats the same design.
A one-replication check of the bias, RMSE, and coverage calculations
lives in `tests/testthat/test-recovery.R` and `python/tests/test_recovery.py`.

