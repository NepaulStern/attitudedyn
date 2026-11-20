# GSS exhibits

The tables compare this package with the replication log
`validation/fortran/src/gss_log.txt`. The point is the replication estimates (the X vector printed before iteration 8).
Draws are Philox, seed 1, with 10 antithetic pairs. Standard errors are
the centred outer product of the person scores. The persistence matrix
is on (-1, 1), and the impulse responses are its powers.

| Quantity | Median absolute difference | Maximum absolute difference |
| --- | ---: | ---: |
| Coefficients | 0 | 0 |
| Standard errors | 0.00051589 | 0.00515438 |
| Persistence matrix | 0 | 0 |
| Variance decomposition | 0.000272752 | 0.000489453 |
| Impulse responses | 0 | 0 |

Log-likelihood -111570.9702226093. People in the likelihood: 11252. Free coefficients: 164.
Spectral radius 0.636678 (replication coefficients 0.636678).
Threads: 4. Elapsed: 7.7 seconds.

