# Roadmap

The package is the C++ and R reimplementation of the replication program, checked on the GSS application and on small simulated panels.

## British election panel

Not in this version. Still to do: read that panel, use its time scale, and evaluate its likelihood with the same core. No download, conversion, or estimation of that dataset belongs in the package until that work is done.

## Other languages

The C++ core has no R headers, so the same library can sit under another front end.

- Python: a small extension module that accepts the same arrays the R wrapper builds and returns the log-likelihood, the gradient, and the post-estimation objects.
- Stata: a plugin with the same inputs, writing `e(b)` and `e(V)` from the unconstrained coefficients and one of the covariance estimators.

## Bayesian estimation

The likelihood and its gradient are in place. A later version can put a prior on the unconstrained coefficients and sample the posterior, using the same data layout and the same recursion. That sampler is not part of this package.
