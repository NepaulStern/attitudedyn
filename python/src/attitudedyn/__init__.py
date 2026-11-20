"""Dynamic ordered-probit models for ordinal panel attitudes.

The numerical work is the C++ library in ``src/core/``. This package is a
wrapper: it prepares a long panel, calls that library for the likelihood
and its analytic gradient, and maximises the likelihood. BFGS follows
the same iteration as the R package. L-BFGS-B and Nelder-Mead use
``scipy.optimize``.
"""

from __future__ import annotations

from attitudedyn.fit import fit
from attitudedyn.result import AttitudeDynResult, ImpulseResponse
from attitudedyn.simulate import simulate_attitudes

__all__ = [
    "AttitudeDynResult",
    "ImpulseResponse",
    "fit",
    "simulate_attitudes",
]

__version__ = "0.1.0"
