"""The Nepaul-Stern Model of Attitude Stability and Constraint.

Estimates a model of several Likert-scale attitudes measured on the same
people across survey waves. The latent attitude carries forward from one
wave to the next, and the recorded answer is the Likert category that
latent attitude falls into. The same model is available in the R package
``attitudedyn``.
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
