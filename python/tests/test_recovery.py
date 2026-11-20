"""A fit on simulated data should recover the data-generating parameters."""

from __future__ import annotations

import numpy as np
import pytest

import attitudedyn


def test_simulated_panel_recovers_persistence_and_slope():
    panel = attitudedyn.simulate_attitudes(
        n_persons=150, n_waves=4, n_attitudes=2, n_factors=1, n_categories=3, seed=7
    )
    truth = panel.attrs["theta"]
    common = dict(
        formula="cbind(y1, y2) ~ x",
        id_col="id",
        wave_col="wave",
        time_col="time",
        n_factors=1,
        draws=6,
        threads=1,
        vcov="none",
        keep_scores=False,
        seed=1,
        include_wave_effects=False,
    )
    at_truth = attitudedyn.fit(panel, maxit=0, start=truth, **common)
    fitted = attitudedyn.fit(panel, maxit=40, reltol=1e-8, **common)
    assert np.isfinite(at_truth.loglik)
    assert fitted.loglik + 1e-3 >= at_truth.loglik
    rho = fitted.persistence_matrix()
    assert rho.loc["y1", "y1"] == pytest.approx(0.45, abs=0.15)
    assert rho.loc["y2", "y2"] == pytest.approx(0.45, abs=0.15)
    assert fitted.params["beta[y1,x]"] == pytest.approx(0.5, abs=0.25)
    assert fitted.params["beta[y2,x]"] == pytest.approx(0.5, abs=0.25)
    assert np.all(np.isfinite(fitted.params.to_numpy()))
