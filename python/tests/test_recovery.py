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


def _fit_tiny(panel, maxit: int, draws: int):
    return attitudedyn.fit(
        panel,
        formula="cbind(y1, y2) ~ x",
        id_col="id",
        wave_col="wave",
        time_col="time",
        n_factors=1,
        draws=draws,
        maxit=maxit,
        reltol=1e-6,
        threads=1,
        vcov="opg",
        include_wave_effects=False,
        seed=1,
    )


def test_tiny_panel_reports_bias_rmse_and_coverage():
    """One short fit still produces finite bias, RMSE, and a coverage flag."""
    panel = attitudedyn.simulate_attitudes(
        n_persons=36, n_waves=3, n_attitudes=2, n_factors=1, n_categories=3, seed=15
    )
    truth = panel.attrs["theta"]
    result = _fit_tiny(panel, maxit=4, draws=3)
    estimate = result.params
    aligned_truth = truth.reindex(estimate.index).to_numpy()
    error = estimate.to_numpy() - aligned_truth
    critical = float(np.abs(norm_critical()))
    covered = np.abs(error) <= critical * result.bse.to_numpy()
    assert np.all(np.isfinite(error))
    assert np.sqrt(np.mean(error**2)) >= 0
    assert 0.0 <= float(covered.mean()) <= 1.0

    probabilities = result.predict(type="prob", n_draws=2, seed=3)
    totals = probabilities.groupby(["id", "wave", "attitude"], sort=False)["probability"].sum()
    assert np.allclose(totals.to_numpy(), 1.0)
    classes = result.predict(type="class", n_draws=2, seed=3)
    assert set(classes["y1"]).issubset({1, 2, 3})


def test_slope_bias_at_the_default_start_is_minus_one_half():
    panel = attitudedyn.simulate_attitudes(
        n_persons=20, n_waves=3, n_attitudes=2, n_factors=1, n_categories=3, seed=6
    )
    result = _fit_tiny(panel, maxit=0, draws=2)
    truth = float(panel.attrs["theta"]["beta[y1,x]"])
    estimate = float(result.params["beta[y1,x]"])
    assert truth == pytest.approx(0.5)
    assert estimate == pytest.approx(0.0)
    assert estimate - truth == pytest.approx(-0.5)


def norm_critical() -> float:
    """Two-sided 95 percent normal critical value, without importing scipy here."""
    return 1.959963984540054
