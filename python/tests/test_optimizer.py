"""Trial steps outside the parameter domain must not abort a fit."""

from __future__ import annotations

import math

import numpy as np
import pytest

import attitudedyn
from attitudedyn.bfgs import minimize_bfgs


def _small():
    panel = attitudedyn.simulate_attitudes(
        n_persons=40, n_waves=3, n_attitudes=2, n_factors=1, n_categories=3, seed=4
    )
    result = attitudedyn.fit(
        panel,
        formula="cbind(y1, y2) ~ x",
        id_col="id",
        wave_col="wave",
        time_col="time",
        n_factors=1,
        draws=3,
        maxit=0,
        threads=1,
        vcov="none",
        keep_scores=False,
    )
    return panel, result


def test_overflowing_threshold_is_a_nonfinite_objective():
    _panel, result = _small()
    threshold = [name for name in result.params.index if name.startswith("threshold[")]
    assert threshold
    extreme = result.params.to_numpy().copy()
    extreme[list(result.params.index).index(threshold[0])] = 80
    evaluation = result._evaluate(extreme, want_gradient=True, want_scores=False)
    assert not np.isfinite(evaluation["objective"])
    assert len(evaluation["gradient"]) == len(extreme)
    assert np.all(np.isfinite(evaluation["gradient"]))

    nonfinite = result.params.to_numpy().copy()
    nonfinite[list(result.params.index).index(threshold[0])] = np.inf
    inf_eval = result._evaluate(nonfinite, want_gradient=True, want_scores=False)
    assert not np.isfinite(inf_eval["objective"])
    assert np.all(np.isfinite(inf_eval["gradient"]))

    with pytest.raises(ValueError, match="starting values"):
        attitudedyn.fit(
            _panel,
            formula="cbind(y1, y2) ~ x",
            id_col="id",
            wave_col="wave",
            time_col="time",
            n_factors=1,
            draws=3,
            maxit=0,
            threads=1,
            vcov="none",
            start=extreme,
        )


def _overflow_start():
    panel = attitudedyn.simulate_attitudes(
        n_persons=120, n_waves=3, n_attitudes=2, n_factors=1, n_categories=3, seed=4
    )
    base = attitudedyn.fit(
        panel,
        formula="cbind(y1, y2) ~ x",
        id_col="id",
        wave_col="wave",
        time_col="time",
        n_factors=1,
        draws=2,
        maxit=0,
        threads=1,
        vcov="none",
        keep_scores=False,
    )
    start = base.params.copy()
    start.loc[start.index.str.startswith("threshold[")] = -8
    return panel, start


@pytest.mark.parametrize("optimizer", ["bfgs", "l-bfgs-b"])
def test_optimizers_reject_an_overflowing_threshold_trial(optimizer):
    panel, start = _overflow_start()
    result = attitudedyn.fit(
        panel,
        formula="cbind(y1, y2) ~ x",
        id_col="id",
        wave_col="wave",
        time_col="time",
        n_factors=1,
        draws=2,
        maxit=4,
        threads=1,
        vcov="none",
        keep_scores=False,
        optimizer=optimizer,
        start=start.to_numpy(),
    )
    assert np.isfinite(result.loglik)
    assert result.loglik > -1e6
    assert np.all(np.isfinite(result.params.to_numpy()))
    thresholds = result.params.loc[result.params.index.str.startswith("threshold[")]
    assert np.all(thresholds.to_numpy() <= 50)


def test_bfgs_matches_r_optim_on_a_smooth_function():
    # stats::optim(c(0, 0, 0), fn, gr, method = "BFGS",
    #              control = list(maxit = 40, reltol = 1e-8))
    # value 0.2822355573785619, convergence 0, 11 function and 8 gradient calls.
    def objective(point):
        return (
            (point[0] - 1.0) ** 2
            + 3.0 * (point[1] + 2.0) ** 2
            + 0.1 * (point[0] * point[1]) ** 2
            + math.sin(point[2]) ** 2
        )

    def slope(point):
        return np.array(
            [
                2.0 * (point[0] - 1.0) + 0.2 * point[0] * point[1] ** 2,
                6.0 * (point[1] + 2.0) + 0.2 * point[1] * point[0] ** 2,
                2.0 * math.sin(point[2]) * math.cos(point[2]),
            ]
        )

    fitted = minimize_bfgs(objective, slope, np.zeros(3), maxit=40, reltol=1e-8)
    assert fitted.status == 0
    assert fitted.nfev == 11
    assert fitted.njev == 8
    assert fitted.fun == pytest.approx(0.2822355573785619, abs=0.0, rel=0.0)
    assert fitted.x == pytest.approx(
        [7.2125051478962587e-01, -1.9659099529748947e00, 0.0],
        abs=0.0,
        rel=0.0,
    )
