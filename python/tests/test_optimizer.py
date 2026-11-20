"""Trial steps outside the parameter domain must not abort a fit."""

from __future__ import annotations

import numpy as np
import pytest

import attitudedyn


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
