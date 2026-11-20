"""Small-panel checks for the public fit API."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import attitudedyn
from attitudedyn.prepare import parse_formula


def _panel():
    return attitudedyn.simulate_attitudes(
        n_persons=40, n_waves=3, n_attitudes=2, n_factors=1, n_categories=3, seed=4
    )


def test_formula_parser_drops_intercept_and_keeps_backticks():
    outcomes, covariates, add_intercept = parse_formula(
        "cbind(ideology, military) ~ 0 + `(Intercept)` + age + female"
    )
    assert outcomes == ["ideology", "military"]
    assert covariates == ["(Intercept)", "age", "female"]
    assert add_intercept is False
    outcomes, covariates, add_intercept = parse_formula("cbind(y1, y2) ~ x")
    assert outcomes == ["y1", "y2"]
    assert covariates == ["x"]
    assert add_intercept is True


def test_bad_inputs_are_rejected():
    panel = _panel()
    with pytest.raises(ValueError, match="not in data"):
        attitudedyn.fit(panel, outcomes=["y1", "y2"], covariates=["x"], id_col="missing", wave_col="wave")
    shifted = panel.copy()
    shifted["y1"] = shifted["y1"] - 1
    with pytest.raises(ValueError, match="1, 2"):
        attitudedyn.fit(
            shifted, outcomes=["y1", "y2"], covariates=["x"], id_col="id", wave_col="wave", time_col="time"
        )
    with pytest.raises(ValueError, match="draws"):
        attitudedyn.fit(
            panel,
            outcomes=["y1", "y2"],
            covariates=["x"],
            id_col="id",
            wave_col="wave",
            time_col="time",
            draws=0,
            maxit=0,
            vcov="none",
        )
    broken = panel.copy()
    broken.loc[0, "x"] = np.inf
    with pytest.raises(ValueError, match="covariate"):
        attitudedyn.fit(
            broken, outcomes=["y1", "y2"], covariates=["x"], id_col="id", wave_col="wave", time_col="time"
        )


def test_thread_count_does_not_change_the_fit():
    panel = _panel()
    kwargs = dict(
        outcomes=["y1", "y2"],
        covariates=["x"],
        id_col="id",
        wave_col="wave",
        time_col="time",
        n_factors=1,
        draws=4,
        maxit=0,
        vcov="none",
        keep_scores=False,
        seed=1,
    )
    one = attitudedyn.fit(panel, threads=1, **kwargs)
    four = attitudedyn.fit(panel, threads=4, **kwargs)
    assert one.loglik == four.loglik
    assert np.array_equal(one.params.to_numpy(), four.params.to_numpy())
    assert one.gradient is not None and four.gradient is not None
    assert np.array_equal(one.gradient, four.gradient)


def test_analytic_gradient_matches_a_central_difference():
    panel = _panel()
    result = attitudedyn.fit(
        panel,
        formula="cbind(y1, y2) ~ x",
        id_col="id",
        wave_col="wave",
        time_col="time",
        n_factors=1,
        draws=4,
        maxit=0,
        threads=1,
        vcov="none",
        keep_scores=False,
    )
    assert result.gradient is not None
    step = 1e-5
    up = result.params.to_numpy().copy()
    down = up.copy()
    up[0] += step
    down[0] -= step
    plus = result._evaluate(up, want_gradient=False, want_scores=False)["objective"]
    minus = result._evaluate(down, want_gradient=False, want_scores=False)["objective"]
    numeric = (plus - minus) / (2 * step)
    assert result.gradient[0] == pytest.approx(numeric, abs=1e-4, rel=1e-4)


def test_result_methods_on_a_short_fit():
    panel = _panel()
    result = attitudedyn.fit(
        panel,
        formula="cbind(y1, y2) ~ x",
        id_col="id",
        wave_col="wave",
        time_col="time",
        n_factors=1,
        draws=4,
        maxit=2,
        threads=1,
        vcov="opg",
    )
    assert np.isfinite(result.loglik)
    assert result.nobs == panel["id"].nunique()
    assert np.all(np.isfinite(result.params.to_numpy()))
    variance = result.vcov()
    assert variance.shape == (len(result.params), len(result.params))
    assert result.bse.shape == (len(result.params),)
    assert result.pvalues.shape == (len(result.params),)
    interval = result.conf_int()
    assert interval.shape == (len(result.params), 2)
    table = result.to_dataframe()
    assert {"estimate", "std_error", "pvalue"}.issubset(table.columns)
    text = result.summary()
    assert "Log-likelihood" in text
    rho = result.persistence_matrix()
    assert rho.shape == (2, 2)
    assert np.all(np.abs(rho.to_numpy()) < 1)
    response = result.impulse_response(horizon=2)
    assert response.response.shape == (3, 2, 2)
    assert np.allclose(response.response[0], np.eye(2))
    parts = result.variance_decomposition()
    assert len(parts) == 2
    predicted = result.predict(type="class", n_draws=2)
    assert set(predicted["y1"]).issubset({1, 2, 3})
    latent = result.predict(type="latent", n_draws=2)
    assert {"y1", "y2"}.issubset(latent.columns)
    probabilities = result.predict(type="prob", n_draws=2)
    assert {"attitude", "category", "probability"}.issubset(probabilities.columns)
    drawn = result.simulate(nsim=1, seed=3)
    assert isinstance(drawn, pd.DataFrame)
    assert {"id", "wave", "y1", "y2"}.issubset(drawn.columns)
    assert np.isfinite(result.aic) and np.isfinite(result.bic)
