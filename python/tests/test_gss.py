"""GSS evaluation against the replication log and the Philox reference likelihood."""

from __future__ import annotations

import numpy as np
import pytest

import attitudedyn
from attitudedyn.gss import (
    COVARIATES,
    EXPECTED_COUNTS,
    OUTCOMES,
    PHILOX_LOGLIK,
    TIME_GRID,
    category_counts,
    fortran_names,
    gss_paths,
    read_gss,
    read_gss_estimate_table,
    read_gss_estimates,
    read_gss_variance,
    repository_root,
)

ROOT = repository_root()
DATA, LOG = gss_paths(ROOT)


@pytest.fixture(scope="module")
def gss_panel():
    panel = read_gss(DATA)
    return panel


def test_replication_log_parses():
    estimates = read_gss_estimates(LOG)
    assert len(estimates) == 164
    assert estimates.index[0] == "fixef101"
    assert estimates.iloc[0] == pytest.approx(2.1941423491)
    table = read_gss_estimate_table(LOG)
    assert len(table) == 164
    assert table.loc[0, "fortran"] == "fixef101"
    assert table.loc[0, "published_se"] == pytest.approx(0.052)
    variance = read_gss_variance(LOG)
    assert len(variance) == 30
    assert variance.loc[0, "component"] == "covariate"
    assert variance.loc[0, "published"] == pytest.approx(0.055)


def test_gss_category_counts_and_sample_size(gss_panel):
    assert gss_panel.attrs["n_people"] == 11282
    assert gss_panel.attrs["n_included"] == 11252
    counts = category_counts(gss_panel)
    for name, expected in EXPECTED_COUNTS.items():
        assert counts[name] == expected


def test_free_coefficient_order_matches_the_log(gss_panel):
    from attitudedyn import _core
    from attitudedyn.prepare import layout_frame, prepare_panel

    estimates = read_gss_estimates(LOG)
    prepared = prepare_panel(
        gss_panel,
        outcomes=OUTCOMES,
        covariates=COVARIATES,
        id_col="id",
        wave_col="wave",
        time_grid=TIME_GRID,
        weights=gss_panel["weight"].to_numpy(),
        add_intercept=False,
        n_factors=2,
    )
    layout = layout_frame(_core.layout(prepared.spec))
    names = fortran_names(layout)
    assert names == list(estimates.index)
    assert "fixef401" not in names
    assert len(names) == 164


def test_philox_loglik_and_gradient_match_across_threads(gss_panel):
    estimates = read_gss_estimates(LOG)
    kwargs = dict(
        outcomes=OUTCOMES,
        covariates=COVARIATES,
        id_col="id",
        wave_col="wave",
        time_grid=TIME_GRID,
        weights="weight",
        n_factors=2,
        draws=10,
        add_intercept=False,
        start=estimates.to_numpy(),
        maxit=0,
        vcov="none",
        keep_scores=False,
        use_gradient=True,
        seed=1,
        antithetic=True,
        draw_type="philox",
    )
    one = attitudedyn.fit(gss_panel, threads=1, **kwargs)
    four = attitudedyn.fit(gss_panel, threads=4, **kwargs)
    assert one.nobs == 11252
    # The reference value is one double. Its unit in the last place is about
    # 1.5e-11, so an absolute gap below 1e-12 is the same double. glibc on
    # the machine that produced the reference matches it. Another libm moves
    # the sum by a few ulps (about 1.5e-11 on a newer glibc, about 1.8e-10
    # with the Windows CRT) because erfc, exp, and log are not bit-identical.
    # 1e-9 still rejects a wrong likelihood and covers that library gap.
    assert one.loglik == pytest.approx(PHILOX_LOGLIK, abs=1e-9, rel=0)
    assert four.loglik == pytest.approx(PHILOX_LOGLIK, abs=1e-9, rel=0)
    assert one.loglik == four.loglik
    assert one.gradient is not None and four.gradient is not None
    assert one.gradient.shape == (164,)
    assert np.all(np.isfinite(one.gradient))
    assert np.array_equal(one.gradient, four.gradient)
