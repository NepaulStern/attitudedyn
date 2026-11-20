"""Draw a panel from the structural model."""

from __future__ import annotations

import numpy as np
import pandas as pd

from attitudedyn import _core
from attitudedyn.coefficients import default_theta
from attitudedyn.prepare import layout_frame
from attitudedyn.result import _panel_frame


def simulate_attitudes(
    n_persons: int = 100,
    n_waves: int = 4,
    n_attitudes: int = 2,
    n_factors: int = 1,
    n_categories: int | list[int] = 3,
    seed: int = 1,
    wave_time: list[float] | None = None,
    include_wave_effects: bool = False,
) -> pd.DataFrame:
    """Simulate a panel from the structural model.

    Covariates are an intercept and one standard normal regressor. Outcomes
    are named ``y1``, ``y2``, and so on. Category 1 is the lowest category.
    The returned frame is ready for :func:`attitudedyn.fit` with
    ``formula="cbind(y1, y2) ~ x"``. The data-generating persistence on the
    diagonal is 0.45, the slope on ``x`` is 0.5, and the intercept is 0.2.

    Parameters
    ----------
    n_persons, n_waves, n_attitudes, n_factors :
        Dimensions of the panel.
    n_categories :
        Categories per attitude. One integer is recycled.
    seed :
        Seed for the counter-based generator.
    wave_time :
        Wave dates. The default is ``1, 2, ..., n_waves``.
    include_wave_effects :
        If true, the data-generating process includes wave effects. The
        default matches the examples, which do not.

    Returns
    -------
    pandas.DataFrame
        Columns ``id``, ``wave``, ``time``, ``x``, and the outcomes.
        ``attrs["theta"]`` is the unconstrained parameter vector.
    """
    if n_persons < 1 or n_waves < 1 or n_attitudes < 1:
        raise ValueError("n_persons, n_waves, and n_attitudes must be positive.")
    if n_factors < 0:
        raise ValueError("n_factors must be nonnegative.")
    if isinstance(n_categories, int):
        categories = [int(n_categories)] * n_attitudes
    else:
        categories = [int(value) for value in n_categories]
        if len(categories) != n_attitudes:
            raise ValueError("n_categories must be one integer or one integer per attitude.")
    if wave_time is None:
        grid = [float(t) for t in range(1, n_waves + 1)]
    else:
        grid = [float(value) for value in wave_time]
        if len(grid) != n_waves:
            raise ValueError("wave_time must have one date per wave.")
    attitude_names = [f"y{j + 1}" for j in range(n_attitudes)]
    covariate_names = ["(Intercept)", "x"]
    spec = {
        "n_waves": int(n_waves),
        "n_attitudes": int(n_attitudes),
        "n_covariates": 2,
        "n_factors": int(n_factors),
        "n_categories": np.asarray(categories, dtype=np.int32),
        "wave_time": np.asarray(grid, dtype=np.float64),
        "include_wave_effects": bool(include_wave_effects),
        "free_contemporaneous": False,
        "fix_last_loading": n_factors >= 2,
        "fix_intercepts": False,
        "intercept_columns": np.asarray([0], dtype=np.int32),
        "empty_wave_effect": np.zeros(0, dtype=np.int32),
        "attitude_names": attitude_names,
        "covariate_names": covariate_names,
        "wave_names": [str(t) for t in range(1, n_waves + 1)],
    }
    layout = layout_frame(_core.layout(spec))
    theta = default_theta(layout)
    slope = (layout["block"].to_numpy() == "beta") & (layout["index2"].to_numpy() == 1)
    intercept = (layout["block"].to_numpy() == "beta") & (layout["index2"].to_numpy() == 0)
    diagonal = (layout["block"].to_numpy() == "persistence") & (
        layout["attitude"].to_numpy() == layout["index2"].to_numpy()
    )
    theta[slope] = 0.5
    theta[intercept] = 0.2
    theta[diagonal] = np.log((1.0 + 0.45) / (1.0 - 0.45))
    panel = _core.simulate_panel(spec, theta, int(n_persons), int(seed))
    frame = _panel_frame(panel, attitude_names, covariate_names)
    frame.attrs["theta"] = pd.Series(theta, index=layout["name"].tolist(), name="theta")
    return frame
