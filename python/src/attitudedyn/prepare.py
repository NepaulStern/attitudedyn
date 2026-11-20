"""Turn a long panel into the arrays the likelihood uses.

One row is one person in one wave. A missing outcome is a missing item and
is integrated out. A wave that does not appear for a person, but sits
between that person's first and last observed waves, is a gap: the latent
state is simulated through it. Covariates on that wave are used if the row
was read and are zero if the wave is absent.

Categories are the integers ``1, ..., K``. The cut between category 1 and
category 2 is normalised at zero, so the numbering is part of the model.
People stay in first-appearance order, which is the order of the
random-number streams.
"""

from __future__ import annotations

import re
import warnings
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

_INTERCEPT_NAMES = {"(Intercept)", "Intercept"}


@dataclass
class PreparedPanel:
    """Arrays passed to the C++ problem builder.

    Attributes
    ----------
    spec, data, draws :
        Dictionaries consumed by :mod:`attitudedyn._core`.
    layout inputs :
        Names, category counts, and the person index needed to label results.
    """

    spec: dict
    data: dict
    person_ids: np.ndarray
    wave_values: np.ndarray
    attitude_names: list[str]
    covariate_names: list[str]
    n_categories: np.ndarray
    n_persons: int
    n_waves: int
    n_attitudes: int
    n_covariates: int


def parse_formula(formula: str) -> tuple[list[str], list[str], bool]:
    """Split a small formula into outcome names, covariate names, and an intercept flag.

    The left-hand side is one name or ``cbind(y1, y2)``. The right-hand side
    is a sum of column names. ``0``, ``+ 0``, or ``- 1`` drops the intercept.
    Names with characters such as parentheses may be wrapped in backticks.
    Functions of columns are not expanded; build those columns first.

    Parameters
    ----------
    formula :
        For example ``"cbind(y1, y2) ~ x"`` or ``"cbind(y1, y2) ~ 0 + `(Intercept)` + age"``.

    Returns
    -------
    outcomes, covariates, add_intercept :
        Column names, and whether a column of ones named ``(Intercept)``
        should be added when the covariate list does not already contain one.
    """
    if not isinstance(formula, str) or "~" not in formula:
        raise ValueError("formula must be a string containing '~', for example 'cbind(y1, y2) ~ x'.")
    lhs, rhs = formula.split("~", 1)
    lhs = lhs.strip()
    rhs = rhs.strip()
    if not lhs:
        raise ValueError("formula is missing a left-hand side.")
    if lhs.startswith("cbind(") and lhs.endswith(")"):
        outcomes = [_unbacktick(part) for part in _split_top_level(lhs[len("cbind(") : -1], ",")]
    else:
        outcomes = [_unbacktick(lhs)]
    outcomes = [name for name in outcomes if name]
    if not outcomes:
        raise ValueError("The left-hand side must name at least one outcome.")
    add_intercept = True
    covariates: list[str] = []
    for part in _split_top_level(rhs, "+"):
        token = part.strip()
        if not token:
            continue
        collapsed = re.sub(r"\s+", "", token)
        if collapsed in {"0", "-1"}:
            add_intercept = False
            continue
        if collapsed == "1":
            add_intercept = True
            continue
        if collapsed.endswith("-1") and not token.startswith("`"):
            add_intercept = False
            token = token[: token.rfind("-")].strip()
            if not token:
                continue
        if len(token) >= 2 and token[0] == "`" and token[-1] == "`":
            name = token[1:-1]
            if not name:
                raise ValueError("A backticked formula name is empty.")
        else:
            name = token
            if any(character in name for character in "()+~"):
                raise ValueError(
                    f"Formula term {token!r} is not a bare column name. "
                    "Build transformed columns first, or pass covariates as a list."
                )
        covariates.append(name)
    return outcomes, covariates, add_intercept


def _split_top_level(text: str, separator: str) -> list[str]:
    parts: list[str] = []
    buffer: list[str] = []
    tick = False
    depth = 0
    for character in text:
        if character == "`":
            tick = not tick
            buffer.append(character)
            continue
        if not tick and character == "(":
            depth += 1
        elif not tick and character == ")":
            depth = max(0, depth - 1)
        if character == separator and not tick and depth == 0:
            parts.append("".join(buffer))
            buffer = []
        else:
            buffer.append(character)
    parts.append("".join(buffer))
    return parts


def _unbacktick(token: str) -> str:
    token = token.strip()
    if len(token) >= 2 and token[0] == "`" and token[-1] == "`":
        return token[1:-1]
    return token


def _wave_name(value: object) -> str:
    number = float(value)  # type: ignore[arg-type]
    if np.isfinite(number) and number == int(number):
        return str(int(number))
    return str(value)


def prepare_panel(
    data: pd.DataFrame,
    outcomes: Sequence[str],
    covariates: Sequence[str],
    id_col: str,
    wave_col: str,
    *,
    time_col: str | None = None,
    time_grid: Sequence[float] | None = None,
    weights: np.ndarray | None = None,
    add_intercept: bool = True,
    include_wave_effects: bool = True,
    n_factors: int = 2,
    free_contemporaneous: bool = False,
    fix_last_loading: bool = True,
    fix_intercept: bool = True,
    na_action: str = "pass",
) -> PreparedPanel:
    """Pack a long data frame into the core's person arrays.

    Parameters
    ----------
    data :
        Long-format frame, one row per person-wave.
    outcomes, covariates :
        Column names. Outcomes are ordinal and coded ``1, ..., K``.
    id_col, wave_col :
        Person identifier and wave. Both are required and cannot be missing.
    time_col, time_grid :
        Wave dates. ``time_grid`` is one date per distinct wave, in sorted
        wave order, and takes precedence over ``time_col``. Dates must be
        finite and strictly increasing. With neither argument the dates are
        ``1, 2, ..., T`` and a warning is issued.
    weights :
        Non-negative person weights already aligned with ``data``, constant
        within each person. ``None`` uses weight 1.
    add_intercept :
        Add a column of ones named ``(Intercept)`` when the covariate list
        does not already name an intercept.
    include_wave_effects, n_factors, free_contemporaneous, fix_last_loading, fix_intercept :
        Specification flags passed through to the core. An intercept is held
        at zero only when wave effects are included and ``fix_intercept`` is
        true.
    na_action :
        ``"pass"`` rejects missing covariates. ``"omit"`` drops person-waves
        whose covariates are missing. Missing outcomes are never dropped.

    Returns
    -------
    PreparedPanel
        Spec, flat arrays, and the labels needed to report a fit.
    """
    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas DataFrame.")
    if na_action not in {"pass", "omit"}:
        raise ValueError("na_action must be 'pass' or 'omit'.")
    frame = data.reset_index(drop=True)
    for column, label in ((id_col, "id_col"), (wave_col, "wave_col")):
        if column not in frame.columns:
            raise ValueError(f"Column {column!r} ({label}) is not in data.")
    outcomes = list(outcomes)
    covariates = list(covariates)
    if not outcomes:
        raise ValueError("Name at least one outcome.")
    missing_outcomes = [name for name in outcomes if name not in frame.columns]
    if missing_outcomes:
        raise ValueError(f"Outcome columns are not in data: {missing_outcomes}.")
    if add_intercept and not any(name in _INTERCEPT_NAMES for name in covariates):
        covariates = ["(Intercept)", *covariates]
        frame["(Intercept)"] = 1.0
    if not covariates:
        raise ValueError("Pass at least one covariate, or set add_intercept=True.")
    missing_covariates = [name for name in covariates if name not in frame.columns]
    if missing_covariates:
        raise ValueError(f"Covariate columns are not in data: {missing_covariates}.")

    design = frame.loc[:, covariates].apply(pd.to_numeric, errors="coerce")
    covariate_missing = ~np.isfinite(design.to_numpy(dtype=np.float64))
    if covariate_missing.any():
        if na_action != "omit":
            raise ValueError(
                "Covariates contain missing or non-finite values. "
                "Missing outcomes are integrated out and can stay missing. "
                "Use na_action='omit' to drop person-waves with missing covariates, or impute them."
            )
        keep = ~covariate_missing.any(axis=1)
        frame = frame.loc[keep].reset_index(drop=True)
        design = design.loc[keep].reset_index(drop=True)
        if weights is not None:
            weights = np.asarray(weights, dtype=np.float64)[np.asarray(keep)]

    if frame[id_col].isna().any() or frame[wave_col].isna().any():
        raise ValueError("id and wave cannot contain missing values.")
    if not pd.api.types.is_numeric_dtype(frame[wave_col]):
        raise ValueError("wave must be numeric.")

    design_values = design.to_numpy(dtype=np.float64)
    id_vector = frame[id_col].to_numpy()
    wave_vector = frame[wave_col].to_numpy()
    person_ids = pd.unique(id_vector)
    n_persons = int(len(person_ids))
    id_index = pd.Series(np.arange(n_persons), index=person_ids)
    person_index = id_index.loc[id_vector].to_numpy(dtype=np.int32)

    wave_values = np.unique(wave_vector)
    n_waves = int(wave_values.size)
    wave_index_map = pd.Series(np.arange(n_waves), index=wave_values)
    wave_index = wave_index_map.loc[wave_vector].to_numpy(dtype=np.int32)
    if pd.MultiIndex.from_arrays([person_index, wave_index]).duplicated().any():
        raise ValueError("Some person-wave pairs appear more than once.")

    if time_grid is not None:
        grid = np.asarray(list(time_grid), dtype=np.float64)
        if grid.size != n_waves:
            raise ValueError(f"time_grid must have one date for each distinct wave ({n_waves}).")
    elif time_col is not None:
        if time_col not in frame.columns:
            raise ValueError(f"Column {time_col!r} (time_col) is not in data.")
        time_vector = pd.to_numeric(frame[time_col], errors="coerce").to_numpy(dtype=np.float64)
        grid = np.empty(n_waves, dtype=np.float64)
        for wave in range(n_waves):
            slice_ = time_vector[wave_index == wave]
            if slice_.size == 0 or not np.all(np.isfinite(slice_)) or np.ptp(slice_) > 1e-8:
                raise ValueError(f"Time is missing or not constant within wave {wave_values[wave]!r}.")
            grid[wave] = float(np.mean(slice_))
    else:
        grid = np.arange(1, n_waves + 1, dtype=np.float64)
        warnings.warn(
            "No wave dates were supplied, so gaps are 1 between successive wave codes. "
            "The GSS estimates use time_grid = [2, 3, 4, 5, 6, 7, 8, 10].",
            UserWarning,
            stacklevel=2,
        )
    if not np.all(np.isfinite(grid)) or np.any(np.diff(grid) <= 0):
        raise ValueError("Wave dates must be finite and strictly increasing.")

    n_attitudes = len(outcomes)
    category_matrix = np.zeros((len(frame), n_attitudes), dtype=np.int32)
    n_categories = np.empty(n_attitudes, dtype=np.int32)
    for j, name in enumerate(outcomes):
        column = frame[name]
        _fill_outcome(column, name, category_matrix[:, j], n_categories, j)

    n_covariates = len(covariates)
    interviewed = np.zeros((n_persons, n_waves), dtype=np.uint8)
    categories = np.zeros((n_persons, n_attitudes, n_waves), dtype=np.int32)
    covariate_cube = np.zeros((n_persons, n_waves, n_covariates), dtype=np.float64)
    interviewed[person_index, wave_index] = 1
    for j in range(n_attitudes):
        categories[person_index, j, wave_index] = category_matrix[:, j]
    covariate_cube[person_index, wave_index, :] = design_values

    if weights is None:
        person_weight = np.ones(n_persons, dtype=np.float64)
    else:
        weight_vector = np.asarray(weights, dtype=np.float64)
        if weight_vector.shape != (len(frame),):
            raise ValueError("weights must have one value per row of data.")
        if not np.all(np.isfinite(weight_vector)) or np.any(weight_vector < 0):
            raise ValueError("Weights must be finite and non-negative, and constant across a person's waves.")
        lower = np.full(n_persons, np.inf)
        upper = np.full(n_persons, -np.inf)
        np.minimum.at(lower, person_index, weight_vector)
        np.maximum.at(upper, person_index, weight_vector)
        if np.any(upper - lower > 1e-10):
            bad = int(np.flatnonzero(upper - lower > 1e-10)[0])
            raise ValueError(
                f"Weights for person {person_ids[bad]!r} must be finite, non-negative, "
                "and constant across waves."
            )
        person_weight = lower

    first_wave = np.full(n_persons, n_waves, dtype=np.int32)
    last_wave = np.full(n_persons, -1, dtype=np.int32)
    np.minimum.at(first_wave, person_index, wave_index)
    np.maximum.at(last_wave, person_index, wave_index)

    observed = (categories >= 1).any(axis=0)
    empty_wave = np.where(observed.reshape(-1), 0, 1).astype(np.int32)
    intercept_columns = [index for index, name in enumerate(covariates) if name in _INTERCEPT_NAMES]

    if int(np.sum(person_weight > 0)) < 2:
        raise ValueError("At least two people with positive weight are required.")
    if n_waves < 2:
        warnings.warn(
            "The panel has a single wave, so persistence is not informed by within-person change.",
            UserWarning,
            stacklevel=2,
        )

    spec = {
        "n_waves": n_waves,
        "n_attitudes": n_attitudes,
        "n_covariates": n_covariates,
        "n_factors": int(n_factors),
        "n_categories": n_categories,
        "wave_time": np.asarray(grid, dtype=np.float64),
        "include_wave_effects": bool(include_wave_effects),
        "free_contemporaneous": bool(free_contemporaneous),
        "fix_last_loading": bool(fix_last_loading),
        "fix_intercepts": bool(fix_intercept) and bool(include_wave_effects),
        "intercept_columns": np.asarray(intercept_columns, dtype=np.int32),
        "empty_wave_effect": empty_wave,
        "attitude_names": list(outcomes),
        "covariate_names": list(covariates),
        "wave_names": [_wave_name(value) for value in wave_values],
    }
    packed = {
        "first_wave": first_wave,
        "last_wave": last_wave,
        "interviewed": interviewed.reshape(-1),
        "category": np.ascontiguousarray(categories).reshape(-1),
        "covariates": np.ascontiguousarray(covariate_cube).reshape(-1),
        "weight": np.asarray(person_weight, dtype=np.float64),
    }
    return PreparedPanel(
        spec=spec,
        data=packed,
        person_ids=np.asarray(person_ids),
        wave_values=np.asarray(wave_values),
        attitude_names=list(outcomes),
        covariate_names=list(covariates),
        n_categories=n_categories.copy(),
        n_persons=n_persons,
        n_waves=n_waves,
        n_attitudes=n_attitudes,
        n_covariates=n_covariates,
    )


def _fill_outcome(
    column: pd.Series,
    name: str,
    destination: np.ndarray,
    n_categories: np.ndarray,
    attitude: int,
) -> None:
    observed = column.notna().to_numpy()
    if not observed.any():
        raise ValueError(f"Outcome {name!r} is entirely missing.")
    values = column.to_numpy()
    if isinstance(column.dtype, pd.CategoricalDtype) or str(column.dtype) == "category":
        codes = column.cat.codes.to_numpy()
        if np.any(codes[observed] < 0):
            raise ValueError(f"Factor outcome {name!r} has an empty level.")
        destination[observed] = codes[observed] + 1
        n_categories[attitude] = int(len(column.cat.categories))
        return
    numeric = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=np.float64)
    seen = numeric[observed]
    if not np.all(np.isfinite(seen)) or np.any(np.abs(seen - np.round(seen)) > 1e-8) or np.any(seen < 1):
        raise ValueError(
            f"Outcome {name!r} must be coded as integers 1, 2, ..., K. "
            "The cut between category 1 and category 2 is normalised at zero."
        )
    coded = np.round(seen).astype(np.int32)
    k = int(coded.max())
    present = set(np.unique(coded).tolist())
    if present != set(range(1, k + 1)):
        raise ValueError(
            f"Outcome {name!r} skips a category code between 1 and {k}. "
            f"Observed codes: {sorted(present)}."
        )
    destination[observed] = coded
    n_categories[attitude] = k


def layout_frame(raw: dict) -> pd.DataFrame:
    """Build a layout table from the dictionary returned by the core."""
    return pd.DataFrame(
        {
            "name": list(raw["name"]),
            "block": list(raw["block"]),
            "attitude": np.asarray(raw["attitude"], dtype=np.int32),
            "index2": np.asarray(raw["index2"], dtype=np.int32),
            "is_free": np.asarray(raw["is_free"], dtype=np.int32),
        }
    )
