"""Unconstrained starting values and the free/full coefficient maps."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd


def default_theta(layout: pd.DataFrame) -> np.ndarray:
    """Starting vector used by the replication package.

    Diagonal persistence starts at the unconstrained value whose natural
    image is 0.3, free loadings start at 0.5, and threshold increments
    start at ``log(0.6)``. Every other coefficient starts at 0.

    Parameters
    ----------
    layout :
        One row per coefficient, with columns ``block``, ``attitude``,
        ``index2``, and ``is_free``.

    Returns
    -------
    numpy.ndarray
        Unconstrained coefficients, in layout order.
    """
    theta = np.zeros(len(layout), dtype=np.float64)
    diagonal = (layout["block"].to_numpy() == "persistence") & (
        layout["attitude"].to_numpy() == layout["index2"].to_numpy()
    )
    theta[diagonal] = np.log((1.0 + 0.3) / (1.0 - 0.3))
    loading = (layout["block"].to_numpy() == "loading") & (layout["is_free"].to_numpy() == 1)
    theta[loading] = 0.5
    theta[layout["block"].to_numpy() == "threshold"] = np.log(0.6)
    return theta


def apply_fix(
    layout: pd.DataFrame, theta: np.ndarray, fix: Mapping[str, float] | None
) -> tuple[pd.DataFrame, np.ndarray]:
    """Hold named coefficients at the given unconstrained values.

    Parameters
    ----------
    layout, theta :
        Layout and unconstrained vector, modified on copies.
    fix :
        Map from coefficient name to the value to hold fixed. ``None``
        changes nothing.

    Returns
    -------
    layout, theta :
        Copies. Fixed rows have ``is_free`` 0.
    """
    layout = layout.copy()
    theta = np.array(theta, dtype=np.float64, copy=True)
    if not fix:
        return layout, theta
    unknown = [name for name in fix if name not in set(layout["name"])]
    if unknown:
        raise ValueError(f"Unknown coefficient in fix: {unknown}.")
    names = layout["name"].tolist()
    for name, value in fix.items():
        row = names.index(name)
        number = float(value)
        if not np.isfinite(number):
            raise ValueError(f"Fixed value for {name} must be finite.")
        theta[row] = number
        layout.loc[layout.index[row], "is_free"] = 0
    return layout, theta


def apply_start(
    layout: pd.DataFrame,
    theta: np.ndarray,
    start: Mapping[str, float] | Sequence[float] | np.ndarray | pd.Series | None,
) -> np.ndarray:
    """Insert user starting values into the unconstrained vector.

    A sequence is read in free-coefficient order. A mapping or a named
    Series may set any coefficient, free or fixed, by name.

    Parameters
    ----------
    layout, theta :
        Layout and the current unconstrained vector.
    start :
        Starting values, or ``None`` to keep ``theta``.

    Returns
    -------
    numpy.ndarray
        Copy of ``theta`` with the starting values written in.
    """
    theta = np.array(theta, dtype=np.float64, copy=True)
    if start is None:
        return theta
    names = layout["name"].tolist()
    free = layout["is_free"].to_numpy() == 1
    if isinstance(start, pd.Series):
        if start.index.dtype == object or any(isinstance(key, str) for key in start.index):
            start = {str(key): float(value) for key, value in start.items()}
        else:
            start = np.asarray(start, dtype=np.float64)
    if isinstance(start, Mapping):
        unknown = [name for name in start if name not in names]
        if unknown:
            raise ValueError(f"Unknown names in start: {unknown}.")
        for name, value in start.items():
            number = float(value)
            if not np.isfinite(number):
                raise ValueError(f"Starting value for {name} must be finite.")
            theta[names.index(name)] = number
        return theta
    values = np.asarray(list(start) if not isinstance(start, np.ndarray) else start, dtype=np.float64)
    if values.ndim != 1:
        raise ValueError("start must be a one-dimensional sequence or a map of coefficient names.")
    n_free = int(free.sum())
    if values.size != n_free:
        raise ValueError(f"start has length {values.size} but there are {n_free} free coefficients.")
    if not np.all(np.isfinite(values)):
        raise ValueError("start must be finite.")
    theta[free] = values
    return theta


def free_from_full(theta_full: np.ndarray, layout: pd.DataFrame) -> np.ndarray:
    """Return the free entries of an unconstrained vector, in layout order."""
    return np.asarray(theta_full, dtype=np.float64)[layout["is_free"].to_numpy() == 1].copy()


def full_from_free(theta_free: np.ndarray, theta_template: np.ndarray, layout: pd.DataFrame) -> np.ndarray:
    """Write free coefficients into a copy of the full template."""
    theta = np.array(theta_template, dtype=np.float64, copy=True)
    theta[layout["is_free"].to_numpy() == 1] = np.asarray(theta_free, dtype=np.float64)
    return theta
