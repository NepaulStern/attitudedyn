"""GSS panel used by the replication program.

People stay in first-appearance order, including people the likelihood
skips, because the random-number streams are keyed by that index. A wave
inside a person's window that was not in the file keeps a zero regressor,
except for age, which is filled from the first wave whose age is at least
16: ``age_k = age_j + 2 * (k - j)``, with ``k`` and ``j`` the wave indices
``1..8``.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

NSEG = np.array([7, 3, 3, 5, 7, 3], dtype=np.int32)
TIME_GRID = np.array([2, 3, 4, 5, 6, 7, 8, 10], dtype=np.float64)
OUTCOMES = ["ideology", "military", "welfare", "immigrants", "party", "environment"]
COVARIATES = [
    "(Intercept)",
    "age",
    "hs",
    "college",
    "adv",
    "female",
    "white",
    "black",
    "inc_mid",
    "inc_high",
]
EXPECTED_COUNTS = {
    "ideology": [865, 2668, 2479, 7762, 2933, 3333, 870],
    "military": [3566, 3893, 2998],
    "welfare": [4180, 3746, 2461],
    "immigrants": [2907, 2963, 4898, 1354, 639],
    "party": [3771, 3585, 2767, 4145, 2067, 2816, 2300],
    "environment": [890, 2954, 6655],
}
PHILOX_LOGLIK = -111570.97022260929

_RAW_COLUMNS = [
    "id",
    "year",
    "age",
    "educ",
    "ideology",
    "race",
    "sex",
    "panel",
    "idate",
    "month",
    "trvdate",
    "military",
    "junk1",
    "environment",
    "welfare",
    "junk3",
    "junk4",
    "immigrants",
    "party",
    "income",
]
_BETA_STEM = ["constnt", "age", "hs", "college", "advdegr", "female", "white", "black"]
_ESTIMATE_ROW = re.compile(
    r"^\s*\d+\s+([A-Za-z0-9]+)\s+(-?[0-9.]+)\**\s+([0-9.]+)\s+(-?[0-9.]+)"
)
_NUMBER = re.compile(r"-?[0-9]+\.[0-9]+")


def repository_root() -> Path:
    """Return the checkout that contains ``core/`` and the GSS file.

    ``ATTDYN_REPOSITORY_ROOT`` overrides the search. Wheel tests set it
    because pytest's working directory is outside the checkout.
    """
    override = os.environ.get("ATTDYN_REPOSITORY_ROOT")
    if override:
        candidate = Path(override)
        data = candidate / "validation" / "fortran" / "data" / "gss2.prn"
        core = candidate / "core" / "include" / "attitudedyn" / "api.hpp"
        if not (data.exists() and core.exists()):
            raise FileNotFoundError(
                f"ATTDYN_REPOSITORY_ROOT={override} does not contain the GSS file and the C++ core."
            )
        return candidate
    starts = [Path.cwd(), Path(__file__).resolve()]
    for origin in starts:
        for candidate in [origin, *origin.parents]:
            data = candidate / "validation" / "fortran" / "data" / "gss2.prn"
            core = candidate / "core" / "include" / "attitudedyn" / "api.hpp"
            if data.exists() and core.exists():
                return candidate
    raise FileNotFoundError(
        "Could not find the attitudedyn repository (expected validation/fortran/data/gss2.prn). "
        "Run from a checkout of the repository."
    )


def gss_paths(root: Path | None = None) -> tuple[Path, Path]:
    """Return ``(gss2.prn, gss_log.txt)``."""
    base = repository_root() if root is None else root
    return (
        base / "validation" / "fortran" / "data" / "gss2.prn",
        base / "validation" / "fortran" / "src" / "gss_log.txt",
    )


def read_gss(path: Path | None = None, n_rows: int = 21648) -> pd.DataFrame:
    """Read ``gss2.prn`` the way the replication loader does.

    Parameters
    ----------
    path :
        File path. The default is the copy shipped under ``validation/fortran/data``.
    n_rows :
        Number of person-wave rows the program reads.

    Returns
    -------
    pandas.DataFrame
        Long panel. ``attrs["n_people"]`` is everyone in the file and
        ``attrs["n_included"]`` is everyone with positive weight.
    """
    if path is None:
        path = gss_paths()[0]
    lines: list[str] = []
    with path.open() as handle:
        for line in handle:
            stripped = line.lstrip()
            if stripped[:1].isdigit() or (stripped.startswith("-") and len(stripped) > 1 and stripped[1].isdigit()):
                lines.append(line)
            if len(lines) >= n_rows:
                break
    if len(lines) < n_rows:
        raise ValueError(f"gss2.prn has {len(lines)} data rows; the program reads {n_rows}.")
    raw = np.loadtxt(lines[:n_rows])
    if raw.shape[1] != len(_RAW_COLUMNS):
        raise ValueError(f"Expected {len(_RAW_COLUMNS)} columns in gss2.prn, found {raw.shape[1]}.")
    year = raw[:, 1].astype(np.int64)
    wave = ((year - 2006) // 2) + 1
    if np.any(wave < 1) or np.any(wave > 8):
        raise ValueError("A GSS year falls outside waves 1..8.")
    # First-appearance order. np.unique would sort, and the random-number
    # streams are keyed by this index.
    order: dict[float, int] = {}
    for value in raw[:, 0]:
        if value not in order:
            order[value] = len(order)
    ids = np.asarray(list(order.keys()))
    if np.all(ids == np.round(ids)):
        ids = np.round(ids).astype(np.int64)
    person_index = np.asarray([order[value] for value in raw[:, 0]], dtype=np.int32)
    n_people = int(ids.size)
    n_waves = 8
    y = np.zeros((n_people, 6, n_waves), dtype=np.int32)
    x = np.zeros((n_people, n_waves, 10), dtype=np.float64)
    for row in range(raw.shape[0]):
        person = int(person_index[row])
        t = int(wave[row]) - 1
        y[person, 0, t] = int(raw[row, 4])
        y[person, 1, t] = int(raw[row, 11])
        y[person, 2, t] = int(raw[row, 14])
        y[person, 3, t] = int(raw[row, 17])
        y[person, 4, t] = int(raw[row, 18]) + 1
        y[person, 5, t] = int(raw[row, 13])
        x[person, t, 0] = 1.0
        x[person, t, 1] = raw[row, 2]
        educ = int(raw[row, 3])
        if educ in (1, 2):
            x[person, t, 2] = 1.0
        if educ == 3:
            x[person, t, 3] = 1.0
        if educ == 4:
            x[person, t, 4] = 1.0
        x[person, t, 5] = raw[row, 6] - 1.0
        if int(raw[row, 5]) == 1:
            x[person, t, 6] = 1.0
        if int(raw[row, 5]) == 2:
            x[person, t, 7] = 1.0
        if int(raw[row, 19]) == 2:
            x[person, t, 8] = 1.0
        if int(raw[row, 19]) == 3:
            x[person, t, 9] = 1.0

    included = np.zeros(n_people, dtype=bool)
    for person in range(n_people):
        for t in range(n_waves):
            for j in range(6):
                value = int(y[person, j, t])
                if 1 <= value <= int(NSEG[j]):
                    included[person] = True

    for person in range(n_people):
        donor = -1
        for t in range(n_waves):
            if x[person, t, 1] >= 16:
                donor = t
                break
        if donor < 0:
            included[person] = False
        else:
            age0 = x[person, donor, 1]
            for k in range(n_waves):
                if x[person, k, 1] < 16:
                    x[person, k, 1] = age0 + 2 * (k - donor)
    for person in range(n_people):
        donor = -1
        for t in range(n_waves):
            female = x[person, t, 5]
            if 0.0 <= female <= 1.0:
                donor = t
                break
        if donor < 0:
            included[person] = False
        else:
            female = x[person, donor, 5]
            for k in range(n_waves):
                if x[person, k, 5] < 0:
                    x[person, k, 5] = female

    start = np.zeros(n_people, dtype=np.int32)
    end = np.zeros(n_people, dtype=np.int32)
    for person in range(n_people):
        if not included[person]:
            continue
        for t in range(n_waves):
            if np.any(y[person, :, t] > 0):
                if start[person] == 0:
                    start[person] = t + 1
                end[person] = t + 1
        if start[person] == 0:
            included[person] = False

    rows: list[pd.DataFrame] = []
    for person in range(n_people):
        if not included[person]:
            rows.append(
                pd.DataFrame(
                    {
                        "id": [ids[person]],
                        "wave": [1],
                        "weight": [0.0],
                        "(Intercept)": [0.0],
                        "age": [0.0],
                        "hs": [0.0],
                        "college": [0.0],
                        "adv": [0.0],
                        "female": [0.0],
                        "white": [0.0],
                        "black": [0.0],
                        "inc_mid": [0.0],
                        "inc_high": [0.0],
                        "ideology": [np.nan],
                        "military": [np.nan],
                        "welfare": [np.nan],
                        "immigrants": [np.nan],
                        "party": [np.nan],
                        "environment": [np.nan],
                    }
                )
            )
            continue
        span = np.arange(int(start[person]) - 1, int(end[person]))
        block = {
            "id": np.full(span.size, ids[person]),
            "wave": span + 1,
            "weight": np.ones(span.size),
            "(Intercept)": x[person, span, 0],
            "age": x[person, span, 1],
            "hs": x[person, span, 2],
            "college": x[person, span, 3],
            "adv": x[person, span, 4],
            "female": x[person, span, 5],
            "white": x[person, span, 6],
            "black": x[person, span, 7],
            "inc_mid": x[person, span, 8],
            "inc_high": x[person, span, 9],
        }
        for j, name in enumerate(OUTCOMES):
            values = y[person, j, span].astype(np.float64)
            invalid = (values < 1) | (values > int(NSEG[j]))
            values = values.copy()
            values[invalid] = np.nan
            block[name] = values
        rows.append(pd.DataFrame(block))
    panel = pd.concat(rows, ignore_index=True)
    panel.attrs["n_people"] = n_people
    panel.attrs["n_included"] = int(included.sum())
    return panel


def category_counts(panel: pd.DataFrame) -> dict[str, list[int]]:
    """Category counts among people with positive weight."""
    keep = panel["weight"].to_numpy() > 0
    counts: dict[str, list[int]] = {}
    for name, nseg in zip(OUTCOMES, NSEG):
        values = panel.loc[keep, name].dropna().to_numpy()
        values = np.round(values).astype(np.int32)
        counts[name] = np.bincount(values, minlength=int(nseg) + 1)[1 : int(nseg) + 1].tolist()
    return counts


def read_gss_estimates(log_path: Path | None = None) -> pd.Series:
    """Free coefficients from the last ``X VECTOR`` before ``iteration = 8``."""
    if log_path is None:
        log_path = gss_paths()[1]
    lines = log_path.read_text(errors="replace").splitlines()
    marks = [index for index, line in enumerate(lines) if "iteration = 8" in line]
    if not marks:
        raise ValueError("The log has no iteration 8 line.")
    mark = marks[0]
    headers = [index for index, line in enumerate(lines) if "X VECTOR" in line and index < mark]
    if not headers:
        raise ValueError("The log has no X VECTOR before iteration 8.")
    names: list[str] = []
    values: list[float] = []
    for line in lines[headers[-1] + 1 : mark]:
        text = line.strip()
        if not text:
            continue
        parts = text.split()
        names.append(parts[0])
        values.append(float(parts[1].replace("D", "E").replace("d", "e")))
    return pd.Series(values, index=names, name="published")


def read_gss_estimate_table(log_path: Path | None = None) -> pd.DataFrame:
    """Point estimates and standard errors from the replication estimates table."""
    if log_path is None:
        log_path = gss_paths()[1]
    lines = log_path.read_text(errors="replace").splitlines()
    starts = [index for index, line in enumerate(lines) if "estimates and significance" in line]
    if not starts:
        raise ValueError("The log has no estimates table.")
    rows = []
    for line in lines[starts[-1] + 1 :]:
        found = _ESTIMATE_ROW.match(line)
        if found is None:
            if rows:
                break
            continue
        rows.append(
            {
                "fortran": found.group(1),
                "published": float(found.group(2)),
                "published_se": float(found.group(3)),
                "published_t": float(found.group(4)),
            }
        )
    if not rows:
        raise ValueError("The estimates table was empty.")
    return pd.DataFrame(rows)


def read_gss_variance(log_path: Path | None = None) -> pd.DataFrame:
    """Six variance blocks, in attitude order, from the end of the replication log."""
    if log_path is None:
        log_path = gss_paths()[1]
    lines = log_path.read_text(errors="replace").splitlines()
    starts = [index for index, line in enumerate(lines) if "variance decomposition for each dependent variable" in line]
    if not starts:
        raise ValueError("The log has no variance decomposition.")
    sizes: list[float] = []
    for line in lines[starts[-1] + 1 :]:
        found = _NUMBER.findall(line)
        if len(found) == 2:
            sizes.append(float(found[0]))
        elif len(found) == 1 and len(sizes) % 5 == 4:
            sizes.append(float(found[0]))
    if len(sizes) < 30:
        raise ValueError(f"Expected 6 variance blocks, found {len(sizes)} numbers.")
    sizes = sizes[:30]
    return pd.DataFrame(
        {
            "attitude_index": np.repeat(np.arange(1, 7), 5),
            "component": ["covariate", "wave_effect", "factor", "shock", "total"] * 6,
            "published": sizes,
        }
    )


def fortran_name(name: str, block: str, attitude: int, index2: int) -> str:
    """Replication-program name of one layout row.

    ``attitude`` and ``index2`` are the 0-based values stored by the core,
    except threshold cuts: those are read from the coefficient name, where
    ``tau_3`` is the first free cut.
    """
    attitude_1 = int(attitude) + 1
    index_1 = int(index2) + 1
    if block == "wave_effect":
        return f"fixef{attitude_1}{index_1:02d}"
    if block == "beta":
        if index_1 <= len(_BETA_STEM):
            return f"{_BETA_STEM[index_1 - 1]}{attitude_1}"
        if index_1 == 9:
            return f"income{attitude_1}1"
        if index_1 == 10:
            return f"income{attitude_1}2"
        raise ValueError(f"No Fortran name for covariate index {index_1}.")
    if block == "persistence":
        return f"rho{attitude_1}{index_1}"
    if block == "loading":
        return f"alamb{index_1}{attitude_1}"
    if block == "threshold":
        cut = int(name[name.rfind(",") + 1 : -1])
        return f"tau{attitude_1}{cut}"
    if block == "chol":
        return f"chol{attitude_1}{index_1}"
    raise ValueError(f"No Fortran name for block {block}.")


def fortran_names(layout: pd.DataFrame) -> list[str]:
    """Fortran names of the free coefficients, in layout order."""
    free = layout.loc[layout["is_free"].to_numpy() == 1]
    return [
        fortran_name(name, block, int(attitude), int(index2))
        for name, block, attitude, index2 in zip(free["name"], free["block"], free["attitude"], free["index2"])
    ]


def natural_persistence(unconstrained: np.ndarray) -> np.ndarray:
    """Map unconstrained persistence through ``2 exp(c) / (1 + exp(c)) - 1``."""
    unconstrained = np.asarray(unconstrained, dtype=np.float64)
    return 2.0 * np.exp(unconstrained) / (1.0 + np.exp(unconstrained)) - 1.0
