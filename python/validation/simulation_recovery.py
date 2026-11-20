#!/usr/bin/env python3
"""Parameter recovery for a panel with the GSS wave dates.

The design matches ``validation/simulation_recovery.R``. From the repository
root, with the package installed::

    python python/validation/simulation_recovery.py --threads 4
    python python/validation/simulation_recovery.py --replications 1

``--notes`` rewrites ``notes/SIMULATION_RECOVERY.md``. The committed note is
the one written by the R script.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd

import attitudedyn


# Eight GSS wave dates, three ordinal attitudes, one factor, one covariate.
RECOVERY_DESIGN = {
    "number_of_persons": 400,
    "number_of_waves": 8,
    "number_of_attitudes": 3,
    "number_of_factors": 1,
    "number_of_categories": [5, 3, 3],
    "wave_dates": [2, 3, 4, 5, 6, 7, 8, 10],
    "number_of_draws": 8,
    "maximum_iterations": 100,
    "relative_tolerance": 1e-8,
    "starting_seed": 20261120,
}
NORMAL_CRITICAL_VALUE = 1.959963984540054


def natural_persistence(unconstrained_coefficient: np.ndarray) -> np.ndarray:
    """Map an unconstrained persistence coefficient onto (-1, 1)."""
    return 2.0 * np.exp(unconstrained_coefficient) / (1.0 + np.exp(unconstrained_coefficient)) - 1.0


def persistence_standard_error(natural_value: np.ndarray, unconstrained_standard_error: np.ndarray) -> np.ndarray:
    """Delta-method standard error. The derivative is 0.5 * (1 - rho^2)."""
    return np.abs(0.5 * (1.0 - natural_value**2)) * unconstrained_standard_error


def outcome_formula(number_of_attitudes: int) -> str:
    names = ", ".join(f"y{index}" for index in range(1, number_of_attitudes + 1))
    return f"cbind({names}) ~ x"


def simulate_recovery_panel(design: dict, replication_index: int) -> pd.DataFrame:
    return attitudedyn.simulate_attitudes(
        n_persons=design["number_of_persons"],
        n_waves=design["number_of_waves"],
        n_attitudes=design["number_of_attitudes"],
        n_factors=design["number_of_factors"],
        n_categories=design["number_of_categories"],
        seed=design["starting_seed"] + replication_index,
        wave_time=design["wave_dates"],
        include_wave_effects=False,
    )


def fit_recovery_replication(panel: pd.DataFrame, design: dict, number_of_threads: int):
    return attitudedyn.fit(
        panel,
        formula=outcome_formula(design["number_of_attitudes"]),
        id_col="id",
        wave_col="wave",
        time_col="time",
        n_factors=design["number_of_factors"],
        draws=design["number_of_draws"],
        include_wave_effects=False,
        maxit=design["maximum_iterations"],
        reltol=design["relative_tolerance"],
        threads=number_of_threads,
        vcov="opg",
        seed=1,
    )


def replication_rows(result, truth: pd.Series, replication_index: int) -> pd.DataFrame:
    """One row per free coefficient, with Wald coverage on both scales."""
    estimate = result.params
    standard_error = result.bse.reindex(estimate.index)
    truth_free = truth.reindex(estimate.index)
    lower = estimate.to_numpy() - NORMAL_CRITICAL_VALUE * standard_error.to_numpy()
    upper = estimate.to_numpy() + NORMAL_CRITICAL_VALUE * standard_error.to_numpy()
    persistence = estimate.index.str.startswith("persistence[")
    natural_truth = natural_persistence(truth_free.to_numpy())
    natural_estimate = natural_persistence(estimate.to_numpy())
    natural_se = persistence_standard_error(natural_estimate, standard_error.to_numpy())
    natural_lower = natural_estimate - NORMAL_CRITICAL_VALUE * natural_se
    natural_upper = natural_estimate + NORMAL_CRITICAL_VALUE * natural_se
    covered = (truth_free.to_numpy() >= lower) & (truth_free.to_numpy() <= upper)
    natural_covered = (natural_truth >= natural_lower) & (natural_truth <= natural_upper)
    return pd.DataFrame(
        {
            "replication": replication_index,
            "parameter": estimate.index.to_numpy(),
            "truth": truth_free.to_numpy(),
            "estimate": estimate.to_numpy(),
            "standard_error": standard_error.to_numpy(),
            "covered": covered,
            "natural_truth": np.where(persistence, natural_truth, np.nan),
            "natural_estimate": np.where(persistence, natural_estimate, np.nan),
            "natural_covered": np.where(persistence, natural_covered, np.nan),
            "converged": result.convergence == 0,
        }
    )


def summarize_recovery(detail: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for parameter_name, block in detail.groupby("parameter", sort=False):
        error = block["estimate"].to_numpy() - block["truth"].to_numpy()
        rows.append(
            {
                "parameter": parameter_name,
                "truth": float(block["truth"].iloc[0]),
                "mean_estimate": float(block["estimate"].mean()),
                "bias": float(error.mean()),
                "rmse": float(np.sqrt(np.mean(error**2))),
                "coverage": float(block["covered"].mean()),
            }
        )
    return pd.DataFrame(rows)


def summarize_natural_persistence(detail: pd.DataFrame) -> pd.DataFrame:
    block = detail[detail["parameter"].str.startswith("persistence[")]
    rows = []
    for parameter_name, slice_ in block.groupby("parameter", sort=False):
        error = slice_["natural_estimate"].to_numpy() - slice_["natural_truth"].to_numpy()
        rows.append(
            {
                "parameter": parameter_name,
                "truth": float(slice_["natural_truth"].iloc[0]),
                "mean_estimate": float(slice_["natural_estimate"].mean()),
                "bias": float(error.mean()),
                "rmse": float(np.sqrt(np.mean(error**2))),
                "coverage": float(slice_["natural_covered"].mean()),
            }
        )
    return pd.DataFrame(rows)


def format_recovery_table(summary: pd.DataFrame) -> str:
    lines = [
        "| Parameter | Truth | Mean estimate | Bias | RMSE | 95% coverage |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summary.itertuples(index=False):
        lines.append(
            f"| {row.parameter} | {row.truth:.4f} | {row.mean_estimate:.4f} | "
            f"{row.bias:.4f} | {row.rmse:.4f} | {row.coverage:.3f} |"
        )
    return "\n".join(lines)


def write_recovery_note(path: Path, design: dict, detail: pd.DataFrame, elapsed_seconds: float, number_of_threads: int) -> None:
    unconstrained = summarize_recovery(detail)
    natural = summarize_natural_persistence(detail)
    error = detail["estimate"].to_numpy() - detail["truth"].to_numpy()
    natural_rows = detail[detail["parameter"].str.startswith("persistence[")]
    natural_error = natural_rows["natural_estimate"].to_numpy() - natural_rows["natural_truth"].to_numpy()
    converged = int(detail.groupby("replication")["converged"].first().sum())
    note = "\n".join(
        [
            "# Parameter recovery",
            "",
            "Panels are drawn from `simulate_attitudes()` and refit.",
            "The wave dates are the GSS grid `c(2, 3, 4, 5, 6, 7, 8, 10)`.",
            "This file was written by `python/validation/simulation_recovery.py`.",
            "",
            "## Summary",
            "",
            f"- Replications: {design['number_of_replications']}",
            f"- Converged: {converged} of {design['number_of_replications']}",
            f"- Threads: {number_of_threads}",
            f"- Mean signed bias, unconstrained: {error.mean():.4f}",
            f"- RMSE, unconstrained: {np.sqrt(np.mean(error**2)):.4f}",
            f"- Coverage of unconstrained 95% intervals: {detail['covered'].mean():.3f}",
            f"- Mean signed bias, natural persistence: {natural_error.mean():.4f}",
            f"- RMSE, natural persistence: {np.sqrt(np.mean(natural_error**2)):.4f}",
            f"- Coverage of natural-persistence 95% intervals: {natural_rows['natural_covered'].mean():.3f}",
            f"- Elapsed seconds: {elapsed_seconds:.1f}",
            "",
            "## Unconstrained coefficients",
            "",
            format_recovery_table(unconstrained),
            "",
            "## Persistence on (-1, 1)",
            "",
            format_recovery_table(natural),
            "",
        ]
    )
    path.write_text(note + "\n", encoding="utf-8")


def run_recovery(number_of_replications: int, number_of_threads: int) -> tuple[pd.DataFrame, float]:
    design = dict(RECOVERY_DESIGN)
    design["number_of_replications"] = number_of_replications
    started = time.perf_counter()
    pieces = []
    for replication_index in range(1, number_of_replications + 1):
        panel = simulate_recovery_panel(design, replication_index)
        result = fit_recovery_replication(panel, design, number_of_threads)
        pieces.append(replication_rows(result, panel.attrs["theta"], replication_index))
        print(
            f"  {replication_index}/{number_of_replications}  "
            f"converged={int(result.convergence == 0)}  loglik={result.loglik:.2f}",
            flush=True,
        )
    detail = pd.concat(pieces, ignore_index=True)
    return detail, time.perf_counter() - started


def main() -> None:
    parser = argparse.ArgumentParser(description="Recover known attitude-dynamics parameters.")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--replications", type=int, default=40)
    parser.add_argument(
        "--notes",
        action="store_true",
        help="Write notes/SIMULATION_RECOVERY.md from this Python run.",
    )
    parsed = parser.parse_args()
    print(
        f"Recovery: {parsed.replications} replications, "
        f"{RECOVERY_DESIGN['number_of_persons']} persons, {parsed.threads} threads.",
        flush=True,
    )
    detail, elapsed_seconds = run_recovery(parsed.replications, parsed.threads)
    error = detail["estimate"].to_numpy() - detail["truth"].to_numpy()
    print(
        f"Bias {error.mean():.4f}  RMSE {np.sqrt(np.mean(error**2)):.4f}  "
        f"coverage {detail['covered'].mean():.3f}  elapsed {elapsed_seconds:.1f}s",
        flush=True,
    )
    if parsed.notes:
        root = Path(__file__).resolve().parents[2]
        design = dict(RECOVERY_DESIGN)
        design["number_of_replications"] = parsed.replications
        note_path = root / "notes" / "SIMULATION_RECOVERY.md"
        write_recovery_note(note_path, design, detail, elapsed_seconds, parsed.threads)
        print(f"Wrote {note_path}", flush=True)


if __name__ == "__main__":
    main()
