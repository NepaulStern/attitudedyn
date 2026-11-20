#!/usr/bin/env python3
"""GSS tables and the impulse-response figure, without a refit.

From the repository root, with the package installed::

    python python/validation/replicate_gss_exhibits.py --threads 4

The default call evaluates the replication estimates (Philox, seed 1, 10
antithetic draws) and computes centred outer-product standard errors. It
does not start the optimiser. ``--fit`` displays a coefficient vector saved
by ``run_gss.py --refit`` instead of optimising it again.

Tables and ``gss_impulse_response.png`` go to ``validation/output/``. The
figure is also copied to ``man/figures/`` for the README.
``notes/GSS_EXHIBITS.md`` records the evaluation.
"""

from __future__ import annotations

import argparse
import shutil
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import attitudedyn
from attitudedyn import _core
from attitudedyn.gss import (
    COVARIATES,
    EXPECTED_COUNTS,
    OUTCOMES,
    TIME_GRID,
    category_counts,
    fortran_names,
    ensure_replication_files,
    natural_persistence,
    read_gss,
    read_gss_estimate_table,
    read_gss_estimates,
    read_gss_variance,
    repository_root,
)


def significance_stars(estimate: float, standard_error: float) -> str:
    """Normal significance marks at 10, 5, and 1 percent."""
    if not np.isfinite(standard_error) or standard_error == 0.0 or not np.isfinite(estimate):
        return ""
    ratio = abs(estimate / standard_error)
    if ratio > 2.576:
        return "***"
    if ratio > 1.960:
        return "**"
    if ratio > 1.645:
        return "*"
    return ""


def fit_at_coefficients(panel: pd.DataFrame, coefficients: np.ndarray, threads: int):
    """Likelihood, scores, and post-estimation at a fixed coefficient vector."""
    return attitudedyn.fit(
        panel,
        outcomes=OUTCOMES,
        covariates=COVARIATES,
        id_col="id",
        wave_col="wave",
        time_grid=TIME_GRID,
        weights="weight",
        n_factors=2,
        draws=10,
        add_intercept=False,
        start=np.asarray(coefficients, dtype=np.float64),
        threads=threads,
        maxit=0,
        vcov="opg_centered",
        use_gradient=True,
        keep_scores=True,
        seed=1,
    )


def impulse_table(result, published_persistence: pd.DataFrame, horizon: int = 8) -> tuple[pd.DataFrame, float, float]:
    """Long impulse-response table and the two spectral radii."""
    fitted = result.impulse_response(horizon=horizon)
    published = _core.impulse_response(
        np.ascontiguousarray(published_persistence.to_numpy(dtype=np.float64)),
        published_persistence.shape[0],
        horizon,
    )
    published_response = np.asarray(published["response"], dtype=np.float64).reshape(
        horizon + 1, published_persistence.shape[0], published_persistence.shape[0]
    )
    rows = []
    for to_index, receiver in enumerate(fitted.attitudes):
        for from_index, source in enumerate(fitted.attitudes):
            for step in range(horizon + 1):
                estimate = float(fitted.response[step, to_index, from_index])
                reference = float(published_response[step, to_index, from_index])
                rows.append(
                    {
                        "to": receiver,
                        "from": source,
                        "horizon": step,
                        "estimate": estimate,
                        "published": reference,
                        "difference": estimate - reference,
                    }
                )
    return pd.DataFrame(rows), float(fitted.spectral_radius), float(published["spectral_radius"])


def exhibit_frames(result, log_path: Path) -> dict:
    """Coefficient, persistence, variance, and impulse-response comparisons."""
    layout_names = fortran_names(result.layout)
    shipped = read_gss_estimates(log_path)
    if layout_names != list(shipped.index):
        raise RuntimeError("Free-coefficient order does not match the replication log.")
    table = read_gss_estimate_table(log_path)
    if list(table["fortran"]) != list(shipped.index):
        raise RuntimeError("The estimates table is not in the same order as the X VECTOR.")
    estimates = result.params.to_numpy()
    try:
        standard_error = np.sqrt(np.clip(np.diag(result.vcov(type="opg_centered").to_numpy()), 0.0, np.inf))
    except (ValueError, np.linalg.LinAlgError):
        standard_error = np.full(estimates.shape, np.nan)
    coefficients = pd.DataFrame(
        {
            "fortran": list(shipped.index),
            "term": list(result.params.index),
            "estimate": estimates,
            "published": shipped.to_numpy(),
            "difference": estimates - shipped.to_numpy(),
            "std.error": standard_error,
            "published_se": table["published_se"].to_numpy(),
            "published_t": table["published_t"].to_numpy(),
        }
    )
    coefficients["se_difference"] = coefficients["std.error"] - coefficients["published_se"]
    with np.errstate(divide="ignore", invalid="ignore"):
        coefficients["statistic"] = coefficients["estimate"] / coefficients["std.error"]
    coefficients["stars"] = [
        significance_stars(estimate, error)
        for estimate, error in zip(coefficients["estimate"], coefficients["std.error"])
    ]

    published_matrix = natural_persistence(shipped.filter(regex=r"^rho").to_numpy()).reshape(6, 6)
    persistence = result.persistence_matrix()
    published_persistence = pd.DataFrame(
        published_matrix, index=persistence.index, columns=persistence.columns
    )
    persistence_out = persistence.copy()
    persistence_out.insert(0, "attitude", persistence_out.index)
    published_out = published_persistence.copy()
    published_out.insert(0, "attitude", published_out.index)

    published_variance = read_gss_variance(log_path)
    fitted_variance = result.variance_decomposition()
    rows = []
    for attitude_row in fitted_variance.itertuples(index=False):
        for component in ("covariate", "wave_effect", "factor", "shock", "total"):
            published = published_variance.loc[
                (published_variance["component"] == component)
                & (published_variance["attitude_index"] == OUTCOMES.index(attitude_row.attitude) + 1),
                "published",
            ].iloc[0]
            rows.append(
                {
                    "attitude": attitude_row.attitude,
                    "component": component,
                    "estimate": getattr(attitude_row, component),
                    "published": float(published),
                }
            )
    variance = pd.DataFrame(rows)
    variance["difference"] = variance["estimate"] - variance["published"]
    impulse, radius, published_radius = impulse_table(result, published_persistence)
    return {
        "coefficients": coefficients,
        "persistence": persistence_out,
        "published_persistence": published_out,
        "variance": variance,
        "impulse": impulse,
        "spectral_radius": radius,
        "published_spectral_radius": published_radius,
    }


def comparison_markdown(frames: dict, *, loglik: float, n_persons: int, n_parameters: int,
                        elapsed: float, threads: int, point: str) -> str:
    """Record of the GSS evaluation written to notes/GSS_EXHIBITS.md."""
    body = [
        "# GSS exhibits",
        "",
        f"Evaluation at the {point}.",
        "Draws are Philox, seed 1, with 10 antithetic pairs. The persistence matrix",
        "is on (-1, 1), and the impulse responses are its powers.",
        "",
        f"Log-likelihood {loglik:.10f}. People in the likelihood: {n_persons}. "
        f"Free coefficients: {n_parameters}.",
        f"Spectral radius {frames['spectral_radius']:.6f}.",
        f"Threads: {threads}. Elapsed: {elapsed:.1f} seconds.",
        "",
        "Regenerate with `Rscript validation/replicate_gss_exhibits.R` or",
        "`python python/validation/replicate_gss_exhibits.py`.",
        "The General Social Survey file comes from the CC0 Harvard Dataverse deposit",
        "doi:10.7910/DVN/DUCHYS and is downloaded into `validation/fortran/vendor/`.",
        "",
    ]
    return "\n".join(body)


def summary_lines(frames: dict, loglik: float) -> str:
    """Machine-readable companion to the markdown table."""
    coefficient_gap = np.abs(frames["coefficients"]["difference"].to_numpy())
    error_gap = np.abs(frames["coefficients"]["se_difference"].to_numpy())
    persistence_gap = np.abs(
        frames["persistence"].drop(columns="attitude").to_numpy()
        - frames["published_persistence"].drop(columns="attitude").to_numpy()
    )
    variance_gap = np.abs(frames["variance"]["difference"].to_numpy())
    impulse_gap = np.abs(frames["impulse"]["difference"].to_numpy())
    lines = [
        f"loglik {loglik:.16f}",
        f"coefficient_abs_diff_max {np.max(coefficient_gap):.6g}",
        f"coefficient_abs_diff_median {np.median(coefficient_gap):.6g}",
        f"se_abs_diff_max {np.nanmax(error_gap):.6g}",
        f"se_abs_diff_median {np.nanmedian(error_gap):.6g}",
        f"persistence_abs_diff_max {np.max(persistence_gap):.6g}",
        f"persistence_abs_diff_median {np.median(persistence_gap):.6g}",
        f"variance_abs_diff_max {np.nanmax(variance_gap):.6g}",
        f"variance_abs_diff_median {np.nanmedian(variance_gap):.6g}",
        f"impulse_abs_diff_max {np.nanmax(impulse_gap):.6g}",
        f"impulse_abs_diff_median {np.nanmedian(impulse_gap):.6g}",
        f"spectral_radius {frames['spectral_radius']:.6f}",
        f"spectral_radius_published {frames['published_spectral_radius']:.6f}",
    ]
    return "\n".join(lines) + "\n"


def plot_impulse_responses(frame: pd.DataFrame, path: Path) -> None:
    """Six-by-six grid. Rows receive the shock named on the left. Columns are the source."""
    attitudes = list(dict.fromkeys(frame["to"]))
    horizon = np.sort(frame["horizon"].unique())
    largest_gap = float(np.nanmax(np.abs(frame["difference"].to_numpy())))
    show_replication = largest_gap > 1e-3
    figure, axes = plt.subplots(len(attitudes), len(attitudes), figsize=(11.2, 9.4), sharex=True)
    for row, receiver in enumerate(attitudes):
        for column, source in enumerate(attitudes):
            panel = frame[(frame["to"] == receiver) & (frame["from"] == source)].sort_values("horizon")
            axis = axes[row, column]
            axis.axhline(0.0, color="#b0b0b0", linewidth=0.6, zorder=0)
            axis.plot(panel["horizon"], panel["estimate"], color="#1f4e79", linewidth=1.6, label="package")
            if show_replication:
                axis.plot(
                    panel["horizon"],
                    panel["published"],
                    color="#b85c38",
                    linewidth=1.15,
                    linestyle="--",
                    label="replication log",
                )
            axis.tick_params(labelsize=7, length=2)
            axis.set_xlim(float(horizon[0]), float(horizon[-1]))
            if row == len(attitudes) - 1:
                axis.set_xlabel(source, fontsize=8)
            else:
                axis.tick_params(labelbottom=False)
            if column == 0:
                axis.set_ylabel(receiver, fontsize=8)
            else:
                axis.tick_params(labelleft=False)
    figure.suptitle("GSS impulse responses", fontsize=14)
    figure.supxlabel("Shock to", fontsize=11)
    figure.supylabel("Response of", fontsize=11)
    if show_replication:
        handles, labels = axes[0, 0].get_legend_handles_labels()
        figure.legend(handles, labels, loc="upper right", frameon=False, fontsize=8)
    figure.tight_layout(rect=(0.02, 0.02, 1.0, 0.96))
    figure.savefig(path, dpi=140)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description="Write the GSS exhibits from the replication estimates.")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument(
        "--fit",
        type=Path,
        default=None,
        help="Optional .npz with a 'coefficients' array. The vector is evaluated, not re-optimised.",
    )
    args = parser.parse_args()
    root = repository_root()
    data_path, log_path = ensure_replication_files(root)
    panel = read_gss(data_path)
    counts = category_counts(panel)
    for name, expected in EXPECTED_COUNTS.items():
        if counts[name] != expected:
            raise SystemExit(f"Category counts for {name} do not match the replication log.")
    if args.fit is None:
        coefficients = read_gss_estimates(log_path).to_numpy()
        point = "replication estimates (the X vector printed before iteration 8)"
    else:
        coefficients = np.asarray(np.load(args.fit)["coefficients"], dtype=np.float64)
        point = f"saved coefficients in {args.fit.name}"
    started = time.perf_counter()
    result = fit_at_coefficients(panel, coefficients, args.threads)
    elapsed = time.perf_counter() - started
    frames = exhibit_frames(result, log_path)
    out_dir = root / "validation" / "output"
    out_dir.mkdir(parents=True, exist_ok=True)
    frames["coefficients"].to_csv(out_dir / "gss_coefficients.csv", index=False)
    frames["persistence"].to_csv(out_dir / "gss_persistence.csv", index=False)
    frames["published_persistence"].to_csv(out_dir / "gss_persistence_published.csv", index=False)
    frames["variance"].to_csv(out_dir / "gss_variance.csv", index=False)
    frames["impulse"].to_csv(out_dir / "gss_impulse_response.csv", index=False)
    figure_path = out_dir / "gss_impulse_response.png"
    plot_impulse_responses(frames["impulse"], figure_path)
    # The README embeds this copy. validation/output/ is not part of the package.
    figure_dir = root / "man" / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(figure_path, figure_dir / "gss_impulse_response.png")
    (out_dir / "gss_summary.txt").write_text(result.summary(type="opg_centered") + "\n")
    text = summary_lines(frames, result.loglik)
    (out_dir / "gss_exhibit_comparison.txt").write_text(text)
    note = comparison_markdown(
        frames,
        loglik=result.loglik,
        n_persons=result.nobs,
        n_parameters=len(result.params),
        elapsed=elapsed,
        threads=args.threads,
        point=point,
    )
    note_path = root / "notes" / "GSS_EXHIBITS.md"
    note_path.write_text(note)
    print(text, end="")
    print(f"elapsed_seconds {elapsed:.3f}")
    print(f"threads {args.threads}")
    print(f"Wrote exhibits to {out_dir}")
    print(f"Wrote {note_path}")


if __name__ == "__main__":
    main()
