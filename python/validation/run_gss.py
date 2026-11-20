#!/usr/bin/env python3
"""Evaluate the GSS model at the replication estimates.

From the repository root, with the Python package installed::

    python python/validation/run_gss.py
    python python/validation/run_gss.py --threads 4
    python python/validation/run_gss.py --refit --threads 4 --maxit 200

The default call evaluates the shipped free coefficients with Philox draws
and does not start the optimiser. ``--refit`` continues from those values.
Category counts are checked against the replication log before the
likelihood is called.

The coefficient vector and the variance decomposition are read from
``validation/fortran/src/gss_log.txt`` (the replication log). That log is
the file the R script compares against; the repository does not ship a
separate ``gss_coefficients.txt`` or ``gss_variance.txt``.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd

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
    natural_persistence,
    read_gss,
    read_gss_estimate_table,
    read_gss_estimates,
    read_gss_variance,
    repository_root,
)


def comparison(result: attitudedyn.AttitudeDynResult, log_path: Path) -> dict[str, pd.DataFrame]:
    """Point estimates, persistence, and variance against the replication log."""
    layout_names = fortran_names(result.layout)
    shipped = read_gss_estimates(log_path)
    if layout_names != list(shipped.index):
        raise RuntimeError(
            "Free-coefficient order does not match the replication log.\n"
            f"First layout name: {layout_names[0]}\n"
            f"First log name: {shipped.index[0]}"
        )
    table = read_gss_estimate_table(log_path)
    if list(table["fortran"]) != list(shipped.index):
        raise RuntimeError("The estimates table is not in the same order as the X VECTOR.")
    estimates = result.params.to_numpy()
    try:
        se = np.sqrt(np.clip(np.diag(result.vcov(type="opg_centered").to_numpy()), 0.0, np.inf))
    except (ValueError, np.linalg.LinAlgError):
        se = np.full(estimates.shape, np.nan)
    coefficients = pd.DataFrame(
        {
            "fortran": list(shipped.index),
            "term": list(result.params.index),
            "estimate": estimates,
            "published": shipped.to_numpy(),
            "difference": estimates - shipped.to_numpy(),
            "std.error": se,
            "published_se": table["published_se"].to_numpy(),
        }
    )
    coefficients["se_difference"] = coefficients["std.error"] - coefficients["published_se"]
    published_rho = natural_persistence(shipped.filter(regex=r"^rho").to_numpy()).reshape(6, 6)
    fitted_rho = result.persistence_matrix()
    published_rho = pd.DataFrame(published_rho, index=fitted_rho.index, columns=fitted_rho.columns)
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
    return {
        "coefficients": coefficients,
        "persistence": fitted_rho,
        "published_persistence": published_rho,
        "variance": variance,
    }


def print_comparison(compared: dict[str, pd.DataFrame], loglik: float, convergence: int, elapsed: float) -> None:
    coefficients = compared["coefficients"]
    abs_diff = np.abs(coefficients["difference"].to_numpy())
    se_diff = np.abs(coefficients["se_difference"].to_numpy())
    rho_diff = np.abs(compared["persistence"].to_numpy() - compared["published_persistence"].to_numpy())
    var_diff = np.abs(compared["variance"]["difference"].to_numpy())
    print(f"convergence {convergence}")
    print(f"elapsed_seconds {elapsed:.3f}")
    print(f"loglik {loglik:.16f}")
    print("published_target -111389.88")
    print(f"coefficient_abs_diff_max {np.max(abs_diff):.6g}")
    print(f"coefficient_abs_diff_median {np.median(abs_diff):.6g}")
    print(f"se_abs_diff_max {np.nanmax(se_diff):.6g}")
    print(f"se_abs_diff_median {np.nanmedian(se_diff):.6g}")
    print(f"persistence_abs_diff_max {np.max(rho_diff):.6g}")
    print(f"persistence_abs_diff_median {np.median(rho_diff):.6g}")
    print(f"variance_abs_diff_max {np.nanmax(var_diff):.6g}")
    print(f"variance_abs_diff_median {np.nanmedian(var_diff):.6g}")
    print("\nCoefficients (unconstrained), fit minus the replication X VECTOR:")
    show = coefficients.copy()
    for column in ("estimate", "published", "difference"):
        show[column] = show[column].map(lambda value: f"{value:.6f}")
    for column in ("std.error", "published_se"):
        show[column] = show[column].map(lambda value: f"{value:.4f}" if np.isfinite(value) else "nan")
    print(show.to_string(index=False))
    print("\nPersistence on (-1, 1), fit minus the transformed replication rho:")
    print((compared["persistence"] - compared["published_persistence"]).round(4).to_string())
    print("\nVariance decomposition, fit minus the replication log:")
    variance = compared["variance"].copy()
    for column in ("estimate", "published", "difference"):
        variance[column] = variance[column].map(lambda value: f"{value:.4f}")
    print(variance.to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate or refit the GSS attitude model.")
    parser.add_argument("--refit", action="store_true", help="Start BFGS from the shipped coefficients.")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--maxit", type=int, default=200)
    args = parser.parse_args()
    root = repository_root()
    data_path, log_path = gss_paths(root)
    panel = read_gss(data_path)
    counts = category_counts(panel)
    for name, expected in EXPECTED_COUNTS.items():
        if counts[name] != expected:
            raise SystemExit(
                f"Category counts for {name} do not match the replication log: {counts[name]}"
            )
    print(f"People in the file: {panel.attrs['n_people']}. People in the likelihood: {panel.attrs['n_included']}.")
    if panel.attrs["n_included"] != 11252:
        raise SystemExit("Expected 11252 people in the likelihood.")
    estimates = read_gss_estimates(log_path)
    if len(estimates) != 164:
        raise SystemExit(f"Expected 164 free coefficients, found {len(estimates)}.")
    started = time.perf_counter()
    result = attitudedyn.fit(
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
        start=estimates.to_numpy(),
        threads=args.threads,
        maxit=args.maxit if args.refit else 0,
        vcov="opg_centered" if args.refit else "none",
        use_gradient=bool(args.refit),
        keep_scores=bool(args.refit),
        seed=1,
    )
    elapsed = time.perf_counter() - started
    order = fortran_names(result.layout)
    if order != list(estimates.index):
        raise SystemExit(
            "Free-coefficient order does not match the replication log. "
            f"First layout name: {order[0]}. First log name: {estimates.index[0]}."
        )
    print(f"free_coefficients {len(result.params)}")
    print(f"persons {result.nobs}")
    print(f"loglik {result.loglik:.16f}")
    print(f"philox_reference {PHILOX_LOGLIK:.16f}")
    print(f"loglik_abs_diff {abs(result.loglik - PHILOX_LOGLIK):.3e}")
    print("published_target -111389.88")
    print(f"elapsed_seconds {elapsed:.3f}")
    persistence = result.persistence_matrix()
    variance = result.variance_decomposition()
    response = result.impulse_response(horizon=8)
    print("\nPersistence on (-1, 1):")
    print(persistence.round(4).to_string())
    print(f"\nSpectral radius: {response.spectral_radius:.6f}")
    print("\nImpulse response at horizon 1 (rows receive, columns are the source):")
    print(pd.DataFrame(response.response[1], index=response.attitudes, columns=response.attitudes).round(4).to_string())
    print("\nVariance decomposition:")
    print(variance.round(4).to_string(index=False))
    if variance.attrs.get("note"):
        print(variance.attrs["note"])
    out_dir = root / "validation" / "output"
    out_dir.mkdir(parents=True, exist_ok=True)
    persistence_out = persistence.copy()
    persistence_out.insert(0, "attitude", persistence_out.index)
    persistence_out.to_csv(out_dir / "gss_persistence.csv", index=False)
    variance.to_csv(out_dir / "gss_variance.csv", index=False)
    response.to_dataframe().to_csv(out_dir / "gss_impulse_response.csv", index=False)
    if args.refit:
        compared = comparison(result, log_path)
        print_comparison(compared, result.loglik, result.convergence, elapsed)
        compared["coefficients"].to_csv(out_dir / "gss_coefficient_comparison.csv", index=False)
        compared["variance"].to_csv(out_dir / "gss_variance_comparison.csv", index=False)
        print(f"Saved the comparison under {out_dir}")
    else:
        print(f"\nSaved persistence, variance, and impulse responses under {out_dir}")


if __name__ == "__main__":
    main()
