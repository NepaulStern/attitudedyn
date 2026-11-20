"""Parameter-recovery coverage for the three covariance estimators.

The data-generating process is ``simulate_attitudes``: two attitudes, four
waves, one factor, three categories, no wave effects. Estimation starts at
the data-generating coefficients and uses BFGS with the analytic gradient.
Each replication draws a new panel (seed ``1000 + replication``) and reuses
estimation seed 1.

Coverage is the fraction of free coefficients, among replications whose
optimiser returned code 0, whose 95 percent Wald interval contains the
data-generating value. The interval is estimate ± ``scipy.stats.norm.ppf(0.975)``
times the standard error. A coefficient whose covariance could not be
computed is counted as not covered.

Run from the repository root:

    python validation/simulation_recovery.py --threads 4
"""

from __future__ import annotations

import argparse
import json
import time
import warnings
from pathlib import Path

import numpy as np
from scipy.stats import norm

import attitudedyn

COVARIANCE_TYPES = ("opg_centered", "hessian", "sandwich")
WALD_Z = float(norm.ppf(0.975))


def one_replication(n_persons: int, replication: int, draws: int, threads: int, maxit: int) -> dict:
    """Fit one panel and record, for each covariance, which coefficients are covered."""
    panel = attitudedyn.simulate_attitudes(
        n_persons=n_persons,
        n_waves=4,
        n_attitudes=2,
        n_factors=1,
        n_categories=3,
        seed=1000 + replication,
    )
    truth_full = panel.attrs["theta"]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = attitudedyn.fit(
            panel,
            formula="cbind(y1, y2) ~ x",
            id_col="id",
            wave_col="wave",
            time_col="time",
            n_factors=1,
            draws=draws,
            threads=threads,
            vcov="none",
            keep_scores=True,
            seed=1,
            include_wave_effects=False,
            maxit=maxit,
            reltol=1e-8,
            start=truth_full,
        )
    names = list(result.params.index)
    truth = truth_full.loc[names].to_numpy(dtype=np.float64)
    estimate = result.params.to_numpy(dtype=np.float64)
    covered = {}
    failures = {}
    for covariance_type in COVARIANCE_TYPES:
        try:
            variance = result.vcov(type=covariance_type).to_numpy(dtype=np.float64)
            standard_error = np.sqrt(np.clip(np.diag(variance), 0.0, np.inf))
            inside = (
                np.isfinite(standard_error)
                & (estimate - WALD_Z * standard_error <= truth)
                & (truth <= estimate + WALD_Z * standard_error)
            )
            covered[covariance_type] = [bool(flag) for flag in inside]
            failures[covariance_type] = 0
        except (ValueError, np.linalg.LinAlgError):
            covered[covariance_type] = [False] * len(names)
            failures[covariance_type] = 1
    return {
        "convergence": int(result.convergence),
        "names": names,
        "truth": [float(value) for value in truth],
        "covered": covered,
        "covariance_failures": failures,
        "warning_count": len(caught),
    }


def summarise_sample_size(n_persons: int, replications: list[dict], elapsed_seconds: float) -> dict:
    """Collapse replication records into counts."""
    converged = [record for record in replications if record["convergence"] == 0]
    names = converged[0]["names"] if converged else []
    per_coefficient = {
        name: {covariance_type: 0 for covariance_type in COVARIANCE_TYPES} for name in names
    }
    covered_counts = {covariance_type: 0 for covariance_type in COVARIANCE_TYPES}
    failure_counts = {covariance_type: 0 for covariance_type in COVARIANCE_TYPES}
    for record in converged:
        for covariance_type in COVARIANCE_TYPES:
            flags = record["covered"][covariance_type]
            covered_counts[covariance_type] += int(sum(flags))
            failure_counts[covariance_type] += int(record["covariance_failures"][covariance_type])
            for name, flag in zip(record["names"], flags):
                per_coefficient[name][covariance_type] += int(flag)
    intervals = len(converged) * len(names)
    coverage = {
        covariance_type: (covered_counts[covariance_type] / intervals if intervals else None)
        for covariance_type in COVARIANCE_TYPES
    }
    return {
        "n_persons": n_persons,
        "replications_requested": len(replications),
        "replications_converged": len(converged),
        "free_coefficients": names,
        "truth": converged[0]["truth"] if converged else [],
        "intervals": intervals,
        "covered_counts": covered_counts,
        "coverage": coverage,
        "covariance_failures": failure_counts,
        "warning_count": int(sum(record["warning_count"] for record in replications)),
        "per_coefficient_covered": per_coefficient,
        "elapsed_seconds": elapsed_seconds,
    }


def run_study(sample_sizes: list[int], replications: int, draws: int, threads: int, maxit: int) -> dict:
    """Run every sample size and return the JSON-ready study."""
    started = time.perf_counter()
    runs = []
    for n_persons in sample_sizes:
        sample_started = time.perf_counter()
        records = [
            one_replication(n_persons, replication, draws, threads, maxit)
            for replication in range(replications)
        ]
        runs.append(summarise_sample_size(n_persons, records, time.perf_counter() - sample_started))
        print(
            f"n={n_persons} converged={runs[-1]['replications_converged']}/{replications} "
            f"elapsed={runs[-1]['elapsed_seconds']:.1f}s "
            + " ".join(f"{name}={runs[-1]['coverage'][name]:.4f}" for name in COVARIANCE_TYPES),
            flush=True,
        )
    return {
        "n_waves": 4,
        "n_attitudes": 2,
        "n_factors": 1,
        "n_categories": 3,
        "include_wave_effects": False,
        "draws": draws,
        "antithetic": True,
        "draw_type": "philox",
        "estimation_seed": 1,
        "data_seed_first": 1000,
        "optimizer": "BFGS",
        "reltol": 1e-8,
        "maxit": maxit,
        "threads": threads,
        "start": "data-generating coefficients",
        "wald_z": WALD_Z,
        "confidence_level": 0.95,
        "hessian_step": "1e-5 * max(1, abs(coefficient)), central difference of the analytic gradient, then symmetrised",
        "runs": runs,
        "elapsed_seconds": time.perf_counter() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Coverage of the three covariance estimators.")
    parser.add_argument("--replications", type=int, default=80)
    parser.add_argument("--draws", type=int, default=8)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--maxit", type=int, default=50)
    parser.add_argument("--n-persons", type=int, nargs="+", default=[200, 400, 1600])
    parser.add_argument("--output", type=Path, default=Path("validation/output/simulation_recovery.json"))
    args = parser.parse_args()
    study = run_study(args.n_persons, args.replications, args.draws, args.threads, args.maxit)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(study, indent=2) + "\n")
    print(f"wrote {args.output} in {study['elapsed_seconds']:.1f}s", flush=True)


if __name__ == "__main__":
    main()
