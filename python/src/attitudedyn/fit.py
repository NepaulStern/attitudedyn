"""Simulated maximum likelihood for The Nepaul-Stern Model."""

from __future__ import annotations

import warnings
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from attitudedyn import _core
from attitudedyn.bfgs import minimize_bfgs
from attitudedyn.coefficients import apply_fix, apply_start, default_theta, free_from_full
from attitudedyn.prepare import layout_frame, parse_formula, prepare_panel
from attitudedyn.result import AttitudeDynResult

_DRAW_KIND = {"philox": 0, "halton": 1, "sobol": 2}
_OPTIMIZERS = {"bfgs": "BFGS", "l-bfgs-b": "L-BFGS-B", "lbfgsb": "L-BFGS-B", "nelder-mead": "Nelder-Mead"}


def fit(
    data: pd.DataFrame,
    outcomes: Sequence[str] | None = None,
    covariates: Sequence[str] | None = None,
    id_col: str = "id",
    wave_col: str = "wave",
    time_grid: Sequence[float] | None = None,
    time_col: str | None = None,
    formula: str | None = None,
    weights: str | np.ndarray | Sequence[float] | None = None,
    subset: Sequence[bool] | np.ndarray | None = None,
    n_factors: int = 2,
    draws: int = 10,
    include_wave_effects: bool = True,
    contemporaneous: str = "identity",
    start: Mapping[str, float] | Sequence[float] | np.ndarray | pd.Series | None = None,
    fix: Mapping[str, float] | None = None,
    threads: int = 1,
    optimizer: str = "bfgs",
    maxit: int = 200,
    reltol: float = 1e-8,
    use_gradient: bool = True,
    vcov: str = "opg_centered",
    numerics: str = "stable",
    fortran_ghk: bool = False,
    fix_last_loading: bool = True,
    fix_intercept: bool = True,
    antithetic: bool = True,
    draw_type: str = "philox",
    seed: int = 1,
    keep_scores: bool = True,
    add_intercept: bool | None = None,
    na_action: str = "pass",
    trace: bool = False,
) -> AttitudeDynResult:
    """Fit The Nepaul-Stern Model of Attitude Stability and Constraint.

    Several Likert items, measured on the same people across survey waves,
    are linked by a matrix of persistence coefficients and by person-level
    factors. The latent attitude carries forward from one wave to the next.
    The recorded answer is the Likert category that latent attitude falls
    into. The likelihood is a Geweke-Hajivassiliou-Keane (GHK) simulated
    likelihood.

    The default search follows the same steps as
    ``stats::optim(method = "BFGS")`` in the R package, and it uses the
    slope of the likelihood. ``"l-bfgs-b"`` and ``"nelder-mead"`` use
    ``scipy.optimize.minimize``. A trial step that is not finite, or whose
    threshold increment would overflow, is rejected: the objective reported
    to the optimiser is a large finite penalty and the gradient is zero, so
    the line search continues.

    Parameters
    ----------
    data :
        Long-format frame, one row per person-wave.
    outcomes, covariates :
        Likert items, coded as integers ``1, ..., K``, and covariate
        columns. Omit both when ``formula`` is set.
    id_col, wave_col :
        Person identifier and wave column.
    time_grid :
        Numeric wave dates, in order of the sorted unique values of
        ``wave_col``. Gaps in the persistence exponent are differences of
        these dates. The GSS application uses ``[2, 3, 4, 5, 6, 7, 8, 10]``.
    time_col :
        Optional time column. The mean within each wave is the common wave
        date, and it must not vary within a wave. Ignored when ``time_grid``
        is set.
    formula :
        Optional ``cbind(y1, y2) ~ x1 + x2``. A leading ``0`` or ``- 1``
        drops the intercept. Only bare column names are accepted.
    weights :
        Column name or array of non-negative weights, constant within each
        person. The sample log-likelihood is the weighted sum of the person
        contributions.
    subset :
        Boolean mask selecting rows of ``data`` before the model is built.
    n_factors :
        Number of standard normal person factors. Use 0 for none.
    draws :
        GHK draws per person. With antithetic sampling the average uses
        ``2 * draws`` paths.
    include_wave_effects :
        Include a wave effect for every attitude.
    contemporaneous :
        ``"identity"`` keeps the innovations uncorrelated within a wave,
        which is the published specification. ``"free"`` estimates a
        correlation factor. That specification has no analytic gradient.
    start :
        Starting values for the free coefficients, named or in
        free-coefficient order, on the unconstrained scale.
    fix :
        Named unconstrained coefficients to hold fixed.
    threads :
        Number of threads used to loop over people. Results do not depend
        on this value. The default is 1.
    optimizer :
        ``"bfgs"`` or ``"l-bfgs-b"``. ``"nelder-mead"`` is derivative-free.
        ``"bfgs"`` is the iteration ``stats::optim`` uses for BFGS.
    maxit :
        Maximum optimiser iterations. ``0`` evaluates ``start`` and returns.
    reltol :
        Relative tolerance on the objective. For BFGS this is ``optim``'s
        ``reltol``: a step that does not reduce the objective by
        ``reltol * (|f| + reltol)`` is not an improvement. L-BFGS-B passes
        the same number to SciPy as ``ftol`` and ``gtol``.
    use_gradient :
        Pass the analytic gradient to BFGS and L-BFGS-B. Turn this off only
        to compare with a numeric gradient.
    vcov :
        Covariance computed with the fit and returned by
        :meth:`AttitudeDynResult.vcov`. ``"opg_centered"``, the default,
        centres the person scores and inverts their outer product.
        ``"opg"`` is that inverse without centring. ``"hessian"`` inverts
        the negative numerical Hessian of the analytic gradient.
        ``"sandwich"`` is that inverse, times the uncentred outer product
        of the scores, times the inverse. ``"none"`` skips it.
    numerics :
        ``"stable"``, the default, uses a double-precision log of a normal
        interval probability and supplies analytic gradients. ``"shared"``
        reproduces the original program's normal-probability routine for
        comparisons.
    fortran_ghk :
        Use the original within-wave conditioning when the contemporaneous
        correlation is not the identity. With an identity correlation the
        two conditionings agree.
    fix_last_loading :
        Fix the loading of the last attitude on the last factor at zero when
        there are two or more factors.
    fix_intercept :
        Fix each attitude's intercept at zero when wave effects are included.
    antithetic :
        Pair each draw with its sign-flipped twin.
    draw_type :
        ``"philox"``, ``"halton"``, or ``"sobol"``. Streams are keyed by
        person, so people do not share one sequence and the result does not
        depend on ``threads``.
    seed :
        Non-negative seed. It does not change NumPy's random state.
    keep_scores :
        Retain the person-score matrix when a covariance is requested.
    add_intercept :
        Add a column of ones named ``(Intercept)`` when the covariate list
        does not already contain an intercept. ``None`` means true for a
        column list, and means whatever the formula requested otherwise.
    na_action :
        ``"pass"`` rejects missing covariates. ``"omit"`` drops those
        person-waves. Missing outcomes are kept.
    trace :
        Print the log-likelihood at the starting values.

    Returns
    -------
    AttitudeDynResult
        Estimates, likelihood, and the post-estimation methods.

    Examples
    --------
    >>> import attitudedyn
    >>> panel = attitudedyn.simulate_attitudes(n_persons=30, n_waves=3, n_factors=1, seed=1)
    >>> result = attitudedyn.fit(
    ...     panel, formula="cbind(y1, y2) ~ x", id_col="id", wave_col="wave",
    ...     time_col="time", n_factors=1, draws=4, maxit=0, threads=1, vcov="none",
    ... )
    >>> result.params.head()
    """
    _check_controls(
        n_factors=n_factors,
        draws=draws,
        threads=threads,
        maxit=maxit,
        reltol=reltol,
        seed=seed,
        contemporaneous=contemporaneous,
        numerics=numerics,
        draw_type=draw_type,
        vcov=vcov,
        optimizer=optimizer,
    )
    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas DataFrame.")
    if subset is not None:
        mask = np.asarray(subset, dtype=bool)
        if mask.shape != (len(data),):
            raise ValueError("subset must be a boolean vector with one entry per row.")
        data = data.loc[mask].copy()
    if formula is not None and (outcomes is not None or covariates is not None):
        raise ValueError("Pass either formula or outcomes and covariates, not both.")
    if formula is not None:
        outcomes, covariates, formula_intercept = parse_formula(formula)
        if add_intercept is None:
            add_intercept = formula_intercept
    elif outcomes is None or covariates is None:
        raise ValueError("Pass outcomes and covariates, or a formula.")
    elif add_intercept is None:
        add_intercept = True
    weight_vector = _weight_vector(data, weights)

    method = _OPTIMIZERS[optimizer.lower().replace("_", "-")]
    analytic = bool(use_gradient) and numerics == "stable" and contemporaneous == "identity"
    if use_gradient and not analytic:
        warnings.warn(
            "Analytic gradients need stable numerics and an identity contemporaneous correlation. "
            "Using numeric derivatives.",
            UserWarning,
            stacklevel=2,
        )
    if method == "Nelder-Mead":
        analytic = False

    prepared = prepare_panel(
        data,
        outcomes=list(outcomes),
        covariates=list(covariates),
        id_col=id_col,
        wave_col=wave_col,
        time_col=time_col,
        time_grid=time_grid,
        weights=weight_vector,
        add_intercept=bool(add_intercept),
        include_wave_effects=include_wave_effects,
        n_factors=int(n_factors),
        free_contemporaneous=contemporaneous == "free",
        fix_last_loading=fix_last_loading,
        fix_intercept=fix_intercept,
        na_action=na_action,
    )
    layout = layout_frame(_core.layout(prepared.spec))
    theta_full = default_theta(layout)
    layout, theta_full = apply_fix(layout, theta_full, fix)
    theta_full = apply_start(layout, theta_full, start)
    free_names = layout.loc[layout["is_free"].to_numpy() == 1, "name"].tolist()
    if not free_names:
        raise ValueError("Every coefficient is fixed, so there is nothing to estimate.")
    empty_fixed = (layout["block"].to_numpy() == "wave_effect") & (layout["is_free"].to_numpy() == 0)
    if empty_fixed.any():
        fixed_names = layout.loc[empty_fixed, "name"].tolist()
        warnings.warn(
            f"Wave effects with no observed responses are fixed at zero ({len(fixed_names)}): {fixed_names}.",
            UserWarning,
            stacklevel=2,
        )
    draws_dict = {
        "n_draws": int(draws),
        "antithetic": bool(antithetic),
        "kind": _DRAW_KIND[draw_type],
        "seed": int(seed),
    }
    problem = _core.build(
        prepared.spec,
        prepared.data,
        np.asarray(theta_full, dtype=np.float64),
        layout["is_free"].to_numpy(),
        draws_dict,
    )
    free_start = free_from_full(theta_full, layout)

    cache: dict[tuple, dict] = {}

    def evaluate(theta: np.ndarray, *, want_gradient: bool, want_scores: bool) -> dict:
        point = np.ascontiguousarray(theta, dtype=np.float64)
        key = (point.tobytes(), want_gradient, want_scores)
        found = cache.get(key)
        if found is not None:
            return found
        if want_gradient and not want_scores:
            richer = cache.get((point.tobytes(), True, True))
            if richer is not None:
                return richer
        evaluation = problem.evaluate(
            point,
            True,
            1 if numerics == "shared" else 0,
            bool(fortran_ghk),
            int(threads),
            bool(want_gradient),
            bool(want_scores),
        )
        cache[key] = evaluation
        return evaluation

    best = {"loglik": -np.inf, "theta": free_start.copy()}

    def remember(theta: np.ndarray, value: float) -> None:
        if np.isfinite(value) and value > best["loglik"]:
            best["loglik"] = float(value)
            best["theta"] = np.array(theta, dtype=np.float64, copy=True)

    def objective(theta: np.ndarray) -> float:
        evaluation = evaluate(theta, want_gradient=analytic, want_scores=False)
        value = float(evaluation["objective"])
        remember(np.asarray(theta, dtype=np.float64), value)
        if not np.isfinite(value):
            return 1e12
        return -value

    def jacobian(theta: np.ndarray) -> np.ndarray:
        evaluation = evaluate(theta, want_gradient=True, want_scores=False)
        gradient = np.asarray(evaluation["gradient"], dtype=np.float64)
        if (
            gradient.shape != np.shape(theta)
            or not np.all(np.isfinite(gradient))
            or not np.isfinite(evaluation["objective"])
        ):
            return np.zeros(np.shape(theta), dtype=np.float64)
        return -gradient

    initial = evaluate(free_start, want_gradient=analytic, want_scores=False)
    if not np.isfinite(initial["objective"]):
        raise ValueError("The log-likelihood is not finite at the starting values.")
    remember(free_start, float(initial["objective"]))
    if trace:
        print(f"initial loglik {float(initial['objective']):.6f}")

    convergence = 0
    coefficients = free_start.copy()
    if maxit > 0:
        # Analytic BFGS uses the same Nash iteration as stats::optim, so a
        # refit from the same coefficients stops at the same point as R.
        # The other methods stay with SciPy.
        use_nash_bfgs = method == "BFGS" and analytic

        def report_iteration(iteration: int, neg_loglik: float) -> None:
            if trace:
                print(f"iter {iteration:4d} loglik {-neg_loglik:.6f}")

        options: dict = {"maxiter": int(maxit)}
        if method == "L-BFGS-B":
            options["ftol"] = float(reltol)
            options["gtol"] = float(reltol)
        elif method != "BFGS":
            options["xatol"] = float(reltol)
            options["fatol"] = float(reltol)
        try:
            if use_nash_bfgs:
                fitted = minimize_bfgs(
                    objective,
                    jacobian,
                    free_start,
                    maxit=int(maxit),
                    reltol=float(reltol),
                    report=report_iteration,
                )
            else:
                if method == "BFGS":
                    options["gtol"] = float(reltol)
                fitted = minimize(
                    objective,
                    free_start,
                    method=method,
                    jac=jacobian if analytic and method != "Nelder-Mead" else None,
                    options=options,
                )
        except (FloatingPointError, np.linalg.LinAlgError) as exc:
            warnings.warn(
                f"The optimiser stopped after a failed trial step ({exc}). The best finite point is reported.",
                UserWarning,
                stacklevel=2,
            )
            fitted = None
        if fitted is None:
            coefficients = best["theta"].copy()
            convergence = 1
        else:
            coefficients = np.asarray(fitted.x, dtype=np.float64)
            reported = -float(fitted.fun) if np.isfinite(fitted.fun) else np.nan
            if (
                not np.all(np.isfinite(coefficients))
                or not np.isfinite(reported)
                or (np.isfinite(best["loglik"]) and best["loglik"] > reported + 1e-8)
            ):
                coefficients = best["theta"].copy()
            convergence = 0 if fitted.success else int(fitted.status)

    need_scores = bool(keep_scores) or vcov in {"opg", "opg_centered", "sandwich"}
    final = evaluate(coefficients, want_gradient=analytic or need_scores, want_scores=need_scores)
    if not np.isfinite(final["objective"]) and np.isfinite(best["loglik"]):
        coefficients = best["theta"].copy()
        final = evaluate(coefficients, want_gradient=analytic or need_scores, want_scores=need_scores)
    gradient = np.asarray(final["gradient"], dtype=np.float64) if analytic or need_scores else None
    scores = np.asarray(final["person_score"], dtype=np.float64) if need_scores and "person_score" in final else None
    params = pd.Series(np.asarray(coefficients, dtype=np.float64), index=free_names, name="params")
    result = AttitudeDynResult(
        problem=problem,
        layout=layout,
        theta_full=theta_full,
        params=params,
        loglik=float(final["objective"]),
        penalty=float(final["penalty"]),
        nobs=int(final["n_persons_used"]),
        convergence=convergence,
        gradient=gradient,
        scores=scores,
        prepared=prepared,
        n_factors=int(n_factors),
        draws=int(draws),
        threads=int(threads),
        seed=int(seed),
        antithetic=bool(antithetic),
        draw_type=draw_type,
        numerics=numerics,
        fortran_ghk=bool(fortran_ghk),
        contemporaneous=contemporaneous,
        vcov_type=vcov,
        use_gradient=analytic,
        include_wave_effects=bool(include_wave_effects),
        optimizer=method,
        maxit=int(maxit),
    )
    _warn_after_fit(result, analytic)
    if vcov != "none":
        try:
            result.vcov(type=vcov)
        except (ValueError, np.linalg.LinAlgError) as exc:
            warnings.warn(f"The covariance matrix was not computed: {exc}", UserWarning, stacklevel=2)
    return result


def _weight_vector(data: pd.DataFrame, weights: str | np.ndarray | Sequence[float] | None) -> np.ndarray | None:
    if weights is None:
        return None
    if isinstance(weights, str):
        if weights not in data.columns:
            raise ValueError(f"Column {weights!r} (weights) is not in data.")
        return pd.to_numeric(data[weights], errors="coerce").to_numpy(dtype=np.float64)
    values = np.asarray(weights, dtype=np.float64)
    if values.shape != (len(data),):
        raise ValueError("weights must have one value per row of data (after subset).")
    return values


def _check_controls(**flags: object) -> None:
    if int(flags["n_factors"]) < 0:  # type: ignore[arg-type]
        raise ValueError("n_factors must be a nonnegative integer.")
    if int(flags["draws"]) < 1:  # type: ignore[arg-type]
        raise ValueError("draws must be a positive integer.")
    if int(flags["threads"]) < 1:  # type: ignore[arg-type]
        raise ValueError("threads must be a positive integer.")
    if int(flags["maxit"]) < 0:  # type: ignore[arg-type]
        raise ValueError("maxit must be a non-negative integer.")
    if not np.isfinite(flags["reltol"]) or float(flags["reltol"]) <= 0:  # type: ignore[arg-type]
        raise ValueError("reltol must be a positive number.")
    if not np.isfinite(flags["seed"]) or float(flags["seed"]) < 0:  # type: ignore[arg-type]
        raise ValueError("seed must be a finite, non-negative number.")
    if flags["contemporaneous"] not in {"identity", "free"}:
        raise ValueError("contemporaneous must be 'identity' or 'free'.")
    if flags["numerics"] not in {"stable", "shared"}:
        raise ValueError("numerics must be 'stable' or 'shared'.")
    if flags["draw_type"] not in _DRAW_KIND:
        raise ValueError("draw_type must be 'philox', 'halton', or 'sobol'.")
    if flags["vcov"] not in {"opg", "opg_centered", "hessian", "sandwich", "none"}:
        raise ValueError("vcov must be 'opg', 'opg_centered', 'hessian', 'sandwich', or 'none'.")
    key = str(flags["optimizer"]).lower().replace("_", "-")
    if key not in _OPTIMIZERS:
        raise ValueError("optimizer must be 'bfgs', 'l-bfgs-b', or 'nelder-mead'.")


def _warn_after_fit(result: AttitudeDynResult, analytic: bool) -> None:
    asked = result.maxit > 0
    if asked and result.convergence != 0:
        warnings.warn(
            f"The optimiser did not report convergence (code {result.convergence}).",
            UserWarning,
            stacklevel=3,
        )
    radius = result.impulse_response(horizon=0).spectral_radius
    if np.isfinite(radius) and radius >= 0.999:
        warnings.warn(
            f"The persistence matrix has spectral radius {radius:.3f}. "
            "A radius of 1 or more means the latent recursion is not stationary.",
            UserWarning,
            stacklevel=3,
        )
    if asked and result.convergence == 0 and analytic and result._scores is not None:
        keep = result._prepared.data["weight"] > 0
        scores = result._scores[keep]
        total = scores.sum(axis=0)
        scale = np.sqrt(np.sum(scores * scores, axis=0))
        scale = np.where(scale == 0.0, 1.0, scale)
        worst = float(np.max(np.abs(total) / scale))
        if np.isfinite(worst) and worst > 0.05:
            warnings.warn(
                "The largest score, scaled by its outer-product standard deviation, "
                f"is {worst:.3f}. The estimates may not be at a maximum.",
                UserWarning,
                stacklevel=3,
            )
