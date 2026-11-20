"""Fitted dynamic ordinal attitude model."""

from __future__ import annotations

import warnings
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.stats import norm

from attitudedyn import _core
from attitudedyn.coefficients import full_from_free


class ImpulseResponse:
    """Powers of the natural persistence matrix.

    ``response[h, to, from]`` is the effect on attitude ``to``, ``h`` steps of
    one time unit later, of a unit innovation in attitude ``from``. Horizon 0
    is the identity.

    Attributes
    ----------
    response :
        Array of shape ``(horizon + 1, n_attitudes, n_attitudes)``.
    spectral_radius :
        Approximate spectral radius of the persistence matrix.
    attitudes :
        Attitude names, in row and column order.
    """

    def __init__(self, response: np.ndarray, spectral_radius: float, attitudes: Sequence[str]) -> None:
        self.response = response
        self.spectral_radius = float(spectral_radius)
        self.attitudes = list(attitudes)

    def to_dataframe(self) -> pd.DataFrame:
        """Return one row per horizon, receiving attitude, and source attitude."""
        horizon, n_attitudes, _ = self.response.shape
        rows = []
        for step in range(horizon):
            for to in range(n_attitudes):
                for source in range(n_attitudes):
                    rows.append(
                        {
                            "horizon": step,
                            "to": self.attitudes[to],
                            "from": self.attitudes[source],
                            "response": self.response[step, to, source],
                        }
                    )
        return pd.DataFrame(rows)


class AttitudeDynResult:
    """Maximum simulated-likelihood estimates of a dynamic ordinal attitude model.

    Coefficients are on the unconstrained scale used by the optimiser.
    :meth:`persistence_matrix` reports persistence on (-1, 1).

    Attributes
    ----------
    params :
        Free coefficients, unconstrained.
    loglik :
        Sample log-likelihood at ``params``, the weighted sum of the person
        contributions minus the contemporaneous-correlation penalty.
    aic, bic, nobs :
        Information criteria and the number of people with positive weight.
    convergence :
        0 when the optimiser reports success, or when no iterations were asked
        for. A positive code is the optimiser's status.
    gradient :
        Analytic gradient of ``loglik`` with respect to ``params``, when it
        was requested.
    """

    def __init__(
        self,
        *,
        problem: _core.Problem,
        layout: pd.DataFrame,
        theta_full: np.ndarray,
        params: pd.Series,
        loglik: float,
        penalty: float,
        nobs: int,
        convergence: int,
        gradient: np.ndarray | None,
        scores: np.ndarray | None,
        prepared,
        n_factors: int,
        draws: int,
        threads: int,
        seed: int,
        antithetic: bool,
        draw_type: str,
        numerics: str,
        fortran_ghk: bool,
        contemporaneous: str,
        vcov_type: str,
        use_gradient: bool,
        include_wave_effects: bool,
        optimizer: str,
        maxit: int,
    ) -> None:
        self._problem = problem
        self.layout = layout
        self._theta_template = np.asarray(theta_full, dtype=np.float64).copy()
        self._theta_full = np.asarray(theta_full, dtype=np.float64).copy()
        self.params = params
        self.loglik = float(loglik)
        self.penalty = float(penalty)
        self.nobs = int(nobs)
        self.convergence = int(convergence)
        self.gradient = None if gradient is None else np.asarray(gradient, dtype=np.float64).copy()
        self._scores = None if scores is None else np.asarray(scores, dtype=np.float64).copy()
        self._prepared = prepared
        self.n_factors = int(n_factors)
        self.draws = int(draws)
        self.threads = int(threads)
        self.seed = int(seed)
        self.antithetic = bool(antithetic)
        self.draw_type = draw_type
        self.numerics = numerics
        self.fortran_ghk = bool(fortran_ghk)
        self.contemporaneous = contemporaneous
        self._vcov_type = vcov_type
        self.use_gradient = bool(use_gradient)
        self.include_wave_effects = bool(include_wave_effects)
        self.optimizer = optimizer
        self.maxit = int(maxit)
        self._vcov_cache: dict[str, pd.DataFrame] = {}
        k = len(self.params)
        self.aic = -2.0 * self.loglik + 2.0 * k
        self.bic = -2.0 * self.loglik + np.log(self.nobs) * k if self.nobs > 0 else np.nan
        self._sync_theta(self.params.to_numpy())

    def _sync_theta(self, theta_free: np.ndarray) -> None:
        self._theta_full = full_from_free(theta_free, self._theta_template, self.layout)

    def _numerics_index(self) -> int:
        return 1 if self.numerics == "shared" else 0

    def _evaluate(self, theta_free: np.ndarray, *, want_gradient: bool, want_scores: bool) -> dict:
        return self._problem.evaluate(
            np.asarray(theta_free, dtype=np.float64),
            True,
            self._numerics_index(),
            self.fortran_ghk,
            self.threads,
            want_gradient,
            want_scores,
        )

    def _person_scores(self, vcov_type: str) -> np.ndarray:
        if vcov_type == "opg_centered" or self._scores is None:
            evaluation = self._evaluate(self.params.to_numpy(), want_gradient=True, want_scores=True)
            scores = np.asarray(evaluation["person_score"], dtype=np.float64)
        else:
            scores = self._scores
        keep = self._prepared.data["weight"] > 0
        if scores.shape[0] != keep.shape[0]:
            raise RuntimeError("The score matrix does not have one row per person.")
        return scores[keep]

    @property
    def bse(self) -> pd.Series:
        """Standard errors of :attr:`params` from :meth:`vcov`."""
        try:
            variance = self.vcov()
            se = np.sqrt(np.clip(np.diag(variance.to_numpy()), 0.0, np.inf))
        except (ValueError, np.linalg.LinAlgError) as exc:
            warnings.warn(f"Standard errors were not computed: {exc}", UserWarning, stacklevel=2)
            se = np.full(len(self.params), np.nan)
        return pd.Series(se, index=self.params.index, name="bse")

    @property
    def pvalues(self) -> pd.Series:
        """Two-sided normal p-values for :attr:`params`."""
        se = self.bse.to_numpy()
        estimate = self.params.to_numpy()
        with np.errstate(divide="ignore", invalid="ignore"):
            statistic = estimate / se
            probability = 2.0 * norm.sf(np.abs(statistic))
        return pd.Series(probability, index=self.params.index, name="pvalue")

    def vcov(self, type: str | None = None) -> pd.DataFrame:
        """Covariance matrix of the unconstrained free coefficients.

        Parameters
        ----------
        type :
            ``"opg"`` (the default stored with the fit) is the inverse outer
            product of the person scores. ``"opg_centered"`` centres those
            scores first. ``"hessian"`` is the inverse negative numerical
            Hessian of the analytic gradient. ``"sandwich"`` is the robust
            covariance. ``"none"`` raises.

        Returns
        -------
        pandas.DataFrame
            Square covariance, indexed by coefficient name.
        """
        chosen = self._vcov_type if type is None else type
        if chosen not in {"opg", "opg_centered", "hessian", "sandwich", "none"}:
            raise ValueError("type must be 'opg', 'opg_centered', 'hessian', 'sandwich', or 'none'.")
        if chosen == "none":
            raise ValueError("This fit was asked not to compute a covariance. Call vcov with type=.")
        cached = self._vcov_cache.get(chosen)
        if cached is not None:
            return cached
        if chosen in {"hessian", "sandwich"} and (self.numerics != "stable" or self.contemporaneous != "identity"):
            raise ValueError(
                "The numerical Hessian needs analytic gradients, hence stable numerics and an identity correlation."
            )
        names = list(self.params.index)
        if chosen in {"opg", "opg_centered", "sandwich"}:
            scores = self._person_scores(chosen)
            if chosen == "opg_centered":
                scores = scores - scores.mean(axis=0, keepdims=True)
            meat = scores.T @ scores
        else:
            meat = None
        if chosen in {"hessian", "sandwich"}:
            hessian = self._numerical_hessian()
            try:
                bread = np.linalg.inv(-hessian)
            except np.linalg.LinAlgError as exc:
                raise ValueError("The Hessian could not be inverted.") from exc
        else:
            bread = None
        if chosen in {"opg", "opg_centered"}:
            estimate = _invert_information(meat)
        elif chosen == "hessian":
            estimate = bread
        else:
            estimate = bread @ meat @ bread
        frame = pd.DataFrame(estimate, index=names, columns=names)
        self._vcov_cache[chosen] = frame
        return frame

    def _numerical_hessian(self) -> np.ndarray:
        theta = self.params.to_numpy()
        p = theta.size
        hessian = np.zeros((p, p), dtype=np.float64)
        for j in range(p):
            step = 1e-5 * max(1.0, abs(theta[j]))
            up = theta.copy()
            down = theta.copy()
            up[j] += step
            down[j] -= step
            gradient_up = np.asarray(
                self._evaluate(up, want_gradient=True, want_scores=False)["gradient"], dtype=np.float64
            )
            gradient_down = np.asarray(
                self._evaluate(down, want_gradient=True, want_scores=False)["gradient"], dtype=np.float64
            )
            hessian[:, j] = (gradient_up - gradient_down) / (2.0 * step)
        return 0.5 * (hessian + hessian.T)

    def conf_int(self, alpha: float = 0.05, type: str | None = None) -> pd.DataFrame:
        """Wald confidence intervals for the unconstrained coefficients.

        Parameters
        ----------
        alpha :
            Tail probability. ``0.05`` is a 95 percent interval.
        type :
            Covariance passed to :meth:`vcov`.

        Returns
        -------
        pandas.DataFrame
            Two columns, the lower and upper normal limits.
        """
        if not 0.0 < alpha < 1.0:
            raise ValueError("alpha must lie strictly between 0 and 1.")
        estimate = self.params.to_numpy()
        se = np.sqrt(np.clip(np.diag(self.vcov(type=type).to_numpy()), 0.0, np.inf))
        z = float(norm.ppf(1.0 - alpha / 2.0))
        lower = estimate - z * se
        upper = estimate + z * se
        low_name = f"{100 * alpha / 2:g}%"
        high_name = f"{100 * (1.0 - alpha / 2.0):g}%"
        return pd.DataFrame({low_name: lower, high_name: upper}, index=self.params.index)

    def to_dataframe(self, alpha: float = 0.05, type: str | None = None) -> pd.DataFrame:
        """Coefficient table with standard errors, z statistics, p-values, and a Wald interval.

        Parameters
        ----------
        alpha :
            Tail probability of the interval.
        type :
            Covariance passed to :meth:`vcov`.

        Returns
        -------
        pandas.DataFrame
            One row per free coefficient.
        """
        estimate = self.params
        try:
            se = pd.Series(np.sqrt(np.clip(np.diag(self.vcov(type=type).to_numpy()), 0.0, np.inf)), index=estimate.index)
        except (ValueError, np.linalg.LinAlgError):
            se = pd.Series(np.nan, index=estimate.index)
        with np.errstate(divide="ignore", invalid="ignore"):
            statistic = estimate.to_numpy() / se.to_numpy()
            probability = 2.0 * norm.sf(np.abs(statistic))
        interval = None
        if np.any(np.isfinite(se.to_numpy())):
            try:
                interval = self.conf_int(alpha=alpha, type=type)
            except (ValueError, np.linalg.LinAlgError):
                interval = None
        frame = pd.DataFrame(
            {
                "estimate": estimate.to_numpy(),
                "std_error": se.to_numpy(),
                "statistic": statistic,
                "pvalue": probability,
            },
            index=estimate.index,
        )
        if interval is not None:
            frame["conf_low"] = interval.iloc[:, 0].to_numpy()
            frame["conf_high"] = interval.iloc[:, 1].to_numpy()
        return frame

    def summary(self, type: str | None = None) -> str:
        """Return a text summary of the fit, the coefficients, and the persistence matrix.

        Parameters
        ----------
        type :
            Covariance passed to :meth:`vcov`. The printed table uses whatever
            that call can compute; missing standard errors are left blank.

        Returns
        -------
        str
        """
        chosen = self._vcov_type if type is None else type
        header = [
            "Dynamic ordinal attitude model",
            f"Log-likelihood: {self.loglik:.6f}",
            (
                f"People: {self.nobs}   Attitudes: {self._prepared.n_attitudes}   "
                f"Waves: {self._prepared.n_waves}"
            ),
            f"Draws: {self.draws}" + (" (plus antithetic twins)" if self.antithetic else ""),
            f"Free coefficients: {len(self.params)}",
            "Converged." if self.convergence == 0 else f"Optimiser code: {self.convergence}",
            "",
            "Coefficients (unconstrained scale):",
        ]
        try:
            table = self.to_dataframe(type=type)
            body = table.to_string(float_format=lambda value: f"{value: .6g}")
        except (ValueError, np.linalg.LinAlgError) as exc:
            body = self.params.to_frame("estimate").to_string(float_format=lambda value: f"{value: .6g}")
            body += f"\n\nStandard errors were not computed: {exc}"
        persistence = self.persistence_matrix()
        lines = header + [
            body,
            "",
            "Signif. codes use the normal distribution. The 10% cutoff is 1.645.",
            "",
            "Persistence on (-1, 1), receiving attitude by source attitude:",
            persistence.round(3).to_string(),
            "",
            f"Covariance type: {chosen}",
        ]
        return "\n".join(lines)

    def persistence_matrix(self) -> pd.DataFrame:
        """Persistence on (-1, 1).

        Row ``to`` and column ``from`` is the coefficient that multiplies
        attitude ``from`` in the recursion for attitude ``to``, before it is
        raised to the gap between waves.

        Returns
        -------
        pandas.DataFrame
        """
        natural = self._problem.natural(self._theta_full)
        matrix = np.asarray(natural["persistence"], dtype=np.float64).reshape(
            self._prepared.n_attitudes, self._prepared.n_attitudes
        )
        names = self._prepared.attitude_names
        return pd.DataFrame(matrix, index=names, columns=names)

    def cutpoints(self) -> dict[str, np.ndarray]:
        """Ordered cutpoints, one vector per attitude.

        Each vector has ``K + 1`` entries. The ends and the cut between
        category 1 and category 2 are the normalised values -100, 0, and 100.
        """
        natural = self._problem.natural(self._theta_full)
        return {
            name: np.asarray(values, dtype=np.float64)
            for name, values in zip(self._prepared.attitude_names, natural["cutpoints"])
        }

    def variance_decomposition(self) -> pd.DataFrame:
        """Split the latent scale into covariate, wave, factor, and shock contributions.

        The shock column is the long-run variance of the recursion
        ``state = R state_lag + innovation`` with an identity innovation
        covariance. It is not finite when ``R`` is not stationary.

        Returns
        -------
        pandas.DataFrame
            One row per attitude.
        """
        component = self._problem.variance(self._theta_full)
        frame = pd.DataFrame(
            {
                "attitude": self._prepared.attitude_names,
                "covariate": np.asarray(component["covariate"], dtype=np.float64),
                "wave_effect": np.asarray(component["wave_effect"], dtype=np.float64),
                "factor": np.asarray(component["factor"], dtype=np.float64),
                "shock": np.asarray(component["shock"], dtype=np.float64),
                "total": np.asarray(component["total"], dtype=np.float64),
            }
        )
        frame.attrs["note"] = component["note"]
        return frame

    def impulse_response(self, horizon: int = 8) -> ImpulseResponse:
        """Impulse responses of the persistence matrix.

        Parameters
        ----------
        horizon :
            Number of steps after the impact. The returned array includes the
            impact at horizon 0, so its first dimension is ``horizon + 1``.

        Returns
        -------
        ImpulseResponse
        """
        if horizon < 0:
            raise ValueError("horizon must be non-negative.")
        matrix = self.persistence_matrix()
        raw = _core.impulse_response(np.ascontiguousarray(matrix.to_numpy(dtype=np.float64)), matrix.shape[0], int(horizon))
        n_attitudes = matrix.shape[0]
        response = np.asarray(raw["response"], dtype=np.float64).reshape(horizon + 1, n_attitudes, n_attitudes)
        return ImpulseResponse(response, raw["spectral_radius"], self._prepared.attitude_names)

    def influence_shares(self) -> pd.DataFrame:
        """Long-run share of each attitude's variance coming from each innovation.

        Returns
        -------
        pandas.DataFrame
            Empty when the persistence matrix is not stationary. The
            ``variance`` attribute holds the long-run variances when they exist.
        """
        matrix = self.persistence_matrix()
        raw = _core.influence_shares(np.ascontiguousarray(matrix.to_numpy(dtype=np.float64)), matrix.shape[0])
        names = self._prepared.attitude_names
        if not raw["stationary"]:
            warnings.warn(
                "The persistence matrix is not stationary, so long-run shares are not defined.",
                UserWarning,
                stacklevel=2,
            )
            frame = pd.DataFrame(index=names, columns=names, dtype=np.float64)
            frame.attrs["stationary"] = False
            return frame
        share = np.asarray(raw["share"], dtype=np.float64).reshape(len(names), len(names))
        frame = pd.DataFrame(share, index=names, columns=names)
        frame.attrs["variance"] = pd.Series(np.asarray(raw["long_run_variance"], dtype=np.float64), index=names)
        frame.attrs["stationary"] = True
        return frame

    def predict(
        self,
        type: str = "prob",
        n_draws: int | None = None,
        seed: int | None = None,
    ) -> pd.DataFrame:
        """Category probabilities, most likely categories, or mean latent attitudes.

        Paths are structural draws, averaged with Rao-Blackwellised interval
        probabilities. Rows are the person-waves inside each person's window
        with positive weight.

        Parameters
        ----------
        type :
            ``"prob"`` returns a long frame with one row per person, wave,
            attitude, and category. ``"class"`` returns the category with the
            highest probability. ``"latent"`` returns the average latent attitude.
        n_draws :
            Structural draws per person. The default is ``max(20, self.draws)``.
        seed :
            Generator seed. The default is ``self.seed + 1`` and does not
            affect the likelihood draws.

        Returns
        -------
        pandas.DataFrame
        """
        if type not in {"prob", "class", "latent"}:
            raise ValueError("type must be 'prob', 'class', or 'latent'.")
        if n_draws is None:
            n_draws = max(20, self.draws)
        if n_draws < 1:
            raise ValueError("n_draws must be a positive integer.")
        if seed is None:
            seed = self.seed + 1
        prediction = self._problem.predict(self._theta_full, int(n_draws), int(seed), self.threads)
        prepared = self._prepared
        n_persons = prepared.n_persons
        n_waves = prepared.n_waves
        n_attitudes = prepared.n_attitudes
        max_categories = int(prediction["max_categories"])
        probability = np.asarray(prediction["probability"], dtype=np.float64).reshape(
            n_persons, n_waves, n_attitudes, max_categories
        )
        latent = np.asarray(prediction["mean_latent"], dtype=np.float64).reshape(n_persons, n_waves, n_attitudes)
        first = prepared.data["first_wave"]
        last = prepared.data["last_wave"]
        weight = prepared.data["weight"]
        rows_id: list = []
        rows_wave: list = []
        rows_person: list[int] = []
        rows_wave_index: list[int] = []
        for person in range(n_persons):
            if weight[person] <= 0 or last[person] < first[person]:
                continue
            for wave in range(int(first[person]), int(last[person]) + 1):
                rows_id.append(prepared.person_ids[person])
                rows_wave.append(prepared.wave_values[wave])
                rows_person.append(person)
                rows_wave_index.append(wave)
        person_index = np.asarray(rows_person, dtype=np.int32)
        wave_index = np.asarray(rows_wave_index, dtype=np.int32)
        base = {"id": rows_id, "wave": rows_wave}
        names = prepared.attitude_names
        if type == "latent":
            columns = {names[j]: latent[person_index, wave_index, j] for j in range(n_attitudes)}
            return pd.DataFrame({**base, **columns})
        if type == "class":
            classes = {}
            for j, name in enumerate(names):
                k = int(prepared.n_categories[j])
                probs = probability[person_index, wave_index, j, :k]
                classes[name] = np.argmax(probs, axis=1) + 1
            return pd.DataFrame({**base, **classes})
        records = []
        for row in range(len(rows_id)):
            for j, name in enumerate(names):
                k = int(prepared.n_categories[j])
                for category in range(k):
                    records.append(
                        {
                            "id": rows_id[row],
                            "wave": rows_wave[row],
                            "attitude": name,
                            "category": category + 1,
                            "probability": probability[person_index[row], wave_index[row], j, category],
                        }
                    )
        return pd.DataFrame(records)

    def simulate(self, nsim: int = 1, seed: int | None = None) -> pd.DataFrame | list[pd.DataFrame]:
        """Draw panels from the fitted structural model.

        Covariates are standard normal, except intercept columns, which are 1.
        This is a draw from the model, not a resample of the estimation sample.

        Parameters
        ----------
        nsim :
            Number of panels.
        seed :
            Seed of the first panel. Later panels use ``seed + 1``, ``seed + 2``,
            and so on. The default is ``self.seed + 17``.

        Returns
        -------
        pandas.DataFrame or list of pandas.DataFrame
            A single frame when ``nsim`` is 1, otherwise a list. The intercept
            column is omitted so the frame can be passed back to :func:`attitudedyn.fit`
            with ``add_intercept=True``.
        """
        if nsim < 1:
            raise ValueError("nsim must be a positive integer.")
        if seed is None:
            seed = self.seed + 17
        frames = []
        for draw in range(nsim):
            panel = self._problem.simulate(self._theta_full, int(self.nobs), int(seed) + draw)
            frames.append(_panel_frame(panel, self._prepared.attitude_names, self._prepared.covariate_names))
        if nsim == 1:
            return frames[0]
        return frames

    def __repr__(self) -> str:
        state = "converged" if self.convergence == 0 else f"optimiser code {self.convergence}"
        return (
            f"AttitudeDynResult(loglik={self.loglik:.6f}, nobs={self.nobs}, "
            f"n_free={len(self.params)}, {state})"
        )


def _invert_information(information: np.ndarray) -> np.ndarray:
    try:
        return np.linalg.inv(information)
    except np.linalg.LinAlgError:
        diagonal = np.diag(information)
        ridge = 1e-8 * float(np.mean(diagonal)) if diagonal.size else 1e-8
        if not np.isfinite(ridge) or ridge == 0.0:
            ridge = 1e-6
        warnings.warn(
            "The outer product of the scores is singular. A ridge was added.",
            UserWarning,
            stacklevel=3,
        )
        return np.linalg.inv(information + np.eye(information.shape[0]) * ridge)


def _panel_frame(panel: Mapping, attitude_names: Sequence[str], covariate_names: Sequence[str]) -> pd.DataFrame:
    n_rows = int(len(panel["id"]))
    n_attitudes = int(panel["n_attitudes"])
    n_covariates = int(panel["n_covariates"])
    outcomes = np.asarray(panel["outcome"], dtype=np.int32).reshape(n_rows, n_attitudes)
    covariates = np.asarray(panel["covariates"], dtype=np.float64).reshape(n_rows, n_covariates)
    keep = [name != "(Intercept)" and name != "Intercept" for name in covariate_names]
    kept_names = [name for name, flag in zip(covariate_names, keep) if flag]
    frame = pd.DataFrame(
        {
            "id": np.asarray(panel["id"]),
            "wave": np.asarray(panel["wave"]),
            "time": np.asarray(panel["time"], dtype=np.float64),
        }
    )
    for column, name in enumerate(kept_names):
        source = int(np.flatnonzero(keep)[column])
        frame[name] = covariates[:, source]
    for column, name in enumerate(attitude_names):
        frame[name] = outcomes[:, column]
    return frame
