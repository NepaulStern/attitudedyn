"""BFGS iteration used by the R package.

``stats::optim(method = "BFGS")`` minimises with the variable-metric
method in Nash, *Compact Numerical Methods for Computers*. A trial step
is shortened by a factor of 0.2 until it passes the acceptance test. The
inverse-Hessian factor starts as the identity and is rebuilt after a run
of steps that do not improve the objective. The search stops when that
rebuild also fails to reduce the objective by ``reltol * (|f| + reltol)``,
which is the ``reltol`` rule documented for ``optim``. The arithmetic
below follows that iteration, so a fit with the same start, analytic
gradient, and ``reltol`` takes the same steps as ``attitudedyn()`` in R.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

# Step reduction and acceptance tolerance from the Nash iteration.
_STEP_REDUCTION = 0.2
_ACCEPTANCE = 0.0001
# A step that does not change ``10 + x`` in double precision is a zero step.
_NO_CHANGE_SHIFT = 10.0


@dataclass
class BFGSResult:
    """Point returned by :func:`minimize_bfgs`."""

    x: np.ndarray
    fun: float
    success: bool
    status: int
    nit: int
    nfev: int
    njev: int


def minimize_bfgs(
    objective,
    gradient,
    x0: np.ndarray,
    *,
    maxit: int,
    reltol: float,
    abstol: float = -math.inf,
    trace: bool = False,
    report=None,
) -> BFGSResult:
    """Minimise ``objective`` from ``x0`` with the Nash BFGS iteration.

    Parameters
    ----------
    objective, gradient :
        Function and analytic gradient. Both receive a one-dimensional
        ``float64`` vector. A non-finite objective value is a rejected trial.
    x0 :
        Starting point.
    maxit :
        Maximum number of iterations, in the same units ``optim`` uses.
        ``0`` evaluates ``x0`` and returns.
    reltol :
        Relative tolerance on the objective. The search treats a reduction
        smaller than ``reltol * (|f| + reltol)`` as no improvement.
    abstol :
        Absolute tolerance on the objective. The default, ``-inf``, does not
        stop the search, which is ``optim``'s default for BFGS.
    trace :
        Print the objective at each iteration.
    report :
        Optional ``report(iteration, objective)`` called once per iteration
        with the best objective so far.

    Returns
    -------
    BFGSResult
        ``status`` is 0 when the relative tolerance is met and 1 when
        ``maxit`` is reached. ``fun`` is the objective at the returned point.
    """
    start = np.asarray(x0, dtype=np.float64)
    if start.ndim != 1:
        raise ValueError("BFGS minimises a one-dimensional parameter vector.")
    dimension = int(start.size)
    point = [float(value) for value in start]
    if maxit <= 0:
        value = float(objective(np.asarray(point, dtype=np.float64)))
        return BFGSResult(np.asarray(point, dtype=np.float64), value, True, 0, 0, 1, 0)

    def function_at(current: list[float]) -> float:
        return float(objective(np.asarray(current, dtype=np.float64)))

    def gradient_at(current: list[float]) -> list[float]:
        slope = np.asarray(gradient(np.asarray(current, dtype=np.float64)), dtype=np.float64)
        if slope.shape != (dimension,):
            raise ValueError("The BFGS gradient must have one entry per parameter.")
        return [float(value) for value in slope]

    # Lower triangle of the inverse-Hessian factor. Row i stores columns 0..i.
    inverse_hessian = [[0.0 for _ in range(index + 1)] for index in range(dimension)]
    value = function_at(point)
    if not math.isfinite(value):
        raise ValueError("The BFGS objective is not finite at the starting values.")
    if trace:
        print(f"initial  value {value:.6f}")
    best_value = value
    slope = gradient_at(point)
    function_count = 1
    gradient_count = 1
    iteration = 1
    last_reset = gradient_count
    step_count = 0

    while True:
        # A reset puts the inverse-Hessian factor back to the identity.
        if last_reset == gradient_count:
            for row in range(dimension):
                for column in range(row):
                    inverse_hessian[row][column] = 0.0
                inverse_hessian[row][row] = 1.0
        location = point.copy()
        slope_then = slope.copy()
        directional_derivative = 0.0
        direction = [0.0] * dimension
        for row in range(dimension):
            curvature = 0.0
            for column in range(row + 1):
                curvature -= inverse_hessian[row][column] * slope[column]
            for column in range(row + 1, dimension):
                curvature -= inverse_hessian[column][row] * slope[column]
            direction[row] = curvature
            directional_derivative += curvature * slope[row]

        if directional_derivative < 0.0:
            step_length = 1.0
            accepted = False
            while True:
                step_count = 0
                for row in range(dimension):
                    point[row] = location[row] + step_length * direction[row]
                    if (_NO_CHANGE_SHIFT + location[row]) == (_NO_CHANGE_SHIFT + point[row]):
                        step_count += 1
                if step_count < dimension:
                    value = function_at(point)
                    function_count += 1
                    accepted = math.isfinite(value) and (
                        value <= best_value + directional_derivative * step_length * _ACCEPTANCE
                    )
                    if not accepted:
                        step_length *= _STEP_REDUCTION
                if step_count == dimension or accepted:
                    break
            improved_enough = (value > abstol) and (
                abs(value - best_value) > reltol * (abs(best_value) + reltol)
            )
            # A tiny improvement, or a value below abstol, ends this attempt.
            if not improved_enough:
                step_count = dimension
                best_value = value
            if step_count < dimension:
                best_value = value
                slope = gradient_at(point)
                gradient_count += 1
                iteration += 1
                slope_change_norm = 0.0
                for row in range(dimension):
                    direction[row] = step_length * direction[row]
                    slope_then[row] = slope[row] - slope_then[row]
                    slope_change_norm += direction[row] * slope_then[row]
                if slope_change_norm > 0.0:
                    curvature_norm = 0.0
                    curvature_direction = [0.0] * dimension
                    for row in range(dimension):
                        weighted = 0.0
                        for column in range(row + 1):
                            weighted += inverse_hessian[row][column] * slope_then[column]
                        for column in range(row + 1, dimension):
                            weighted += inverse_hessian[column][row] * slope_then[column]
                        curvature_direction[row] = weighted
                        curvature_norm += weighted * slope_then[row]
                    curvature_scale = 1.0 + curvature_norm / slope_change_norm
                    for row in range(dimension):
                        for column in range(row + 1):
                            inverse_hessian[row][column] += (
                                curvature_scale * direction[row] * direction[column]
                                - curvature_direction[row] * direction[column]
                                - direction[row] * curvature_direction[column]
                            ) / slope_change_norm
                else:
                    last_reset = gradient_count
            elif last_reset < gradient_count:
                # The factor is stale. Rebuild it and try one more search.
                step_count = 0
                last_reset = gradient_count
        else:
            # The factor does not give a downhill direction.
            step_count = 0
            if last_reset == gradient_count:
                step_count = dimension
            else:
                last_reset = gradient_count
        if trace:
            print(f"iter{iteration:4d} value {value:.6f}")
        if report is not None:
            report(iteration, best_value)
        if iteration >= maxit:
            break
        if gradient_count - last_reset > 2 * dimension:
            last_reset = gradient_count
        if not (step_count != dimension or last_reset != gradient_count):
            break

    if trace:
        print(f"final  value {best_value:.6f}")
        if iteration < maxit:
            print("converged")
        else:
            print(f"stopped after {iteration} iterations")
    status = 0 if iteration < maxit else 1
    return BFGSResult(
        x=np.asarray(point, dtype=np.float64),
        fun=float(best_value),
        success=status == 0,
        status=status,
        nit=int(iteration),
        nfev=int(function_count),
        njev=int(gradient_count),
    )
