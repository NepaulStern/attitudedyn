#ifndef ATTITUDYN_NORMAL_HPP
#define ATTITUDYN_NORMAL_HPP

// Standard normal distribution functions shared by the C++ likelihood and,
// through validation/fortran/math_bridge.cpp, by the Fortran reference
// oracle. One implementation so a comparison on identical draws is not
// confounded by two different approximations to Phi and Phi^{-1}.
//
// The quantile is Peter J. Acklam's rational approximation (public domain)
// with one Halley refinement step that uses the cdf below.

#include <cmath>
#include <limits>

namespace attitudedyn {

inline double standard_normal_cdf(double x) {
  // Phi(x) = 0.5 * erfc(-x / sqrt(2)). erfc is the correctly rounded
  // (or correctly directed) libm routine and does not underflow to zero
  // nearly as early as a naive series for the left tail.
  return 0.5 * std::erfc(-x * 0.707106781186547524400844362104849039);
}

inline double standard_normal_pdf(double x) {
  return 0.398942280401432677939946059934381868 *
         std::exp(-0.5 * x * x);
}

inline double standard_normal_log_pdf(double x) {
  return -0.918938533204672741780329736405617639 -
         0.5 * x * x;
}

// log Phi(x) for x < 0, via the asymptotic Mills expansion of the left
// tail, stopped before the series diverges. Used only when erfc underflows.
inline double standard_normal_log_cdf_left_tail(double x) {
  const double ax = -x;
  double series = 1.0;
  double term = 1.0;
  for (int k = 1; k <= 20; ++k) {
    const double next = term * (-(2.0 * k - 1.0) / (ax * ax));
    if (std::fabs(next) > std::fabs(term)) break;
    term = next;
    const double trial = series + term;
    if (trial <= 0.0) break;
    series = trial;
    if (std::fabs(term) < 1e-18 * series) break;
  }
  return standard_normal_log_pdf(x) - std::log(ax) + std::log(series);
}

inline double standard_normal_log_cdf(double x) {
  if (x > 8.3) {
    // 1 - Phi(x) is far below 1e-16; log Phi(x) = log1p(-(1-Phi)) ~ 0.
    const double sf = standard_normal_cdf(-x);
    if (sf == 0.0) return 0.0;
    return std::log1p(-sf);
  }
  if (x < -12.0) return standard_normal_log_cdf_left_tail(x);
  const double p = standard_normal_cdf(x);
  if (p <= 0.0) return -std::numeric_limits<double>::infinity();
  return std::log(p);
}

inline double standard_normal_log_sf(double x) {
  return standard_normal_log_cdf(-x);
}

// log(Phi(b) - Phi(a)), stable when both limits sit in the same tail.
// This replaces the Fortran tail blend (probmake1), which is discontinuous
// at -4 and -5. See notes/DISCREPANCIES.md.
inline double log_normal_interval(double lower, double upper) {
  if (!(upper > lower) || !std::isfinite(lower) || !std::isfinite(upper)) {
    return -std::numeric_limits<double>::infinity();
  }
  if (lower > 0.0) {
    const double log_sf_lo = standard_normal_log_sf(lower);
    const double log_sf_hi = standard_normal_log_sf(upper);
    // sf(lower) - sf(upper), and sf(lower) >= sf(upper).
    if (!std::isfinite(log_sf_lo)) return -std::numeric_limits<double>::infinity();
    return log_sf_lo + std::log1p(-std::exp(log_sf_hi - log_sf_lo));
  }
  if (upper < 0.0) {
    const double log_cdf_hi = standard_normal_log_cdf(upper);
    const double log_cdf_lo = standard_normal_log_cdf(lower);
    if (!std::isfinite(log_cdf_hi)) return -std::numeric_limits<double>::infinity();
    return log_cdf_hi + std::log1p(-std::exp(log_cdf_lo - log_cdf_hi));
  }
  const double p = standard_normal_cdf(upper) - standard_normal_cdf(lower);
  if (!(p > 0.0)) return -std::numeric_limits<double>::infinity();
  return std::log(p);
}

// Acklam's quantile, refined once by Halley's method. Accurate to nearly
// a double ulp on (1e-12, 1-1e-12) and still accurate far into the tails.
inline double standard_normal_quantile(double p) {
  if (!(p > 0.0) || !(p < 1.0) || !std::isfinite(p)) {
    if (p == 0.0) return -std::numeric_limits<double>::infinity();
    if (p == 1.0) return std::numeric_limits<double>::infinity();
    return std::numeric_limits<double>::quiet_NaN();
  }
  // Rational coefficients published by Peter J. Acklam.
  const double a1 = -3.969683028665376e+01;
  const double a2 = 2.209460984245205e+02;
  const double a3 = -2.759285104469687e+02;
  const double a4 = 1.383577518672690e+02;
  const double a5 = -3.066479806614716e+01;
  const double a6 = 2.506628277459239e+00;
  const double b1 = -5.447609879822406e+01;
  const double b2 = 1.615858368580409e+02;
  const double b3 = -1.556989798598866e+02;
  const double b4 = 6.680131188771972e+01;
  const double b5 = -1.328068155288572e+01;
  const double c1 = -7.784894002430293e-03;
  const double c2 = -3.223964580411365e-01;
  const double c3 = -2.400758277161838e+00;
  const double c4 = -2.549732539343734e+00;
  const double c5 = 4.374664141464968e+00;
  const double c6 = 2.938163982698783e+00;
  const double d1 = 7.784695709041462e-03;
  const double d2 = 3.224671290700398e-01;
  const double d3 = 2.445134137142996e+00;
  const double d4 = 3.754408661907416e+00;
  const double plow = 0.02425;
  const double phigh = 1.0 - plow;

  double x;
  // The lower-tail rational is already negative. The upper tail is its
  // reflection, so only that branch is negated. Swapping the signs puts a
  // jump at p = 0.02425 and 0.97575.
  if (p < plow) {
    const double q = std::sqrt(-2.0 * std::log(p));
    x = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
        ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
  } else if (p > phigh) {
    const double q = std::sqrt(-2.0 * std::log(1.0 - p));
    x = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
        ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
  } else {
    const double q = p - 0.5;
    const double r = q * q;
    x = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
        (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0);
  }

  // One Halley step: x <- x - u / (1 + x u / 2), u = (Phi(x) - p) / phi(x).
  const double err = standard_normal_cdf(x) - p;
  const double phi = standard_normal_pdf(x);
  if (phi > 0.0 && std::isfinite(err)) {
    const double u = err / phi;
    x -= u / (1.0 + 0.5 * x * u);
  }
  return x;
}

// Fortran probmake1, including the tail blend and pi = 3.1415.
// `cdf` is injected so the shared-math oracle and this function call the
// same Phi. The blend weight is the original one: it is 1 at the inner
// edge (-4) and 0 at -5, which is the reverse of a smooth handover and
// makes the approximation discontinuous at -4. Kept only for the
// Fortran comparison (numerics = "shared").
template <typename Cdf>
inline double fortran_log_interval(double arg1, double arg2, Cdf cdf) {
  const double big1 = 4.0;
  const double pi = 3.1415;
  const double small1 = -4.0;
  const double small2 = -5.0;
  if (arg2 >= small1 && arg1 <= big1) {
    const double diff = cdf(arg2) - cdf(arg1);
    if (!(diff > 0.0)) return -std::numeric_limits<double>::infinity();
    return std::log(diff);
  }
  auto tail = [&](double hi, double lo) {
    const double delta = hi - lo;
    const double fac1 = delta * ((2.0 * hi) - delta) / 2.0;
    const double fac2 = std::log(-hi);
    const double fac3 = hi * hi / 2.0;
    const double fac4 = std::log(2.0 * pi) / 2.0;
    const double approx = -(std::exp(fac1) + fac2 + fac3 + fac4);
    if (hi >= small2) {
      const double diff = cdf(hi) - cdf(lo);
      const double direct = (diff > 0.0)
                                ? std::log(diff)
                                : -std::numeric_limits<double>::infinity();
      const double weight = std::fabs(hi - small2);
      return weight * approx + (1.0 - weight) * direct;
    }
    return approx;
  };
  if (arg2 < small1) return tail(arg2, arg1);
  if (arg1 > big1) return tail(-arg1, -arg2);
  return -std::numeric_limits<double>::infinity();
}

// rho = 2 exp(c) / (1 + exp(c)) - 1, written as in parset212.
// For large |c| the algebraically equivalent tanh form is used so the
// optimiser can step far from the start without producing NaN. On the
// range of the published GSS coefficients the two expressions agree to
// well under 1e-15, and the Fortran comparison uses moderate c.
inline double persistence_from_unconstrained(double unconstrained) {
  if (unconstrained > 40.0 || unconstrained < -40.0) {
    const double e = std::exp(unconstrained > 0.0 ? -unconstrained : unconstrained);
    if (unconstrained > 0.0) return (1.0 - e) / (1.0 + e);
    return (e - 1.0) / (e + 1.0);
  }
  const double e = std::exp(unconstrained);
  return (2.0 * e / (1.0 + e)) - 1.0;
}

// d(rho)/dc = 2 exp(c) / (1 + exp(c))^2 = (1 - rho^2) / 2.
inline double persistence_unconstrained_derivative(double rho) {
  return 0.5 * (1.0 - rho * rho);
}

}  // namespace attitudedyn

#endif
