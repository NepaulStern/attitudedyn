#include "internal.hpp"

#include "attitudedyn/normal.hpp"
#include "draws.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <thread>
#include <utility>

namespace attitudedyn {
namespace detail {
namespace {

struct DrawSource {
  std::uint64_t seed = 1;
  DrawKind kind = DrawKind::kPhilox;
  int person = 0;
  int draw = 0;
  bool negate = false;
  int n_factors = 0;
  int n_waves = 0;
  int n_draws = 0;
  const float* factor_draws = nullptr;
  const float* shock_draws = nullptr;

  double factor(int factor) const {
    if (factor_draws != nullptr) {
      const float value =
          factor_draws[(static_cast<std::size_t>(person) * n_draws + draw) * n_factors + factor];
      return negate ? -static_cast<double>(value) : static_cast<double>(value);
    }
    return standard_normal_draw(seed, person, draw, factor, negate, kind);
  }

  // Generated shock. External arrays are read by external_shock(), which also
  // knows how many attitudes there are. This method is used only when the
  // external array is absent.
  double shock(int attitude, int wave) const {
    const int stream = n_factors + attitude * n_waves + wave;
    return standard_normal_draw(seed, person, draw, stream, negate, kind);
  }
};

double external_shock(const DrawSource& source, int n_attitudes, int attitude, int wave) {
  const float value = source.shock_draws[(((static_cast<std::size_t>(source.person) * source.n_draws +
                                             source.draw) *
                                            n_attitudes) +
                                           attitude) *
                                              source.n_waves +
                                          wave];
  return source.negate ? -static_cast<double>(value) : static_cast<double>(value);
}

void scaled_persistence(double rho, double drho_dc, double gap, double& rhop, double& d_rhop_dc) {
  // sign(rho) * |rho|^gap, as in the Fortran preplik12. The derivative is
  // with respect to the unconstrained coefficient c, with
  // rho = 2e^c/(1+e^c)-1.
  if (rho == 0.0) {
    rhop = 0.0;
    d_rhop_dc = (gap == 1.0) ? drho_dc : 0.0;
    return;
  }
  const double magnitude = std::max(std::fabs(rho), 1e-12);
  rhop = std::copysign(std::pow(std::fabs(rho), gap), rho);
  const double d_rho = gap * std::pow(magnitude, gap - 1.0);
  d_rhop_dc = d_rho * drho_dc;
}

double systematic_part(const Unpacked& theta, const Person& person, int attitude, int wave,
                       bool include_wave_effects) {
  double value = 0.0;
  if (include_wave_effects) {
    value += theta.wave_effect[attitude * theta.n_waves + wave];
  }
  // preplik12 multiplies x beta on every wave in the person's window, including
  // a gap. An unread wave is stored as a zero regressor, so it contributes only
  // the wave effect. A wave that was read but has no usable outcome still
  // contributes x beta.
  const double* x = person.covariates.data() + static_cast<std::size_t>(wave) * theta.n_covariates;
  const double* beta = theta.beta.data() + static_cast<std::size_t>(attitude) * theta.n_covariates;
  for (int column = 0; column < theta.n_covariates; ++column) {
    value += beta[column] * x[column];
  }
  return value;
}

void add_systematic_derivative(const Unpacked& theta, const Person& person, int attitude, int wave,
                               bool include_wave_effects, double* dy) {
  if (include_wave_effects) {
    const int slot = theta.wave_grad[attitude * theta.n_waves + wave];
    if (slot >= 0) dy[slot] += 1.0;
  }
  const double* x = person.covariates.data() + static_cast<std::size_t>(wave) * theta.n_covariates;
  for (int column = 0; column < theta.n_covariates; ++column) {
    const int slot = theta.beta_grad[attitude * theta.n_covariates + column];
    if (slot >= 0) dy[slot] += x[column];
  }
}

void add_cut_derivative(const Unpacked& theta, int attitude, int cut, double sign, double* dest) {
  // Cuts 0 and 1 are the normalised ends (-100 and 0). The top cut is +100.
  // Interior cuts are cumulative sums of exp(c).
  for (const ThresholdIncrement& piece : theta.threshold_increments[attitude]) {
    if (piece.grad >= 0 && cut >= piece.first_cut) dest[piece.grad] += sign * piece.exp_c;
  }
}

struct Truncation {
  double shock = 0.0;
  bool clamped = false;
  double phi_lower = 0.0;
  double phi_upper = 0.0;
  double phi_shock = 0.0;
  double uniform_of_draw = 0.0;
};

Truncation truncate_standard_normal(double raw_normal, double lower_limit, double upper_limit,
                                    Numerics numerics) {
  // Inverse-transform draw from the standard normal truncated to
  // (lower_limit, upper_limit), using Phi(raw_normal) as the uniform.
  // This is normadj1. The Fortran code then casts the uniform to single
  // precision; the shared and stable paths stay in double precision.
  Truncation out;
  const double cdf_lower = standard_normal_cdf(lower_limit);
  const double cdf_upper = standard_normal_cdf(upper_limit);
  const double cdf_draw = standard_normal_cdf(raw_normal);
  double uniform = cdf_lower + cdf_draw * (cdf_upper - cdf_lower);
  const double lowest = (numerics == Numerics::kShared) ? 1e-4 : 1e-12;
  const double highest = (numerics == Numerics::kShared) ? 0.9999 : 1.0 - 1e-12;
  if (uniform < lowest) {
    uniform = lowest;
    out.clamped = true;
  } else if (uniform > highest) {
    uniform = highest;
    out.clamped = true;
  }
  out.shock = standard_normal_quantile(uniform);
  out.phi_lower = standard_normal_pdf(lower_limit);
  out.phi_upper = standard_normal_pdf(upper_limit);
  out.phi_shock = standard_normal_pdf(out.shock);
  out.uniform_of_draw = cdf_draw;
  return out;
}

double log_probability(double lower_limit, double upper_limit, Numerics numerics) {
  if (numerics == Numerics::kShared) {
    return fortran_log_interval(lower_limit, upper_limit, standard_normal_cdf);
  }
  return log_normal_interval(lower_limit, upper_limit);
}

// One GHK path. `gradient` has length n_grad when with_gradient is true and
// receives the derivative of this path's log-likelihood only.
double one_path(const ModelSpec& spec, const Unpacked& theta, const Person& person,
                const DrawSource& source, Numerics numerics, bool with_gradient,
                double* path_gradient) {
  const int J = theta.n_attitudes;
  const int G = theta.n_grad;
  std::vector<double> factors(static_cast<std::size_t>(theta.n_factors));
  for (int factor = 0; factor < theta.n_factors; ++factor) factors[factor] = source.factor(factor);

  std::vector<double> state_prev(static_cast<std::size_t>(J), 0.0);
  std::vector<double> state_cur(static_cast<std::size_t>(J), 0.0);
  std::vector<double> dstate_prev(with_gradient ? static_cast<std::size_t>(J * G) : 0);
  std::vector<double> dstate_cur(with_gradient ? static_cast<std::size_t>(J * G) : 0);
  std::vector<double> d_pre(with_gradient ? static_cast<std::size_t>(G) : 0);
  std::vector<double> d_arg(with_gradient ? static_cast<std::size_t>(G) : 0);

  double path_loglik = 0.0;
  if (with_gradient) std::fill(path_gradient, path_gradient + G, 0.0);

  for (int wave = person.first_wave; wave <= person.last_wave; ++wave) {
    const bool has_lag = wave > person.first_wave;
    const double gap = has_lag ? (spec.wave_time[wave] - spec.wave_time[wave - 1]) : 0.0;
    for (int attitude = 0; attitude < J; ++attitude) {
      const int category = person.category[attitude * spec.n_waves + wave];
      const bool observed = person.interviewed[wave] && category >= 1 &&
                            category <= spec.n_categories[attitude];

      // Lagged contribution: the persistence matrix times the previous latent
      // attitude, as in likobs12. The covariate index, the wave effect and the
      // factor are inside that latent attitude, so they carry forward too.
      double lagged = 0.0;
      if (with_gradient) std::fill(d_pre.begin(), d_pre.end(), 0.0);
      if (has_lag) {
        for (int source_attitude = 0; source_attitude < J; ++source_attitude) {
          const int slot = attitude * J + source_attitude;
          double rhop = 0.0;
          double d_rhop = 0.0;
          scaled_persistence(theta.rho[slot], theta.drho_dc[slot], gap, rhop, d_rhop);
          lagged += rhop * state_prev[source_attitude];
          if (with_gradient) {
            const int g = theta.rho_grad[slot];
            if (g >= 0) d_pre[g] += d_rhop * state_prev[source_attitude];
            const double* previous = dstate_prev.data() + static_cast<std::size_t>(source_attitude) * G;
            for (int column = 0; column < G; ++column) d_pre[column] += rhop * previous[column];
          }
        }
      }

      const double xb = systematic_part(theta, person, attitude, wave, spec.include_wave_effects);
      double factor_mean = 0.0;
      for (int factor = 0; factor < theta.n_factors; ++factor) {
        factor_mean += theta.loading[factor * J + attitude] * factors[factor];
      }
      const double y_pre = xb + factor_mean + lagged;
      if (with_gradient) {
        add_systematic_derivative(theta, person, attitude, wave, spec.include_wave_effects,
                                  d_pre.data());
        for (int factor = 0; factor < theta.n_factors; ++factor) {
          const int g = theta.loading_grad[factor * J + attitude];
          if (g >= 0) d_pre[g] += factors[factor];
        }
      }

      const double raw = (source.shock_draws == nullptr)
                             ? source.shock(attitude, wave)
                             : external_shock(source, J, attitude, wave);

      double state_new = 0.0;
      if (!observed) {
        // A missing item, or a wave the person skipped, is integrated by
        // simulating it forward. It does not contribute to the likelihood.
        // The shock is the untruncated normal draw (likobs12, the second branch).
        state_new = y_pre + raw;
        if (with_gradient) {
          double* dest = dstate_cur.data() + static_cast<std::size_t>(attitude) * G;
          std::copy(d_pre.begin(), d_pre.end(), dest);
        }
      } else {
        const int base = attitude * theta.max_cuts;
        const double lower = theta.threshold[base + (category - 1)];
        const double upper = theta.threshold[base + category];
        // Identity contemporaneous correlation: the innovation is N(0, 1)
        // and is not conditioned on the other attitudes in this wave.
        // Equation (1) plus the ordered-probit interval.
        const double arg_lower = lower - y_pre;
        const double arg_upper = upper - y_pre;
        path_loglik += log_probability(arg_lower, arg_upper, numerics);

        const Truncation draw = truncate_standard_normal(raw, arg_lower, arg_upper, numerics);
        state_new = y_pre + draw.shock;

        if (with_gradient) {
          const double lp = log_probability(arg_lower, arg_upper, Numerics::kStable);
          const double weight_lower =
              std::exp(std::min(700.0, standard_normal_log_pdf(arg_lower) - lp));
          const double weight_upper =
              std::exp(std::min(700.0, standard_normal_log_pdf(arg_upper) - lp));
          // d log P = -w_lo * d(arg_lo) + w_hi * d(arg_hi)
          // arg = cut - y_pre, so the y_pre part is (w_lo - w_hi) * dy_pre
          // and the cut part is added from the threshold increments.
          for (int column = 0; column < G; ++column) {
            path_gradient[column] += (weight_lower - weight_upper) * d_pre[column];
          }
          add_cut_derivative(theta, attitude, category - 1, -weight_lower, path_gradient);
          add_cut_derivative(theta, attitude, category, weight_upper, path_gradient);

          double* dest = dstate_cur.data() + static_cast<std::size_t>(attitude) * G;
          if (draw.clamped || !(draw.phi_shock > 0.0)) {
            std::copy(d_pre.begin(), d_pre.end(), dest);
          } else {
            // dutc/da = phi(limit) * (weight of that limit) / phi(utc).
            const double scale_lower =
                draw.phi_lower * (1.0 - draw.uniform_of_draw) / draw.phi_shock;
            const double scale_upper = draw.phi_upper * draw.uniform_of_draw / draw.phi_shock;
            // utc = qnorm(F(lo) + Phi(raw) * (F(hi) - F(lo)))
            // d utc = scale_lo * d(lo - y_pre) + scale_hi * d(hi - y_pre)
            //       = -(scale_lo + scale_hi) * dy_pre + scale_lo * d lo + scale_hi * d hi
            const double scale_y = scale_lower + scale_upper;
            std::fill(d_arg.begin(), d_arg.end(), 0.0);
            add_cut_derivative(theta, attitude, category - 1, scale_lower, d_arg.data());
            add_cut_derivative(theta, attitude, category, scale_upper, d_arg.data());
            for (int column = 0; column < G; ++column) {
              dest[column] = d_pre[column] * (1.0 - scale_y) + d_arg[column];
            }
          }
        }
      }
      state_cur[attitude] = state_new;
    }
    state_prev.swap(state_cur);
    if (with_gradient) dstate_prev.swap(dstate_cur);
  }
  return path_loglik;
}

bool cholesky_factor(const std::vector<double>& omega, int n, std::vector<double>& factor) {
  factor.assign(static_cast<std::size_t>(n * n), 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      double sum = omega[i * n + j];
      for (int k = 0; k < j; ++k) sum -= factor[i * n + k] * factor[j * n + k];
      if (j == i) {
        if (!(sum > 1e-14)) return false;
        factor[i * n + i] = std::sqrt(sum);
      } else {
        factor[i * n + j] = sum / factor[j * n + j];
      }
    }
  }
  return true;
}

// GHK path when the innovations in a wave are correlated. Gradients are not
// returned: evaluate() rejects that combination. fortran_ghk follows ghk1.
// Each path starts from a zero latent state.
double one_correlated_path(const ModelSpec& spec, const Unpacked& theta, const Person& person,
                           const DrawSource& source, Numerics numerics, bool fortran_ghk) {
  const int J = theta.n_attitudes;
  std::vector<double> factor_loading_draw(static_cast<std::size_t>(theta.n_factors));
  for (int factor = 0; factor < theta.n_factors; ++factor) {
    factor_loading_draw[factor] = source.factor(factor);
  }
  std::vector<double> chol;
  if (!fortran_ghk && !cholesky_factor(theta.omega, J, chol)) {
    return -std::numeric_limits<double>::infinity();
  }
  std::vector<double> state_prev(static_cast<std::size_t>(J), 0.0);
  std::vector<double> state_cur(static_cast<std::size_t>(J), 0.0);
  std::vector<double> eta(static_cast<std::size_t>(J), 0.0);
  double path_loglik = 0.0;

  for (int wave = person.first_wave; wave <= person.last_wave; ++wave) {
    const bool has_lag = wave > person.first_wave;
    const double gap = has_lag ? (spec.wave_time[wave] - spec.wave_time[wave - 1]) : 0.0;
    for (int attitude = 0; attitude < J; ++attitude) {
      double lagged = 0.0;
      if (has_lag) {
        for (int source_attitude = 0; source_attitude < J; ++source_attitude) {
          double rhop = 0.0;
          double unused = 0.0;
          scaled_persistence(theta.rho[attitude * J + source_attitude], 0.0, gap, rhop, unused);
          lagged += rhop * state_prev[source_attitude];
        }
      }
      const double xb = systematic_part(theta, person, attitude, wave, spec.include_wave_effects);
      double factor_mean = 0.0;
      for (int factor = 0; factor < theta.n_factors; ++factor) {
        factor_mean += theta.loading[factor * J + attitude] * factor_loading_draw[factor];
      }
      const double raw = (source.shock_draws == nullptr)
                             ? source.shock(attitude, wave)
                             : external_shock(source, J, attitude, wave);
      const int category = person.category[attitude * spec.n_waves + wave];
      const bool observed = person.interviewed[wave] && category >= 1 &&
                            category <= spec.n_categories[attitude];

      double conditional_mean = 0.0;
      double conditional_sd = 1.0;
      if (!fortran_ghk) {
        for (int earlier = 0; earlier < attitude; ++earlier) {
          conditional_mean += chol[attitude * J + earlier] * eta[earlier];
        }
        conditional_sd = chol[attitude * J + attitude];
      } else if (attitude > 0) {
        // ghk1 on the leading (attitude+1) block, with ystart reversed and
        // aone overwritten on every pass of the inner loop.
        const int n = attitude + 1;
        std::vector<double> block(static_cast<std::size_t>(n * n));
        for (int row = 0; row < n; ++row) {
          for (int col = 0; col < n; ++col) block[row * n + col] = theta.omega[row * J + col];
        }
        // Invert `block` in place by Gauss-Jordan, which for the identity
        // returns the identity and therefore the same moments as ghk1.
        std::vector<double> inv(static_cast<std::size_t>(n * n), 0.0);
        for (int i = 0; i < n; ++i) inv[i * n + i] = 1.0;
        for (int col = 0; col < n; ++col) {
          int pivot = col;
          double best = std::fabs(block[col * n + col]);
          for (int row = col + 1; row < n; ++row) {
            const double value = std::fabs(block[row * n + col]);
            if (value > best) {
              best = value;
              pivot = row;
            }
          }
          if (!(best > 0.0)) return -std::numeric_limits<double>::infinity();
          if (pivot != col) {
            for (int k = 0; k < n; ++k) {
              std::swap(block[col * n + k], block[pivot * n + k]);
              std::swap(inv[col * n + k], inv[pivot * n + k]);
            }
          }
          const double diag = block[col * n + col];
          for (int k = 0; k < n; ++k) {
            block[col * n + k] /= diag;
            inv[col * n + k] /= diag;
          }
          for (int row = 0; row < n; ++row) {
            if (row == col) continue;
            const double f = block[row * n + col];
            for (int k = 0; k < n; ++k) {
              block[row * n + k] -= f * block[col * n + k];
              inv[row * n + k] -= f * inv[col * n + k];
            }
          }
        }
        // Reversed ystart. The current attitude's slot still holds the
        // previous wave's state (or zero on the first wave of the path).
        std::vector<double> ystart(static_cast<std::size_t>(n));
        ystart[0] = has_lag ? state_prev[attitude] : 0.0;
        for (int m = 1; m < n; ++m) ystart[m] = state_cur[attitude - m];
        const double aone = inv[0] * inv[(n - 1) * n + 0];
        double adjmean = 0.0;
        for (int i = 0; i < n - 1; ++i) adjmean -= aone * ystart[i];
        conditional_mean = adjmean;
        conditional_sd = 1.0 / std::sqrt(inv[0]);
      }

      if (!observed) {
        const double innovation = conditional_mean + conditional_sd * raw;
        eta[attitude] = raw;
        state_cur[attitude] = xb + factor_mean + lagged + innovation;
      } else {
        const int base = attitude * theta.max_cuts;
        const double y_without_innovation = xb + factor_mean + lagged;
        const double lower = theta.threshold[base + (category - 1)];
        const double upper = theta.threshold[base + category];
        const double arg_lower = (lower - y_without_innovation - conditional_mean) / conditional_sd;
        const double arg_upper = (upper - y_without_innovation - conditional_mean) / conditional_sd;
        path_loglik += log_probability(arg_lower, arg_upper, numerics);
        const Truncation draw = truncate_standard_normal(raw, arg_lower, arg_upper, numerics);
        eta[attitude] = draw.shock;
        const double innovation = conditional_mean + conditional_sd * draw.shock;
        state_cur[attitude] = y_without_innovation + innovation;
      }
    }
    state_prev.swap(state_cur);
  }
  return path_loglik;
}

}  // namespace

double person_contribution(const ModelSpec& spec, const Unpacked& theta, const Person& person,
                           int person_index, const DrawOptions& draws, const float* factor_draws,
                           const float* shock_draws, Numerics numerics, bool fortran_ghk,
                           double* gradient, double* /*score_slot*/) {
  if (person.last_wave < person.first_wave || !(person.weight > 0.0)) {
    if (gradient != nullptr) std::fill(gradient, gradient + theta.n_grad, 0.0);
    return 0.0;
  }
  const int n_paths = draws.n_draws * (draws.antithetic ? 2 : 1);
  const bool with_gradient = gradient != nullptr && theta.n_grad > 0;
  if (gradient != nullptr && theta.n_grad == 0) {
    // No differentiated coefficients.
  }
  std::vector<double> path_loglik(static_cast<std::size_t>(n_paths));
  std::vector<double> path_gradient(with_gradient ? static_cast<std::size_t>(n_paths * theta.n_grad) : 0);

  const bool correlated = !theta.identity_correlation;
  int path = 0;
  for (int draw = 0; draw < draws.n_draws; ++draw) {
    for (int anti = 0; anti < (draws.antithetic ? 2 : 1); ++anti) {
      DrawSource source;
      source.seed = draws.seed;
      source.kind = draws.kind;
      source.person = person_index;
      source.draw = draw;
      source.negate = anti == 1;
      source.n_factors = theta.n_factors;
      source.n_waves = theta.n_waves;
      source.n_draws = draws.n_draws;
      source.factor_draws = factor_draws;
      source.shock_draws = shock_draws;
      double* path_grad = with_gradient ? path_gradient.data() + static_cast<std::size_t>(path) * theta.n_grad
                                        : nullptr;
      if (correlated) {
        path_loglik[path] =
            one_correlated_path(spec, theta, person, source, numerics, fortran_ghk);
      } else if (with_gradient) {
        path_loglik[path] = one_path(spec, theta, person, source, numerics, true, path_grad);
      } else {
        path_loglik[path] = one_path(spec, theta, person, source, numerics, false, nullptr);
      }
      ++path;
    }
  }

  // log-mean-exp over the antithetic paths. Same accumulation order as
  // likobs12: draw index, then the original sign, then the negated sign.
  double maximum = path_loglik[0];
  for (int i = 1; i < n_paths; ++i) maximum = std::max(maximum, path_loglik[i]);
  if (!std::isfinite(maximum)) {
    if (gradient != nullptr) std::fill(gradient, gradient + theta.n_grad, 0.0);
    return -std::numeric_limits<double>::infinity();
  }
  double sum_exp = 0.0;
  for (int i = 0; i < n_paths; ++i) sum_exp += std::exp(path_loglik[i] - maximum);
  const double loglik = maximum + std::log(sum_exp) - std::log(static_cast<double>(n_paths));

  if (gradient != nullptr) {
    std::fill(gradient, gradient + theta.n_grad, 0.0);
    if (with_gradient && sum_exp > 0.0 && std::isfinite(loglik)) {
      for (int i = 0; i < n_paths; ++i) {
        const double weight = std::exp(path_loglik[i] - maximum) / sum_exp;
        const double* path_grad = path_gradient.data() + static_cast<std::size_t>(i) * theta.n_grad;
        for (int column = 0; column < theta.n_grad; ++column) {
          gradient[column] += weight * path_grad[column];
        }
      }
      for (int column = 0; column < theta.n_grad; ++column) gradient[column] *= person.weight;
    }
  }
  if (!std::isfinite(loglik)) return -std::numeric_limits<double>::infinity();
  return person.weight * loglik;
}

}  // namespace detail

namespace {

void check_problem(const Problem& problem) {
  if (static_cast<int>(problem.theta_full.size()) != problem.layout.n_full) {
    throw std::invalid_argument("The stored coefficient vector does not match the layout.");
  }
  if (problem.draws.n_draws < 1) {
    throw std::invalid_argument("The number of simulation draws must be at least 1.");
  }
  const int J = problem.spec.n_attitudes;
  const int T = problem.spec.n_waves;
  const int X = problem.spec.n_covariates;
  for (int i = 0; i < static_cast<int>(problem.persons.size()); ++i) {
    const Person& person = problem.persons[i];
    if (person.last_wave < person.first_wave) continue;
    if (person.first_wave < 0 || person.last_wave >= T) {
      throw std::invalid_argument("A person's wave window is outside 1..n_waves.");
    }
    if (static_cast<int>(person.interviewed.size()) != T ||
        static_cast<int>(person.category.size()) != J * T ||
        static_cast<int>(person.covariates.size()) != T * X) {
      throw std::invalid_argument(
          "A person's interviewed, category, or covariate array has the wrong length.");
    }
    if (!std::isfinite(person.weight) || person.weight < 0.0) {
      throw std::invalid_argument("A person weight is negative or not finite.");
    }
    for (double value : person.covariates) {
      if (!std::isfinite(value)) {
        throw std::invalid_argument("A covariate is not finite. Missing covariates are not "
                                    "treated as missing outcomes; drop or impute them first.");
      }
    }
  }
}

Evaluation evaluate_impl(const Problem& problem, const std::vector<double>& theta_full,
                         const EvaluateRequest& request) {
  check_problem(problem);
  if (static_cast<int>(theta_full.size()) != problem.layout.n_full) {
    throw std::invalid_argument("theta has the wrong length for this model.");
  }
  const int n_full = problem.layout.n_full;
  std::vector<int> grad_of(static_cast<std::size_t>(n_full), -1);
  for (int column = 0; column < static_cast<int>(request.grad_columns.size()); ++column) {
    const int position = request.grad_columns[column];
    if (position < 0 || position >= n_full) {
      throw std::invalid_argument("A gradient column refers to a coefficient that does not exist.");
    }
    if (grad_of[position] != -1) {
      throw std::invalid_argument("A coefficient was listed twice in the gradient.");
    }
    grad_of[position] = column;
  }
  const int n_grad = static_cast<int>(request.grad_columns.size());
  const detail::Unpacked unpacked = detail::unpack(problem.spec, theta_full, grad_of);
  if (!unpacked.parameters_ok) {
    // An exploratory step walked off the parameter domain. Return a
    // non-finite objective and a finite gradient so a line search can
    // reject the step instead of aborting the fit.
    Evaluation rejected;
    rejected.objective = -std::numeric_limits<double>::infinity();
    rejected.person_loglik.assign(static_cast<std::size_t>(problem.persons.size()), 0.0);
    rejected.gradient.assign(static_cast<std::size_t>(n_grad), 0.0);
    return rejected;
  }
  if (n_grad > 0 && request.numerics != Numerics::kStable) {
    throw std::invalid_argument(
        "Analytic gradients are available for numerics = \"stable\" only. "
        "The Fortran tail blend is not differentiable at its knots.");
  }
  if (n_grad > 0 && !unpacked.identity_correlation) {
    throw std::invalid_argument(
        "Analytic gradients are implemented when the contemporaneous correlation "
        "is the identity, which is the published model. With a free correlation, "
        "ask for the likelihood only and use a numeric gradient.");
  }
  const bool external_factors = !request.factor_draws.empty();
  const bool external_shocks = !request.shock_draws.empty();
  if (problem.spec.n_factors > 0 && external_factors != external_shocks) {
    throw std::invalid_argument("Factor draws and shock draws must be supplied together.");
  }
  if (problem.spec.n_factors == 0 && external_factors) {
    throw std::invalid_argument("This model has no factors, so it does not take factor draws.");
  }
  const int n_persons = static_cast<int>(problem.persons.size());
  if (external_factors || external_shocks) {
    const std::size_t n_factor =
        static_cast<std::size_t>(n_persons) * problem.draws.n_draws * std::max(problem.spec.n_factors, 0);
    const std::size_t n_shock = static_cast<std::size_t>(n_persons) * problem.draws.n_draws *
                                problem.spec.n_attitudes * problem.spec.n_waves;
    if (external_factors && request.factor_draws.size() != n_factor) {
      throw std::invalid_argument("The external factor-draw array has the wrong length.");
    }
    if (!external_shocks || request.shock_draws.size() != n_shock) {
      throw std::invalid_argument("The external shock-draw array has the wrong length.");
    }
  }

  Evaluation out;
  out.penalty = unpacked.penalty;
  out.person_loglik.assign(static_cast<std::size_t>(n_persons), 0.0);
  out.gradient.assign(static_cast<std::size_t>(n_grad), 0.0);
  if (request.want_scores && n_grad > 0) {
    out.person_score.assign(static_cast<std::size_t>(n_persons) * n_grad, 0.0);
  }
  if (n_persons == 0) {
    out.objective = -out.penalty;
    return out;
  }

  const float* factor_ptr = request.factor_draws.empty() ? nullptr : request.factor_draws.data();
  const float* shock_ptr = request.shock_draws.empty() ? nullptr : request.shock_draws.data();
  // Person gradients are stored and then added in person order, so the total
  // does not depend on how the persons are split across threads.
  if (n_grad > 0 && out.person_score.empty()) {
    out.person_score.assign(static_cast<std::size_t>(n_persons) * n_grad, 0.0);
  }
  auto run_range_stored = [&](int begin, int end) {
    for (int person = begin; person < end; ++person) {
      double* score = n_grad > 0 ? out.person_score.data() + static_cast<std::size_t>(person) * n_grad
                                 : nullptr;
      out.person_loglik[person] = detail::person_contribution(
          problem.spec, unpacked, problem.persons[person], person, problem.draws, factor_ptr,
          shock_ptr, request.numerics, request.fortran_ghk, score, score);
    }
  };

  int n_threads = request.n_threads;
  if (n_threads < 1) {
    throw std::invalid_argument("n_threads must be at least 1.");
  }
  if (n_threads > n_persons) n_threads = n_persons;
  if (n_threads == 1) {
    run_range_stored(0, n_persons);
  } else {
    std::vector<std::thread> threads;
    threads.reserve(static_cast<std::size_t>(n_threads));
    const int chunk = (n_persons + n_threads - 1) / n_threads;
    for (int thread_index = 0; thread_index < n_threads; ++thread_index) {
      const int begin = thread_index * chunk;
      const int end = std::min(n_persons, begin + chunk);
      if (begin >= end) break;
      threads.emplace_back(run_range_stored, begin, end);
    }
    for (std::thread& thread : threads) thread.join();
  }

  double total = 0.0;
  for (int person = 0; person < n_persons; ++person) {
    if (problem.persons[person].last_wave >= problem.persons[person].first_wave &&
        problem.persons[person].weight > 0.0) {
      ++out.n_persons_used;
    }
    const double contribution = out.person_loglik[person];
    if (std::isfinite(contribution)) total += contribution;
    else total = -std::numeric_limits<double>::infinity();
  }
  if (n_grad > 0 && std::isfinite(total)) {
    for (int person = 0; person < n_persons; ++person) {
      const double* score = out.person_score.data() + static_cast<std::size_t>(person) * n_grad;
      for (int column = 0; column < n_grad; ++column) out.gradient[column] += score[column];
    }
  }
  // The penalty is subtracted once from the sample log-likelihood, as in
  // likely12. It is not added inside each person (covar does that, which
  // multiplies the penalty's score by the sample size; see the notes).
  out.objective = std::isfinite(total) ? total - out.penalty : total;
  // A non-finite derivative must not reach a quasi-Newton update. Zero it.
  // The objective itself stays put: a finite log-likelihood is still a
  // usable trial, and a non-finite one is already the rejected step.
  if (!std::isfinite(out.objective) || n_grad > 0) {
    for (double& component : out.gradient) {
      if (!std::isfinite(out.objective) || !std::isfinite(component)) component = 0.0;
    }
    if (!std::isfinite(out.objective) && !out.person_score.empty()) {
      std::fill(out.person_score.begin(), out.person_score.end(), 0.0);
    }
  }
  if (!request.want_scores) out.person_score.clear();
  return out;
}

}  // namespace

Evaluation evaluate_full(const Problem& problem, const std::vector<double>& theta_full,
                         const EvaluateRequest& request) {
  return evaluate_impl(problem, theta_full, request);
}

Evaluation evaluate(const Problem& problem, const std::vector<double>& theta_free,
                    const EvaluateRequest& request) {
  if (static_cast<int>(theta_free.size()) != problem.layout.n_free) {
    throw std::invalid_argument("The free coefficient vector has length " +
                                std::to_string(theta_free.size()) + " but the layout has " +
                                std::to_string(problem.layout.n_free) + " free coefficients.");
  }
  std::vector<double> theta = problem.theta_full;
  if (static_cast<int>(theta.size()) != problem.layout.n_full) {
    throw std::invalid_argument("problem.theta_full does not match the layout.");
  }
  int cursor = 0;
  for (int i = 0; i < problem.layout.n_full; ++i) {
    if (problem.layout.entries[i].is_free) theta[i] = theta_free[cursor++];
  }
  EvaluateRequest expanded = request;
  if (expanded.grad_columns.empty() && request.grad_columns.empty()) {
    // The default gradient is the free set, but only when the caller asks for
    // a gradient by putting something in grad_columns. An empty list means
    // "likelihood only", so a fit has to pass the free positions explicitly.
  }
  return evaluate_impl(problem, theta, expanded);
}

}  // namespace attitudedyn
