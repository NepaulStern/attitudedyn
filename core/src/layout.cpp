#include "internal.hpp"

#include "attitudedyn/normal.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace attitudedyn {
namespace {

std::string clean_name(std::string name, const std::string& fallback) {
  if (name.empty()) name = fallback;
  for (char& c : name) {
    if (c == ',' || c == '[' || c == ']' || c == ' ') c = '_';
  }
  return name;
}

std::string attitude_name(const ModelSpec& spec, int attitude) {
  if (attitude >= 0 && attitude < static_cast<int>(spec.attitude_names.size()) &&
      !spec.attitude_names[attitude].empty()) {
    return clean_name(spec.attitude_names[attitude], "");
  }
  return "y" + std::to_string(attitude + 1);
}

std::string covariate_name(const ModelSpec& spec, int covariate) {
  if (covariate >= 0 && covariate < static_cast<int>(spec.covariate_names.size()) &&
      !spec.covariate_names[covariate].empty()) {
    return clean_name(spec.covariate_names[covariate], "");
  }
  return "x" + std::to_string(covariate + 1);
}

std::string wave_name(const ModelSpec& spec, int wave) {
  if (wave >= 0 && wave < static_cast<int>(spec.wave_names.size()) &&
      !spec.wave_names[wave].empty()) {
    return clean_name(spec.wave_names[wave], "");
  }
  return std::to_string(wave + 1);
}

bool is_intercept_column(const ModelSpec& spec, int column) {
  for (int intercept : spec.intercept_columns) {
    if (intercept == column) return true;
  }
  return false;
}

void check_spec(const ModelSpec& spec) {
  if (spec.n_attitudes < 1) {
    throw std::invalid_argument("The model needs at least one ordinal attitude.");
  }
  if (spec.n_waves < 1) {
    throw std::invalid_argument("The model needs at least one wave.");
  }
  if (spec.n_covariates < 1) {
    throw std::invalid_argument("The model needs at least one covariate column.");
  }
  if (spec.n_factors < 0) {
    throw std::invalid_argument("The number of person factors cannot be negative.");
  }
  if (static_cast<int>(spec.n_categories.size()) != spec.n_attitudes) {
    throw std::invalid_argument(
        "n_categories must have one entry for each attitude.");
  }
  if (static_cast<int>(spec.wave_time.size()) != spec.n_waves) {
    throw std::invalid_argument(
        "wave_time must have one strictly increasing value for each wave.");
  }
  for (int t = 0; t < spec.n_waves; ++t) {
    if (!std::isfinite(spec.wave_time[t])) {
      throw std::invalid_argument("A wave time is not finite.");
    }
    if (t > 0 && !(spec.wave_time[t] > spec.wave_time[t - 1])) {
      throw std::invalid_argument(
          "Wave times must be strictly increasing. Gaps are wave_time[t] - wave_time[t-1], "
          "in whatever units the persistence exponent should use.");
    }
  }
  for (int j = 0; j < spec.n_attitudes; ++j) {
    if (spec.n_categories[j] < 2) {
      throw std::invalid_argument(
          "Each attitude needs at least two ordered categories.");
    }
    if (spec.n_categories[j] > 12) {
      throw std::invalid_argument(
          "Attitudes with more than 12 categories are not supported.");
    }
  }
  if (!spec.empty_wave_effect.empty() &&
      static_cast<int>(spec.empty_wave_effect.size()) != spec.n_attitudes * spec.n_waves) {
    throw std::invalid_argument(
        "empty_wave_effect must be empty or have length n_attitudes * n_waves.");
  }
  for (int column : spec.intercept_columns) {
    if (column < 0 || column >= spec.n_covariates) {
      throw std::invalid_argument("An intercept column index is outside the covariate matrix.");
    }
  }
  const int streams = spec.n_factors + spec.n_attitudes * spec.n_waves;
  if (streams > 1025) {
    throw std::invalid_argument(
        "factors + attitudes * waves exceeds 1024 draw streams.");
  }
}

}  // namespace

Layout make_layout(const ModelSpec& spec) {
  check_spec(spec);
  Layout layout;
  auto add = [&](std::string name, std::string block, int attitude, int index2, bool free) {
    LayoutEntry entry;
    entry.name = std::move(name);
    entry.block = std::move(block);
    entry.attitude = attitude;
    entry.index2 = index2;
    entry.is_free = free ? 1 : 0;
    layout.entries.push_back(std::move(entry));
  };

  if (spec.include_wave_effects) {
    for (int attitude = 0; attitude < spec.n_attitudes; ++attitude) {
      for (int wave = 0; wave < spec.n_waves; ++wave) {
        const bool empty = !spec.empty_wave_effect.empty() &&
                           spec.empty_wave_effect[attitude * spec.n_waves + wave];
        add("wave_effect[" + attitude_name(spec, attitude) + "," + wave_name(spec, wave) + "]",
            "wave_effect", attitude, wave, !empty);
      }
    }
  }
  for (int attitude = 0; attitude < spec.n_attitudes; ++attitude) {
    for (int column = 0; column < spec.n_covariates; ++column) {
      const bool intercept = is_intercept_column(spec, column);
      const bool fixed = intercept && spec.fix_intercepts && spec.include_wave_effects;
      add("beta[" + attitude_name(spec, attitude) + "," + covariate_name(spec, column) + "]",
          "beta", attitude, column, !fixed);
    }
  }
  // Persistence is stored and applied as in the Fortran program: row = receiver,
  // column = source. The printed unconstrained coefficient is what the
  // optimiser sees; the recursion uses the transformed value in (-1, 1).
  for (int to = 0; to < spec.n_attitudes; ++to) {
    for (int from = 0; from < spec.n_attitudes; ++from) {
      add("persistence[" + attitude_name(spec, to) + "," + attitude_name(spec, from) + "]",
          "persistence", to, from, true);
    }
  }
  for (int factor = 0; factor < spec.n_factors; ++factor) {
    for (int attitude = 0; attitude < spec.n_attitudes; ++attitude) {
      const bool last = spec.fix_last_loading && spec.n_factors >= 2 &&
                        factor == spec.n_factors - 1 && attitude == spec.n_attitudes - 1;
      add("loading[" + attitude_name(spec, attitude) + ",f" + std::to_string(factor + 1) + "]",
          "loading", attitude, factor, !last);
    }
  }
  for (int to = 1; to < spec.n_attitudes; ++to) {
    for (int from = 0; from < to; ++from) {
      add("chol[" + attitude_name(spec, to) + "," + attitude_name(spec, from) + "]",
          "chol", to, from, spec.free_contemporaneous);
    }
  }
  for (int attitude = 0; attitude < spec.n_attitudes; ++attitude) {
    const int n_increments = spec.n_categories[attitude] - 2;
    for (int increment = 0; increment < n_increments; ++increment) {
      // Cut numbering matches the Fortran labels: tau_3 is the first free
      // cut, the boundary between category 2 and category 3.
      const int cut_number = increment + 3;
      add("threshold[" + attitude_name(spec, attitude) + "," + std::to_string(cut_number) + "]",
          "threshold", attitude, cut_number, true);
    }
  }

  layout.n_full = static_cast<int>(layout.entries.size());
  layout.n_free = 0;
  for (const auto& entry : layout.entries) {
    if (entry.is_free) ++layout.n_free;
  }
  return layout;
}

namespace detail {
namespace {

double build_correlation(int n_attitudes, const std::vector<double>& chol, std::vector<double>& omega,
                         double& penalty) {
  // Contemporaneous correlation of the innovations.
  // The free parameters are the strict lower triangle of a factor, with the
  // factor's own diagonal taken to be 0 (a zero-initialised Fortran local).
  // The correlation diagonal is then forced to 1, and a quadratic penalty is
  // added if that is not possible. This is omegamake1. The published GSS and
  // BESIP runs fix every entry at 0, so the correlation is the identity and
  // the penalty is 0.
  omega.assign(static_cast<std::size_t>(n_attitudes * n_attitudes), 0.0);
  omega[0] = 1.0;
  penalty = 0.0;
  const double tolerance = 0.01;
  const double penalty_rate = 100.0;
  int cursor = 0;
  std::vector<double> factor(static_cast<std::size_t>(n_attitudes * n_attitudes), 0.0);
  for (int i = 1; i < n_attitudes; ++i) {
    for (int j = 0; j < i; ++j) {
      factor[i * n_attitudes + j] = chol[cursor++];
    }
  }
  for (int i = 1; i < n_attitudes; ++i) {
    for (int j = 0; j < i; ++j) {
      double value = 0.0;
      for (int k = 0; k <= j; ++k) {
        value += factor[i * n_attitudes + k] * factor[j * n_attitudes + k];
      }
      omega[i * n_attitudes + j] = value;
      omega[j * n_attitudes + i] = value;
    }
    double residual = 1.0;
    for (int j = 0; j < i; ++j) {
      const double value = omega[i * n_attitudes + j];
      residual -= value * value;
    }
    omega[i * n_attitudes + i] = 1.0;
    if (residual < tolerance) {
      const double deviation = tolerance - residual;
      penalty += penalty_rate * deviation * deviation;
      const double ratio = (1.0 - tolerance) / (1.0 - residual);
      for (int j = 0; j < i; ++j) {
        omega[i * n_attitudes + j] *= ratio;
        omega[j * n_attitudes + i] = omega[i * n_attitudes + j];
      }
    }
  }
  return penalty;
}

}  // namespace

Unpacked unpack(const ModelSpec& spec, const std::vector<double>& theta_full,
                const std::vector<int>& grad_of_full) {
  const Layout layout = make_layout(spec);
  if (static_cast<int>(theta_full.size()) != layout.n_full) {
    throw std::invalid_argument("theta has length " + std::to_string(theta_full.size()) +
                                " but this specification has " + std::to_string(layout.n_full) +
                                " coefficients.");
  }
  if (!grad_of_full.empty() && static_cast<int>(grad_of_full.size()) != layout.n_full) {
    throw std::invalid_argument("The gradient index has a different length from theta.");
  }
  auto grad_at = [&](int index) {
    if (grad_of_full.empty()) return -1;
    return grad_of_full[index];
  };

  Unpacked out;
  out.n_attitudes = spec.n_attitudes;
  out.n_waves = spec.n_waves;
  out.n_covariates = spec.n_covariates;
  out.n_factors = spec.n_factors;
  out.wave_effect.assign(static_cast<std::size_t>(spec.n_attitudes * spec.n_waves), 0.0);
  out.beta.assign(static_cast<std::size_t>(spec.n_attitudes * spec.n_covariates), 0.0);
  out.rho.assign(static_cast<std::size_t>(spec.n_attitudes * spec.n_attitudes), 0.0);
  out.drho_dc.assign(out.rho.size(), 0.0);
  out.loading.assign(static_cast<std::size_t>(spec.n_factors * spec.n_attitudes), 0.0);
  out.wave_grad.assign(out.wave_effect.size(), -1);
  out.beta_grad.assign(out.beta.size(), -1);
  out.rho_grad.assign(out.rho.size(), -1);
  out.loading_grad.assign(out.loading.size(), -1);
  out.threshold_increments.resize(spec.n_attitudes);

  int max_categories = 2;
  for (int k : spec.n_categories) max_categories = std::max(max_categories, k);
  out.max_cuts = max_categories + 1;
  out.threshold.assign(static_cast<std::size_t>(spec.n_attitudes * out.max_cuts), 0.0);

  int n_grad = 0;
  for (int g : grad_of_full) n_grad = std::max(n_grad, g + 1);
  out.n_grad = n_grad;

  int index = 0;
  if (spec.include_wave_effects) {
    for (int attitude = 0; attitude < spec.n_attitudes; ++attitude) {
      for (int wave = 0; wave < spec.n_waves; ++wave) {
        const int slot = attitude * spec.n_waves + wave;
        out.wave_effect[slot] = theta_full[index];
        out.wave_grad[slot] = grad_at(index);
        ++index;
      }
    }
  }
  for (int attitude = 0; attitude < spec.n_attitudes; ++attitude) {
    for (int column = 0; column < spec.n_covariates; ++column) {
      const int slot = attitude * spec.n_covariates + column;
      out.beta[slot] = theta_full[index];
      out.beta_grad[slot] = grad_at(index);
      ++index;
    }
  }
  for (int to = 0; to < spec.n_attitudes; ++to) {
    for (int from = 0; from < spec.n_attitudes; ++from) {
      const int slot = to * spec.n_attitudes + from;
      const double unconstrained = theta_full[index];
      if (!std::isfinite(unconstrained)) {
        throw std::invalid_argument("A persistence coefficient is not finite.");
      }
      out.rho[slot] = persistence_from_unconstrained(unconstrained);
      out.drho_dc[slot] = persistence_unconstrained_derivative(out.rho[slot]);
      out.rho_grad[slot] = grad_at(index);
      ++index;
    }
  }
  for (int factor = 0; factor < spec.n_factors; ++factor) {
    for (int attitude = 0; attitude < spec.n_attitudes; ++attitude) {
      const int slot = factor * spec.n_attitudes + attitude;
      out.loading[slot] = theta_full[index];
      out.loading_grad[slot] = grad_at(index);
      ++index;
    }
  }

  const int n_chol = spec.n_attitudes * (spec.n_attitudes - 1) / 2;
  std::vector<double> chol(static_cast<std::size_t>(n_chol), 0.0);
  bool any_chol = false;
  for (int entry = 0; entry < n_chol; ++entry) {
    chol[entry] = theta_full[index];
    if (chol[entry] != 0.0) any_chol = true;
    if (grad_at(index) >= 0) {
      out.chol_position.push_back(index);
      out.chol_grad.push_back(grad_at(index));
    }
    ++index;
  }
  out.identity_correlation = !any_chol;
  out.omega.assign(static_cast<std::size_t>(spec.n_attitudes * spec.n_attitudes), 0.0);
  for (int attitude = 0; attitude < spec.n_attitudes; ++attitude) {
    out.omega[attitude * spec.n_attitudes + attitude] = 1.0;
  }
  if (!out.identity_correlation) {
    build_correlation(spec.n_attitudes, chol, out.omega, out.penalty);
  }

  for (int attitude = 0; attitude < spec.n_attitudes; ++attitude) {
    const int K = spec.n_categories[attitude];
    const int base = attitude * out.max_cuts;
    out.threshold[base + 0] = -100.0;
    out.threshold[base + 1] = 0.0;
    double accumulated = 0.0;
    const int n_increments = K - 2;
    for (int increment = 0; increment < n_increments; ++increment) {
      const double unconstrained = theta_full[index];
      if (!std::isfinite(unconstrained) || unconstrained > 50.0) {
        throw std::invalid_argument(
            "A threshold increment is not finite, or exp(threshold) overflowed. "
            "Threshold coefficients are the log gaps between interior cuts.");
      }
      const double gap = std::exp(unconstrained);
      accumulated += gap;
      out.threshold[base + increment + 2] = accumulated;
      ThresholdIncrement piece;
      piece.grad = grad_at(index);
      piece.exp_c = gap;
      piece.first_cut = increment + 2;
      out.threshold_increments[attitude].push_back(piece);
      ++index;
    }
    out.threshold[base + K] = 100.0;
  }
  if (index != layout.n_full) {
    throw std::logic_error("unpack did not consume every coefficient.");
  }
  return out;
}

}  // namespace detail
}  // namespace attitudedyn
