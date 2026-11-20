#include "internal.hpp"

#include "draws.hpp"
#include "attitudedyn/normal.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace attitudedyn {
namespace {

int max_categories(const ModelSpec& spec) {
  int m = 2;
  for (int k : spec.n_categories) m = std::max(m, k);
  return m;
}

double infinity_norm(const std::vector<double>& matrix, int n) {
  double norm = 0.0;
  for (int row = 0; row < n; ++row) {
    double sum = 0.0;
    for (int col = 0; col < n; ++col) sum += std::fabs(matrix[row * n + col]);
    norm = std::max(norm, sum);
  }
  return norm;
}

std::vector<double> multiply(const std::vector<double>& left, const std::vector<double>& right, int n) {
  std::vector<double> out(static_cast<std::size_t>(n * n), 0.0);
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < n; ++k) {
      const double lik = left[i * n + k];
      if (lik == 0.0) continue;
      for (int j = 0; j < n; ++j) out[i * n + j] += lik * right[k * n + j];
    }
  }
  return out;
}

}  // namespace

NaturalParameters natural_parameters(const ModelSpec& spec, const std::vector<double>& theta_full) {
  const std::vector<int> no_gradient(theta_full.size(), -1);
  const detail::Unpacked unpacked = detail::unpack(spec, theta_full, no_gradient);
  NaturalParameters out;
  out.persistence = unpacked.rho;
  out.loadings = unpacked.loading;
  out.wave_effects = unpacked.wave_effect;
  out.beta = unpacked.beta;
  out.cutpoints.resize(spec.n_attitudes);
  for (int attitude = 0; attitude < spec.n_attitudes; ++attitude) {
    const int K = spec.n_categories[attitude];
    out.cutpoints[attitude].resize(static_cast<std::size_t>(K + 1));
    for (int cut = 0; cut <= K; ++cut) {
      out.cutpoints[attitude][cut] = unpacked.threshold[attitude * unpacked.max_cuts + cut];
    }
  }
  return out;
}

ImpulseResponse impulse_response(const std::vector<double>& persistence_row_major, int n_attitudes,
                                 int horizon) {
  if (n_attitudes < 1 || horizon < 0) {
    throw std::invalid_argument("The impulse response needs a positive dimension and a non-negative horizon.");
  }
  if (static_cast<int>(persistence_row_major.size()) != n_attitudes * n_attitudes) {
    throw std::invalid_argument("The persistence matrix is not square.");
  }
  ImpulseResponse out;
  out.horizon = horizon;
  out.n_attitudes = n_attitudes;
  out.response.assign(static_cast<std::size_t>(horizon + 1) * n_attitudes * n_attitudes, 0.0);
  std::vector<double> power(static_cast<std::size_t>(n_attitudes * n_attitudes), 0.0);
  for (int i = 0; i < n_attitudes; ++i) power[i * n_attitudes + i] = 1.0;
  for (int step = 0; step <= horizon; ++step) {
    std::copy(power.begin(), power.end(),
              out.response.begin() + static_cast<std::size_t>(step) * n_attitudes * n_attitudes);
    power = multiply(persistence_row_major, power, n_attitudes);
  }
  // Gelfand's formula: the spectral radius is the limit of ||R^n||^{1/n}.
  std::vector<double> grown = persistence_row_major;
  double radius = infinity_norm(grown, n_attitudes);
  for (int n = 1; n < 32; ++n) {
    grown = multiply(persistence_row_major, grown, n_attitudes);
    const double norm = infinity_norm(grown, n_attitudes);
    if (norm == 0.0) {
      radius = 0.0;
      break;
    }
    radius = std::pow(norm, 1.0 / static_cast<double>(n + 1));
  }
  out.spectral_radius = radius;
  return out;
}

InfluenceShares influence_shares(const std::vector<double>& persistence_row_major, int n_attitudes) {
  InfluenceShares out;
  if (n_attitudes < 1 ||
      static_cast<int>(persistence_row_major.size()) != n_attitudes * n_attitudes) {
    throw std::invalid_argument("The persistence matrix is not square.");
  }
  const ImpulseResponse probe = impulse_response(persistence_row_major, n_attitudes, 0);
  out.stationary = probe.spectral_radius < 0.999;
  out.long_run_variance.assign(static_cast<std::size_t>(n_attitudes), 0.0);
  out.share.assign(static_cast<std::size_t>(n_attitudes * n_attitudes), 0.0);
  if (!out.stationary) return out;

  // V = sum_{h=0}^{H} R^h (R^h)'. The (to, from) summand at horizon h is
  // (R^h)[to, from]^2, which is the variance contributed by a unit innovation
  // in attitude `from`. H stops when the added mass is negligible.
  std::vector<double> power(static_cast<std::size_t>(n_attitudes * n_attitudes), 0.0);
  for (int i = 0; i < n_attitudes; ++i) power[i * n_attitudes + i] = 1.0;
  std::vector<double> contribution(out.share.size(), 0.0);
  for (int step = 0; step < 400; ++step) {
    double added = 0.0;
    for (int to = 0; to < n_attitudes; ++to) {
      for (int from = 0; from < n_attitudes; ++from) {
        const double coef = power[to * n_attitudes + from];
        const double piece = coef * coef;
        contribution[to * n_attitudes + from] += piece;
        added += piece;
      }
    }
    if (step > 5 && added < 1e-14) break;
    power = multiply(persistence_row_major, power, n_attitudes);
  }
  for (int to = 0; to < n_attitudes; ++to) {
    double total = 0.0;
    for (int from = 0; from < n_attitudes; ++from) total += contribution[to * n_attitudes + from];
    out.long_run_variance[to] = total;
    if (total > 0.0) {
      for (int from = 0; from < n_attitudes; ++from) {
        out.share[to * n_attitudes + from] = contribution[to * n_attitudes + from] / total;
      }
    }
  }
  return out;
}

namespace {

double sample_variance(const std::vector<double>& values) {
  if (values.size() < 2) return 0.0;
  double mean = 0.0;
  for (double v : values) mean += v;
  mean /= static_cast<double>(values.size());
  double sse = 0.0;
  for (double v : values) {
    const double d = v - mean;
    sse += d * d;
  }
  return sse / static_cast<double>(values.size());
}

}  // namespace

VarianceComponent variance_decomposition(const Problem& problem, const std::vector<double>& theta_full,
                                         Dynamics dynamics) {
  const std::vector<int> no_gradient(theta_full.size(), -1);
  const detail::Unpacked unpacked = detail::unpack(problem.spec, theta_full, no_gradient);
  const int J = problem.spec.n_attitudes;
  VarianceComponent out;
  out.covariate.assign(J, 0.0);
  out.wave_effect.assign(J, 0.0);
  out.factor.assign(J, 0.0);
  out.shock.assign(J, std::numeric_limits<double>::quiet_NaN());
  out.total.assign(J, 0.0);
  std::vector<std::vector<double>> xb(J), delta(J);
  for (const Person& person : problem.persons) {
    if (person.last_wave < person.first_wave) continue;
    for (int wave = person.first_wave; wave <= person.last_wave; ++wave) {
      if (!person.interviewed[wave]) continue;
      for (int attitude = 0; attitude < J; ++attitude) {
        double systematic = 0.0;
        const double* x = person.covariates.data() + static_cast<std::size_t>(wave) * unpacked.n_covariates;
        const double* beta = unpacked.beta.data() + static_cast<std::size_t>(attitude) * unpacked.n_covariates;
        for (int column = 0; column < unpacked.n_covariates; ++column) systematic += beta[column] * x[column];
        xb[attitude].push_back(systematic);
        double wave_effect = 0.0;
        if (problem.spec.include_wave_effects) {
          wave_effect = unpacked.wave_effect[attitude * unpacked.n_waves + wave];
        }
        delta[attitude].push_back(wave_effect);
      }
    }
  }
  for (int attitude = 0; attitude < J; ++attitude) {
    out.covariate[attitude] = sample_variance(xb[attitude]);
    out.wave_effect[attitude] = sample_variance(delta[attitude]);
    double factor = 0.0;
    for (int f = 0; f < unpacked.n_factors; ++f) {
      const double loading = unpacked.loading[f * J + attitude];
      factor += loading * loading;
    }
    out.factor[attitude] = factor;
  }
  const InfluenceShares shares = influence_shares(unpacked.rho, J);
  if (!shares.stationary) {
    out.note = (dynamics == Dynamics::kLatentAttitude)
                   ? "The persistence matrix is not stationary, so the long-run shock variance is not finite."
                   : "The shock persistence matrix is not stationary, so the long-run shock variance is not finite.";
  } else {
    out.shock = shares.long_run_variance;
    out.note = "Shock variance is the long-run variance of the recursion state = R * state_lag + innovation, "
               "with innovation covariance equal to the identity.";
  }
  for (int attitude = 0; attitude < J; ++attitude) {
    const double shock = std::isfinite(out.shock[attitude]) ? out.shock[attitude] : 0.0;
    out.total[attitude] = out.covariate[attitude] + out.wave_effect[attitude] + out.factor[attitude] + shock;
  }
  return out;
}

VarianceComponent variance_decomposition_fortran(const Problem& problem,
                                                  const std::vector<double>& theta_full) {
  // The Fortran vardecomp routine double-counts the wave effect and accumulates
  // products of single persistence entries rather than matrix powers. A
  // line-by-line port is still worth having for the GSS log, but only after it
  // has been checked against that log. Until then this entry point refuses, so
  // a caller cannot mistake an unaudited number for the published table.
  (void)problem;
  (void)theta_full;
  throw std::logic_error(
      "variance_decomposition_fortran is not enabled yet. Use variance_decomposition(), "
      "which splits covariate, wave, factor and long-run shock variance without "
      "the Fortran double count.");
}

SimulatedPanel simulate_panel(const ModelSpec& spec, const std::vector<double>& theta_full,
                              int n_persons, Dynamics dynamics, std::uint64_t seed) {
  if (n_persons < 1) throw std::invalid_argument("Simulate at least one person.");
  const std::vector<int> no_gradient(theta_full.size(), -1);
  const detail::Unpacked unpacked = detail::unpack(spec, theta_full, no_gradient);
  SimulatedPanel panel;
  panel.n_persons = n_persons;
  panel.n_waves = spec.n_waves;
  panel.n_attitudes = spec.n_attitudes;
  panel.n_covariates = spec.n_covariates;
  const int rows = n_persons * spec.n_waves;
  panel.id.resize(rows);
  panel.wave.resize(rows);
  panel.time.resize(rows);
  panel.covariates.resize(static_cast<std::size_t>(rows) * spec.n_covariates);
  panel.outcome.resize(static_cast<std::size_t>(rows) * spec.n_attitudes);

  const bool latent = dynamics == Dynamics::kLatentAttitude;
  const int J = spec.n_attitudes;
  for (int person = 0; person < n_persons; ++person) {
    std::vector<double> factors(static_cast<std::size_t>(spec.n_factors));
    for (int factor = 0; factor < spec.n_factors; ++factor) {
      factors[factor] = detail::standard_normal_draw(seed, person, 0, factor, false, DrawKind::kPhilox);
    }
    std::vector<double> state(static_cast<std::size_t>(J), 0.0);
    for (int wave = 0; wave < spec.n_waves; ++wave) {
      const int row = person * spec.n_waves + wave;
      panel.id[row] = person + 1;
      panel.wave[row] = wave + 1;
      panel.time[row] = spec.wave_time[wave];
      for (int column = 0; column < spec.n_covariates; ++column) {
        double value = 1.0;
        bool intercept = false;
        for (int candidate : spec.intercept_columns) {
          if (candidate == column) intercept = true;
        }
        if (!intercept) {
          const int stream = spec.n_factors + spec.n_attitudes * spec.n_waves + column;
          value = detail::standard_normal_draw(seed, person, wave, stream, false, DrawKind::kPhilox);
        }
        panel.covariates[static_cast<std::size_t>(row) * spec.n_covariates + column] = value;
      }
      std::vector<double> next(static_cast<std::size_t>(J), 0.0);
      const double gap = (wave == 0) ? 0.0 : spec.wave_time[wave] - spec.wave_time[wave - 1];
      for (int attitude = 0; attitude < J; ++attitude) {
        double lagged = 0.0;
        if (wave > 0) {
          for (int source = 0; source < J; ++source) {
            const double rho = unpacked.rho[attitude * J + source];
            const double rhop = (rho == 0.0) ? 0.0 : std::copysign(std::pow(std::fabs(rho), gap), rho);
            lagged += rhop * state[source];
          }
        }
        double xb = 0.0;
        if (spec.include_wave_effects) xb += unpacked.wave_effect[attitude * spec.n_waves + wave];
        for (int column = 0; column < spec.n_covariates; ++column) {
          xb += unpacked.beta[attitude * spec.n_covariates + column] *
                panel.covariates[static_cast<std::size_t>(row) * spec.n_covariates + column];
        }
        double factor_mean = 0.0;
        for (int factor = 0; factor < spec.n_factors; ++factor) {
          factor_mean += unpacked.loading[factor * J + attitude] * factors[factor];
        }
        const int stream = spec.n_factors + attitude * spec.n_waves + wave;
        const double innovation =
            detail::standard_normal_draw(seed ^ 0xA5A5A5A5A5A5A5A5ull, person, 1, stream, false,
                                         DrawKind::kPhilox);
        const double ystar = xb + factor_mean + lagged + innovation;
        next[attitude] = latent ? ystar : lagged + innovation;
        int category = spec.n_categories[attitude];
        for (int k = 1; k <= spec.n_categories[attitude]; ++k) {
          if (ystar <= unpacked.threshold[attitude * unpacked.max_cuts + k]) {
            category = k;
            break;
          }
        }
        panel.outcome[static_cast<std::size_t>(row) * J + attitude] = category;
      }
      state.swap(next);
    }
  }
  return panel;
}

Prediction predict_probabilities(const Problem& problem, const std::vector<double>& theta_full,
                                 Dynamics dynamics, int n_draws, std::uint64_t seed, int n_threads) {
  if (n_draws < 1) throw std::invalid_argument("Prediction needs at least one draw.");
  if (n_threads < 1) throw std::invalid_argument("n_threads must be at least 1.");
  (void)n_threads;  // the likelihood is the threaded hot path; prediction is a post-step
  const std::vector<int> no_gradient(theta_full.size(), -1);
  const detail::Unpacked unpacked = detail::unpack(problem.spec, theta_full, no_gradient);
  if (!unpacked.identity_correlation) {
    Prediction failed;
    failed.error = "Category probabilities with a free contemporaneous correlation are not implemented.";
    return failed;
  }
  const int J = problem.spec.n_attitudes;
  const int T = problem.spec.n_waves;
  const int n_persons = static_cast<int>(problem.persons.size());
  const int Kmax = max_categories(problem.spec);
  Prediction out;
  out.max_categories = Kmax;
  out.probability.assign(static_cast<std::size_t>(n_persons) * T * J * Kmax, 0.0);
  out.mean_latent.assign(static_cast<std::size_t>(n_persons) * T * J, 0.0);
  const bool latent = dynamics == Dynamics::kLatentAttitude;
  const double inv = 1.0 / static_cast<double>(n_draws);

  for (int person_index = 0; person_index < n_persons; ++person_index) {
    const Person& person = problem.persons[person_index];
    if (person.last_wave < person.first_wave) continue;
    for (int draw = 0; draw < n_draws; ++draw) {
      std::vector<double> factors(static_cast<std::size_t>(problem.spec.n_factors));
      for (int factor = 0; factor < problem.spec.n_factors; ++factor) {
        factors[factor] =
            detail::standard_normal_draw(seed, person_index, draw, factor, false, problem.draws.kind);
      }
      std::vector<double> state(static_cast<std::size_t>(J), 0.0);
      for (int wave = person.first_wave; wave <= person.last_wave; ++wave) {
        const double gap = (wave == person.first_wave)
                               ? 0.0
                               : problem.spec.wave_time[wave] - problem.spec.wave_time[wave - 1];
        std::vector<double> next(static_cast<std::size_t>(J), 0.0);
        for (int attitude = 0; attitude < J; ++attitude) {
          double lagged = 0.0;
          if (wave > person.first_wave) {
            for (int source = 0; source < J; ++source) {
              const double rho = unpacked.rho[attitude * J + source];
              const double rhop = (rho == 0.0) ? 0.0 : std::copysign(std::pow(std::fabs(rho), gap), rho);
              lagged += rhop * state[source];
            }
          }
          double xb = 0.0;
          if (problem.spec.include_wave_effects) {
            xb += unpacked.wave_effect[attitude * T + wave];
          }
          if (person.interviewed[wave]) {
            const double* x = person.covariates.data() + static_cast<std::size_t>(wave) * unpacked.n_covariates;
            for (int column = 0; column < unpacked.n_covariates; ++column) {
              xb += unpacked.beta[attitude * unpacked.n_covariates + column] * x[column];
            }
          }
          double factor_mean = 0.0;
          for (int factor = 0; factor < problem.spec.n_factors; ++factor) {
            factor_mean += unpacked.loading[factor * J + attitude] * factors[factor];
          }
          const double y_pre = xb + factor_mean + lagged;
          const int stream = problem.spec.n_factors + attitude * T + wave;
          const double innovation =
              detail::standard_normal_draw(seed, person_index, draw, stream, false, problem.draws.kind);
          next[attitude] = (latent ? y_pre : lagged) + innovation;
          const std::size_t latent_index =
              (static_cast<std::size_t>(person_index) * T + wave) * J + attitude;
          out.mean_latent[latent_index] += y_pre * inv;
          const int K = problem.spec.n_categories[attitude];
          for (int category = 1; category <= K; ++category) {
            const double lower = unpacked.threshold[attitude * unpacked.max_cuts + category - 1];
            const double upper = unpacked.threshold[attitude * unpacked.max_cuts + category];
            const double probability = std::exp(log_normal_interval(lower - y_pre, upper - y_pre));
            const std::size_t slot =
                (((static_cast<std::size_t>(person_index) * T + wave) * J + attitude) * Kmax) +
                (category - 1);
            out.probability[slot] += probability * inv;
          }
        }
        state.swap(next);
      }
    }
  }
  return out;
}

}  // namespace attitudedyn
