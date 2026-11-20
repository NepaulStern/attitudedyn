#include "attitudedyn/api.hpp"
#include "attitudedyn/normal.hpp"
#include "draws.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int g_failures = 0;

void check(bool condition, const std::string& message) {
  if (!condition) {
    std::cerr << "FAIL " << message << "\n";
    ++g_failures;
  }
}

void check_near(double actual, double expected, double tol, const std::string& message) {
  if (!std::isfinite(actual) || std::fabs(actual - expected) > tol) {
    std::cerr << "FAIL " << message << " actual=" << actual << " expected=" << expected << "\n";
    ++g_failures;
  }
}

attitudedyn::Person make_person(int n_waves, int /*n_attitudes*/, int /*n_covariates*/,
                                const std::vector<int>& category,
                                const std::vector<double>& covariates) {
  attitudedyn::Person person;
  person.first_wave = 0;
  person.last_wave = n_waves - 1;
  person.interviewed.assign(n_waves, 1);
  person.category = category;
  person.covariates = covariates;
  person.weight = 1.0;
  return person;
}

}  // namespace

int main() {
  using namespace attitudedyn;
  check_near(standard_normal_cdf(0.0), 0.5, 1e-15, "Phi(0)");
  check_near(standard_normal_quantile(0.5), 0.0, 1e-12, "qnorm(0.5)");
  check_near(standard_normal_quantile(0.975), 1.959963984540054, 1e-12, "qnorm(0.975)");
  check_near(standard_normal_cdf(standard_normal_quantile(1e-6)), 1e-6, 1e-12, "cdf(qnorm) tail");
  check_near(std::exp(log_normal_interval(-0.5, 1.25)),
             standard_normal_cdf(1.25) - standard_normal_cdf(-0.5), 1e-12, "log interval");

  check_near(detail::halton_point(1, 0), 0.5, 1e-15, "halton base2 #1");
  check_near(detail::halton_point(2, 0), 0.25, 1e-15, "halton base2 #2");
  check_near(detail::halton_point(3, 0), 0.75, 1e-15, "halton base2 #3");
  check_near(detail::halton_point(1, 1), 1.0 / 3.0, 1e-15, "halton base3 #1");

  // Joe and Kuo gray-code Sobol points, skipping the origin at index 0.
  const double sobol_ref[7][3] = {
      {0.5, 0.5, 0.5},     {0.75, 0.25, 0.25}, {0.25, 0.75, 0.75}, {0.375, 0.375, 0.625},
      {0.875, 0.875, 0.125}, {0.625, 0.125, 0.875}, {0.125, 0.625, 0.375}};
  for (int i = 0; i < 7; ++i) {
    for (int d = 0; d < 3; ++d) {
      check_near(detail::sobol_point(d, static_cast<std::uint32_t>(i + 1)), sobol_ref[i][d], 1e-12,
                 "sobol point");
    }
  }

  const double z1 = detail::standard_normal_draw(11, 3, 4, 1, false, DrawKind::kPhilox);
  const double z2 = detail::standard_normal_draw(11, 3, 4, 1, false, DrawKind::kPhilox);
  const double z3 = detail::standard_normal_draw(11, 3, 4, 2, false, DrawKind::kPhilox);
  check(z1 == z2, "Philox is deterministic");
  check(z1 != z3, "Philox streams differ");
  check(detail::standard_normal_draw(11, 3, 4, 1, true, DrawKind::kPhilox) == -z1,
        "antithetic negation");

  // Hand likelihood. One binary item, y* = x * beta, category 2, so the
  // interval is (0, 100] and beta = x = 1. The probability is Phi(1).
  ModelSpec binary;
  binary.n_waves = 1;
  binary.n_attitudes = 1;
  binary.n_covariates = 1;
  binary.n_factors = 0;
  binary.n_categories = {2};
  binary.wave_time = {0.0};
  binary.include_wave_effects = false;
  binary.fix_intercepts = false;
  binary.attitude_names = {"item"};
  binary.covariate_names = {"x"};
  Problem binary_problem;
  binary_problem.spec = binary;
  binary_problem.layout = make_layout(binary);
  check(binary_problem.layout.n_full == 2, "binary layout length");
  binary_problem.theta_full = {1.0, 0.0};
  binary_problem.draws.n_draws = 2;
  binary_problem.draws.antithetic = true;
  binary_problem.draws.seed = 7;
  binary_problem.persons.push_back(make_person(1, 1, 1, {2}, {1.0}));
  EvaluateRequest binary_request;
  binary_request.grad_columns = {0, 1};
  const Evaluation binary_eval = evaluate_full(binary_problem, binary_problem.theta_full, binary_request);
  const double expected = std::log(standard_normal_cdf(1.0));
  check_near(binary_eval.objective, expected, 1e-10, "binary loglik");
  const double expected_slope = standard_normal_pdf(1.0) / standard_normal_cdf(1.0);
  check_near(binary_eval.gradient[0], expected_slope, 1e-8, "binary beta gradient");
  check_near(binary_eval.gradient[1], 0.0, 1e-10, "binary persistence gradient is zero in one wave");

  // GSS coefficient count: 6 attitudes, 8 waves, 10 covariates, 2 factors,
  // categories 7,3,3,5,7,3. One empty wave effect, six fixed intercepts,
  // the last loading of the second factor, and the 15 correlation parameters.
  ModelSpec gss;
  gss.n_waves = 8;
  gss.n_attitudes = 6;
  gss.n_covariates = 10;
  gss.n_factors = 2;
  gss.n_categories = {7, 3, 3, 5, 7, 3};
  gss.wave_time = {2, 3, 4, 5, 6, 7, 8, 10};
  gss.include_wave_effects = true;
  gss.fix_intercepts = true;
  gss.fix_last_loading = true;
  gss.free_contemporaneous = false;
  gss.intercept_columns = {0};
  gss.empty_wave_effect.assign(6 * 8, 0);
  gss.empty_wave_effect[3 * 8 + 0] = 1;
  gss.attitude_names = {"ideology", "military", "welfare", "immigrants", "party", "environment"};
  const Layout gss_layout = make_layout(gss);
  check(gss_layout.n_full == 187, "GSS coefficient count");
  check(gss_layout.n_free == 164, "GSS free coefficient count");
  check(gss_layout.entries[3 * 8].name == "wave_effect[immigrants,1]", "empty wave name");
  check(!gss_layout.entries[3 * 8].is_free, "immigrants wave 1 is fixed");

  // Finite-difference check on a small panel, every free coefficient.
  ModelSpec small;
  small.n_waves = 3;
  small.n_attitudes = 2;
  small.n_covariates = 2;
  small.n_factors = 1;
  small.n_categories = {3, 3};
  small.wave_time = {0.0, 1.0, 2.0};
  small.include_wave_effects = true;
  small.fix_intercepts = true;
  small.intercept_columns = {0};
  small.attitude_names = {"a", "b"};
  small.covariate_names = {"(Intercept)", "x"};
  small.wave_names = {"w1", "w2", "w3"};
  Problem small_problem;
  small_problem.spec = small;
  small_problem.layout = make_layout(small);
  small_problem.theta_full.assign(small_problem.layout.n_full, 0.0);
  for (int i = 0; i < small_problem.layout.n_full; ++i) {
    const std::string& block = small_problem.layout.entries[i].block;
    if (!small_problem.layout.entries[i].is_free) continue;
    if (block == "persistence" && small_problem.layout.entries[i].attitude ==
                                       small_problem.layout.entries[i].index2) {
      small_problem.theta_full[i] = 0.4;  // unconstrained; natural rho is about 0.197
    } else if (block == "loading") {
      small_problem.theta_full[i] = 0.3;
    } else if (block == "threshold") {
      small_problem.theta_full[i] = std::log(0.8);
    } else if (block == "beta") {
      small_problem.theta_full[i] = 0.25;
    } else if (block == "wave_effect") {
      small_problem.theta_full[i] = 0.1 * small_problem.layout.entries[i].index2;
    }
  }
  small_problem.draws.n_draws = 3;
  small_problem.draws.antithetic = true;
  small_problem.draws.seed = 19;
  for (int person = 0; person < 6; ++person) {
    std::vector<int> category(2 * 3);
    std::vector<double> x(3 * 2);
    for (int wave = 0; wave < 3; ++wave) {
      x[wave * 2 + 0] = 1.0;
      x[wave * 2 + 1] = ((person + wave) % 5) * 0.3 - 0.6;
      for (int attitude = 0; attitude < 2; ++attitude) {
        category[attitude * 3 + wave] = 1 + ((person + wave + attitude) % 3);
      }
    }
    Person panel_person = make_person(3, 2, 2, category, x);
    if (person == 5) panel_person.interviewed[1] = 0;  // a gap
    small_problem.persons.push_back(panel_person);
  }
  std::vector<int> free_columns;
  for (int i = 0; i < small_problem.layout.n_full; ++i) {
    if (small_problem.layout.entries[i].is_free) free_columns.push_back(i);
  }
  EvaluateRequest gradient_request;
  gradient_request.grad_columns = free_columns;
  gradient_request.n_threads = 1;
  const Evaluation analytic =
      evaluate_full(small_problem, small_problem.theta_full, gradient_request);
  check(std::isfinite(analytic.objective), "small panel loglik is finite");
  const double step = 1e-6;
  for (int column = 0; column < static_cast<int>(free_columns.size()); ++column) {
    std::vector<double> up = small_problem.theta_full;
    std::vector<double> down = small_problem.theta_full;
    up[free_columns[column]] += step;
    down[free_columns[column]] -= step;
    EvaluateRequest value_only;
    value_only.n_threads = 1;
    const double plus = evaluate_full(small_problem, up, value_only).objective;
    const double minus = evaluate_full(small_problem, down, value_only).objective;
    const double numeric = (plus - minus) / (2.0 * step);
    const double scale = std::max(1.0, std::fabs(analytic.gradient[column]));
    check(std::fabs(analytic.gradient[column] - numeric) <= 2e-5 * scale,
          "finite difference " + small_problem.layout.entries[free_columns[column]].name +
              " analytic=" + std::to_string(analytic.gradient[column]) +
              " numeric=" + std::to_string(numeric));
  }

  EvaluateRequest threaded = gradient_request;
  threaded.n_threads = 3;
  const Evaluation threaded_eval = evaluate_full(small_problem, small_problem.theta_full, threaded);
  check(threaded_eval.objective == analytic.objective, "thread count does not change the loglik");
  check(threaded_eval.gradient == analytic.gradient, "thread count does not change the gradient");

  std::vector<double> extreme_theta = small_problem.theta_full;
  bool saw_threshold = false;
  for (int i = 0; i < small_problem.layout.n_full; ++i) {
    if (small_problem.layout.entries[i].block == "threshold") {
      extreme_theta[i] = 80.0;
      saw_threshold = true;
    }
  }
  check(saw_threshold, "small panel has a threshold coefficient");
  const Evaluation extreme_eval = evaluate_full(small_problem, extreme_theta, gradient_request);
  check(!std::isfinite(extreme_eval.objective), "overflowing threshold objective is not finite");
  check(extreme_eval.gradient.size() == free_columns.size(), "overflowing threshold keeps the gradient length");
  bool extreme_gradient_finite = true;
  for (double component : extreme_eval.gradient) {
    if (!std::isfinite(component)) extreme_gradient_finite = false;
  }
  check(extreme_gradient_finite, "overflowing threshold gradient is finite");
  for (int i = 0; i < small_problem.layout.n_full; ++i) {
    if (small_problem.layout.entries[i].block == "persistence") {
      extreme_theta[i] = std::numeric_limits<double>::infinity();
    }
  }
  EvaluateRequest extreme_value;
  extreme_value.n_threads = 1;
  const Evaluation infinite_eval = evaluate_full(small_problem, extreme_theta, extreme_value);
  check(!std::isfinite(infinite_eval.objective), "non-finite persistence is a non-finite objective");

  // One attitude, times 0, 1, 3, persistence 0.5, wave effects 0, 1, 3.
  // The single person contributes waves 2 and 3 to the moments (values 1 and
  // 3, variance 1). The shock contribution is 1.25: the product that uses
  // the never-written first-wave transition is zero, and the one-step
  // product at the second wave is 0.5^2.
  ModelSpec decomposition_spec;
  decomposition_spec.n_waves = 3;
  decomposition_spec.n_attitudes = 1;
  decomposition_spec.n_covariates = 1;
  decomposition_spec.n_factors = 0;
  decomposition_spec.n_categories = {2};
  decomposition_spec.wave_time = {0.0, 1.0, 3.0};
  decomposition_spec.include_wave_effects = true;
  decomposition_spec.fix_intercepts = true;
  decomposition_spec.intercept_columns = {0};
  decomposition_spec.attitude_names = {"item"};
  decomposition_spec.covariate_names = {"(Intercept)"};
  Problem decomposition;
  decomposition.spec = decomposition_spec;
  decomposition.layout = make_layout(decomposition_spec);
  check(decomposition.layout.n_full == 5, "variance-decomposition layout");
  decomposition.theta_full = {0.0, 1.0, 3.0, 0.0, std::log(3.0)};
  Person decomposition_person = make_person(3, 1, 1, {1, 1, 1}, {1.0, 1.0, 1.0});
  decomposition_person.interviewed[1] = 0;
  decomposition.persons.push_back(decomposition_person);
  const VarianceComponent one_person = variance_decomposition(decomposition, decomposition.theta_full);
  check_near(one_person.covariate[0], 1.0, 1e-12, "x-beta variance includes the wave effect");
  check_near(one_person.wave_effect[0], 1.0, 1e-12, "wave-effect variance");
  check_near(one_person.factor[0], 0.0, 1e-12, "no factor");
  check_near(one_person.shock[0], 1.25, 1e-12, "finite-horizon shock contribution");
  check_near(one_person.total[0], 3.25, 1e-12, "variance total");

  // A later-starting person who is stored first leaves the early transition
  // at zero, so that person's shock contribution stays 1. Reversing the
  // order fills the transition before it is read.
  Person later_person = make_person(3, 1, 1, {1, 1, 1}, {1.0, 1.0, 1.0});
  later_person.first_wave = 1;
  Problem later_first = decomposition;
  later_first.persons = {later_person, decomposition_person};
  const VarianceComponent later_first_parts = variance_decomposition(later_first, later_first.theta_full);
  check_near(later_first_parts.shock[0], 1.125, 1e-12, "early transition stays zero until visited");
  check_near(later_first_parts.covariate[0], 8.0 / 9.0, 1e-12, "pooled x-beta variance");
  Problem early_first = decomposition;
  early_first.persons = {decomposition_person, later_person};
  const VarianceComponent early_first_parts = variance_decomposition(early_first, early_first.theta_full);
  check_near(early_first_parts.shock[0], 1.25, 1e-12, "a full window fills the transition for the next person");

  const ImpulseResponse irf = impulse_response(
      natural_parameters(small_problem.spec, small_problem.theta_full).persistence, 2, 4);
  check(irf.response.size() == 5u * 4u, "irf length");
  check_near(irf.response[0], 1.0, 1e-12, "irf horizon 0 is the identity");
  check(std::isfinite(irf.spectral_radius), "spectral radius is finite");

  const double again = detail::standard_normal_draw(11, 3, 4, 1, false, DrawKind::kPhilox);
  check(again == z1, "draws stay fixed after other evaluations");

  if (g_failures == 0) {
    std::cout << "All core tests passed.\n";
    return 0;
  }
  std::cerr << g_failures << " core test(s) failed.\n";
  return 1;
}
