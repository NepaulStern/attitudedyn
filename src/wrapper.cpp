#include <Rcpp.h>
#include "attitudedyn/api.hpp"

#include <cstdint>
#include <string>
#include <vector>

using namespace Rcpp;

namespace {

attitudedyn::DrawKind draw_kind_from_int(int kind) {
  if (kind == 1) return attitudedyn::DrawKind::kHalton;
  if (kind == 2) return attitudedyn::DrawKind::kSobol;
  return attitudedyn::DrawKind::kPhilox;
}

attitudedyn::Numerics numerics_from_int(int numerics) {
  return numerics == 1 ? attitudedyn::Numerics::kShared : attitudedyn::Numerics::kStable;
}

std::vector<std::string> as_strings(const CharacterVector& values) {
  std::vector<std::string> out;
  out.reserve(values.size());
  for (R_xlen_t i = 0; i < values.size(); ++i) {
    if (CharacterVector::is_na(values[i])) out.emplace_back("");
    else out.emplace_back(CHAR(values[i]));
  }
  return out;
}

attitudedyn::Problem* problem_from(SEXP xp) {
  XPtr<attitudedyn::Problem> pointer(xp);
  if (!pointer) stop("The model pointer is no longer valid. Refit, or call the rebuild helper.");
  return pointer.get();
}

std::uint64_t seed_from_numeric(double seed) {
  if (!R_finite(seed) || seed < 0) stop("`seed` must be a finite, non-negative number.");
  return static_cast<std::uint64_t>(seed);
}

}  // namespace

// [[Rcpp::export]]
List attitudedyn_layout_cpp(int n_waves, int n_attitudes, int n_covariates, int n_factors,
                            IntegerVector n_categories, NumericVector wave_time,
                            bool include_wave_effects, bool free_contemporaneous,
                            bool fix_last_loading, bool fix_intercepts,
                            IntegerVector intercept_columns, IntegerVector empty_wave_effect,
                            CharacterVector attitude_names, CharacterVector covariate_names,
                            CharacterVector wave_names) {
  attitudedyn::ModelSpec spec;
  spec.n_waves = n_waves;
  spec.n_attitudes = n_attitudes;
  spec.n_covariates = n_covariates;
  spec.n_factors = n_factors;
  spec.n_categories.assign(n_categories.begin(), n_categories.end());
  spec.wave_time.assign(wave_time.begin(), wave_time.end());
  spec.include_wave_effects = include_wave_effects;
  spec.free_contemporaneous = free_contemporaneous;
  spec.fix_last_loading = fix_last_loading;
  spec.fix_intercepts = fix_intercepts;
  spec.intercept_columns.assign(intercept_columns.begin(), intercept_columns.end());
  spec.empty_wave_effect.assign(empty_wave_effect.begin(), empty_wave_effect.end());
  spec.attitude_names = as_strings(attitude_names);
  spec.covariate_names = as_strings(covariate_names);
  spec.wave_names = as_strings(wave_names);
  const attitudedyn::Layout layout = attitudedyn::make_layout(spec);
  CharacterVector name(layout.n_full), block(layout.n_full);
  IntegerVector attitude(layout.n_full), index2(layout.n_full), is_free(layout.n_full);
  for (int i = 0; i < layout.n_full; ++i) {
    name[i] = layout.entries[i].name;
    block[i] = layout.entries[i].block;
    attitude[i] = layout.entries[i].attitude + 1;
    index2[i] = layout.entries[i].index2 + 1;
    is_free[i] = layout.entries[i].is_free ? 1 : 0;
  }
  return List::create(Named("n_full") = layout.n_full, Named("n_free") = layout.n_free,
                      Named("name") = name, Named("block") = block, Named("attitude") = attitude,
                      Named("index2") = index2, Named("is_free") = is_free);
}

// [[Rcpp::export]]
SEXP attitudedyn_build_cpp(
    int n_waves, int n_attitudes, int n_covariates, int n_factors, IntegerVector n_categories,
    NumericVector wave_time, bool include_wave_effects, bool free_contemporaneous,
    bool fix_last_loading, bool fix_intercepts, IntegerVector intercept_columns,
    IntegerVector empty_wave_effect, CharacterVector attitude_names, CharacterVector covariate_names,
    CharacterVector wave_names, IntegerVector first_wave, IntegerVector last_wave,
    IntegerVector interviewed, IntegerVector category, NumericVector covariates, NumericVector weight,
    NumericVector theta_full, IntegerVector is_free, int n_draws, bool antithetic, int draw_kind,
    double seed) {
  auto* problem = new attitudedyn::Problem();
  attitudedyn::ModelSpec& spec = problem->spec;
  spec.n_waves = n_waves;
  spec.n_attitudes = n_attitudes;
  spec.n_covariates = n_covariates;
  spec.n_factors = n_factors;
  spec.n_categories.assign(n_categories.begin(), n_categories.end());
  spec.wave_time.assign(wave_time.begin(), wave_time.end());
  spec.include_wave_effects = include_wave_effects;
  spec.free_contemporaneous = free_contemporaneous;
  spec.fix_last_loading = fix_last_loading;
  spec.fix_intercepts = fix_intercepts;
  spec.intercept_columns.assign(intercept_columns.begin(), intercept_columns.end());
  spec.empty_wave_effect.assign(empty_wave_effect.begin(), empty_wave_effect.end());
  spec.attitude_names = as_strings(attitude_names);
  spec.covariate_names = as_strings(covariate_names);
  spec.wave_names = as_strings(wave_names);
  problem->layout = attitudedyn::make_layout(spec);
  if (is_free.size() != problem->layout.n_full) {
    stop("is_free must have one flag per coefficient.");
  }
  problem->layout.n_free = 0;
  for (int i = 0; i < problem->layout.n_full; ++i) {
    problem->layout.entries[i].is_free = is_free[i] ? 1 : 0;
    if (problem->layout.entries[i].is_free) ++problem->layout.n_free;
  }
  problem->theta_full.assign(theta_full.begin(), theta_full.end());
  problem->draws.n_draws = n_draws;
  problem->draws.antithetic = antithetic;
  problem->draws.kind = draw_kind_from_int(draw_kind);
  problem->draws.seed = seed_from_numeric(seed);

  const int n_persons = first_wave.size();
  if (last_wave.size() != n_persons || weight.size() != n_persons) {
    stop("first_wave, last_wave, and weight must have one entry per person.");
  }
  if (interviewed.size() != static_cast<R_xlen_t>(n_persons) * n_waves) {
    stop("interviewed must have length n_persons * n_waves.");
  }
  if (category.size() != static_cast<R_xlen_t>(n_persons) * n_attitudes * n_waves) {
    stop("category must have length n_persons * n_attitudes * n_waves.");
  }
  if (covariates.size() != static_cast<R_xlen_t>(n_persons) * n_waves * n_covariates) {
    stop("covariates must have length n_persons * n_waves * n_covariates.");
  }
  problem->persons.resize(n_persons);
  for (int person = 0; person < n_persons; ++person) {
    attitudedyn::Person& slot = problem->persons[person];
    slot.first_wave = first_wave[person];
    slot.last_wave = last_wave[person];
    slot.weight = weight[person];
    slot.interviewed.resize(n_waves);
    slot.category.resize(static_cast<std::size_t>(n_attitudes) * n_waves);
    slot.covariates.resize(static_cast<std::size_t>(n_waves) * n_covariates);
    for (int wave = 0; wave < n_waves; ++wave) {
      slot.interviewed[wave] = static_cast<unsigned char>(interviewed[person * n_waves + wave]);
      for (int column = 0; column < n_covariates; ++column) {
        slot.covariates[static_cast<std::size_t>(wave) * n_covariates + column] =
            covariates[(static_cast<R_xlen_t>(person) * n_waves + wave) * n_covariates + column];
      }
    }
    for (int attitude = 0; attitude < n_attitudes; ++attitude) {
      for (int wave = 0; wave < n_waves; ++wave) {
        const R_xlen_t index =
            (static_cast<R_xlen_t>(person) * n_attitudes + attitude) * n_waves + wave;
        slot.category[static_cast<std::size_t>(attitude) * n_waves + wave] = category[index];
      }
    }
  }
  return XPtr<attitudedyn::Problem>(problem, true);
}

// [[Rcpp::export]]
List attitudedyn_eval_cpp(SEXP xp, NumericVector theta, bool theta_is_free, int numerics,
                          bool fortran_ghk, int n_threads, bool want_gradient, bool want_scores,
                          bool return_person, Nullable<NumericVector> factor_draws,
                          Nullable<NumericVector> shock_draws) {
  attitudedyn::Problem* problem = problem_from(xp);
  attitudedyn::EvaluateRequest request;
  request.numerics = numerics_from_int(numerics);
  request.fortran_ghk = fortran_ghk;
  request.n_threads = n_threads;
  request.want_scores = want_scores;
  if (want_gradient || want_scores) {
    if (theta_is_free) {
      for (int i = 0; i < problem->layout.n_full; ++i) {
        if (problem->layout.entries[i].is_free) request.grad_columns.push_back(i);
      }
    } else {
      request.grad_columns.resize(problem->layout.n_full);
      for (int i = 0; i < problem->layout.n_full; ++i) request.grad_columns[i] = i;
    }
  }
  if (factor_draws.isNotNull()) {
    NumericVector values(factor_draws);
    request.factor_draws.assign(values.begin(), values.end());
  }
  if (shock_draws.isNotNull()) {
    NumericVector values(shock_draws);
    request.shock_draws.assign(values.begin(), values.end());
  }
  // NumericVector iterators are double. The core stores external draws as
  // float so a file of Fortran real*4 values round-trips exactly.
  if (!request.factor_draws.empty() || !request.shock_draws.empty()) {
    // assigned above as double via the vector<float> constructor, which casts.
  }
  std::vector<double> theta_cpp(theta.begin(), theta.end());
  attitudedyn::Evaluation evaluation =
      theta_is_free ? attitudedyn::evaluate(*problem, theta_cpp, request)
                    : attitudedyn::evaluate_full(*problem, theta_cpp, request);
  NumericVector gradient(evaluation.gradient.begin(), evaluation.gradient.end());
  List out = List::create(Named("objective") = evaluation.objective,
                          Named("penalty") = evaluation.penalty,
                          Named("n_persons_used") = evaluation.n_persons_used,
                          Named("gradient") = gradient);
  if (return_person) {
    out["person_loglik"] = NumericVector(evaluation.person_loglik.begin(), evaluation.person_loglik.end());
  }
  if (want_scores && !evaluation.person_score.empty()) {
    out["person_score"] =
        NumericVector(evaluation.person_score.begin(), evaluation.person_score.end());
    out["n_gradient"] = static_cast<int>(evaluation.gradient.size());
  }
  return out;
}

// [[Rcpp::export]]
List attitudedyn_natural_cpp(SEXP xp, NumericVector theta_full) {
  attitudedyn::Problem* problem = problem_from(xp);
  std::vector<double> theta(theta_full.begin(), theta_full.end());
  const attitudedyn::NaturalParameters natural =
      attitudedyn::natural_parameters(problem->spec, theta);
  return List::create(Named("persistence") = NumericVector(natural.persistence.begin(), natural.persistence.end()),
                      Named("loadings") = NumericVector(natural.loadings.begin(), natural.loadings.end()),
                      Named("wave_effects") = NumericVector(natural.wave_effects.begin(), natural.wave_effects.end()),
                      Named("beta") = NumericVector(natural.beta.begin(), natural.beta.end()),
                      Named("cutpoints") = [&] {
                        List cuts(problem->spec.n_attitudes);
                        for (int j = 0; j < problem->spec.n_attitudes; ++j) {
                          cuts[j] = NumericVector(natural.cutpoints[j].begin(), natural.cutpoints[j].end());
                        }
                        return cuts;
                      }());
}

// [[Rcpp::export]]
List attitudedyn_irf_cpp(NumericVector persistence, int n_attitudes, int horizon) {
  std::vector<double> matrix(persistence.begin(), persistence.end());
  const attitudedyn::ImpulseResponse irf =
      attitudedyn::impulse_response(matrix, n_attitudes, horizon);
  return List::create(Named("response") = NumericVector(irf.response.begin(), irf.response.end()),
                      Named("horizon") = irf.horizon, Named("n_attitudes") = irf.n_attitudes,
                      Named("spectral_radius") = irf.spectral_radius);
}

// [[Rcpp::export]]
List attitudedyn_influence_cpp(NumericVector persistence, int n_attitudes) {
  std::vector<double> matrix(persistence.begin(), persistence.end());
  const attitudedyn::InfluenceShares shares = attitudedyn::influence_shares(matrix, n_attitudes);
  return List::create(Named("share") = NumericVector(shares.share.begin(), shares.share.end()),
                      Named("long_run_variance") =
                          NumericVector(shares.long_run_variance.begin(), shares.long_run_variance.end()),
                      Named("stationary") = shares.stationary);
}

// [[Rcpp::export]]
List attitudedyn_variance_cpp(SEXP xp, NumericVector theta_full) {
  attitudedyn::Problem* problem = problem_from(xp);
  std::vector<double> theta(theta_full.begin(), theta_full.end());
  const attitudedyn::VarianceComponent component =
      attitudedyn::variance_decomposition(*problem, theta);
  return List::create(Named("covariate") = NumericVector(component.covariate.begin(), component.covariate.end()),
                      Named("wave_effect") = NumericVector(component.wave_effect.begin(), component.wave_effect.end()),
                      Named("factor") = NumericVector(component.factor.begin(), component.factor.end()),
                      Named("shock") = NumericVector(component.shock.begin(), component.shock.end()),
                      Named("total") = NumericVector(component.total.begin(), component.total.end()),
                      Named("note") = component.note);
}

// [[Rcpp::export]]
List attitudedyn_simulate_cpp(int n_waves, int n_attitudes, int n_covariates, int n_factors,
                              IntegerVector n_categories, NumericVector wave_time,
                              bool include_wave_effects, bool free_contemporaneous, bool fix_last_loading,
                              bool fix_intercepts, IntegerVector intercept_columns,
                              IntegerVector empty_wave_effect, CharacterVector attitude_names,
                              CharacterVector covariate_names, CharacterVector wave_names,
                              NumericVector theta_full, int n_persons, double seed) {
  attitudedyn::ModelSpec spec;
  spec.n_waves = n_waves;
  spec.n_attitudes = n_attitudes;
  spec.n_covariates = n_covariates;
  spec.n_factors = n_factors;
  spec.n_categories.assign(n_categories.begin(), n_categories.end());
  spec.wave_time.assign(wave_time.begin(), wave_time.end());
  spec.include_wave_effects = include_wave_effects;
  spec.free_contemporaneous = free_contemporaneous;
  spec.fix_last_loading = fix_last_loading;
  spec.fix_intercepts = fix_intercepts;
  spec.intercept_columns.assign(intercept_columns.begin(), intercept_columns.end());
  spec.empty_wave_effect.assign(empty_wave_effect.begin(), empty_wave_effect.end());
  spec.attitude_names = as_strings(attitude_names);
  spec.covariate_names = as_strings(covariate_names);
  spec.wave_names = as_strings(wave_names);
  std::vector<double> theta(theta_full.begin(), theta_full.end());
  const attitudedyn::SimulatedPanel panel =
      attitudedyn::simulate_panel(spec, theta, n_persons, seed_from_numeric(seed));
  return List::create(Named("id") = IntegerVector(panel.id.begin(), panel.id.end()),
                      Named("wave") = IntegerVector(panel.wave.begin(), panel.wave.end()),
                      Named("time") = NumericVector(panel.time.begin(), panel.time.end()),
                      Named("covariates") = NumericVector(panel.covariates.begin(), panel.covariates.end()),
                      Named("outcome") = IntegerVector(panel.outcome.begin(), panel.outcome.end()),
                      Named("n_covariates") = panel.n_covariates, Named("n_attitudes") = panel.n_attitudes);
}

// [[Rcpp::export]]
List attitudedyn_predict_cpp(SEXP xp, NumericVector theta_full, int n_draws, double seed,
                             int n_threads) {
  attitudedyn::Problem* problem = problem_from(xp);
  std::vector<double> theta(theta_full.begin(), theta_full.end());
  const attitudedyn::Prediction prediction = attitudedyn::predict_probabilities(
      *problem, theta, n_draws, seed_from_numeric(seed), n_threads);
  if (!prediction.error.empty()) stop("%s", prediction.error.c_str());
  return List::create(
      Named("probability") = NumericVector(prediction.probability.begin(), prediction.probability.end()),
      Named("mean_latent") = NumericVector(prediction.mean_latent.begin(), prediction.mean_latent.end()),
      Named("max_categories") = prediction.max_categories);
}
