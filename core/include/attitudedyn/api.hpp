#ifndef ATTITUDYN_API_HPP
#define ATTITUDYN_API_HPP

// Public C++ API of the attitude-dynamics core.
//
// This header includes only the standard library. The R package, and any
// later Python or Stata wrapper, should talk to the library through the
// types and functions declared here. Nothing in this interface mentions R.

#include <cstdint>
#include <string>
#include <vector>

namespace attitudedyn {

// How the lagged latent state is formed. The default matches the Fortran
// program (the whole latent attitude persists). Shock-only is equation (2)
// of Nepaul & Stern (2026): only the innovation is autoregressive.
enum class Dynamics {
  kLatentAttitude = 0,
  kShockOnly = 1
};

// kStable: double-precision log(Phi(b)-Phi(a)), analytic gradient.
// kShared: Fortran probmake1 tail blend, but with the same double-precision
//          Phi and Phi^{-1} the patched Fortran oracle calls. No gradient.
enum class Numerics {
  kStable = 0,
  kShared = 1
};

// Quasi-random streams are Owen-scrambled and keyed by the person index,
// so two people do not share one sequence and the result does not depend
// on the thread count.
enum class DrawKind {
  kPhilox = 0,
  kHalton = 1,
  kSobol = 2
};

struct DrawOptions {
  DrawKind kind = DrawKind::kPhilox;
  int n_draws = 10;
  bool antithetic = true;
  std::uint64_t seed = 1;
};

struct ModelSpec {
  int n_waves = 0;
  int n_attitudes = 0;
  int n_covariates = 0;
  int n_factors = 0;
  std::vector<int> n_categories;     // length n_attitudes, each >= 2
  std::vector<double> wave_time;     // length n_waves, strictly increasing
  bool include_wave_effects = true;
  bool free_contemporaneous = false;
  bool fix_last_loading = true;      // lambda[J, M] = 0 when M >= 2
  bool fix_intercepts = true;        // collinear with a full set of wave effects
  std::vector<int> intercept_columns;
  // 1 if attitude a has no observed responses in wave t (index a * n_waves + t).
  std::vector<char> empty_wave_effect;
  std::vector<std::string> attitude_names;
  std::vector<std::string> covariate_names;
  std::vector<std::string> wave_names;
};

struct LayoutEntry {
  std::string name;
  std::string block;  // wave_effect, beta, persistence, loading, chol, threshold
  int attitude = -1;  // 0-based, or -1
  int index2 = -1;    // wave, covariate, source attitude, factor, or cut
  char is_free = 1;
};

struct Layout {
  int n_full = 0;
  int n_free = 0;
  std::vector<LayoutEntry> entries;
};

Layout make_layout(const ModelSpec& spec);

// One person. Waves outside [first_wave, last_wave] are not in the
// recursion. A wave inside the window with interviewed[t] == 0 is a gap:
// the latent state is simulated forward and covariates are ignored
// (they must still be finite; the contribution uses a zero regressor
// vector, matching an uninitialised Fortran slot that the loader zeroed).
struct Person {
  int first_wave = 0;  // 0-based, inclusive
  int last_wave = -1;  // inclusive; last < first means "drop this person"
  std::vector<unsigned char> interviewed;  // length n_waves
  std::vector<int> category;               // (attitude * n_waves + wave), 0 = missing
  std::vector<double> covariates;          // (wave * n_covariates + k)
  double weight = 1.0;
};

struct Problem {
  ModelSpec spec;
  Layout layout;
  std::vector<Person> persons;
  DrawOptions draws;
  // Values of every coefficient, including the fixed ones. Free entries
  // are overwritten by the argument of evaluate().
  std::vector<double> theta_full;
};

struct Evaluation {
  double objective = 0.0;   // sum of weighted person log-likelihoods, minus penalty
  double penalty = 0.0;     // non-positive-definite contemporaneous correlation
  int n_persons_used = 0;
  std::vector<double> person_loglik;  // length n_persons; 0 if the person was dropped
  std::vector<double> gradient;       // d(objective) / d(theta_full[grad_columns])
  std::vector<double> person_score;   // row-major n_persons * gradient.size(), if requested
  std::string error;
};

struct EvaluateRequest {
  Dynamics dynamics = Dynamics::kLatentAttitude;
  Numerics numerics = Numerics::kStable;
  // Reproduce the ghk1 overwrite when the contemporaneous correlation is
  // not the identity. Ignored when that correlation is the identity, which
  // is the published specification: both versions then return the same
  // conditional moments.
  bool fortran_ghk = false;
  int n_threads = 1;
  bool want_scores = false;
  // Which full-coefficient positions to differentiate, in gradient order.
  // Empty means "no gradient".
  std::vector<int> grad_columns;
  // Optional standard normal draws. When non-empty, DrawOptions::kind is
  // ignored and these values are used (the Fortran comparison does this).
  // factor: ((person * n_draws + draw) * n_factors + factor)
  // shock:  (((person * n_draws + draw) * n_attitudes + attitude) * n_waves + wave)
  std::vector<float> factor_draws;
  std::vector<float> shock_draws;
};

// theta_free is laid out in the order of the free coefficients. It replaces
// problem.theta_full at those positions for this evaluation only.
Evaluation evaluate(const Problem& problem, const std::vector<double>& theta_free,
                    const EvaluateRequest& request);

// Same evaluation against an explicit full coefficient vector.
Evaluation evaluate_full(const Problem& problem, const std::vector<double>& theta_full,
                         const EvaluateRequest& request);

struct NaturalParameters {
  // Row-major n_attitudes * n_attitudes. Row is the receiving attitude,
  // column is the source, which is Fortran rho(to, from) and the matrix
  // that multiplies the lagged state.
  std::vector<double> persistence;
  // For each attitude, the K+1 cutpoints, with the normalised ends
  // (-100, 0, ..., +100).
  std::vector<std::vector<double>> cutpoints;
  std::vector<double> loadings;       // attitude-major within factor, factor-major: (f * J + j)
  std::vector<double> wave_effects;   // (attitude * n_waves + wave), 0 if omitted
  std::vector<double> beta;           // (attitude * n_covariates + k)
};

NaturalParameters natural_parameters(const ModelSpec& spec, const std::vector<double>& theta_full);

// Impulse response of the persistence matrix. response[(h * J + to) * J + from]
// is the effect, h steps of one time unit later, of a unit innovation in
// attitude `from` on attitude `to`. Horizon 0 is the identity.
struct ImpulseResponse {
  int horizon = 0;
  int n_attitudes = 0;
  std::vector<double> response;
  double spectral_radius = 0.0;
};

ImpulseResponse impulse_response(const std::vector<double>& persistence_row_major,
                                 int n_attitudes, int horizon);

// Share of the long-run variance of attitude `to` that is attributable to
// an innovation in attitude `from`, under the equal-spacing recursion
// y = R y_lag + innovation (or the same recursion for the shock alone).
// shares[to * J + from]. Empty if the spectral radius is not strictly below 1.
struct InfluenceShares {
  std::vector<double> share;          // row-major to, from
  std::vector<double> long_run_variance;
  bool stationary = false;
};

InfluenceShares influence_shares(const std::vector<double>& persistence_row_major,
                                 int n_attitudes);

struct VarianceComponent {
  std::vector<double> covariate;    // per attitude
  std::vector<double> wave_effect;
  std::vector<double> factor;
  std::vector<double> shock;
  std::vector<double> total;
  std::string note;
};

// Sample variances of x beta and of the wave effects, plus ||lambda||^2 and
// the long-run shock variance. Wave effects are not also folded into the
// covariate piece.
VarianceComponent variance_decomposition(const Problem& problem,
                                         const std::vector<double>& theta_full,
                                         Dynamics dynamics);

// The decomposition implemented by the Fortran vardecomp subroutine,
// including the double appearance of the wave effect (once inside x beta
// and once on its own) and the product-of-scalars shock accumulation.
VarianceComponent variance_decomposition_fortran(const Problem& problem,
                                                  const std::vector<double>& theta_full);

struct SimulatedPanel {
  int n_persons = 0;
  int n_waves = 0;
  int n_attitudes = 0;
  int n_covariates = 0;
  std::vector<int> id;          // length n_persons * n_waves, person-major
  std::vector<int> wave;        // 1-based wave index
  std::vector<double> time;
  std::vector<double> covariates;
  std::vector<int> outcome;     // (row * n_attitudes + attitude), categories 1..K
};

// Draw a panel from the structural model (not from the GHK truncation).
// Covariates are N(0, 1), except intercept columns which are 1.
SimulatedPanel simulate_panel(const ModelSpec& spec, const std::vector<double>& theta_full,
                              int n_persons, Dynamics dynamics, std::uint64_t seed);

// Category probabilities for every person-wave in `persons`, averaged over
// `n_draws` structural paths (Rao-Blackwellised interval probabilities).
// probabilities layout:
//   ((person * n_waves + wave) * n_attitudes + attitude) * max_categories + (k-1)
// max_categories is the maximum of n_categories. Unused category slots are 0.
// `mean_latent` is the average simulated latent attitude, same person/wave/attitude
// order without the category dimension.
struct Prediction {
  int max_categories = 0;
  std::vector<double> probability;
  std::vector<double> mean_latent;
  std::string error;
};

Prediction predict_probabilities(const Problem& problem, const std::vector<double>& theta_full,
                                 Dynamics dynamics, int n_draws, std::uint64_t seed,
                                 int n_threads);

}  // namespace attitudedyn

#endif
