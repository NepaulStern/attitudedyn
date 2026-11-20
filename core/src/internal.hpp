#ifndef ATTITUDYN_INTERNAL_HPP
#define ATTITUDYN_INTERNAL_HPP

#include "attitudedyn/api.hpp"

#include <vector>

namespace attitudedyn {
namespace detail {

struct ThresholdIncrement {
  int grad = -1;       // column of the requested gradient, or -1
  double exp_c = 0.0;  // the positive gap between successive interior cuts
  int first_cut = 0;   // smallest 0-based cut index that includes this increment
};

struct Unpacked {
  int n_attitudes = 0;
  int n_waves = 0;
  int n_covariates = 0;
  int n_factors = 0;
  int n_grad = 0;
  bool identity_correlation = true;
  double penalty = 0.0;

  std::vector<double> wave_effect;  // J * T
  std::vector<double> beta;         // J * X
  std::vector<double> rho;          // J * J, natural scale
  std::vector<double> drho_dc;      // J * J
  std::vector<double> loading;      // (factor * J + attitude)
  std::vector<double> threshold;    // (attitude * max_cuts + cut), max_cuts = max(K)+1
  int max_cuts = 0;
  std::vector<double> omega;        // J * J row-major correlation

  // Gradient column of each structural element, or -1.
  std::vector<int> wave_grad;
  std::vector<int> beta_grad;
  std::vector<int> rho_grad;
  std::vector<int> loading_grad;
  std::vector<std::vector<ThresholdIncrement>> threshold_increments;

  // chol gradient columns in the order the coefficients appear, for the
  // finite-difference penalty gradient. Empty if none of them are differentiated.
  std::vector<int> chol_position;  // index into theta_full
  std::vector<int> chol_grad;
};

// grad_of_full[i] is the gradient column of coefficient i, or -1.
Unpacked unpack(const ModelSpec& spec, const std::vector<double>& theta_full,
                const std::vector<int>& grad_of_full);

double person_contribution(const ModelSpec& spec, const Unpacked& theta,
                           const Person& person, int person_index,
                           const DrawOptions& draws, const float* factor_draws,
                           const float* shock_draws, Numerics numerics, bool fortran_ghk,
                           double* gradient, double* score_if_wanted);

}  // namespace detail
}  // namespace attitudedyn

#endif
