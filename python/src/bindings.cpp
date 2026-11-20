// Python binding for the attitude-dynamics core.
//
// This file translates NumPy arrays into the types declared in
// attitudedyn/api.hpp and copies the results back. The likelihood, the
// analytic gradient, the draws, and the post-estimation routines are the
// functions in core/. Nothing here reimplements them.

#include "attitudedyn/api.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace py = pybind11;

namespace {

template <typename F>
auto guard(F&& function) -> decltype(function()) {
  try {
    return function();
  } catch (const std::invalid_argument& exc) {
    throw py::value_error(exc.what());
  }
}

std::vector<double> as_f64(const py::handle& object) {
  py::array_t<double, py::array::c_style | py::array::forcecast> array(
      py::reinterpret_borrow<py::object>(object));
  const auto buffer = array.request();
  const auto* pointer = static_cast<const double*>(buffer.ptr);
  return std::vector<double>(pointer, pointer + buffer.size);
}

std::vector<int> as_i32(const py::handle& object) {
  py::array_t<std::int64_t, py::array::c_style | py::array::forcecast> array(
      py::reinterpret_borrow<py::object>(object));
  const auto buffer = array.request();
  const auto* pointer = static_cast<const std::int64_t*>(buffer.ptr);
  std::vector<int> out(static_cast<std::size_t>(buffer.size));
  for (py::ssize_t i = 0; i < buffer.size; ++i) out[static_cast<std::size_t>(i)] = static_cast<int>(pointer[i]);
  return out;
}

std::vector<std::string> as_strings(const py::handle& object) {
  std::vector<std::string> out;
  for (const py::handle item : object) out.push_back(py::cast<std::string>(item));
  return out;
}

std::uint64_t as_seed(const py::handle& object) {
  const double value = py::cast<double>(object);
  if (!std::isfinite(value) || value < 0.0) {
    throw std::invalid_argument("`seed` must be a finite, non-negative number.");
  }
  return static_cast<std::uint64_t>(value);
}

py::array_t<double> as_array(const std::vector<double>& values) {
  py::array_t<double> out(static_cast<py::ssize_t>(values.size()));
  if (!values.empty()) {
    std::memcpy(out.mutable_data(), values.data(), values.size() * sizeof(double));
  }
  return out;
}

py::array_t<int> as_array_int(const std::vector<int>& values) {
  py::array_t<int> out(static_cast<py::ssize_t>(values.size()));
  if (!values.empty()) {
    std::memcpy(out.mutable_data(), values.data(), values.size() * sizeof(int));
  }
  return out;
}

py::dict layout_dict(const attitudedyn::Layout& layout) {
  py::list name, block;
  py::array_t<int> attitude(layout.n_full);
  py::array_t<int> index2(layout.n_full);
  py::array_t<int> is_free(layout.n_full);
  auto attitude_ptr = attitude.mutable_data();
  auto index_ptr = index2.mutable_data();
  auto free_ptr = is_free.mutable_data();
  for (int i = 0; i < layout.n_full; ++i) {
    name.append(layout.entries[i].name);
    block.append(layout.entries[i].block);
    attitude_ptr[i] = layout.entries[i].attitude;
    index_ptr[i] = layout.entries[i].index2;
    free_ptr[i] = layout.entries[i].is_free ? 1 : 0;
  }
  py::dict out;
  out["n_full"] = layout.n_full;
  out["n_free"] = layout.n_free;
  out["name"] = name;
  out["block"] = block;
  out["attitude"] = attitude;
  out["index2"] = index2;
  out["is_free"] = is_free;
  return out;
}

attitudedyn::ModelSpec parse_spec(const py::dict& spec) {
  attitudedyn::ModelSpec out;
  out.n_waves = spec["n_waves"].cast<int>();
  out.n_attitudes = spec["n_attitudes"].cast<int>();
  out.n_covariates = spec["n_covariates"].cast<int>();
  out.n_factors = spec["n_factors"].cast<int>();
  out.n_categories = as_i32(spec["n_categories"]);
  out.wave_time = as_f64(spec["wave_time"]);
  out.include_wave_effects = spec["include_wave_effects"].cast<bool>();
  out.free_contemporaneous = spec["free_contemporaneous"].cast<bool>();
  out.fix_last_loading = spec["fix_last_loading"].cast<bool>();
  out.fix_intercepts = spec["fix_intercepts"].cast<bool>();
  out.intercept_columns = as_i32(spec["intercept_columns"]);
  const std::vector<int> empty = as_i32(spec["empty_wave_effect"]);
  out.empty_wave_effect.assign(empty.begin(), empty.end());
  out.attitude_names = as_strings(spec["attitude_names"]);
  out.covariate_names = as_strings(spec["covariate_names"]);
  out.wave_names = as_strings(spec["wave_names"]);
  return out;
}

attitudedyn::DrawKind draw_kind_from_int(int kind) {
  if (kind == 1) return attitudedyn::DrawKind::kHalton;
  if (kind == 2) return attitudedyn::DrawKind::kSobol;
  return attitudedyn::DrawKind::kPhilox;
}

attitudedyn::Numerics numerics_from_int(int numerics) {
  return numerics == 1 ? attitudedyn::Numerics::kShared : attitudedyn::Numerics::kStable;
}

class Problem {
 public:
  attitudedyn::Problem impl;

  py::dict layout() const { return layout_dict(impl.layout); }

  int n_persons() const { return static_cast<int>(impl.persons.size()); }
};

void fill_persons(attitudedyn::Problem& problem, const py::dict& data) {
  const std::vector<int> first_wave = as_i32(data["first_wave"]);
  const std::vector<int> last_wave = as_i32(data["last_wave"]);
  const std::vector<double> weight = as_f64(data["weight"]);
  const std::vector<int> interviewed = as_i32(data["interviewed"]);
  const std::vector<int> category = as_i32(data["category"]);
  const std::vector<double> covariates = as_f64(data["covariates"]);
  const int n_persons = static_cast<int>(first_wave.size());
  const int n_waves = problem.spec.n_waves;
  const int n_attitudes = problem.spec.n_attitudes;
  const int n_covariates = problem.spec.n_covariates;
  if (static_cast<int>(last_wave.size()) != n_persons || static_cast<int>(weight.size()) != n_persons) {
    throw std::invalid_argument("first_wave, last_wave, and weight must have one entry per person.");
  }
  if (static_cast<int>(interviewed.size()) != n_persons * n_waves) {
    throw std::invalid_argument("interviewed must have length n_persons * n_waves.");
  }
  if (static_cast<int>(category.size()) != n_persons * n_attitudes * n_waves) {
    throw std::invalid_argument("category must have length n_persons * n_attitudes * n_waves.");
  }
  if (static_cast<int>(covariates.size()) != n_persons * n_waves * n_covariates) {
    throw std::invalid_argument("covariates must have length n_persons * n_waves * n_covariates.");
  }
  problem.persons.resize(static_cast<std::size_t>(n_persons));
  for (int person = 0; person < n_persons; ++person) {
    attitudedyn::Person& slot = problem.persons[person];
    slot.first_wave = first_wave[person];
    slot.last_wave = last_wave[person];
    slot.weight = weight[person];
    slot.interviewed.resize(static_cast<std::size_t>(n_waves));
    slot.category.resize(static_cast<std::size_t>(n_attitudes) * n_waves);
    slot.covariates.resize(static_cast<std::size_t>(n_waves) * n_covariates);
    for (int wave = 0; wave < n_waves; ++wave) {
      slot.interviewed[wave] = static_cast<unsigned char>(interviewed[person * n_waves + wave]);
      for (int column = 0; column < n_covariates; ++column) {
        slot.covariates[static_cast<std::size_t>(wave) * n_covariates + column] =
            covariates[(static_cast<std::size_t>(person) * n_waves + wave) * n_covariates + column];
      }
    }
    for (int attitude = 0; attitude < n_attitudes; ++attitude) {
      for (int wave = 0; wave < n_waves; ++wave) {
        const std::size_t index =
            (static_cast<std::size_t>(person) * n_attitudes + attitude) * n_waves + wave;
        slot.category[static_cast<std::size_t>(attitude) * n_waves + wave] = category[index];
      }
    }
  }
}

py::dict pack_evaluation(const attitudedyn::Evaluation& evaluation, bool want_scores) {
  py::dict out;
  out["objective"] = evaluation.objective;
  out["penalty"] = evaluation.penalty;
  out["n_persons_used"] = evaluation.n_persons_used;
  out["gradient"] = as_array(evaluation.gradient);
  out["error"] = evaluation.error;
  if (!evaluation.person_loglik.empty()) out["person_loglik"] = as_array(evaluation.person_loglik);
  if (want_scores && !evaluation.person_score.empty() && !evaluation.gradient.empty()) {
    const int n_grad = static_cast<int>(evaluation.gradient.size());
    const int n_persons = static_cast<int>(evaluation.person_score.size() / static_cast<std::size_t>(n_grad));
    py::array_t<double> scores({n_persons, n_grad});
    std::memcpy(scores.mutable_data(), evaluation.person_score.data(),
                evaluation.person_score.size() * sizeof(double));
    out["person_score"] = scores;
  }
  return out;
}

py::dict pack_panel(const attitudedyn::SimulatedPanel& panel) {
  py::dict out;
  out["n_persons"] = panel.n_persons;
  out["n_waves"] = panel.n_waves;
  out["n_attitudes"] = panel.n_attitudes;
  out["n_covariates"] = panel.n_covariates;
  out["id"] = as_array_int(panel.id);
  out["wave"] = as_array_int(panel.wave);
  out["time"] = as_array(panel.time);
  out["covariates"] = as_array(panel.covariates);
  out["outcome"] = as_array_int(panel.outcome);
  return out;
}

}  // namespace

PYBIND11_MODULE(_core, module) {
  module.doc() =
      "Binding to the shared attitudedyn C++ core. "
      "The numerical work is implemented in core/; this module only moves arrays across the boundary.";

  module.def(
      "layout",
      [](const py::dict& spec) {
        return guard([&] {
          const attitudedyn::Layout layout = attitudedyn::make_layout(parse_spec(spec));
          return layout_dict(layout);
        });
      },
      py::arg("spec"));

  py::class_<Problem>(module, "Problem")
      .def("layout", &Problem::layout)
      .def("n_persons", &Problem::n_persons)
      .def(
          "evaluate",
          [](Problem& problem, const py::array& theta, bool theta_is_free, int numerics, bool fortran_ghk,
             int n_threads, bool want_gradient, bool want_scores) {
            return guard([&] {
              attitudedyn::EvaluateRequest request;
              request.numerics = numerics_from_int(numerics);
              request.fortran_ghk = fortran_ghk;
              request.n_threads = n_threads;
              request.want_scores = want_scores;
              if (want_gradient || want_scores) {
                if (theta_is_free) {
                  for (int i = 0; i < problem.impl.layout.n_full; ++i) {
                    if (problem.impl.layout.entries[i].is_free) request.grad_columns.push_back(i);
                  }
                } else {
                  request.grad_columns.resize(static_cast<std::size_t>(problem.impl.layout.n_full));
                  for (int i = 0; i < problem.impl.layout.n_full; ++i) request.grad_columns[i] = i;
                }
              }
              const std::vector<double> theta_cpp = as_f64(theta);
              const attitudedyn::Evaluation evaluation =
                  theta_is_free ? attitudedyn::evaluate(problem.impl, theta_cpp, request)
                                : attitudedyn::evaluate_full(problem.impl, theta_cpp, request);
              return pack_evaluation(evaluation, want_scores);
            });
          },
          py::arg("theta"), py::arg("theta_is_free"), py::arg("numerics"), py::arg("fortran_ghk"),
          py::arg("n_threads"), py::arg("want_gradient"), py::arg("want_scores"))
      .def(
          "natural",
          [](Problem& problem, const py::array& theta_full) {
            return guard([&] {
              const attitudedyn::NaturalParameters natural =
                  attitudedyn::natural_parameters(problem.impl.spec, as_f64(theta_full));
              py::dict out;
              out["persistence"] = as_array(natural.persistence);
              out["loadings"] = as_array(natural.loadings);
              out["wave_effects"] = as_array(natural.wave_effects);
              out["beta"] = as_array(natural.beta);
              py::list cuts;
              for (const std::vector<double>& cut : natural.cutpoints) cuts.append(as_array(cut));
              out["cutpoints"] = cuts;
              return out;
            });
          },
          py::arg("theta_full"))
      .def(
          "variance",
          [](Problem& problem, const py::array& theta_full) {
            return guard([&] {
              const attitudedyn::VarianceComponent component =
                  attitudedyn::variance_decomposition(problem.impl, as_f64(theta_full));
              py::dict out;
              out["covariate"] = as_array(component.covariate);
              out["wave_effect"] = as_array(component.wave_effect);
              out["factor"] = as_array(component.factor);
              out["shock"] = as_array(component.shock);
              out["total"] = as_array(component.total);
              out["note"] = component.note;
              return out;
            });
          },
          py::arg("theta_full"))
      .def(
          "predict",
          [](Problem& problem, const py::array& theta_full, int n_draws, const py::object& seed, int n_threads) {
            return guard([&] {
              const attitudedyn::Prediction prediction = attitudedyn::predict_probabilities(
                  problem.impl, as_f64(theta_full), n_draws, as_seed(seed), n_threads);
              if (!prediction.error.empty()) throw std::invalid_argument(prediction.error);
              py::dict out;
              out["probability"] = as_array(prediction.probability);
              out["mean_latent"] = as_array(prediction.mean_latent);
              out["max_categories"] = prediction.max_categories;
              return out;
            });
          },
          py::arg("theta_full"), py::arg("n_draws"), py::arg("seed"), py::arg("n_threads"))
      .def(
          "simulate",
          [](Problem& problem, const py::array& theta_full, int n_persons, const py::object& seed) {
            return guard([&] {
              return pack_panel(
                  attitudedyn::simulate_panel(problem.impl.spec, as_f64(theta_full), n_persons, as_seed(seed)));
            });
          },
          py::arg("theta_full"), py::arg("n_persons"), py::arg("seed"));

  module.def(
      "build",
      [](const py::dict& spec, const py::dict& data, const py::array& theta_full, const py::array& is_free,
         const py::dict& draws) {
        return guard([&] {
          auto problem = std::make_unique<Problem>();
          problem->impl.spec = parse_spec(spec);
          problem->impl.layout = attitudedyn::make_layout(problem->impl.spec);
          const std::vector<int> free = as_i32(is_free);
          if (static_cast<int>(free.size()) != problem->impl.layout.n_full) {
            throw std::invalid_argument("is_free must have one flag per coefficient.");
          }
          problem->impl.layout.n_free = 0;
          for (int i = 0; i < problem->impl.layout.n_full; ++i) {
            problem->impl.layout.entries[i].is_free = free[i] ? 1 : 0;
            if (problem->impl.layout.entries[i].is_free) ++problem->impl.layout.n_free;
          }
          problem->impl.theta_full = as_f64(theta_full);
          if (static_cast<int>(problem->impl.theta_full.size()) != problem->impl.layout.n_full) {
            throw std::invalid_argument("theta_full must have one value per coefficient.");
          }
          problem->impl.draws.n_draws = draws["n_draws"].cast<int>();
          problem->impl.draws.antithetic = draws["antithetic"].cast<bool>();
          problem->impl.draws.kind = draw_kind_from_int(draws["kind"].cast<int>());
          problem->impl.draws.seed = as_seed(draws["seed"]);
          fill_persons(problem->impl, data);
          return problem;
        });
      },
      py::arg("spec"), py::arg("data"), py::arg("theta_full"), py::arg("is_free"), py::arg("draws"));

  module.def(
      "simulate_panel",
      [](const py::dict& spec, const py::array& theta_full, int n_persons, const py::object& seed) {
        return guard([&] {
          const attitudedyn::ModelSpec model = parse_spec(spec);
          return pack_panel(attitudedyn::simulate_panel(model, as_f64(theta_full), n_persons, as_seed(seed)));
        });
      },
      py::arg("spec"), py::arg("theta_full"), py::arg("n_persons"), py::arg("seed"));

  module.def(
      "impulse_response",
      [](const py::array& persistence, int n_attitudes, int horizon) {
        return guard([&] {
          const attitudedyn::ImpulseResponse response =
              attitudedyn::impulse_response(as_f64(persistence), n_attitudes, horizon);
          py::dict out;
          out["response"] = as_array(response.response);
          out["horizon"] = response.horizon;
          out["n_attitudes"] = response.n_attitudes;
          out["spectral_radius"] = response.spectral_radius;
          return out;
        });
      },
      py::arg("persistence"), py::arg("n_attitudes"), py::arg("horizon"));

  module.def(
      "influence_shares",
      [](const py::array& persistence, int n_attitudes) {
        return guard([&] {
          const attitudedyn::InfluenceShares shares =
              attitudedyn::influence_shares(as_f64(persistence), n_attitudes);
          py::dict out;
          out["share"] = as_array(shares.share);
          out["long_run_variance"] = as_array(shares.long_run_variance);
          out["stationary"] = shares.stationary;
          return out;
        });
      },
      py::arg("persistence"), py::arg("n_attitudes"));
}
