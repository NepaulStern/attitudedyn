// Compare the C++ likelihood with the Fortran oracle dump.
//
// The dump is written by validation/fortran (method = 0) at the shipped
// GSS coefficients. Draws are the normals Fortran stored, in this order:
//   factor: ((person * n_draws + draw) * n_factors + factor)
//   shock:  (((person * n_draws + draw) * n_attitudes + attitude) * n_waves + wave)
// Person 0 is Fortran person 1. An excluded person has last_wave < first_wave.

#include "attitudedyn/api.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

template <typename T>
T read_pod(std::istream& in) {
  T value{};
  in.read(reinterpret_cast<char*>(&value), sizeof(T));
  if (!in) throw std::runtime_error("The oracle dump ended early.");
  return value;
}

struct Dump {
  int n_persons = 0;
  int n_attitudes = 0;
  int n_covariates = 0;
  int n_waves = 0;
  int n_factors = 0;
  int n_draws = 0;
  int n_coef = 0;
  std::vector<int> n_categories;
  std::vector<double> time;
  std::vector<double> coef;
  std::vector<int> icoef;
  std::vector<double> person_loglik;
  std::vector<attitudedyn::Person> persons;
  std::vector<float> factor_draws;
  std::vector<float> shock_draws;
};

Dump read_dump(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("Cannot open " + path);
  Dump dump;
  const std::int32_t magic = read_pod<std::int32_t>(in);
  if (magic != 20261120) throw std::runtime_error("oracle_body.bin has the wrong magic number.");
  auto i32 = [&]() { return static_cast<int>(read_pod<std::int32_t>(in)); };
  dump.n_persons = i32();
  dump.n_attitudes = i32();
  dump.n_covariates = i32();
  dump.n_waves = i32();
  dump.n_factors = i32();
  dump.n_draws = i32();
  dump.n_coef = i32();
  dump.n_categories.resize(static_cast<std::size_t>(dump.n_attitudes));
  for (int& k : dump.n_categories) k = i32();
  dump.time.resize(static_cast<std::size_t>(dump.n_waves));
  for (double& t : dump.time) t = read_pod<double>(in);
  dump.coef.resize(static_cast<std::size_t>(dump.n_coef));
  for (double& c : dump.coef) c = read_pod<double>(in);
  dump.icoef.resize(static_cast<std::size_t>(dump.n_coef));
  for (int& flag : dump.icoef) flag = i32();
  dump.person_loglik.resize(static_cast<std::size_t>(dump.n_persons));
  for (double& ll : dump.person_loglik) ll = read_pod<double>(in);

  const int T = dump.n_waves;
  const int J = dump.n_attitudes;
  const int X = dump.n_covariates;
  dump.persons.resize(static_cast<std::size_t>(dump.n_persons));
  for (int person = 0; person < dump.n_persons; ++person) {
    attitudedyn::Person& p = dump.persons[person];
    const int included = i32();
    const int istart = i32();
    const int iend = i32();
    p.interviewed.assign(static_cast<std::size_t>(T), 0);
    for (int wave = 0; wave < T; ++wave) p.interviewed[wave] = static_cast<unsigned char>(i32());
    p.category.assign(static_cast<std::size_t>(J * T), 0);
    for (int attitude = 0; attitude < J; ++attitude) {
      for (int wave = 0; wave < T; ++wave) {
        const int iy = i32();
        if (iy >= 1 && iy <= dump.n_categories[attitude]) {
          p.category[attitude * T + wave] = iy;
        }
      }
    }
    p.covariates.assign(static_cast<std::size_t>(T * X), 0.0);
    for (int wave = 0; wave < T; ++wave) {
      for (int column = 0; column < X; ++column) {
        p.covariates[wave * X + column] = read_pod<float>(in);
      }
    }
    if (included) {
      p.first_wave = istart - 1;
      p.last_wave = iend - 1;
      p.weight = 1.0;
    } else {
      p.first_wave = 0;
      p.last_wave = -1;
      p.weight = 0.0;
    }
  }

  const std::size_t n_factor =
      static_cast<std::size_t>(dump.n_persons) * dump.n_draws * dump.n_factors;
  const std::size_t n_shock =
      static_cast<std::size_t>(dump.n_persons) * dump.n_draws * dump.n_attitudes * dump.n_waves;
  dump.factor_draws.resize(n_factor);
  dump.shock_draws.resize(n_shock);
  in.read(reinterpret_cast<char*>(dump.factor_draws.data()),
          static_cast<std::streamsize>(n_factor * sizeof(float)));
  in.read(reinterpret_cast<char*>(dump.shock_draws.data()),
          static_cast<std::streamsize>(n_shock * sizeof(float)));
  if (!in) throw std::runtime_error("The draw arrays in the oracle dump are short.");
  return dump;
}

void print_frequencies(const Dump& dump) {
  const int J = dump.n_attitudes;
  const int T = dump.n_waves;
  std::cout << "category_counts\n";
  for (int attitude = 0; attitude < J; ++attitude) {
    std::vector<int> count(static_cast<std::size_t>(dump.n_categories[attitude]), 0);
    for (const auto& person : dump.persons) {
      if (person.last_wave < person.first_wave) continue;
      for (int wave = person.first_wave; wave <= person.last_wave; ++wave) {
        const int category = person.category[attitude * T + wave];
        if (category >= 1) ++count[category - 1];
      }
    }
    std::cout << "attitude " << (attitude + 1);
    for (int count_k : count) std::cout << " " << count_k;
    std::cout << "\n";
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const std::string dump_path = argc > 1 ? argv[1] : "oracle_body.bin";
    const std::string obj_path = argc > 2 ? argv[2] : "oracle_obj.txt";
    const int n_threads = argc > 3 ? std::stoi(argv[3]) : 4;
    const Dump dump = read_dump(dump_path);

    std::ifstream obj_in(obj_path);
    std::string obj_text((std::istreambuf_iterator<char>(obj_in)), std::istreambuf_iterator<char>());
    for (char& c : obj_text) {
      if (c == 'D' || c == 'd') c = 'E';
    }
    std::stringstream obj_values(obj_text);
    double fortran_obj = 0.0;
    double fortran_penalty = 0.0;
    obj_values >> fortran_obj >> fortran_penalty;
    if (!obj_values) throw std::runtime_error("Cannot read " + obj_path);

    attitudedyn::ModelSpec spec;
    spec.n_waves = dump.n_waves;
    spec.n_attitudes = dump.n_attitudes;
    spec.n_covariates = dump.n_covariates;
    spec.n_factors = dump.n_factors;
    spec.n_categories = dump.n_categories;
    spec.wave_time = dump.time;
    spec.include_wave_effects = true;
    spec.free_contemporaneous = false;
    spec.fix_last_loading = true;
    spec.fix_intercepts = true;
    spec.intercept_columns = {0};
    const char* attitude_names[] = {"ideology", "military", "welfare",
                                    "immigrants", "party", "environment"};
    const char* covariate_names[] = {"const", "age", "hs", "college", "adv",
                                     "female", "white", "black", "inc_mid", "inc_high"};
    for (int i = 0; i < dump.n_attitudes; ++i) spec.attitude_names.emplace_back(attitude_names[i]);
    for (int i = 0; i < dump.n_covariates; ++i) spec.covariate_names.emplace_back(covariate_names[i]);
    for (int t = 0; t < dump.n_waves; ++t) spec.wave_names.push_back(std::to_string(t + 1));

    // Wave effects fixed in the parameter file are the empty cells.
    spec.empty_wave_effect.assign(static_cast<std::size_t>(dump.n_attitudes * dump.n_waves), 0);
    attitudedyn::Layout probe = attitudedyn::make_layout(spec);
    if (probe.n_full != dump.n_coef) {
      throw std::runtime_error("Layout length " + std::to_string(probe.n_full) +
                               " does not equal Fortran ncoef " + std::to_string(dump.n_coef));
    }
    for (int i = 0; i < probe.n_full; ++i) {
      if (probe.entries[i].block == "wave_effect" && dump.icoef[i] != 0) {
        spec.empty_wave_effect[probe.entries[i].attitude * dump.n_waves + probe.entries[i].index2] = 1;
      }
    }
    attitudedyn::Layout layout = attitudedyn::make_layout(spec);
    int mismatches = 0;
    for (int i = 0; i < layout.n_full; ++i) {
      const bool free = layout.entries[i].is_free != 0;
      const bool fortran_free = dump.icoef[i] == 0;
      if (free != fortran_free) {
        if (mismatches < 8) {
          std::cerr << "free-flag mismatch " << layout.entries[i].name << " cpp " << free
                    << " fortran " << fortran_free << "\n";
        }
        ++mismatches;
      }
    }
    if (mismatches) {
      throw std::runtime_error(std::to_string(mismatches) + " coefficients do not match the Fortran flags.");
    }

    attitudedyn::Problem problem;
    problem.spec = spec;
    problem.layout = layout;
    problem.persons = dump.persons;
    problem.theta_full = dump.coef;
    problem.draws.n_draws = dump.n_draws;
    problem.draws.antithetic = true;
    problem.draws.seed = 1;
    problem.draws.kind = attitudedyn::DrawKind::kPhilox;

    attitudedyn::EvaluateRequest request;
    request.n_threads = n_threads;
    request.factor_draws = dump.factor_draws;
    request.shock_draws = dump.shock_draws;

    request.numerics = attitudedyn::Numerics::kShared;
    const attitudedyn::Evaluation shared =
        attitudedyn::evaluate_full(problem, dump.coef, request);
    request.numerics = attitudedyn::Numerics::kStable;
    const attitudedyn::Evaluation stable =
        attitudedyn::evaluate_full(problem, dump.coef, request);
    request.factor_draws.clear();
    request.shock_draws.clear();
    const attitudedyn::Evaluation philox =
        attitudedyn::evaluate_full(problem, dump.coef, request);

    double person_sum = 0.0;
    int used = 0;
    double max_abs_person = 0.0;
    int worst = 0;
    for (int i = 0; i < dump.n_persons; ++i) {
      person_sum += dump.person_loglik[i];
      if (dump.persons[i].last_wave >= dump.persons[i].first_wave) ++used;
      const double gap = std::abs(shared.person_loglik[i] - dump.person_loglik[i]);
      if (gap >= max_abs_person) {
        max_abs_person = gap;
        worst = i;
      }
    }
    std::cout << "worst_person " << worst << " fortran " << dump.person_loglik[worst]
              << " cpp " << shared.person_loglik[worst] << "\n";
    const auto& wp = dump.persons[worst];
    std::cout << "worst_window " << wp.first_wave << " " << wp.last_wave << "\n";

    std::cout.setf(std::ios::scientific);
    std::cout.precision(16);
    std::cout << "persons_used " << used << "\n";
    std::cout << "fortran_obj " << fortran_obj << "\n";
    std::cout << "fortran_penalty " << fortran_penalty << "\n";
    std::cout << "fortran_person_sum " << person_sum << "\n";
    std::cout << "cpp_shared " << shared.objective << "\n";
    std::cout << "cpp_stable_same_draws " << stable.objective << "\n";
    std::cout << "cpp_stable_philox " << philox.objective << "\n";
    std::cout << "max_abs_person_shared " << max_abs_person << "\n";
    std::cout << "published_target " << -111389.88 << "\n";
    print_frequencies(dump);
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << "\n";
    return 1;
  }
}
