#include "draws.hpp"

#include "attitudedyn/normal.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

// Sobol direction numbers for dimensions 2..1025.
// Copyright (c) 2008, Frances Y. Kuo and Stephen Joe. All rights reserved.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//   * Redistributions of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//   * Redistributions in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in the
//     documentation and/or other materials provided with the distribution.
//   * Neither the names of the copyright holders nor the names of the
//     University of New South Wales and the University of Waikato and its
//     contributors may be used to endorse or promote products derived from
//     this software without specific prior written permission.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
// THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#include "sobol_polynomials.inc"

namespace attitudedyn {
namespace detail {
namespace {

constexpr int kMaxStreams = 1025;

std::uint64_t splitmix64(std::uint64_t x) {
  x += 0x9E3779B97F4A7C15ull;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
  return x ^ (x >> 31);
}

const std::vector<int>& primes() {
  static std::vector<int> table = [] {
    const int limit = 20000;
    std::vector<char> composite(limit, 0);
    std::vector<int> out;
    for (int i = 2; i < limit && static_cast<int>(out.size()) < kMaxStreams; ++i) {
      if (composite[i]) continue;
      out.push_back(i);
      if (i > 1 && i < limit / i) {
        for (int j = i * i; j < limit; j += i) composite[j] = 1;
      }
    }
    if (static_cast<int>(out.size()) < kMaxStreams) {
      throw std::logic_error("prime table is shorter than the Sobol table");
    }
    return out;
  }();
  return table;
}

// Direction numbers V[stream * 33 + bit], bit = 1..32, scaled by 2^32.
// Initialised from Joe and Kuo's recurrence (the same one as their sobol.cc).
const std::vector<std::uint32_t>& sobol_directions() {
  static std::vector<std::uint32_t> V = [] {
    std::vector<std::uint32_t> v(static_cast<std::size_t>(kSobolDimensions) * 33, 0);
    for (int bit = 1; bit <= 32; ++bit) {
      v[bit] = 1u << (32 - bit);
    }
    for (int row = 0; row < 1024; ++row) {
      const int stream = row + 1;
      const int s = kSobolPoly[row].s;
      const unsigned a = kSobolPoly[row].a;
      auto at = [&](int bit) -> std::uint32_t& {
        return v[static_cast<std::size_t>(stream) * 33 + bit];
      };
      for (int bit = 1; bit <= s; ++bit) {
        at(bit) = kSobolPoly[row].m[bit - 1] << (32 - bit);
      }
      for (int bit = s + 1; bit <= 32; ++bit) {
        at(bit) = at(bit - s) ^ (at(bit - s) >> s);
        for (int k = 1; k <= s - 1; ++k) {
          const unsigned use = (a >> (s - 1 - k)) & 1u;
          if (use) at(bit) ^= at(bit - k);
        }
      }
    }
    return v;
  }();
  return V;
}

std::uint32_t philox4x32(std::uint64_t seed, int person, int draw_index, int stream) {
  // Philox-4x32-10 (Salmon, Moraes, Dror and Shaw, SC11).
  // Counter words are the person, the draw and the stream, so the value does
  // not depend on which thread evaluates the person.
  const std::uint32_t mult0 = 0xD2511F53u;
  const std::uint32_t mult1 = 0xCD9E8D57u;
  const std::uint32_t bump0 = 0x9E3779B9u;
  const std::uint32_t bump1 = 0xBB67AE85u;
  auto mulhilo = [](std::uint32_t a, std::uint32_t b, std::uint32_t& hi, std::uint32_t& lo) {
    const std::uint64_t product = static_cast<std::uint64_t>(a) * static_cast<std::uint64_t>(b);
    lo = static_cast<std::uint32_t>(product);
    hi = static_cast<std::uint32_t>(product >> 32);
  };
  std::uint32_t r0 = static_cast<std::uint32_t>(person);
  std::uint32_t r1 = static_cast<std::uint32_t>(static_cast<unsigned>(person) >> 16);
  std::uint32_t r2 = static_cast<std::uint32_t>(draw_index);
  std::uint32_t r3 = static_cast<std::uint32_t>(stream);
  std::uint32_t k0 = static_cast<std::uint32_t>(seed);
  std::uint32_t k1 = static_cast<std::uint32_t>(seed >> 32);
  for (int round = 0; round < 10; ++round) {
    std::uint32_t hi0, lo0, hi1, lo1;
    mulhilo(mult0, r0, hi0, lo0);
    mulhilo(mult1, r2, hi1, lo1);
    const std::uint32_t n0 = hi1 ^ r1 ^ k0;
    const std::uint32_t n1 = lo1;
    const std::uint32_t n2 = hi0 ^ r3 ^ k1;
    const std::uint32_t n3 = lo0;
    r0 = n0;
    r1 = n1;
    r2 = n2;
    r3 = n3;
    k0 += bump0;
    k1 += bump1;
  }
  return r0;
}

double uniform_from_bits(std::uint32_t bits) {
  return (static_cast<double>(bits) + 0.5) / 4294967296.0;
}

}  // namespace

double halton_point(std::uint64_t index, int dimension) {
  if (dimension < 0 || dimension >= kMaxStreams) {
    throw std::invalid_argument("Halton dimension is outside 0..1024");
  }
  if (index == 0) return 0.0;
  const int base = primes()[dimension];
  double inverse = 1.0 / static_cast<double>(base);
  double result = 0.0;
  std::uint64_t n = index;
  while (n > 0) {
    result += static_cast<double>(n % static_cast<std::uint64_t>(base)) * inverse;
    n /= static_cast<std::uint64_t>(base);
    inverse /= static_cast<double>(base);
  }
  return result;
}

double sobol_point(int dimension, std::uint32_t index) {
  if (dimension < 0 || dimension >= kSobolDimensions) {
    throw std::invalid_argument(
        "Sobol dimension exceeds 1024. Reduce the number of attitudes, waves, or factors.");
  }
  const auto& V = sobol_directions();
  const std::uint32_t gray = index ^ (index >> 1);
  std::uint32_t point = 0;
  std::uint32_t bits = gray;
  int bit = 1;
  while (bits != 0) {
    if (bits & 1u) {
      point ^= V[static_cast<std::size_t>(dimension) * 33 + bit];
    }
    bits >>= 1;
    ++bit;
  }
  return static_cast<double>(point) / 4294967296.0;
}

namespace {

// Owen scramble of a base-`base` radical inverse. The permutation of each
// digit is an affine bijection (base is prime) whose coefficients depend on
// the digits already emitted, which is the nested uniform scramble.
double owen_halton(std::uint64_t index, int dimension, std::uint64_t seed) {
  const int base = primes()[dimension];
  int digits[64];
  int n_digits = 0;
  std::uint64_t n = index;
  while (n > 0 && n_digits < 64) {
    digits[n_digits++] = static_cast<int>(n % static_cast<std::uint64_t>(base));
    n /= static_cast<std::uint64_t>(base);
  }
  // Pad so the trailing zeros of a short expansion are scrambled too.
  const int length = n_digits < 16 ? 16 : n_digits;
  double place = 1.0;
  double result = 0.0;
  std::uint64_t state = splitmix64(seed ^ (static_cast<std::uint64_t>(dimension) * 0xD1B54A32D192ED03ull));
  for (int i = length - 1; i >= 0; --i) {
    place /= static_cast<double>(base);
    const int digit = (i < n_digits) ? digits[i] : 0;
    const std::uint64_t h = splitmix64(state);
    const int slope = static_cast<int>(h % static_cast<std::uint64_t>(base - 1)) + 1;
    const int shift = static_cast<int>(splitmix64(h) % static_cast<std::uint64_t>(base));
    const int permuted =
        static_cast<int>((static_cast<long long>(slope) * digit + shift) % base);
    result += static_cast<double>(permuted) * place;
    state = splitmix64(state ^ (static_cast<std::uint64_t>(permuted) << 17) ^
                       static_cast<std::uint64_t>(i));
  }
  // Keep the point inside (0, 1). An exact 0 or 1 would make Phi^{-1} infinite.
  if (result <= 0.0) result = 0.5 / std::pow(static_cast<double>(base), length);
  if (result >= 1.0) result = 1.0 - 0.5 / std::pow(static_cast<double>(base), length);
  return result;
}

std::uint32_t owen_base2(std::uint32_t x, std::uint64_t seed) {
  std::uint32_t out = 0;
  std::uint32_t prefix = 0;
  for (int bit = 31; bit >= 0; --bit) {
    std::uint32_t value = (x >> bit) & 1u;
    const std::uint64_t h = splitmix64(seed ^ (static_cast<std::uint64_t>(prefix) << 8) ^
                                        static_cast<std::uint64_t>(bit));
    value ^= static_cast<std::uint32_t>(h & 1u);
    out |= value << bit;
    prefix = (prefix << 1) | value;
  }
  return out;
}

double owen_sobol_uniform(int dimension, std::uint32_t index, std::uint64_t seed) {
  const auto& V = sobol_directions();
  const std::uint32_t gray = index ^ (index >> 1);
  std::uint32_t point = 0;
  std::uint32_t bits_left = gray;
  int bit = 1;
  while (bits_left != 0) {
    if (bits_left & 1u) point ^= V[static_cast<std::size_t>(dimension) * 33 + bit];
    bits_left >>= 1;
    ++bit;
  }
  const std::uint64_t stream_seed =
      splitmix64(seed ^ (static_cast<std::uint64_t>(dimension) * 0x94D049BB133111EBull));
  return uniform_from_bits(owen_base2(point, stream_seed));
}

}  // namespace

double standard_normal_draw(std::uint64_t seed, int person, int draw_index, int stream,
                            bool negate, DrawKind kind) {
  if (stream < 0 || stream >= kMaxStreams) {
    throw std::invalid_argument(
        "A draw stream index is outside 0..1024. The panel has more latent "
        "dimensions (factors + attitudes x waves) than the quasi-random tables cover.");
  }
  double uniform = 0.5;
  switch (kind) {
    case DrawKind::kPhilox: {
      // Mix the person into the counter. `person` is non-negative.
      const std::uint32_t bits = philox4x32(seed, person, draw_index, stream);
      uniform = uniform_from_bits(bits);
      break;
    }
    case DrawKind::kHalton: {
      // Index starts at 1. The scramble seed depends on the person, so the
      // panels are independent while each person's draws stay a scrambled net.
      const std::uint64_t person_seed =
          splitmix64(seed ^ (static_cast<std::uint64_t>(person) * 0xBF58476D1CE4E5B9ull));
      uniform = owen_halton(static_cast<std::uint64_t>(draw_index) + 1ull, stream, person_seed);
      break;
    }
    case DrawKind::kSobol: {
      const std::uint64_t person_seed =
          splitmix64(seed ^ (static_cast<std::uint64_t>(person) * 0x94D049BB133111EBull) ^
                     0xD1B54A32D192ED03ull);
      // Skip the origin (index 0).
      uniform = owen_sobol_uniform(stream, static_cast<std::uint32_t>(draw_index) + 1u, person_seed);
      break;
    }
  }
  double normal = standard_normal_quantile(uniform);
  if (negate) normal = -normal;
  return normal;
}

}  // namespace detail
}  // namespace attitudedyn
