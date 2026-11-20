#ifndef ATTITUDYN_DRAWS_HPP
#define ATTITUDYN_DRAWS_HPP

#include "attitudedyn/api.hpp"

#include <cstdint>

namespace attitudedyn {
namespace detail {

// Unscrambled sequences, used by the unit tests and by the scrambled
// generators. `index` is the Sobol/Halton index starting at 0. Sobol
// index 0 is the origin. Halton index 0 is defined as 0; the radical
// inverse used for Monte Carlo starts at index 1.
double halton_point(std::uint64_t index, int dimension);
double sobol_point(int dimension, std::uint32_t index);

// Standard normal draw for one stream. Antithetic draws negate the normal,
// which is what errget2 does for the multi-attitude model.
double standard_normal_draw(std::uint64_t seed, int person, int draw_index,
                            int stream, bool negate, DrawKind kind);

}  // namespace detail
}  // namespace attitudedyn

#endif
