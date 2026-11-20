// Double-precision Phi and Phi^{-1} shared with the C++ likelihood.
// The Fortran oracle calls these when ATTDYN_SHARED_MATH=1, so the
// comparison is not mixed up with anormc's single-precision erfc or
// with mdnris. gfortran appends one underscore.

#include "attitudedyn/normal.hpp"

extern "C" void attdyn_cdf_(double* x, double* y) {
  *y = attitudedyn::standard_normal_cdf(*x);
}

extern "C" void attdyn_qnorm_(double* u, double* z) {
  *z = attitudedyn::standard_normal_quantile(*u);
}
