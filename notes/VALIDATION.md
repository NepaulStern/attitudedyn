# Validation log

The comparison is the GSS application at the free coefficients printed in the replication log immediately before `iteration = 8` (164 numbers). The log's objective at that point is **-111389.88**. That run was compiled with ifx. The driver in this repository is gfortran, and the package's own generator is Philox, so those three numbers are not the same and are not expected to be.

People in the file: 11282. People in the likelihood: 11252. Category counts match the log:

| Attitude | Counts |
| --- | --- |
| ideology | 865, 2668, 2479, 7762, 2933, 3333, 870 |
| military | 3566, 3893, 2998 |
| welfare | 4180, 3746, 2461 |
| immigrants | 2907, 2963, 4898, 1354, 639 |
| party | 3771, 3585, 2767, 4145, 2067, 2816, 2300 |
| environment | 890, 2954, 6655 |

The penalty is 0. One wave effect is fixed because it has no observations: immigrants in wave 1.

## Same draws, shared normal cdf

The comparison binary and the C++ core were given the same draws and the same double-precision normal cdf and quantile. The person loop in the comparison binary is serial. `ulimit -s unlimited` is required: the data reader puts a large array on the stack.

| Quantity | Value |
| --- | --- |
| Fortran objective | -111659.08497174930 |
| C++ objective, shared numerics | -111659.08497174931 |
| Largest absolute gap, one person | 8.53e-14 |
| C++ objective, stable numerics, same draws | -111661.53451460011 |

## Original single-precision normal cdf

The same driver, with the program's own single-precision cdf and quantile, produces its own draws and its own objective.

| Quantity | Value |
| --- | --- |
| Fortran objective | -111591.85211179820 |
| C++ shared numerics on those draws | -111591.84769206618 |
| Largest absolute gap, one person | 4.39e-3 |
| C++ stable numerics on those draws | -111593.88618143339 |

The original program evaluates the normal cdf in single precision and blends the two tails between -4 and -5. The package default is a double-precision log of the interval probability. The shared row above puts the same double-precision functions on both sides of the comparison.

## Package generator

Philox, seed 1, 10 antithetic draws, stable numerics, the same coefficients and the same prepared sample:

```text
-111570.97022260929
```

`Rscript validation/run_gss.R` prints that value. The C++ thread count does not change it. A gfortran OpenMP build of the original person loop does not reproduce the serial Fortran sum; leave the comparison driver serial. The ifx log used 96 threads and a different random-number stream, which is why its objective is -111389.88.

## Commands on a 20-core Ubuntu machine

R 4.4 from micromamba, from the repository root:

```bash
micromamba create -y -n attitudedyn -c conda-forge \
  "r-base=4.4" r-rcpp r-cli r-broom r-testthat r-knitr r-rmarkdown \
  r-roxygen2 gfortran gxx make
micromamba activate attitudedyn

R CMD INSTALL .

cd validation/fortran
make
mkdir -p run && cd run
ln -s ../data data
cp ../build/alecia.p alecia.p
ulimit -s unlimited
OMP_NUM_THREADS=1 ../build/alecia_oracle
ATTDYN_SHARED_MATH=1 OMP_NUM_THREADS=1 ../build/alecia_oracle
../build/compare_gss oracle_body.bin oracle_obj.txt 20

cd ../../..
Rscript validation/run_gss.R --threads 20
```

`validation/run_gss.R --refit` starts BFGS from the shipped coefficients. That is optional. The check above does not refit.

The core unit tests and the R tests use small simulated panels:

```bash
make -C core test
Rscript -e 'testthat::test_local(".")'
```
