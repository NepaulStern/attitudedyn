# The R package, in plain language

`attitudedyn` is an R package. R is the language a user types. The package is the directory of files R installs and then loads with `library(attitudedyn)`. The numerical work is a C++ library that does not know R exists. This note is for someone who has not built an R package before.

## What you install

A user needs R, a C++17 compiler, and the packages listed under `Imports` in `DESCRIPTION`: Rcpp (the bridge from R to C++), cli (error messages), stats (the optimiser and a few utilities that ship with R), and broom (the `tidy` and `glance` tables). The suggested packages are only for tests and the vignettes.

From a clone of this repository:

```bash
R CMD INSTALL .
```

After that, `library(attitudedyn)` works in any R session. You do not need to stay inside the source directory.

## What each directory is

- `DESCRIPTION` is the package's identity card: name, version, author, licence, and which other packages it needs. R refuses to install a package whose `DESCRIPTION` is missing a required field.
- `NAMESPACE` is the list of names the package offers to users, and the names it borrows from other packages. It is generated. Do not edit it by hand.
- `R/` is the R code users call: `attitudedyn()`, the print and summary methods, simulation, and the post-estimation helpers.
- `src/` is the thin C++ file that translates R vectors into the core's C++ types and back. `Makevars` tells the compiler to use C++17 and to look in `core/`.
- `core/` is the C++ library. It includes only the C++ standard library. A Python or Stata wrapper can link the same files.
- `man/` is the help pages. `?attitudedyn` reads these. They are generated from the comments above each function.
- `tests/testthat/` is the R test suite.
- `vignettes/` are the long-form articles. They are rebuilt when the package is checked.
- `notes/` and `validation/` are for people working on the source. They are not installed with the package.

## How a help page gets written

Above an R function there is a comment that starts with `#'`. Those lines are the help page. `roxygen2` reads them and writes `man/*.Rd` and `NAMESPACE`. After changing a help comment or an exported name:

```bash
Rscript -e 'Rcpp::compileAttributes(); roxygen2::roxygenise()'
```

`compileAttributes` rewrites `src/RcppExports.cpp` and `R/RcppExports.R` from the `// [[Rcpp::export]]` markers in `src/wrapper.cpp`.

## How a fit runs

`attitudedyn()` is ordinary R. It builds a model matrix, packs the panel into a few long vectors, and hands them to the C++ core. The core returns a log-likelihood and, when asked, a gradient. R then calls `stats::optim` (or a short scoring loop) and stores the result as an object of class `"attitudedyn"`. `summary`, `coef`, `vcov`, and the other methods are functions R finds by that class name.

The default number of threads is 1. The sum over people is taken in person order, so a larger `n_threads` changes the speed and not the value.

## How to check your own change

```bash
Rscript -e 'testthat::test_local(".")'
R CMD build .
R CMD check --as-cran attitudedyn_0.1.0.tar.gz
```

`R CMD check` is the test CRAN runs. It rebuilds the package from the tarball, runs the examples, the tests, and the vignettes, and reports errors, warnings, and notes. Check the tarball, not the source directory: the source directory contains git files and build products that the tarball correctly leaves out.

The C++ core has its own test, separate from R:

```bash
make -C core test
```
