## Test environments

* Local Ubuntu 24.04, R 4.3.3, `R CMD check --as-cran` on the built tarball.
* GitHub Actions (rhub / win-builder style): macOS and Windows on R release, Ubuntu on R release and R devel. Each run is `R CMD check --as-cran`. The PDF manual is built in the local check. The GitHub Actions runs pass `--no-manual` because those runners do not carry a full TeX tree.

## R CMD check results

There were no ERRORs or WARNINGs.

The incoming-feasibility check reports one NOTE: this is a new submission. That NOTE also prints `LICENSE`, which is the usual check for `License: MIT + file LICENSE`.

Two further NOTEs appear only on this Ubuntu build of R, and not from the package:

* `unable to verify current time`. The clock check could not reach a time server from this machine.
* Compilation used `-mno-omit-leaf-frame-pointer`. Ubuntu's `Makeconf` puts that flag in `CXX17FLAGS`. `src/Makevars` does not set it. CRAN's toolchain does not add it, and the GitHub Actions R builds do not either.

The C++ specification check is OK. It prints the standard line that not all R platforms support C++17, which is what `SystemRequirements: C++17` asks for. The HTML version of the manual was checked with HTML Tidy and reported no problems.

Examples for every exported function finish in well under a second. The spelling word list is clean. `DESCRIPTION` cites Nepaul and Stern (2026) as `<doi:10.1017/pan.2026.10034>`. The project site, the article, and the replication-package URLs in `DESCRIPTION` resolve.

## This is a new submission

No previous CRAN submission of `attitudedyn`.
