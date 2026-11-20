# Find the source tree both from testthat::test_local and from R CMD check.
# During a check the sources live in <pkg>.Rcheck/00_pkg_src/attitudedyn,
# while the tests run under <pkg>.Rcheck/tests/testthat.
package_source_root <- function() {
  candidates <- character()
  here <- normalizePath(getwd(), mustWork = FALSE)
  bases <- unique(c(here, dirname(here), dirname(dirname(here)), dirname(dirname(dirname(here)))))
  for (base in bases) {
    candidates <- c(
      candidates,
      base,
      file.path(base, "00_pkg_src", "attitudedyn"),
      file.path(base, "attitudedyn")
    )
  }
  candidates <- unique(candidates)
  for (candidate in candidates) {
    description <- file.path(candidate, "DESCRIPTION")
    if (!file.exists(description)) next
    fields <- read.dcf(description)
    if (identical(unname(fields[1, "Package"]), "attitudedyn")) {
      return(normalizePath(candidate))
    }
  }
  NULL
}

test_that("Rd files, vignettes, and the description pass the spelling check", {
  skip_if_not_installed("spelling")
  package_root <- package_source_root()
  expect_false(is.null(package_root))
  misspelled <- spelling::spell_check_package(package_root)
  expect_equal(nrow(misspelled), 0, info = paste(misspelled$word, collapse = ", "))
})
