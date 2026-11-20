# Read gss2.prn the way datget3 does.
#
# People stay in first-appearance order, including people the likelihood
# skips, because the random-number streams are keyed by that index.
# A wave inside a person's window that was not in the file keeps a zero
# regressor, except for age, which is filled from the first wave whose
# age is at least 16: age_k = age_j + 2 * (k - j), with k and j the wave
# indices 1..8.

# Download gss2.prn, the replication log, and the Fortran comparison sources
# from the CC0 Harvard Dataverse deposit doi:10.7910/DVN/DUCHYS into
# validation/fortran/vendor/, which is git-ignored.
ensure_replication_files <- function(root) {
  vendor <- file.path(root, "validation", "fortran", "vendor")
  needed <- c(
    file.path(vendor, "data", "gss2.prn"),
    file.path(vendor, "gss_log.txt"),
    file.path(vendor, "alecia_parallel.f"),
    file.path(vendor, "eigrs.f"),
    file.path(vendor, "opt.f")
  )
  script <- file.path(root, "validation", "fortran", "fetch_vendor.py")
  status <- system2("python3", script)
  if (!identical(status, 0L) || any(!file.exists(needed))) {
    stop(
      "Could not download the replication files from Harvard Dataverse (doi:10.7910/DVN/DUCHYS).",
      call. = FALSE
    )
  }
  list(
    data = file.path(vendor, "data", "gss2.prn"),
    log = file.path(vendor, "gss_log.txt"),
    vendor = vendor
  )
}

gss_nseg <- function() c(7L, 3L, 3L, 5L, 7L, 3L)

gss_time_grid <- function() c(2, 3, 4, 5, 6, 7, 8, 10)

gss_frequencies <- function() {
  list(
    ideology = c(865L, 2668L, 2479L, 7762L, 2933L, 3333L, 870L),
    military = c(3566L, 3893L, 2998L),
    welfare = c(4180L, 3746L, 2461L),
    immigrants = c(2907L, 2963L, 4898L, 1354L, 639L),
    party = c(3771L, 3585L, 2767L, 4145L, 2067L, 2816L, 2300L),
    environment = c(890L, 2954L, 6655L)
  )
}

read_gss <- function(path, n_rows = 21648L) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[grepl("^[[:space:]]*-?[0-9]", lines)]
  if (length(lines) < n_rows) {
    stop("gss2.prn has ", length(lines), " data rows; the program reads ", n_rows, ".")
  }
  raw <- utils::read.table(
    text = lines[seq_len(n_rows)], header = FALSE,
    col.names = c(
      "id", "year", "age", "educ", "ideology", "race", "sex", "panel",
      "idate", "month", "trvdate", "military", "junk1", "environment",
      "welfare", "junk3", "junk4", "immigrants", "party", "income"
    )
  )
  raw$wave <- ((raw$year - 2006L) %/% 2L) + 1L
  if (any(raw$wave < 1L | raw$wave > 8L)) {
    stop("A GSS year falls outside waves 1..8.")
  }
  ids <- unique(raw$id)
  n <- length(ids)
  n_waves <- 8L
  nseg <- gss_nseg()
  index <- match(raw$id, ids)
  y <- array(0L, c(n, 6L, n_waves))
  x <- array(0, c(n, n_waves, 10L))
  for (r in seq_len(nrow(raw))) {
    i <- index[[r]]
    t <- raw$wave[[r]]
    y[i, 1L, t] <- raw$ideology[[r]]
    y[i, 2L, t] <- raw$military[[r]]
    y[i, 3L, t] <- raw$welfare[[r]]
    y[i, 4L, t] <- raw$immigrants[[r]]
    y[i, 5L, t] <- raw$party[[r]] + 1L
    y[i, 6L, t] <- raw$environment[[r]]
    x[i, t, 1L] <- 1
    x[i, t, 2L] <- raw$age[[r]]
    educ <- raw$educ[[r]]
    if (educ == 1L || educ == 2L) x[i, t, 3L] <- 1
    if (educ == 3L) x[i, t, 4L] <- 1
    if (educ == 4L) x[i, t, 5L] <- 1
    x[i, t, 6L] <- raw$sex[[r]] - 1
    if (raw$race[[r]] == 1L) x[i, t, 7L] <- 1
    if (raw$race[[r]] == 2L) x[i, t, 8L] <- 1
    if (raw$income[[r]] == 2L) x[i, t, 9L] <- 1
    if (raw$income[[r]] == 3L) x[i, t, 10L] <- 1
  }

  included <- rep(FALSE, n)
  for (i in seq_len(n)) {
    for (t in seq_len(n_waves)) {
      for (j in seq_len(6L)) {
        if (y[i, j, t] >= 1L && y[i, j, t] <= nseg[[j]]) included[[i]] <- TRUE
      }
    }
  }
  for (i in seq_len(n)) {
    donor <- 0L
    for (t in seq_len(n_waves)) {
      if (x[i, t, 2L] >= 16) {
        donor <- t
        break
      }
    }
    if (donor == 0L) {
      included[[i]] <- FALSE
    } else {
      age0 <- x[i, donor, 2L]
      for (k in seq_len(n_waves)) {
        if (x[i, k, 2L] < 16) x[i, k, 2L] <- age0 + 2 * (k - donor)
      }
    }
  }
  for (i in seq_len(n)) {
    donor <- 0L
    for (t in seq_len(n_waves)) {
      female <- x[i, t, 6L]
      if (female >= 0 && female <= 1) {
        donor <- t
        break
      }
    }
    if (donor == 0L) {
      included[[i]] <- FALSE
    } else {
      female <- x[i, donor, 6L]
      for (k in seq_len(n_waves)) {
        if (x[i, k, 6L] < 0) x[i, k, 6L] <- female
      }
    }
  }

  start <- integer(n)
  end <- integer(n)
  for (i in seq_len(n)) {
    if (!included[[i]]) next
    for (t in seq_len(n_waves)) {
      if (any(y[i, , t] > 0L)) {
        if (start[[i]] == 0L) start[[i]] <- t
        end[[i]] <- t
      }
    }
    if (start[[i]] == 0L) included[[i]] <- FALSE
  }

  rows <- vector("list", n)
  for (i in seq_len(n)) {
    if (!included[[i]]) {
      rows[[i]] <- data.frame(
        id = ids[[i]], wave = 1L, weight = 0,
        `(Intercept)` = 0, age = 0, hs = 0, college = 0, adv = 0,
        female = 0, white = 0, black = 0, inc_mid = 0, inc_high = 0,
        ideology = NA_integer_, military = NA_integer_, welfare = NA_integer_,
        immigrants = NA_integer_, party = NA_integer_, environment = NA_integer_,
        check.names = FALSE
      )
      next
    }
    span <- start[[i]]:end[[i]]
    block <- data.frame(
      id = ids[[i]], wave = span, weight = 1,
      `(Intercept)` = x[i, span, 1L],
      age = x[i, span, 2L],
      hs = x[i, span, 3L],
      college = x[i, span, 4L],
      adv = x[i, span, 5L],
      female = x[i, span, 6L],
      white = x[i, span, 7L],
      black = x[i, span, 8L],
      inc_mid = x[i, span, 9L],
      inc_high = x[i, span, 10L],
      check.names = FALSE
    )
    outcomes <- c("ideology", "military", "welfare", "immigrants", "party", "environment")
    for (j in seq_along(outcomes)) {
      values <- y[i, j, span]
      values[values < 1L | values > nseg[[j]]] <- NA_integer_
      block[[outcomes[[j]]]] <- values
    }
    rows[[i]] <- block
  }
  panel <- do.call(rbind, rows)
  rownames(panel) <- NULL
  attr(panel, "n_people") <- n
  attr(panel, "n_included") <- sum(included)
  panel
}

gss_category_counts <- function(panel) {
  outcomes <- c("ideology", "military", "welfare", "immigrants", "party", "environment")
  nseg <- gss_nseg()
  counts <- vector("list", length(outcomes))
  names(counts) <- outcomes
  keep <- panel$weight > 0
  for (j in seq_along(outcomes)) {
    values <- panel[[outcomes[[j]]]][keep]
    values <- values[!is.na(values)]
    counts[[j]] <- tabulate(values, nbins = nseg[[j]])
  }
  counts
}

# Fortran name of one layout row. Attitude and index2 are the 1-based
# values stored on the fitted layout, except threshold cuts: those are
# already 1-based in the coefficient name (tau_3 is the first free cut).
gss_fortran_name <- function(name, block, attitude, index2) {
  if (block == "wave_effect") return(sprintf("fixef%d%02d", attitude, index2))
  if (block == "beta") {
    stem <- c(
      "constnt", "age", "hs", "college", "advdegr", "female", "white", "black"
    )
    if (index2 <= length(stem)) return(sprintf("%s%d", stem[[index2]], attitude))
    if (index2 == 9L) return(sprintf("income%d1", attitude))
    if (index2 == 10L) return(sprintf("income%d2", attitude))
  }
  if (block == "persistence") return(sprintf("rho%d%d", attitude, index2))
  if (block == "loading") return(sprintf("alamb%d%d", index2, attitude))
  if (block == "threshold") {
    cut <- as.integer(sub(".*,", "", sub("\\]$", "", name)))
    return(sprintf("tau%d%d", attitude, cut))
  }
  if (block == "chol") return(sprintf("chol%d%d", attitude, index2))
  stop("No Fortran name for block ", block, ".")
}

gss_fortran_names <- function(layout) {
  free <- layout[layout$is_free == 1L, , drop = FALSE]
  vapply(seq_len(nrow(free)), function(i) {
    gss_fortran_name(free$name[[i]], free$block[[i]], free$attitude[[i]], free$index2[[i]])
  }, character(1))
}

# rho = 2 exp(c) / (1 + exp(c)) - 1, the map used for the published GSS
# coefficients, which are well inside the range where this expression and
# the large-|c| form agree.
gss_natural_persistence <- function(unconstrained) {
  2 * exp(unconstrained) / (1 + exp(unconstrained)) - 1
}

# Free coefficients from the last X VECTOR in the replication log, the
# vector printed just before iteration 8.
read_gss_estimates <- function(log_path) {
  lines <- readLines(log_path, warn = FALSE)
  mark <- grep("iteration = 8", lines)
  if (!length(mark)) stop("The log has no iteration 8 line.")
  mark <- mark[[1L]]
  header <- grep("X VECTOR", lines)
  header <- header[header < mark]
  if (!length(header)) stop("The log has no X VECTOR before iteration 8.")
  header <- max(header)
  block <- lines[(header + 1L):(mark - 1L)]
  block <- trimws(block)
  block <- block[nzchar(block)]
  parts <- strsplit(block, "[[:space:]]+")
  name <- vapply(parts, `[[`, character(1), 1L)
  value <- vapply(parts, function(p) as.numeric(gsub("D", "E", p[[2L]], fixed = TRUE)), numeric(1))
  stats::setNames(value, name)
}

# Point estimates and standard errors from the "estimates and significance"
# table printed after the last iteration.
read_gss_estimate_table <- function(log_path) {
  lines <- readLines(log_path, warn = FALSE)
  start <- grep("estimates and significance", lines)
  if (!length(start)) stop("The log has no estimates table.")
  start <- start[[length(start)]]
  rows <- list()
  for (line in lines[(start + 1L):length(lines)]) {
    match <- regexec(
      "^\\s*[0-9]+\\s+([A-Za-z0-9]+)\\s+(-?[0-9.]+)\\**\\s+([0-9.]+)\\s+(-?[0-9.]+)",
      line
    )
    found <- regmatches(line, match)[[1L]]
    if (!length(found)) {
      if (length(rows)) break
      next
    }
    rows[[length(rows) + 1L]] <- data.frame(
      fortran = found[[2L]],
      published = as.numeric(found[[3L]]),
      published_se = as.numeric(found[[4L]]),
      published_t = as.numeric(found[[5L]]),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) stop("The estimates table was empty.")
  do.call(rbind, rows)
}

# Six variance blocks, in attitude order, from the end of the replication log.
read_gss_variance <- function(log_path) {
  lines <- readLines(log_path, warn = FALSE)
  start <- grep("variance decomposition for each dependent variable", lines)
  if (!length(start)) stop("The log has no variance decomposition.")
  start <- start[[length(start)]]
  body <- lines[(start + 1L):length(lines)]
  numbers <- gregexpr("-?[0-9]+\\.[0-9]+", body)
  values <- regmatches(body, numbers)
  # Each block contributes x-beta, delta, lambda, omega (size and share) and total.
  sizes <- list()
  cursor <- 1L
  for (line_values in values) {
    if (!length(line_values)) next
    if (length(line_values) == 2L) {
      sizes[[length(sizes) + 1L]] <- as.numeric(line_values[[1L]])
    } else if (length(line_values) == 1L && length(sizes) %% 5L == 4L) {
      sizes[[length(sizes) + 1L]] <- as.numeric(line_values[[1L]])
    }
  }
  if (length(sizes) < 30L) stop("Expected 6 variance blocks, found ", length(sizes), " numbers.")
  sizes <- unlist(sizes)[seq_len(30L)]
  data.frame(
    attitude_index = rep(seq_len(6L), each = 5L),
    component = rep(c("covariate", "wave_effect", "factor", "shock", "total"), 6L),
    published = sizes,
    stringsAsFactors = FALSE
  )
}

gss_comparison <- function(fit, log_path) {
  layout_names <- gss_fortran_names(fit$layout)
  shipped <- read_gss_estimates(log_path)
  if (!identical(layout_names, names(shipped))) {
    stop(
      "Free-coefficient order does not match the replication log.\n",
      "First layout name: ", layout_names[[1L]], "\n",
      "First log name: ", names(shipped)[[1L]]
    )
  }
  table <- read_gss_estimate_table(log_path)
  if (!identical(table$fortran, names(shipped))) {
    stop("The estimates table is not in the same order as the X VECTOR.")
  }
  estimates <- coef(fit)
  se <- tryCatch(sqrt(diag(vcov(fit, type = "opg_centered"))), error = function(e) {
    rep(NA_real_, length(estimates))
  })
  coefficients <- data.frame(
    fortran = names(shipped),
    term = names(estimates),
    estimate = as.numeric(estimates),
    published = as.numeric(shipped),
    difference = as.numeric(estimates) - as.numeric(shipped),
    std.error = as.numeric(se),
    published_se = table$published_se,
    stringsAsFactors = FALSE
  )
  coefficients$se_difference <- coefficients$std.error - coefficients$published_se

  rho_rows <- fit$layout$block == "persistence" & fit$layout$is_free == 1L
  # Persistence is free for every entry, so the natural matrix is the
  # transform of the unconstrained coefficients in receiver-major order.
  published_rho <- matrix(gss_natural_persistence(shipped[grep("^rho", names(shipped))]), 6L, 6L, byrow = TRUE)
  fitted_rho <- persistence_matrix(fit)
  dimnames(published_rho) <- dimnames(fitted_rho)
  if (!identical(rownames(fitted_rho), fit$attitude_names)) {
    stop("Persistence matrix is not in attitude order.")
  }

  published_variance <- read_gss_variance(log_path)
  fitted_variance <- variance_decomposition(fit)
  variance <- data.frame(
    attitude = rep(fit$attitude_names, each = 5L),
    component = published_variance$component,
    estimate = NA_real_,
    published = published_variance$published,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(fitted_variance))) {
    row <- fitted_variance[i, ]
    for (component in c("covariate", "wave_effect", "factor", "shock", "total")) {
      slot <- variance$attitude == row$attitude & variance$component == component
      variance$estimate[slot] <- row[[component]]
    }
  }
  variance$difference <- variance$estimate - variance$published

  list(
    coefficients = coefficients,
    persistence = fitted_rho,
    published_persistence = published_rho,
    variance = variance,
    loglik = as.numeric(fit$loglik),
    convergence = fit$convergence
  )
}

print_gss_comparison <- function(comparison, elapsed = NA_real_) {
  coefs <- comparison$coefficients
  abs_diff <- abs(coefs$difference)
  se_diff <- abs(coefs$se_difference)
  rho_diff <- abs(comparison$persistence - comparison$published_persistence)
  var_diff <- abs(comparison$variance$difference)
  cat(sprintf("convergence %d\n", comparison$convergence))
  cat(sprintf("elapsed_seconds %.3f\n", elapsed))
  cat(sprintf("loglik %.16f\n", comparison$loglik))
  cat(sprintf("coefficient_abs_diff_max %.6g\n", max(abs_diff)))
  cat(sprintf("coefficient_abs_diff_median %.6g\n", stats::median(abs_diff)))
  cat(sprintf("se_abs_diff_max %.6g\n", max(se_diff, na.rm = TRUE)))
  cat(sprintf("se_abs_diff_median %.6g\n", stats::median(se_diff, na.rm = TRUE)))
  cat(sprintf("persistence_abs_diff_max %.6g\n", max(rho_diff)))
  cat(sprintf("persistence_abs_diff_median %.6g\n", stats::median(rho_diff)))
  cat(sprintf("variance_abs_diff_max %.6g\n", max(var_diff, na.rm = TRUE)))
  cat(sprintf("variance_abs_diff_median %.6g\n", stats::median(var_diff, na.rm = TRUE)))
  cat("\nCoefficients (unconstrained), refit minus the replication X VECTOR:\n")
  show <- data.frame(
    fortran = coefs$fortran,
    term = coefs$term,
    estimate = round(coefs$estimate, 6),
    published = round(coefs$published, 6),
    difference = round(coefs$difference, 6),
    std.error = round(coefs$std.error, 4),
    published_se = round(coefs$published_se, 4),
    stringsAsFactors = FALSE
  )
  print(show, row.names = FALSE, right = FALSE)
  cat("\nPersistence on (-1, 1), refit minus the transformed replication rho:\n")
  print(round(comparison$persistence - comparison$published_persistence, 4))
  cat("\nVariance decomposition, refit minus the replication log:\n")
  print(data.frame(
    attitude = comparison$variance$attitude,
    component = comparison$variance$component,
    estimate = round(comparison$variance$estimate, 4),
    published = round(comparison$variance$published, 4),
    difference = round(comparison$variance$difference, 4),
    stringsAsFactors = FALSE
  ), row.names = FALSE)
  invisible(comparison)
}
