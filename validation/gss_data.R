# Read gss2.prn the way datget3 does.
#
# People stay in first-appearance order, including people the likelihood
# skips, because the random-number streams are keyed by that index.
# A wave inside a person's window that was not in the file keeps a zero
# regressor, except for age, which is filled from the first wave whose
# age is at least 16: age_k = age_j + 2 * (k - j), with k and j the wave
# indices 1..8.

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
