column_arg <- function(value, argument, data) {
  if (inherits(value, "formula")) {
    vars <- all.vars(value)
    if (length(vars) != 1L) {
      cli::cli_abort("{.arg {argument}} must name exactly one column.")
    }
    return(vars[[1L]])
  }
  if (is.character(value) && length(value) == 1L) return(value)
  cli::cli_abort("{.arg {argument}} must be a one-sided formula, such as {.code ~ id}, or a column name.")
}

#' Turn a long panel into the arrays the likelihood uses.
#'
#' One row is one person in one wave. A missing outcome is a missing item and
#' is integrated out. A wave that does not appear for a person, but sits
#' between that person's first and last observed waves, is a gap: the latent
#' state is simulated through it. Covariates on that wave are used if the
#' row was read and are zero if the wave is absent, which is what the
#' Fortran loader does.
#'
#' Categories are the integers `1, ..., K`. The cut between category 1 and
#' category 2 is normalised at zero, so the numbering is part of the model.
#'
#' @param formula A [formula] with [cbind()] on the left-hand side.
#' @param data A data frame in long form.
#' @param id,wave Column holding the person id and the wave, as a formula or a name.
#' @param time Optional numeric column of interview time. Values must be
#'   constant within a wave; the common value is the wave date.
#' @param time_grid Optional numeric vector of wave dates, in order of the
#'   sorted unique wave values. This is the GSS grid
#'   `c(2, 3, 4, 5, 6, 7, 8, 10)`.
#' @param weight_vector Optional numeric weights, already aligned with `data`.
#' @param na.action What to do with missing covariates. Missing outcomes are
#'   never dropped. `na.omit` drops person-waves whose covariates are missing.
#' @param control An [attitudedyn_control()] list. Only the normalisation
#'   flags are used here.
#' @return A list of inputs for the C++ problem builder. Not exported.
#' @keywords internal
prepare_attitudedyn <- function(formula, data, id, wave, time = NULL, time_grid = NULL,
                                weight_vector = NULL, na.action = na.pass,
                                control = attitudedyn_control()) {
  if (!is.data.frame(data)) cli::cli_abort("{.arg data} must be a data frame.")
  if (!inherits(formula, "formula")) cli::cli_abort("{.arg formula} must be a formula.")
  id_name <- column_arg(id, "id", data)
  wave_name <- column_arg(wave, "wave", data)
  if (!id_name %in% names(data)) cli::cli_abort("Column {.val {id_name}} is not in {.arg data}.")
  if (!wave_name %in% names(data)) cli::cli_abort("Column {.val {wave_name}} is not in {.arg data}.")

  if (!is.null(weight_vector) && length(weight_vector) != nrow(data)) {
    cli::cli_abort("{.arg weights} must have one value per row of {.arg data}.")
  }

  time_vector <- NULL
  if (!missing(time) && !is.null(time)) {
    time_name <- column_arg(time, "time", data)
    if (!time_name %in% names(data)) cli::cli_abort("Column {.val {time_name}} is not in {.arg data}.")
    time_vector <- data[[time_name]]
  }

  frame <- stats::model.frame(formula, data = data, na.action = stats::na.pass, drop.unused.levels = TRUE)
  response <- stats::model.response(frame, type = "any")
  design <- stats::model.matrix(formula, frame)
  if (is.null(dim(response))) response <- cbind(response)
  if (!is.matrix(response)) {
    cli::cli_abort("The left-hand side must be one outcome or {.fn cbind} of outcomes.")
  }
  outcome_names <- colnames(response)
  if (is.null(outcome_names)) outcome_names <- paste0("y", seq_len(ncol(response)))

  id_vector <- data[[id_name]]
  wave_vector <- data[[wave_name]]
  if (nrow(frame) != nrow(data)) {
    cli::cli_abort("The model frame dropped rows unexpectedly. Check for missing variables.")
  }
  if (anyNA(id_vector) || anyNA(wave_vector)) {
    cli::cli_abort("{.arg id} and {.arg wave} cannot contain missing values.")
  }

  covariate_missing <- !stats::complete.cases(design)
  if (any(covariate_missing)) {
    omit <- identical(na.action, stats::na.omit) || identical(na.action, stats::na.exclude)
    if (!omit) {
      cli::cli_abort(c(
        "Covariates contain missing values.",
        "i" = "Missing outcomes are integrated out and can stay missing.",
        "i" = "Use {.code na.action = na.omit} to drop person-waves with missing covariates, or impute them."
      ))
    }
    keep <- !covariate_missing
    design <- design[keep, , drop = FALSE]
    response <- response[keep, , drop = FALSE]
    id_vector <- id_vector[keep]
    wave_vector <- wave_vector[keep]
    if (!is.null(weight_vector)) weight_vector <- weight_vector[keep]
    if (!is.null(time_vector)) time_vector <- time_vector[keep]
  }

  n_attitudes <- ncol(response)
  n_categories <- integer(n_attitudes)
  category_matrix <- matrix(0L, nrow(response), n_attitudes)
  for (j in seq_len(n_attitudes)) {
    column <- response[, j]
    observed <- !is.na(column)
    if (!any(observed)) {
      cli::cli_abort("Outcome {.val {outcome_names[j]}} is entirely missing.")
    }
    values <- column[observed]
    if (is.factor(values) || is.ordered(values)) {
      coded <- as.integer(values)
      if (any(coded < 1L)) cli::cli_abort("Factor outcome {.val {outcome_names[j]}} has an empty level.")
      category_matrix[observed, j] <- coded
      n_categories[j] <- nlevels(response[, j])
    } else {
      if (!is.numeric(values) || any(values != round(values)) || any(values < 1)) {
        cli::cli_abort(c(
          "Outcome {.val {outcome_names[j]}} must be coded as integers 1, 2, ..., K.",
          "i" = "The cut between category 1 and category 2 is normalised at zero."
        ))
      }
      coded <- as.integer(round(values))
      k <- max(coded)
      if (!setequal(unique(coded), seq_len(k))) {
        cli::cli_abort(c(
          "Outcome {.val {outcome_names[j]}} skips a category code between 1 and {k}.",
          "i" = "Observed codes: {.val {sort(unique(coded))}}."
        ))
      }
      category_matrix[observed, j] <- coded
      n_categories[j] <- k
    }
  }

  wave_values <- sort(unique(wave_vector))
  wave_index <- match(wave_vector, wave_values)
  n_waves <- length(wave_values)
  if (!is.null(time_grid)) {
    if (length(time_grid) != n_waves) {
      cli::cli_abort("{.arg time_grid} must have one date for each distinct wave ({n_waves}).")
    }
    grid <- as.numeric(time_grid)
  } else if (!is.null(time_vector)) {
    grid <- numeric(n_waves)
    for (t in seq_len(n_waves)) {
      slice <- time_vector[wave_index == t]
      if (anyNA(slice) || diff(range(slice)) > 1e-8) {
        cli::cli_abort("Time is missing or not constant within wave {.val {wave_values[t]}}.")
      }
      grid[t] <- mean(slice)
    }
  } else {
    grid <- as.numeric(seq_len(n_waves))
    cli::cli_warn(c(
      "No wave dates were supplied, so gaps are 1 between successive wave codes.",
      "i" = "The GSS estimates use {.code time_grid = c(2, 3, 4, 5, 6, 7, 8, 10)}."
    ))
  }
  if (any(!is.finite(grid)) || any(diff(grid) <= 0)) {
    cli::cli_abort("Wave dates must be finite and strictly increasing.")
  }

  if (any(duplicated(cbind(id_vector, wave_index)))) {
    cli::cli_abort("Some person-wave pairs appear more than once.")
  }

  person_ids <- unique(id_vector)
  n_persons <- length(person_ids)
  n_covariates <- ncol(design)
  interviewed <- matrix(0L, n_persons, n_waves)
  categories <- array(0L, dim = c(n_persons, n_attitudes, n_waves))
  covariates <- array(0, dim = c(n_persons, n_waves, n_covariates))
  person_weight <- rep(1, n_persons)
  first_wave <- integer(n_persons)
  last_wave <- integer(n_persons)

  id_split <- split(seq_len(nrow(design)), factor(id_vector, levels = person_ids))
  for (i in seq_len(n_persons)) {
    rows <- id_split[[i]]
    waves <- wave_index[rows]
    first_wave[i] <- min(waves) - 1L
    last_wave[i] <- max(waves) - 1L
    interviewed[i, waves] <- 1L
    for (r in seq_along(rows)) {
      categories[i, , waves[r]] <- category_matrix[rows[r], ]
      covariates[i, waves[r], ] <- design[rows[r], ]
    }
    if (!is.null(weight_vector)) {
      w <- weight_vector[rows]
      if (any(!is.finite(w)) || any(w < 0) || diff(range(w)) > 1e-10) {
        cli::cli_abort("Weights for person {.val {person_ids[i]}} must be finite, non-negative, and constant across waves.")
      }
      person_weight[i] <- w[1L]
    }
  }

  empty_wave <- integer(n_attitudes * n_waves)
  for (j in seq_len(n_attitudes)) {
    for (t in seq_len(n_waves)) {
      observed <- any(categories[, j, t] >= 1L)
      if (!observed) empty_wave[(j - 1L) * n_waves + t] <- 1L
    }
  }
  intercept_names <- c("(Intercept)", "`(Intercept)`")
  intercept_columns <- which(colnames(design) %in% intercept_names) - 1L

  # Flat storage: person, attitude, wave with wave fastest; covariates with
  # the covariate index fastest.
  category_flat <- integer(n_persons * n_attitudes * n_waves)
  covariate_flat <- numeric(n_persons * n_waves * n_covariates)
  interviewed_flat <- integer(n_persons * n_waves)
  for (i in seq_len(n_persons)) {
    interviewed_flat[((i - 1L) * n_waves + 1L):(i * n_waves)] <- interviewed[i, ]
    for (j in seq_len(n_attitudes)) {
      base <- ((i - 1L) * n_attitudes + (j - 1L)) * n_waves
      category_flat[base + seq_len(n_waves)] <- categories[i, j, ]
    }
    for (t in seq_len(n_waves)) {
      base <- ((i - 1L) * n_waves + (t - 1L)) * n_covariates
      covariate_flat[base + seq_len(n_covariates)] <- covariates[i, t, ]
    }
  }

  list(
    formula = formula,
    terms = stats::terms(frame),
    n_waves = n_waves,
    n_attitudes = n_attitudes,
    n_covariates = n_covariates,
    n_categories = as.integer(n_categories),
    wave_time = grid,
    wave_values = wave_values,
    wave_names = as.character(wave_values),
    attitude_names = outcome_names,
    covariate_names = colnames(design),
    intercept_columns = as.integer(intercept_columns),
    empty_wave_effect = as.integer(empty_wave),
    include_wave_effects = TRUE, # replaced by the caller
    n_persons = n_persons,
    person_ids = person_ids,
    first_wave = as.integer(first_wave),
    last_wave = as.integer(last_wave),
    interviewed = interviewed_flat,
    category = category_flat,
    covariates = covariate_flat,
    weight = person_weight,
    fix_last_loading = control$fix_last_loading,
    fix_intercept = control$fix_intercept
  )
}
