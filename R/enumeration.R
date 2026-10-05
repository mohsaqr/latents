#' Enumerate numbers of individual profiles and group classes
#'
#' Fits every requested combination and retains errors and warnings alongside
#' successful fits. Information criteria are descriptive: neither the smallest
#' BIC nor high entropy guarantees the correct number of classes. Unconverged,
#' boundary, and unreplicated fits are reported, not silently selected.
#' @param data Data frame.
#' @param vars Continuous indicator names.
#' @param id The column that identifies the groups the observations are
#'   nested in (students in schools, reports in students). `NULL` enumerates
#'   single-level models, with a message; [enumerate_lpa()] and
#'   [enumerate_lca()] fit those by name.
#' @param n_profiles Positive integer profile counts to try.
#' @param n_group_classes Positive integer group-class counts to try. With
#'   `id = NULL` it is 1 and may be left out.
#' @param model The covariance structures to cross with the class counts.
#'   `"basic"`, the default, fits the four structures that combine variances
#'   equal or varying across profiles with covariances absent or present:
#'   `"EEI"` (equal variances, no covariances), `"VVI"` (varying variances, no
#'   covariances), `"EEE"` (one covariance matrix shared by all profiles) and
#'   `"VVV"` (a covariance matrix for each profile). They are the choices that
#'   matter most in practice, since omitting covariances between correlated
#'   indicators makes the grid favour additional profiles. `"all"` fits all
#'   14 structures of Celeux and Govaert (1995). Otherwise name the structures
#'   by their three-letter codes, as [multilpa()]'s `model` takes them: the
#'   letters give the volume, shape and orientation of each profile's
#'   covariance matrix, each equal (`E`) or varying (`V`) across profiles,
#'   with `I` for no covariances. `NULL` crosses no structures and fits
#'   whatever the other arguments ask for. With one continuous indicator the
#'   structures reduce to equal or varying variance, and codes that coincide
#'   are fitted once; with only categorical indicators there is no structure
#'   to cross. A structure set through `variance_model`, `covariance_model`,
#'   `volume`, `shape` or `orientation` replaces the default.
#' @param seed Optional reproducible seed for each fit.
#' @param family `"profiles"` (the default) enumerates the profile model over
#'   `n_profiles`, `n_group_classes` and `model`. One or more of the
#'   group-class families `"additive"`, `"dispersion"` and
#'   `"additive_dispersion"` enumerates those instead, over `family`,
#'   `n_group_classes` and `between_variance`; `n_profiles` and `model` are
#'   then refused, and the result is a `latents_family_enumeration` read with
#'   [get_results.latents_family_enumeration()] and [candidate_fit()].
#' @param between_variance For group-class families: the between-variance
#'   restrictions to cross, `"varying"`, `"equal"` or both (the default). The
#'   dispersion family is always fitted with `"equal"`.
#' @param time The occasion column: enumerate latent transition models
#'   ([lta()]) over `n_profiles` (default `2:4`) and `n_group_classes`
#'   (default `1`), holding every other `lta()` argument given in `...`
#'   fixed. `model`, when given, crosses those covariance structures. The
#'   result is a `latents_transition_enumeration`, read with
#'   [get_results.latents_transition_enumeration()] and [candidate_fit()].
#' @param ... Further arguments to [multilpa()] (or to [lta()] with `time`); for group-class families only
#'   `n_starts`, `max_iter`, `tol` and `min_variance`.
#' @return An object of class `multilpa_enumeration`. Read it with the verbs
#'   that describe it rather than by reaching into it: [as.data.frame()] gives
#'   one row per candidate model with every criterion and diagnostic,
#'   [summary()] gives one row per information criterion naming the candidate
#'   that minimises it, [plot()] draws one criterion across the grid, and
#'   [candidate_fit()] returns the fitted model for one cell of the grid.
#' @seealso [candidate_fit()] to take one fitted model out of the grid,
#'   [summary.multilpa_enumeration()] for the criterion-by-criterion comparison.
#' @section Conditions:
#'   `latents_bad_argument` when `model` names an unknown structure, when it is
#'   given together with another way of setting the structure, or when it is
#'   given although every indicator is categorical.
#' @references
#' Celeux, G., & Govaert, G. (1995). Gaussian parsimonious clustering models.
#' *Pattern Recognition*, 28(5), 781--793.
#' @examples
#' set.seed(1)
#' d <- data.frame(g = rep(1:10, each = 10), y = rnorm(100))
#' candidates <- enumerate_classes(d, "y", "g", n_profiles = 1:2,
#'                               n_group_classes = 1, n_starts = 2, seed = 1)
#' as.data.frame(candidates)
#' summary(candidates)
#'
#' # Single level, crossing class counts with two named structures:
#' single <- enumerate_lpa(iris, c("Sepal.Length", "Sepal.Width",
#'                                 "Petal.Length", "Petal.Width"),
#'                         n_profiles = 1:3, model = c("VVV", "EEE"),
#'                         n_starts = 2, seed = 1)
#' summary(single)
#' @export
enumerate_classes <- function(data, vars, id, n_profiles = 1:4,
                              n_group_classes = 1:3, model = "basic",
                              seed = NULL, family = "profiles",
                              between_variance = c("varying", "equal"),
                              time = NULL, ...) {
  if (missing(id)) {
    .multilpa_missing_id("enumerate_classes",
                         "`enumerate_lpa()` or `enumerate_lca()`")
  }
  if (!is.null(time)) {
    # Transition models: profiles x group classes (x structures when `model`
    # names codes); every other lta() argument is held fixed through `...`.
    if (!identical(family, "profiles") || !missing(between_variance)) {
      stop(errorCondition(
        "Transition models (`time`) are enumerated for the profile family only.",
        class = "latents_bad_argument", call = NULL))
    }
    return(.lta_enumerate(data, vars, id, time,
                          if (missing(n_profiles)) 2:4 else n_profiles,
                          if (missing(n_group_classes)) 1L else n_group_classes,
                          if (missing(model)) NULL else model, seed, list(...),
                          match.call()))
  }
  if (!identical(family, "profiles")) {
    # Group-class families have no individual profiles and no covariance
    # structure codes; they are enumerated over family and class count.
    if (!missing(n_profiles) || !missing(model) || is.null(id)) {
      stop(errorCondition(paste(
        "Group-class families take `id`, `family`, `n_group_classes` and",
        "`between_variance`; `n_profiles` and `model` belong to the profile model."),
        class = "latents_bad_argument", call = NULL))
    }
    return(.additive_enumerate(data, vars, id, family, n_group_classes,
                               between_variance, seed, list(...), match.call()))
  }
  if (!missing(between_variance)) {
    stop(errorCondition(
      "`between_variance` belongs to the group-class families (`family = `).",
      class = "latents_bad_argument", call = NULL))
  }
  model_given <- !missing(model)
  single_level <- is.null(id)
  if (single_level) {
    # One observation per unit leaves no composition for a group class to
    # differ in, so the only group-class count is 1.
    if (!missing(n_group_classes) &&
        !isTRUE(all(as.integer(n_group_classes) == 1L))) {
      stop(errorCondition(paste(
        "`id = NULL` fits single-level models, which have one group class;",
        "leave `n_group_classes` out, or pass `id` to enumerate group classes."),
        class = "latents_bad_argument", call = NULL))
    }
    n_group_classes <- 1L
  }
  stopifnot(is.data.frame(data), is.character(vars),
            "`id` must be a single column name, or NULL for single-level models" =
              single_level || (is.character(id) && length(id) == 1L),
            is.numeric(n_profiles), length(n_profiles) > 0L,
            all(is.finite(n_profiles)), all(n_profiles >= 1), all(n_profiles == as.integer(n_profiles)),
            is.numeric(n_group_classes), length(n_group_classes) > 0L,
            all(is.finite(n_group_classes)), all(n_group_classes >= 1),
            all(n_group_classes == as.integer(n_group_classes)),
            "`model` must be \"basic\", \"all\", covariance model codes, or NULL" =
              is.null(model) || (is.character(model) &&
                                   length(model) > 0L &&
                                   !anyNA(model)))
  # Every name this guard used to reject is now a formal, so R refuses the call
  # with "matched by multiple actual arguments" before `...` is assembled.
  extra <- list(...)
  # Anything `multilpa()` does not take would otherwise fail inside every
  # candidate and be recorded as a grid of failures; refuse it here instead.
  unknown_arguments <- setdiff(names(extra), names(formals(multilpa)))
  if (length(unknown_arguments) > 0L) {
    hint <- if ("structure" %in% unknown_arguments)
      " Covariance models are named with `model`." else ""
    stop(errorCondition(sprintf(
      "%s %s not an argument of multilpa().%s",
      paste(sprintf("`%s`", unknown_arguments), collapse = ", "),
      if (length(unknown_arguments) == 1L) "is" else "are", hint),
      class = "latents_bad_argument", call = NULL))
  }
  # `NA` stands for "whatever the other arguments already said": no structure
  # is crossed, and each candidate is fitted as `multilpa()` would fit it.
  grid <- expand.grid(n_profiles = unique(n_profiles),
                      n_group_classes = unique(n_group_classes),
                      model = .multilpa_enumeration_models(
                        model, model_given, vars, extra),
                      stringsAsFactors = FALSE)
  runs <- lapply(seq_len(nrow(grid)), function(i) {
    warnings <- character()
    error_text <- NA_character_
    requested <- .multilpa_structure_arguments(grid$model[i])
    fit <- tryCatch(withCallingHandlers(do.call(multilpa,
      c(list(data = data, vars = vars, id = id,
             n_profiles = grid$n_profiles[i], n_group_classes = grid$n_group_classes[i], seed = seed),
        requested, extra)),
      # The single-level notice is given once for the whole grid below,
      # not once per candidate.
      latents_single_level = function(notice) invokeRestart("muffleMessage"),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
      }), error = function(error) {
        error_text <<- conditionMessage(error)
        NULL
      })
    row <- cbind(
      data.frame(n_profiles = grid$n_profiles[i], n_group_classes = grid$n_group_classes[i],
        # A latent class model has no covariance structure; naming one (the
        # default the fit records for its absent Gaussian block) misleads.
        model = if (is.null(fit)) grid$model[i] else
          if (length(.multilpa_continuous_names(fit)) == 0L) NA_character_ else
            fit$covariance_structure,
        log_likelihood = if (is.null(fit)) NA_real_ else fit$log_likelihood,
        n_parameters = if (is.null(fit)) NA_integer_ else fit$n_parameters,
        stringsAsFactors = FALSE),
      .multilpa_enumeration_indices(fit),
      data.frame(
        profile_entropy = if (is.null(fit)) NA_real_ else
          .multilpa_relative_entropy(fit$subject_posteriors),
        group_entropy = if (is.null(fit)) NA_real_ else
          .multilpa_relative_entropy(fit$group_posteriors),
        converged = !is.null(fit) && fit$converged,
        boundary = if (is.null(fit)) NA else fit$boundary,
        n_best_replicated = if (is.null(fit)) NA_integer_ else fit$n_best_replicated,
        warnings = paste(unique(warnings), collapse = "; "), error = error_text))
    list(fit = fit, row = row)
  })
  if (single_level) {
    .multilpa_single_level_notice(
      latent_class = all(vars %in% extra$categorical))
  }
  result <- list(table = do.call(rbind, lapply(runs, `[[`, "row")),
       fits = lapply(runs, `[[`, "fit"), call = match.call())
  class(result) <- "multilpa_enumeration"
  result
}

#' The covariance structures an enumeration crosses
#'
#' Expands the named sets and drops the structures the data cannot tell apart:
#' with one continuous indicator there are no covariances and no orientation,
#' so only the volume (equal or varying variance) distinguishes two codes, and
#' with none there is no covariance structure at all.
#'
#' @param model `enumerate_classes()`'s `model`.
#' @param model_given Whether the caller passed `model` rather than taking the
#'   default.
#' @param vars The indicator names.
#' @param extra The further arguments for `multilpa()`, as `list(...)`.
#' @return A character vector of structure codes, or `NA_character_` when no
#'   structure is to be crossed.
#' @noRd
.multilpa_enumeration_models <- function(model, model_given, vars, extra) {
  switches <- intersect(names(extra), c("variance_model", "covariance_model",
                                        "volume", "shape", "orientation"))
  if (length(switches) > 0L) {
    # The structure was set the other way; the default set gives way to it,
    # an explicit `model` contradicts it.
    if (model_given && !is.null(model)) {
      stop(errorCondition(sprintf(paste(
        "`model` and %s both set the covariance structure;",
        "name the structures with `model` alone."),
        paste(sprintf("`%s`", switches), collapse = ", ")),
        class = "latents_bad_argument", call = NULL))
    }
    return(NA_character_)
  }
  if (is.null(model)) return(NA_character_)
  sets <- list(basic = c("EEI", "VVI", "EEE", "VVV"),
               all = .multilpa_structures())
  codes <- unique(unlist(lapply(model, \(name) sets[[name]] %||% name),
                         use.names = FALSE))
  unknown <- setdiff(codes, .multilpa_structures())
  if (length(unknown) > 0L) {
    stop(errorCondition(sprintf(paste(
      "`model` names %s, which this package does not fit. Use \"basic\",",
      "\"all\", or codes from %s."),
      paste(sprintf("`%s`", unknown), collapse = ", "),
      paste(.multilpa_structures(), collapse = ", ")),
      class = "latents_bad_argument", call = NULL))
  }
  continuous <- setdiff(vars, c(extra$categorical, extra$ordinal, extra$count))
  if (length(continuous) == 0L) {
    if (model_given) {
      stop(errorCondition(paste(
        "No indicator is continuous, so there is no covariance structure",
        "for `model` to choose; leave `model` out."),
        class = "latents_bad_argument", call = NULL))
    }
    return(NA_character_)
  }
  if (length(continuous) == 1L) {
    # One variable: every code reduces to its volume letter.
    codes <- codes[!duplicated(substr(codes, 1L, 1L))]
  }
  codes
}

#' Take one fitted model out of an enumeration grid
#'
#' An enumeration keeps every candidate it fitted. This returns one of them, so
#' that a caller who wants to plot, summarise or test a particular candidate
#' names it by its class counts instead of indexing the grid by position.
#'
#' @param x An `multilpa_enumeration` result from [enumerate_classes()].
#' @param n_profiles Number of individual profiles identifying the candidate.
#' @param n_group_classes Number of group classes identifying the candidate.
#' @param model The covariance model identifying the candidate, needed when
#'   [enumerate_classes()] was given several `model`s and the class counts
#'   alone name several candidates. Naming counts that match more than one
#'   candidate raises `latents_unknown_candidate` listing the models it could
#'   have meant.
#' @return The fitted model for that cell of the grid: an object of class
#'   `multilpa`, exactly as [multilpa()] returned it, with every verb of this
#'   package available on it.
#' @section Conditions:
#'   `latents_unknown_candidate` when the requested class counts are not in
#'   the grid, and `latents_failed_candidate` when they are in the grid but
#'   that fit raised an error, so no model exists to return. Neither situation
#'   returns `NULL`, because a `NULL` would flow silently into whatever the
#'   caller did next.
#' @examples
#' set.seed(1)
#' d <- data.frame(g = rep(1:10, each = 10), y = rnorm(100))
#' candidates <- enumerate_classes(d, "y", "g", n_profiles = 1:2,
#'                               n_group_classes = 1, n_starts = 2, seed = 1)
#' candidate_fit(candidates, n_profiles = 2, n_group_classes = 1, model = "VVI")
#' @export
candidate_fit <- function(x, n_profiles, n_group_classes = 1L,
                          model = NULL) {
  if (inherits(x, "latents_transition_enumeration")) {
    at <- which(x$table$n_profiles == n_profiles &
                  x$table$n_group_classes == n_group_classes &
                  (is.null(model) | x$table$model %in% (model %||% "")))
    if (length(at) != 1L) {
      stop(errorCondition(sprintf(
        "%d candidates have %d profiles and %d group classes%s.", length(at),
        as.integer(n_profiles), as.integer(n_group_classes),
        if (length(at) > 1L) "; name one with `model`" else ""),
        class = "latents_unknown_candidate", call = NULL))
    }
    if (is.null(x$fits[[at]])) {
      stop(errorCondition(sprintf("That candidate could not be fitted: %s",
                                  x$table$error[at]),
                          class = "latents_failed_candidate", call = NULL))
    }
    return(x$fits[[at]])
  }
  if (inherits(x, "latents_family_enumeration")) {
    at <- which(x$table$n_group_classes == n_group_classes &
                  (is.null(model) | x$table$model %in% (model %||% "")))
    if (length(at) != 1L) {
      stop(errorCondition(sprintf(paste(
        "%d candidates have %d group classes%s; name one with `model`, one of %s."),
        length(at), as.integer(n_group_classes),
        if (is.null(model)) "" else sprintf(" and model %s", model),
        paste(sprintf("\"%s\"", unique(x$table$model)), collapse = ", ")),
        class = "latents_unknown_candidate", call = NULL))
    }
    if (is.null(x$fits[[at]])) {
      stop(errorCondition(sprintf("That candidate could not be fitted: %s",
                                  x$table$error[at]),
                          class = "latents_failed_candidate", call = NULL))
    }
    return(x$fits[[at]])
  }
  stopifnot(
    "`x` must be an `multilpa_enumeration` result" =
      inherits(x, "multilpa_enumeration"),
    "`n_profiles` must be a single positive integer" =
      is.numeric(n_profiles) && length(n_profiles) == 1L &&
      is.finite(n_profiles) && n_profiles >= 1 &&
      n_profiles == floor(n_profiles),
    "`n_group_classes` must be a single positive integer" =
      is.numeric(n_group_classes) && length(n_group_classes) == 1L &&
      is.finite(n_group_classes) && n_group_classes >= 1 &&
      n_group_classes == floor(n_group_classes),
    "`model` must be a single model code, or NULL" =
      is.null(model) || (is.character(model) &&
                           length(model) == 1L && !is.na(model)))
  grid <- x$table
  at <- which(grid$n_profiles == n_profiles &
                grid$n_group_classes == n_group_classes)
  described <- sprintf("%d profiles and %d group classes",
                       as.integer(n_profiles), as.integer(n_group_classes))
  if (!is.null(model)) {
    at <- at[grid$model[at] %in% model]
    described <- sprintf("%s under %s", described, model)
  }
  if (length(at) > 1L) {
    # A grid crossed with covariance structures has several candidates per
    # pair of class counts, so the counts alone no longer name one.
    stop(errorCondition(sprintf(paste(
      "%d candidates have %s: %s. Name one with `model`."),
      length(at), described, paste(grid$model[at], collapse = ", ")),
      class = "latents_unknown_candidate", call = NULL))
  }
  if (length(at) != 1L) {
    stop(errorCondition(sprintf("No candidate with %s was enumerated.",
                                described),
                        class = "latents_unknown_candidate", call = NULL))
  }
  fit <- x$fits[[at]]
  if (is.null(fit)) {
    stop(errorCondition(sprintf(
      "The candidate with %s could not be fitted: %s", described,
      grid$error[at]),
      class = "latents_failed_candidate", call = NULL))
  }
  fit
}

#' @rdname latents-summary
#' @export
summary.multilpa_enumeration <- function(object, ...) {
  stopifnot("`object` must be an `multilpa_enumeration` result" =
              inherits(object, "multilpa_enumeration"))
  grid <- object$table
  result <- list(criteria = .multilpa_enumeration_minima(grid),
                 candidates = grid,
                 n_candidates = nrow(grid),
                 n_converged = sum(grid$converged),
                 n_failed = sum(is.na(grid$log_likelihood)),
                 n_boundary = sum(grid$boundary %in% TRUE),
                 call = object$call)
  # Every table the fit can produce, built once here, so `get_results()` on the
  # summary serves the same tables the fit would and `print()` can show them
  # all without recomputing anything.
  result$tables <- get_results(object, "all")
  class(result) <- "summary_multilpa_enumeration"
  result
}

#' Which candidate minimises each information criterion
#' @param grid The enumeration table.
#' @return One row per criterion carried by the grid.
#' @noRd
.multilpa_enumeration_minima <- function(grid) {
  present <- intersect(.multilpa_enumeration_criteria(), names(grid))
  eligible <- grid$converged %in% TRUE
  rows <- lapply(present, function(name) {
    values <- grid[[name]]
    values[!eligible] <- NA_real_
    best <- if (all(is.na(values))) NA_integer_ else which.min(values)
    data.frame(
      criterion = sub("_(groups|individual)$", "", name),
      convention = if (grepl("_groups$", name)) "groups" else
        if (grepl("_individual$", name)) "individuals" else NA_character_,
      n_profiles = if (is.na(best)) NA_integer_ else
        as.integer(grid$n_profiles[best]),
      n_group_classes = if (is.na(best)) NA_integer_ else
        as.integer(grid$n_group_classes[best]),
      model = if (is.na(best)) NA_character_ else
        as.character(grid$model[best]),
      value = if (is.na(best)) NA_real_ else values[best],
      row.names = NULL, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' @rdname latents-print
#' @export
print.summary_multilpa_enumeration <- function(x, digits = 4L, rows = 10L, ...) {
  stopifnot("`x` must be a `summary_multilpa_enumeration` object" =
              inherits(x, "summary_multilpa_enumeration"))
  .multilpa_check_print_arguments(digits, rows)
  cat(sprintf("Class enumeration: %d candidates, %d converged, %d failed to fit\n",
              x$n_candidates, x$n_converged, x$n_failed))
  identified <- c("n_profiles", "n_group_classes", "model")
  present <- x$criteria[!is.na(x$criteria$n_profiles), identified, drop = FALSE]
  disagreement <- unique(present)
  cat(sprintf("%d distinct candidate(s) are minimal under some criterion.\n",
              nrow(disagreement)))
  cat("No candidate is selected automatically. Choose one convention and keep it.\n")
  tables <- x$tables
  if (!is.null(tables$candidates)) {
    full_width <- ncol(tables$candidates)
    tables$candidates <- .multilpa_compact_candidates(tables$candidates)
    cat(sprintf(paste(
      "The candidates table shows %d of its %d columns;",
      "get_results(x, what = \"candidates\") returns all of them.\n"),
      ncol(tables$candidates), full_width))
  }
  .multilpa_print_tables(tables, rows = rows, digits = digits)
  .multilpa_print_table_footer(tables)
  invisible(x)
}

#' The columns of an enumeration grid worth printing
#'
#' The grid carries every criterion at both levels, which wraps a printed
#' table across several screens. The print keeps the columns a reader compares
#' candidates on, drops identifiers that do not vary, and shows a criterion
#' once where its group-level and individual-level versions are the same
#' numbers, as they are in a single-level grid. The full table is unchanged.
#'
#' @param grid The enumeration table.
#' @return A data frame with fewer columns.
#' @noRd
.multilpa_compact_candidates <- function(grid) {
  varies <- function(name) name %in% names(grid) &&
    length(unique(grid[[name]])) > 1L
  same_levels <- function(stem) {
    both <- paste0(stem, c("_groups", "_individual"))
    all(both %in% names(grid)) &&
      isTRUE(all.equal(grid[[both[1L]]], grid[[both[2L]]]))
  }
  bic <- if (same_levels("bic")) "bic_groups" else c("bic_groups", "bic_individual")
  keep <- c("n_profiles",
            if (varies("n_group_classes")) "n_group_classes",
            if (varies("model")) "model",
            "log_likelihood", "n_parameters", "aic", bic, "icl_individual",
            "profile_entropy",
            if (!all(is.na(grid$group_entropy))) "group_entropy",
            "converged", "boundary")
  compact <- grid[, intersect(keep, names(grid)), drop = FALSE]
  if (same_levels("bic")) names(compact)[names(compact) == "bic_groups"] <- "bic"
  compact
}

#' @rdname latents-as-data-frame
#' @export
as.data.frame.summary_multilpa_enumeration <- function(x, row.names = NULL, optional = FALSE, ...) {
  stopifnot("`x` must be an object of class `summary_multilpa_enumeration`" = inherits(x, "summary_multilpa_enumeration"))
  .multilpa_coerce(x, row.names, list(...))
}

#' Generate observations from a fitted discrete multilevel model
#' @param object Fitted model.
#' @return Simulated data frame with the original group layout.
#' @noRd
.multilpa_simulate <- function(object) {
  stopifnot(inherits(object, "multilpa"), !inherits(object, "multilpa_covariates"))
  group_class <- sample.int(object$n_group_classes, object$n_groups,
                            replace = TRUE, prob = object$group_probabilities)
  profile <- .latents_draw_rows(
    object$profile_probabilities[group_class[object$group_index], , drop = FALSE])
  result <- .multilpa_draw_indicators(object, profile)
  result[[object$id]] <- object$group_values[object$group_index]
  result
}

#' Generate observations from a fitted covariate model, holding the covariates
#'
#' The model conditions on its covariates, so the simulation does too: every
#' group keeps its observed group covariates and every row its profile
#' covariates. The group class is drawn from the group logits, the profile from
#' the row's profile logits in that class, and the indicators from the profile.
#'
#' @param object A fitted `multilpa_covariates` model.
#' @param frame The fitting data; its indicator columns are replaced, and every
#'   other column (identifier and covariates) is kept as observed.
#' @return `frame` with simulated indicators.
#' @noRd
.multilpa_cov_simulate <- function(object, frame) {
  stopifnot(inherits(object, "multilpa_covariates"))
  group_prior <- .multilpa_softmax(object$group_design, object$group_coefficients)
  group_class <- .latents_draw_rows(group_prior)[object$group_index]
  priors <- lapply(object$profile_design, function(design) {
    .multilpa_softmax(design, object$profile_coefficients)
  })
  # Each row's profile distribution is its own group class's.
  profile_prior <- Reduce(`+`, lapply(seq_along(priors), function(class) {
    priors[[class]] * (group_class == class)
  }))
  profile <- .latents_draw_rows(profile_prior)
  simulated <- .multilpa_draw_indicators(object, profile)
  frame[object$vars] <- simulated[object$vars]
  frame
}

#' Give simulated data the observed data's missing values
#'
#' A bootstrap of a `missing = "fiml"` fit must refit replicates that lose the
#' same information the observed data lost, or the reference distribution is
#' that of a complete-data statistic. The observed pattern is carried over
#' cell for cell, which treats it as fixed: independent of the profiles and of
#' the values that went missing (missing completely at random given the
#' pattern).
#'
#' @param simulated The simulated data frame.
#' @param observed The fitting data, in the same row order.
#' @param vars The indicator columns.
#' @return `simulated` with `NA` wherever `observed` has one among `vars`.
#' @noRd
.multilpa_carry_missingness <- function(simulated, observed, vars) {
  simulated[vars] <- Map(function(values, reference) {
    values[is.na(reference)] <- NA
    values
  }, simulated[vars], observed[vars])
  simulated
}

#' Draw indicators for given profiles from a fitted measurement model
#' @param object A fitted model carrying the measurement blocks.
#' @param profile Integer profile per observation.
#' @return A data frame with one column per indicator, in `object$vars` order.
#' @noRd
.multilpa_draw_indicators <- function(object, profile) {
  continuous <- .multilpa_continuous_names(object)
  dimension <- length(continuous)
  result <- data.frame(row.names = seq_len(object$n_observations))
  if (dimension > 0L) {
    noise <- matrix(stats::rnorm(object$n_observations * dimension),
                    object$n_observations)
    values <- if (identical(object$covariance_model, "full")) {
      factors <- lapply(seq_len(object$n_profiles),
                        function(k) chol(object$covariances[, , k]))
      t(matrix(vapply(seq_len(object$n_observations), function(i) {
        as.vector(noise[i, , drop = FALSE] %*% factors[[profile[i]]]) +
          object$means[profile[i], ]
      }, numeric(dimension)), nrow = dimension, ncol = object$n_observations))
    } else {
      object$means[profile, , drop = FALSE] +
        noise * sqrt(object$variances[profile, , drop = FALSE])
    }
    values <- as.data.frame(values)
    names(values) <- continuous
    result <- cbind(result, values)
  }
  ## Categorical indicators return their original value type and factor level
  ## order, so refitting re-encodes the same categories.
  blocks <- object$response_probabilities
  if (!is.null(blocks) && length(blocks) > 0L) {
    drawn <- lapply(names(blocks), function(indicator) {
      probabilities <- blocks[[indicator]][profile, , drop = FALSE]
      values <- object$categorical_values[[indicator]]
      if (is.null(values)) {
        # Older fits retained labels but not their input types. Keep their
        # previous numeric-label fallback for saved objects.
        values <- object$categorical_levels[[indicator]]
        # Coercion warns for labels such as "a"; the NA result immediately
        # below decides whether to retain those labels unchanged.
        numeric_values <- suppressWarnings(as.numeric(values)) # nolint: undesirable_function_linter.
        if (!anyNA(numeric_values)) values <- numeric_values
      }
      values[.latents_draw_rows(probabilities)]
    })
    names(drawn) <- names(blocks)
    result <- cbind(result, as.data.frame(drawn, stringsAsFactors = FALSE))
  }
  extra <- .latents_draw_extra(object, profile)
  if (length(extra) > 0L) result[names(extra)] <- extra
  result[, object$vars, drop = FALSE]
}

#' The measurement constraint both bootstrap refits must carry
#'
#' A fit made with `fixed` treats its measurement blocks as known, so its
#' likelihood is conditional on those values and its free parameters are only
#' the ones the constraint left. A bootstrap that refits the replicates without
#' the constraint compares unrestricted models, while the observed statistic
#' compares constrained ones, and the resulting p-value answers a question
#' nobody asked. This decides, once, what every refit must be given, and
#' refuses the pairs whose constrained nesting cannot be established.
#'
#' Nesting holds when both models hold the same blocks at the same values and
#' differ only in the number of group classes: the extra class is a mixing
#' parameter, and the null is the alternative with that class emptied. It does
#' not hold when the models differ in the number of profiles, because the held
#' measurement then has a different shape in each model and the smaller one is
#' not a restriction of the larger one.
#'
#' @param null_model The smaller fitted `multilpa` model.
#' @param alternative_model The larger fitted `multilpa` model.
#' @return A list with `start`, the `multilpa_start` carrying the held
#'   measurement values (`NULL` when nothing is held), and `fixed`, the
#'   character vector of held block names (`character(0)` when nothing is
#'   held). Raises `latents_bad_nesting` when the two constrained models are
#'   not nested.
#' @noRd
.multilpa_bootstrap_constraint <- function(null_model, alternative_model) {
  stopifnot("`null_model` must be a fitted `multilpa` model" =
              inherits(null_model, "multilpa"),
            "`alternative_model` must be a fitted `multilpa` model" =
              inherits(alternative_model, "multilpa"))
  held_null <- null_model$fixed %||% character()
  held_alternative <- alternative_model$fixed %||% character()
  if (length(held_null) == 0L && length(held_alternative) == 0L) {
    return(list(start = NULL, fixed = character()))
  }
  if (length(held_null) == 0L || length(held_alternative) == 0L) {
    stop(errorCondition(paste(
      "One model holds measurement blocks fixed and the other estimates them,",
      "so the smaller is not nested in the larger. Fit both with the same",
      "`fixed` specification, or neither."),
      class = "latents_bad_nesting", call = NULL))
  }
  if (!identical(sort(held_null), sort(held_alternative))) {
    stop(errorCondition(sprintf(paste(
      "The models hold different measurement blocks fixed (%s against %s), so",
      "they are not nested. Give both the same `fixed` specification."),
      paste(sprintf("\"%s\"", sort(held_null)), collapse = ", "),
      paste(sprintf("\"%s\"", sort(held_alternative)), collapse = ", ")),
      class = "latents_bad_nesting", call = NULL))
  }
  if (!identical(null_model$n_profiles, alternative_model$n_profiles)) {
    stop(errorCondition(sprintf(paste(
      "A held measurement cannot be compared across %d and %d profiles: the",
      "held blocks have different shapes, so the null is not a restriction of",
      "the alternative. Compare models that differ by one group class, or",
      "refit both without `fixed`."),
      null_model$n_profiles, alternative_model$n_profiles),
      class = "latents_bad_nesting", call = NULL))
  }
  covariance_model <- null_model$covariance_model %||% "diagonal"
  start <- starting_values(null_model, what = "measurement")
  values_null <- .multilpa_held_parameters(start, held_null, covariance_model)
  values_alternative <- .multilpa_held_parameters(
    starting_values(alternative_model, what = "measurement"), held_alternative,
    alternative_model$covariance_model %||% "diagonal")
  agreement <- all.equal(values_null, values_alternative, tolerance = 1e-8)
  if (!isTRUE(agreement)) {
    stop(errorCondition(sprintf(paste(
      "The models hold the same blocks at different values, so they are not",
      "nested: %s. Hold both at one measurement solution, as fit_staged()",
      "does."), paste(agreement, collapse = "; ")),
      class = "latents_bad_nesting", call = NULL))
  }
  if (alternative_model$n_parameters <= null_model$n_parameters) {
    stop(errorCondition(sprintf(paste(
      "The alternative must estimate more free parameters than the null; with",
      "this constraint it estimates %d against %d."),
      alternative_model$n_parameters, null_model$n_parameters),
      class = "latents_bad_nesting", call = NULL))
  }
  list(start = start, fixed = held_null)
}

#' Are two covariate models nested for a bootstrap comparison
#'
#' The alternative must be the null with one more profile or group class and
#' the same membership regressions. A null with one group class cannot carry
#' group covariates or slopes by group class (both need a second class), and it
#' is still the alternative with one class emptied, so there the alternative
#' may add them.
#'
#' @param null_model,alternative_model Fitted `multilpa_covariates` models.
#' @return `NULL`, invisibly; raises `latents_bad_nesting`.
#' @noRd
.multilpa_check_covariate_nesting <- function(null_model, alternative_model) {
  refuse <- function(message) {
    stop(errorCondition(message, class = "latents_bad_nesting", call = NULL))
  }
  if (!identical(null_model$profile_covariates, alternative_model$profile_covariates)) {
    refuse("The models use different profile covariates, so they are not nested.")
  }
  one_class_null <- null_model$n_group_classes == 1L
  if (!one_class_null &&
      !identical(null_model$group_covariates, alternative_model$group_covariates)) {
    refuse("The models use different group covariates, so they are not nested.")
  }
  if (one_class_null && length(null_model$group_covariates) > 0L) {
    refuse("A one-group-class null cannot carry group covariates.")
  }
  slopes <- function(model) model$profile_slopes %||% "shared"
  if (!one_class_null && !identical(slopes(null_model), slopes(alternative_model))) {
    refuse(paste(
      "One model lets profile slopes differ by group class and the other does",
      "not, so they differ by more than one class."))
  }
  invisible(NULL)
}

#' Check that data reproduce a covariate fit's likelihood
#' @param model A fitted `multilpa_covariates` model.
#' @param data The fitting data.
#' @return `NULL`, invisibly; raises `latents_bad_inference_data`.
#' @noRd
.multilpa_cov_check_likelihood <- function(model, data) {
  broken <- function() {
    stop(errorCondition("data do not reproduce the fitted model likelihood.",
                        class = "latents_bad_inference_data", call = NULL))
  }
  needed <- unique(c(model$vars, model$id, model$profile_covariates,
                     model$group_covariates))
  if (!is.data.frame(data) || !all(needed %in% names(data)) ||
      nrow(data) != model$n_observations ||
      !identical(match(data[[model$id]], model$group_values), model$group_index)) {
    broken()
  }
  first_rows <- match(seq_len(model$n_groups), model$group_index)
  designs <- tryCatch(
    .multilpa_cov_designs(data, model$vars, model$profile_covariates,
                          model$group_covariates, first_rows,
                          model$n_group_classes, model$categorical %||% character(),
                          model$min_probability %||% 1e-10,
                          model$missing %||% "error",
                          model$profile_slopes %||% "shared",
                          model$ordinal %||% character(), model$count %||% character(),
                          model$extra_data$count_model %||% "poisson",
                          model$extra_data$count_dispersion %||% "varying"),
    latents_bad_data = function(condition) broken())
  parameters <- model[intersect(c("variances", "covariances", "response_probabilities",
                                  "ordinal_intercepts", "ordinal_locations",
                                  "count_means", "count_dispersion"), names(model))]
  parameters$means <- if (ncol(designs$x) > 0L) {
    sweep(model$means, 2L, designs$center, "-")
  } else model$means
  likelihood <- .multilpa_cov_expectation(
    designs$x, model$group_index, parameters, designs$profile_design,
    designs$w, model$profile_coefficients, model$group_coefficients,
    designs$codes, extra = designs$extra)$log_likelihood
  if (abs(likelihood - model$log_likelihood) > 1e-7 * (1 + abs(likelihood))) broken()
  invisible(NULL)
}

#' Parametric bootstrap likelihood-ratio comparison
#'
#' Simulates indicators under the null model while preserving observed group
#' sizes, refits both models, and compares their likelihood differences.
#' Models must differ by exactly one individual profile or one group class, with
#' the other count, covariance structure and centering mode fixed. Grand-mean
#' centering is repeated in every simulated refit. Person-centred fits are
#' refused because this model does not specify a generative distribution for
#' the group baselines removed by that transformation. This is a native
#' parametric bootstrap, not an implementation of Mplus TECH14. It does not use
#' a chi-square reference distribution. Any failed/nonconverged or reversed
#' replicate makes the p-value NA, avoiding silent deletion of difficult fits.
#'
#' Fits made with `missing = "fiml"` are supported when both models were: every
#' replicate is given the observed data's missing cells before it is refitted,
#' preserving the pattern of available measurements. Calibration requires
#' missingness independent of the indicators and latent classes, conditional
#' on any fixed covariates. Carrying a fixed mask does not reproduce a general
#' MAR mechanism that depends on observed indicators, or an MNAR mechanism.
#' FIML estimation under MAR does not by itself validate this bootstrap under
#' MAR; the retained size study covers MCAR only.
#'
#' Membership-covariate fits (`multilpa(profile_covariates = ,
#' group_covariates = )`) are supported: the covariates are held at their
#' observed values, group classes are drawn from the fitted group logits and
#' profiles from each row's profile logits, as the model conditions on the
#' covariates. Both models must use the same profile covariates, group
#' covariates and slope specification, except that a one-group-class null,
#' which cannot carry group covariates or slopes by group class, may be
#' compared with an alternative that adds them. When `data` is omitted it is
#' rebuilt from the alternative, whose columns include the null's. A
#' covariate replicate counts as valid when its likelihood has converged even
#' if a membership logit is still drifting, as it does for an empty or
#' separated class in an over-fitted alternative: the statistic reads only the
#' maximized likelihood, which such a fit has reached. The replicate table's
#' `logits_settled` column records whether both refits' logits had also
#' converged. The original models must be fully converged.
#'
#' Models fitted with `fixed`, including those from [fit_staged()], are
#' supported: every replicate is refitted with the same blocks held at the same
#' values, so the simulated statistics compare the same two constrained models
#' as the observed statistic does. The p-value is then conditional on that
#' measurement solution, which is treated as known and whose own uncertainty is
#' not propagated. A constrained pair whose nesting cannot be established is
#' refused rather than given a p-value.
#' @section Conditions:
#'   `latents_bad_nesting` when only one of the two models holds measurement
#'   blocks fixed, when they hold different blocks, when they hold the same
#'   blocks at different values, when they hold a measurement fixed while
#'   differing in the number of profiles (the held blocks then have different
#'   shapes, so the null is not a restriction of the alternative), or when the
#'   alternative does not estimate more free parameters than the null.
#'   `latents_failed_replicates` is warned when some replicate fails
#'   validation, and the p-value is `NA`.
#'   `latents_unsupported_bootstrap` refuses person-centred fits, for which
#'   this simulator has no group-baseline distribution.
#'   `latents_unsupported_prior` refuses fits estimated with a prior: their
#'   likelihoods are evaluated at posterior modes, not ML estimates.
#'   `latents_incomparable_models` when the two models differ in their
#'   observations, grouping, covariance structure, centering or missing-data
#'   handling; `latents_bad_data` when the data have missing indicators and
#'   the models were not fitted with `missing = "fiml"`; `latents_bad_nesting`
#'   also when two covariate models use different covariates or slope
#'   specifications (other than a one-group-class null).
#' @param null_model Smaller, converged [multilpa()] model, with or without
#'   membership covariates.
#' @param alternative_model Larger model fitted to exactly the same data.
#' @param data Optional. The data both models were fitted to, used to verify
#'   both fitted likelihoods; when omitted it is rebuilt from what the null
#'   model stores.
#' @param iter Number of simulated datasets (at least two; use many for inference).
#' @param n_starts Number of starts for each simulated fit.
#' @param max_iter Maximum EM iterations for each simulated fit.
#' @param tol Relative likelihood convergence tolerance.
#' @param seed Optional seed: any whole number `set.seed()` accepts. The
#'   caller's random-number state is restored on exit.
#' @return An object of class `multilpa_bootstrap_lrt`, carrying the observed
#'   statistic, the finite-simulation corrected p-value, its Monte Carlo
#'   standard error, the measurement blocks held fixed in every fit, and one
#'   record per replicate. Read it with the verbs that describe it:
#'   [as.data.frame()] gives the single-row test result,
#'   `get_results(what = "replicates")` one row per simulated replicate,
#'   [summary()] the
#'   test beside the replicate diagnostics, and [plot()] the simulated null
#'   distribution with the observed statistic marked. Inspect failed starts,
#'   boundary flags and likelihood replication before interpreting results.
#' @examples
#' set.seed(1)
#' example_data <- data.frame(
#'   school = rep(seq_len(10), each = 10),
#'   score = rnorm(100, rep(c(-2, 2), each = 50))
#' )
#' smaller <- multilpa(example_data, "score", "school", n_profiles = 1,
#'                     n_group_classes = 1, n_starts = 2, seed = 1)
#' larger <- multilpa(example_data, "score", "school", n_profiles = 2,
#'                    n_group_classes = 1, n_starts = 2, seed = 1)
#' # `iter` is small so the example runs quickly; use many more for inference.
#' comparison <- bootstrap_lrt(smaller, larger, iter = 9, n_starts = 1,
#'                             max_iter = 2000, tol = 1e-6, seed = 1)
#' comparison
#' @export
bootstrap_lrt <- function(null_model, alternative_model, data = NULL,
                                 iter = 199L, n_starts = 10L, max_iter = 1000L,
                                 tol = 1e-8, seed = NULL) {
  if (.latents_is_weighted(null_model) || .latents_is_weighted(alternative_model)) {
    .latents_refuse_weights(paste(
      "`bootstrap_lrt()`: the parametric bootstrap draws samples from the",
      "fitted model, not from the sampling design, so its reference",
      "distribution ignores the weights. Compare weighted fits on their",
      "information criteria"))
  }
  transition_classes <- c("multilpa_transitions", "multilpa_lta")
  if (inherits(null_model, transition_classes) ||
      inherits(alternative_model, transition_classes)) {
    if (!inherits(null_model, transition_classes) ||
        !inherits(alternative_model, transition_classes)) {
      stop(errorCondition("bootstrap_lrt() compares two transition fits, or two of another family.",
                          class = "latents_bad_argument", call = NULL))
    }
    return(.lta_bootstrap_lrt(null_model, alternative_model, data, iter, n_starts,
                              max_iter, tol, seed, match.call()))
  }
  if (inherits(null_model, "multilpa_additive") ||
      inherits(alternative_model, "multilpa_additive")) {
    return(.additive_bootstrap_lrt(null_model, alternative_model, data, iter,
                                   n_starts, max_iter, tol, seed, match.call()))
  }
  covariate_family <- inherits(null_model, "multilpa_covariates")
  stopifnot(
    "both models must be multilpa() fits of one family, with or without membership covariates" =
      (inherits(null_model, "multilpa") && inherits(alternative_model, "multilpa")) ||
      (covariate_family && inherits(alternative_model, "multilpa_covariates")))
  .multilpa_refuse_noise(null_model, "bootstrap_lrt()")
  .multilpa_refuse_noise(alternative_model, "bootstrap_lrt()")
  if (!is.null(null_model$prior) || !is.null(alternative_model$prior)) {
    stop(errorCondition(paste(
      "bootstrap_lrt() requires maximum-likelihood fits; a fit with `prior`",
      "is estimated at a posterior mode. Refit both models without `prior`."),
      class = "latents_unsupported_prior", call = NULL))
  }
  # The null model is the one being simulated from, so its own columns are the
  # ones that matter when the caller does not supply data.
  # A covariate alternative carries every covariate the null does (the nesting
  # rules below require it) and possibly group covariates the one-class null
  # cannot, so its stored columns are the ones both refits need.
  data <- if (covariate_family && is.null(data)) {
    .multilpa_cov_stored_data(alternative_model)
  } else .multilpa_number_single_level(null_model,
                                       .multilpa_resolve_data(null_model, data))
  stopifnot(is.data.frame(data), is.numeric(iter), length(iter) == 1L,
            is.finite(iter), iter >= 2L, iter == as.integer(iter),
            is.numeric(n_starts), length(n_starts) == 1L, is.finite(n_starts),
            n_starts >= 1, n_starts == as.integer(n_starts),
            is.numeric(max_iter), length(max_iter) == 1L, is.finite(max_iter),
            max_iter >= 1, max_iter == as.integer(max_iter),
            is.numeric(tol), length(tol) == 1L, is.finite(tol), tol > 0)
  .multilpa_check_seed(seed)
  .latents_local_seed(seed)
  fields <- c("vars", "id", "group_values", "group_index", "variance_model",
              "min_variance", "covariance_model")
  fields <- c(fields, "categorical", "categorical_levels", "min_probability",
              "ordinal", "count")
  extra_spec <- function(model) {
    list(levels = model$extra_data$ordinal_levels,
         count_model = model$extra_data$count_model %||% "poisson",
         count_dispersion = model$extra_data$count_dispersion %||% "varying")
  }
  missing_for <- function(model) model$missing %||% "error"
  structure_for <- function(model) model$covariance_structure %||%
    .multilpa_resolve_structure(model$variance_model,
                                model$covariance_model %||% "diagonal")
  centering_for <- function(model) model$centering %||% "none"
  same_fields <- all(vapply(fields, function(field)
    identical(null_model[[field]], alternative_model[[field]]), logical(1)))
  if (!same_fields || !identical(extra_spec(null_model), extra_spec(alternative_model)) ||
      !identical(structure_for(null_model), structure_for(alternative_model)) ||
      !identical(centering_for(null_model), centering_for(alternative_model)) ||
      !identical(missing_for(null_model), missing_for(alternative_model))) {
    stop(errorCondition(paste(
      "Models must use the same observations, group layout, covariance",
      "structure, centering mode and missing-data handling."),
        class = "latents_incomparable_models", call = NULL))
  }
  if (covariate_family) .multilpa_check_covariate_nesting(null_model, alternative_model)
  if (identical(centering_for(null_model), "person")) {
    stop(errorCondition(paste(
      "A person-centred fit removes each group's baseline, but this bootstrap",
      "does not model those baselines for simulation. Compare fits without",
      "person centering."),
      class = "latents_unsupported_bootstrap", call = NULL))
  }
  delta <- c(alternative_model$n_profiles - null_model$n_profiles,
             alternative_model$n_group_classes - null_model$n_group_classes)
  if (!all(delta >= 0) || sum(delta) != 1) stop(errorCondition("Models must differ by exactly one class at one level.",
        class = "latents_bad_nesting", call = NULL))
  # The observed statistic compares two constrained fits whenever the models
  # hold measurement blocks fixed, so every replicate must be refitted under
  # the same constraint or the reference distribution belongs to a different
  # pair of models than the statistic it is being compared against. Nesting is
  # a property of the two models, so it is settled before the data are read.
  constraint <- if (covariate_family) list(start = NULL, fixed = character()) else
    .multilpa_bootstrap_constraint(null_model, alternative_model)
  fiml <- identical(missing_for(null_model), "fiml")
  if (covariate_family) {
    invisible(lapply(list(null_model, alternative_model), function(model) {
      if (!isTRUE(model$converged) || isTRUE(model$boundary)) {
        stop(errorCondition(
          "Original models must be converged with inactive variance bounds.",
          class = "latents_no_converge", call = NULL))
      }
      .multilpa_cov_check_likelihood(model, data)
    }))
  } else invisible(lapply(list(null_model, alternative_model), function(model) {
    if (!isTRUE(model$converged) || isTRUE(model$boundary)) {
      stop(errorCondition(
        "Original models must be converged with inactive variance bounds.",
        class = "latents_no_converge", call = NULL))
    }
    if (!all(c(model$vars, model$id) %in% names(data)) ||
        nrow(data) != model$n_observations ||
        !identical(data[[model$id]], model$group_values[model$group_index])) {
      stop(errorCondition("data must preserve the original observations and group ordering.",
        class = "latents_bad_inference_data", call = NULL))
    }
    continuous <- .multilpa_continuous_names(model)
    x <- if (length(continuous) == 0L) matrix(numeric(0), nrow(data), 0L) else
      .multilpa_center_like(null_model, as.matrix(data[continuous]))
    # Missing values are the ones the fit integrated out; anything else that
    # is not finite is a broken contract.
    if (!is.numeric(x) || any(is.infinite(x)) || any(is.nan(x)) ||
        (!fiml && anyNA(x))) {
      stop(errorCondition(paste(
        "Bootstrap needs finite indicators; missing values are allowed only",
        "when both models were fitted with `missing = \"fiml\"`."),
        class = "latents_bad_data", call = NULL))
    }
    if (!is.null(model$indicator_data) && length(continuous) > 0L &&
        !identical(x, model$indicator_data)) {
      stop(errorCondition("data must reproduce the original indicator data and row order.",
        class = "latents_bad_inference_data", call = NULL))
    }
    encoded <- if (length(model$categorical %||% character()) == 0L) NULL else
      .multilpa_encode_categorical(data[, model$categorical, drop = FALSE])
    codes <- encoded$codes
    if (!fiml && anyNA(codes)) {
      stop(errorCondition(paste(
        "Bootstrap needs complete indicators unless both models were fitted",
        "with `missing = \"fiml\"`."), class = "latents_bad_data", call = NULL))
    }
    if (!is.null(encoded) && !identical(encoded$levels, model$categorical_levels)) {
      stop(errorCondition("data must reproduce the original categorical levels and coding.",
        class = "latents_bad_inference_data", call = NULL))
    }
    if (!is.null(codes) && !identical(unname(codes), unname(model$categorical_data))) {
      stop(errorCondition("data must reproduce the original categorical indicators and row order.",
        class = "latents_bad_inference_data", call = NULL))
    }
    if (!is.null(model$extra_data)) .multilpa_check_alignment(model, data)
    extra <- if (is.null(model$extra_data)) NULL else
      .latents_prepare_extra_like(model, data)
    likelihood <- .multilpa_expectation(x, model$group_index, model,
                                        codes, extra = extra)$log_likelihood
    if (abs(likelihood - model$log_likelihood) > 1e-7 * (1 + abs(likelihood))) {
      stop(errorCondition("data do not reproduce the fitted model likelihood.",
        class = "latents_bad_inference_data", call = NULL))
    }
  }))
  # How far below zero a statistic may dip and still be convergence noise rather
  # than a real reversal. EM stops on a RELATIVE change, so each fit sits within
  # about `tol * |log_likelihood|` of its optimum and the statistic doubles that.
  # A fixed absolute window instead rejected healthy replicates on any data whose
  # likelihood is large: at `tol = 1e-8` on a log likelihood of -3813 the
  # boundary noise is 7.6e-05, seven times the old fixed -1e-5. The multiplier
  # allows for the several EM steps that noise accumulates over, and the window
  # tightens automatically when the caller tightens `tol`.
  reversal_window <- 100 * tol * (1 + abs(null_model$log_likelihood))
  observed <- 2 * (alternative_model$log_likelihood - null_model$log_likelihood)
  if (observed < -reversal_window) {
    stop(errorCondition(sprintf(
      paste("Alternative has lower likelihood by %.3g, beyond the %.3g that",
            "`tol = %.3g` can explain; improve its optimization first."),
      -observed, reversal_window, tol),
      class = "latents_reversed_likelihood", call = NULL))
  }
  observed <- max(0, observed)
  replicates <- .latents_blrt_replicates(iter, function(i) {
    simulated <- if (covariate_family) .multilpa_cov_simulate(null_model, data) else
      .multilpa_simulate(null_model)
    if (fiml) {
      simulated <- .multilpa_carry_missingness(simulated, data, null_model$vars)
    }
    models <- lapply(list(null_model, alternative_model), function(model) {
      ## `variance_model` and `covariance_model` cannot express a structure
      ## that constrains the volume, the shape or the orientation: refitting
      ## from them alone returned a VEI replicate as VVI, two parameters
      ## wider, so the reference distribution belonged to a different pair of
      ## models than the statistic compared against it.
      family <- .multilpa_structure_arguments(
        model$covariance_structure %||% NA_character_)
      if (length(family) == 0L) {
        family <- list(variance_model = model$variance_model,
                       covariance_model = model$covariance_model %||% "diagonal")
      }
      do.call(multilpa, c(
        list(data = simulated, vars = model$vars, id = model$id,
             n_profiles = model$n_profiles,
             n_group_classes = model$n_group_classes),
        family,
        list(n_starts = n_starts, max_iter = max_iter, tol = tol,
             min_variance = model$min_variance,
             centering = model$centering %||% "none",
             categorical = model$categorical %||% character(),
             ordinal = model$ordinal %||% character(),
             count = model$count %||% character(),
             count_model = model$extra_data$count_model %||% "poisson",
             count_dispersion = model$extra_data$count_dispersion %||% "varying",
             min_probability = model$min_probability %||% 1e-10,
             missing = missing_for(model),
             start = constraint$start, fixed = constraint$fixed),
        if (covariate_family) list(
          profile_covariates = model$profile_covariates,
          group_covariates = model$group_covariates,
          profile_slopes = model$profile_slopes %||% "shared")))
    })
    statistic <- 2 * (models[[2L]]$log_likelihood - models[[1L]]$log_likelihood)
    # The statistic reads only the maximized likelihoods. A covariate refit
    # whose likelihood has settled but whose membership logits still drift
    # (an empty or separated class in an over-fitted alternative) sits at
    # the likelihood's supremum, so it is a valid replicate; the logits'
    # state is recorded rather than silently accepted.
    at_maximum <- function(model) {
      isTRUE(model$converged) ||
        (covariate_family && isTRUE(model$likelihood_converged))
    }
    logits_settled <- all(vapply(models, function(model) isTRUE(model$converged),
                                 logical(1)))
    valid <- all(vapply(models, at_maximum, logical(1))) &&
      statistic >= -reversal_window
    list(statistic = statistic, valid = valid,
         boundary = any(vapply(models, `[[`, logical(1), "boundary")),
         logits_settled = logits_settled,
         # Covariate fits do not count replicated maxima.
         null_replications = models[[1L]]$n_best_replicated %||% NA_integer_,
         alternative_replications = models[[2L]]$n_best_replicated %||% NA_integer_)
  }, muffle_warnings = FALSE)
  tally <- .latents_blrt_p_value(replicates, observed, iter)
  result <- list(statistic = observed, p_value = tally$p_value,
       monte_carlo_se = tally$monte_carlo_se,
       iter = iter, n_valid = tally$n_valid, replicates = replicates,
       fixed = constraint$fixed,
       null_profiles = null_model$n_profiles,
       null_group_classes = null_model$n_group_classes,
       alternative_profiles = alternative_model$n_profiles,
       alternative_group_classes = alternative_model$n_group_classes,
       call = match.call())
  class(result) <- "multilpa_bootstrap_lrt"
  result
}

#' @rdname latents-print
#' @export
print.multilpa_bootstrap_lrt <- function(x, ...) {
  stopifnot("`x` must be an `multilpa_bootstrap_lrt` result" =
              inherits(x, "multilpa_bootstrap_lrt"))
  cat(sprintf("Parametric bootstrap likelihood-ratio comparison\n"))
  count <- function(n, word) sprintf("%d %s%s", n, word, if (n == 1L) "" else
    if (endsWith(word, "s")) "es" else "s")
  model <- function(profiles, group_classes, family) {
    # A group-class family has no profiles; name the family instead.
    if (!is.null(family) && !identical(family, "profiles")) {
      return(paste0(family, ", ", count(group_classes, "group class")))
    }
    paste0(count(profiles, "profile"),
           if (x$null_group_classes == 1L && x$alternative_group_classes == 1L) ""
           else paste0(", ", count(group_classes, "group class")))
  }
  cat(sprintf("Null: %s; alternative: %s\n",
              model(x$null_profiles, x$null_group_classes, x$null_family),
              model(x$alternative_profiles, x$alternative_group_classes,
                    x$alternative_family)))
  cat(sprintf("Observed statistic: %.6f\n", x$statistic))
  cat(sprintf("p-value: %s (Monte Carlo SE %s) from %d of %d valid replicates\n",
              format(x$p_value, digits = 4L), format(x$monte_carlo_se, digits = 3L),
              x$n_valid, x$iter))
  if (length(x$fixed %||% character()) > 0L) {
    cat(sprintf(paste("Every fit held %s fixed at one measurement solution,",
                      "so the p-value is conditional on it.\n"),
                paste(x$fixed, collapse = ", ")))
  }
  if (is.na(x$p_value)) {
    cat(paste("The p-value is withheld because not every replicate was valid;",
              "get_results(x, \"replicates\") lists them.\n"))
  }
  invisible(x)
}

#' @rdname latents-summary
#' @export
summary.multilpa_bootstrap_lrt <- function(object, ...) {
  stopifnot("`object` must be an `multilpa_bootstrap_lrt` result" =
              inherits(object, "multilpa_bootstrap_lrt"))
  result <- list(test = .multilpa_bootstrap_test_frame(object),
                 replicates = object$replicates,
                 n_boundary = sum(object$replicates$boundary %in% TRUE),
                 n_errors = sum(!is.na(object$replicates$error)),
                 call = object$call)
  # Every table the fit can produce, built once here, so `get_results()` on the
  # summary serves the same tables the fit would and `print()` can show them
  # all without recomputing anything.
  result$tables <- get_results(object, "all")
  class(result) <- "summary_multilpa_bootstrap_lrt"
  result
}

#' The bootstrap test as a single tidy row
#' @param x An `multilpa_bootstrap_lrt` result.
#' @return A one-row `data.frame` describing the comparison.
#' @noRd
.multilpa_bootstrap_test_frame <- function(x) {
  data.frame(null_profiles = x$null_profiles,
             null_group_classes = x$null_group_classes,
             alternative_profiles = x$alternative_profiles,
             alternative_group_classes = x$alternative_group_classes,
             statistic = x$statistic, p_value = x$p_value,
             monte_carlo_se = x$monte_carlo_se,
             iter = x$iter, n_valid = x$n_valid,
             # Named rather than flagged, because which blocks were held is what
             # the conditional p-value is conditional on.
             fixed = if (length(x$fixed %||% character()) == 0L) NA_character_
               else paste(x$fixed, collapse = ", "),
             null_family = x$null_family %||% "profiles",
             alternative_family = x$alternative_family %||% "profiles",
             row.names = NULL, stringsAsFactors = FALSE)
}

#' @rdname latents-print
#' @export
print.summary_multilpa_bootstrap_lrt <- function(x, digits = 4L, rows = 10L,
                                                 ...) {
  stopifnot("`x` must be a `summary_multilpa_bootstrap_lrt` object" =
              inherits(x, "summary_multilpa_bootstrap_lrt"))
  .multilpa_check_print_arguments(digits, rows)
  cat("Parametric bootstrap likelihood-ratio comparison\n")
  cat(sprintf("%d replicate(s) reached a parameter boundary; %d raised an error.\n",
              x$n_boundary, x$n_errors))
  cat("The reference distribution is simulated, not chi-square.\n")
  .multilpa_print_tables(x$tables, rows = rows, digits = digits)
  .multilpa_print_table_footer(x$tables)
  invisible(x)
}

#' @rdname latents-as-data-frame
#' @export
as.data.frame.multilpa_bootstrap_lrt <- function(x, row.names = NULL, optional = FALSE, ...) {
  stopifnot("`x` must be an object of class `multilpa_bootstrap_lrt`" = inherits(x, "multilpa_bootstrap_lrt"))
  .multilpa_coerce(x, row.names, list(...))
}

#' @rdname latents-as-data-frame
#' @export
as.data.frame.summary_multilpa_bootstrap_lrt <- function(x, row.names = NULL, optional = FALSE, ...) {
  stopifnot("`x` must be an object of class `summary_multilpa_bootstrap_lrt`" = inherits(x, "summary_multilpa_bootstrap_lrt"))
  .multilpa_coerce(x, row.names, list(...))
}

#' Plot a simulated bootstrap null distribution
#'
#' Draws the replicate likelihood-ratio statistics as a histogram with the
#' observed statistic marked, so that the p-value can be read as a tail area of
#' the distribution that was actually simulated. Invalid replicates carry no
#' statistic and are counted in the subtitle rather than dropped silently.
#'
#' @param x An `multilpa_bootstrap_lrt` result.
#' @param main,subtitle Panel title and secondary line, or `NULL` for defaults.
#' @param palette A vector of at least two colours: the histogram fill and the
#'   observed-statistic marker. `NULL` uses the package palette.
#' @param style A list of visual constants, as built by `.multilpa_style()`.
#' @param ... Further named visual constants, merged into `style`.
#' @return The input, invisibly. Called for the side effect of drawing.
#' @section Conditions:
#'   `latents_nothing_to_plot` when no replicate produced a finite statistic.
#' @examples
#' set.seed(1)
#' example_data <- data.frame(
#'   school = rep(seq_len(10), each = 10),
#'   score = rnorm(100, rep(c(-2, 2), each = 50))
#' )
#' smaller <- multilpa(example_data, "score", "school", n_profiles = 1,
#'                     n_group_classes = 1, n_starts = 2, seed = 1)
#' larger <- multilpa(example_data, "score", "school", n_profiles = 2,
#'                    n_group_classes = 1, n_starts = 2, seed = 1)
#' # `iter` is small so the example runs quickly; use many more for inference.
#' comparison <- bootstrap_lrt(smaller, larger, iter = 9, n_starts = 1,
#'                             max_iter = 2000, tol = 1e-6, seed = 1)
#' plot(comparison)
#' @export
plot.multilpa_bootstrap_lrt <- function(x, main = NULL, subtitle = NULL,
                                        palette = NULL,
                                        style = .multilpa_style(), ...) {
  stopifnot("`x` must be an `multilpa_bootstrap_lrt` result" =
              inherits(x, "multilpa_bootstrap_lrt"))
  values <- x$replicates$statistic
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    stop(errorCondition(
      "No replicate produced a finite statistic; there is nothing to plot.",
      class = "latents_nothing_to_plot", call = NULL))
  }
  style <- utils::modifyList(style, list(...))
  colours <- if (is.null(palette)) .multilpa_palette(2L) else
    rep(palette, length.out = 2L)
  previous <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(previous), add = TRUE, after = FALSE)
  graphics::par(xpd = NA, mar = style$margins)
  breaks <- pretty(range(c(values, x$statistic, 0)), n = 12L)
  counts <- graphics::hist(values, breaks = breaks, plot = FALSE)
  .multilpa_panel(xlim = range(breaks), ylim = c(0, max(counts$counts) * 1.12),
    xlab = "Simulated likelihood-ratio statistic", ylab = "Replicates",
    main = if (is.null(main)) "Simulated null distribution" else main,
    subtitle = if (is.null(subtitle)) sprintf(
      "%d of %d replicates valid; observed %.3f, p = %s", x$n_valid, x$iter,
      x$statistic, format(x$p_value, digits = 3L)) else subtitle,
    style = style)
  graphics::rect(counts$breaks[-length(counts$breaks)], 0,
                 counts$breaks[-1L], counts$counts, col = colours[1L],
                 border = style$panel_fill, lwd = 1.2)
  graphics::abline(v = x$statistic, col = colours[2L], lwd = style$line_width,
                   lty = 2L)
  graphics::text(x$statistic, max(counts$counts) * 1.06, " observed",
                 adj = c(0, 0.5), col = colours[2L],
                 cex = style$label_text_size, font = 2L)
  invisible(x)
}

#' Spread the information criteria of one candidate into a single row
#' @param fit A fitted `multilpa` model, or `NULL` for a failed candidate.
#' @return A one-row `data.frame` of criteria, all `NA_real_` when `fit` is `NULL`.
#' @noRd
.multilpa_enumeration_indices <- function(fit) {
  names_wanted <- .multilpa_enumeration_criteria()
  if (is.null(fit)) {
    return(as.data.frame(stats::setNames(
      rep(list(NA_real_), length(names_wanted)), names_wanted)))
  }
  ## The same pivot `get_results(what = "information_criteria")` performs in its
  ## wide form, so the
  ## grid and a single fit cannot name the same quantity differently. The wide
  ## form leads with the likelihood and parameter count, which the grid already
  ## carries in its own columns.
  wide <- .multilpa_information_criteria(fit, format = "wide")
  as.data.frame(as.list(stats::setNames(
    lapply(names_wanted, function(name) wide[[name]]), names_wanted)))
}

#' The information criteria an enumeration grid carries
#'
#' Named once so that the failed-candidate row and the fitted row cannot drift
#' apart, and so that a renamed criterion is caught rather than silently
#' turned into a missing column.
#'
#' @return A character vector of grid column names, in grid order.
#' @noRd
.multilpa_enumeration_criteria <- function() {
  c("aic", "kic", "bic_groups", "bic_individual", "sabic_groups",
    "sabic_individual", "caic_groups", "caic_individual", "awe_groups",
    "awe_individual", "icl_groups", "icl_individual", "clc_groups",
    "clc_individual")
}
