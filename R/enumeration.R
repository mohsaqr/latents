utils::globalVariables("statistic")

#' Class enumeration and comparison of mixture models
#'
#' @description
#' `enumerate_classes()` fits every candidate in a grid of class counts
#' (and, for continuous indicators, covariance structures) and collects
#' their information criteria, classification quality and estimation
#' diagnostics in one object. It covers latent profile, latent class and
#' latent transition models, single-level or two-level, and the additive,
#' dispersion and additive-dispersion families of [multilpa()]. The fitted
#' candidates stay in the result (see *Value*).
#'
#' Information criteria summarise relative support; entropy and average
#' posteriors describe classification certainty within each candidate.
#' Neither picks a solution: also check that the best likelihood was
#' replicated across starts, that no estimate sits at a bound, and that an
#' added class brings interpretable differences in the response patterns.
#' More profiles can stand in for covariances a diagonal structure cannot
#' hold, so compare class counts alongside plausible structures. The grid
#' is the analyst's: when every criterion is lowest at its edge, extend it.
#'
#' @param data A data frame with one row per observation: the indicators,
#'   the group identifier and, when used, the occasion column and
#'   covariates. Use the same rows and scale for every candidate.
#' @param vars Names of the indicator columns.
#' @param id The column identifying the groups the rows are nested in.
#'   `NULL` enumerates single-level models, with a message; [enumerate_lpa()]
#'   and [enumerate_lca()] do that directly.
#' @param n_profiles Profile (class) counts to try (default `1:4`; `2:4` for
#'   transition models).
#' @param n_group_classes Group-class counts to try (default `1:3`; default
#'   `1` with `time`, and only `1` with `id = NULL`).
#' @param model Covariance structures to cross with the class counts.
#'   `"basic"` (the default for profile grids) fits EEI (equal variances, no
#'   covariances), VVI (varying variances), EEE (one covariance matrix
#'   shared by all profiles) and VVV (a covariance matrix per profile).
#'   `"all"` fits all fourteen structures.
#'   Otherwise mclust's three-letter codes (volume, shape and orientation,
#'   each equal `E` or varying `V`, `I` for no covariances), or `NULL` for
#'   none, which is the default for a transition grid. With one continuous
#'   indicator, codes that coincide are fitted once; with only categorical
#'   indicators no structure is crossed. Giving `variance_model`,
#'   `covariance_model`, `volume`, `shape` or `orientation` in `...` sets
#'   one structure instead.
#' @param seed Optional integer seed, reproducing the starting values for
#'   the same data and controls.
#' @param family `"profiles"` (default) enumerates the profile model over
#'   `n_profiles`, `n_group_classes` and `model`. One or more of
#'   `"additive"`, `"dispersion"` and `"additive_dispersion"` enumerates
#'   those group-class families over `family`, `n_group_classes` and
#'   `between_variance` instead; they need `id`, refuse `n_profiles` and
#'   `model`, and return a `latents_family_enumeration`.
#' @param between_variance For the group-class families, the restrictions to
#'   cross: `"varying"`, `"equal"` or both (default). The dispersion family
#'   is fitted with `"equal"` only.
#' @param time The occasion column. Enumerates latent transition models
#'   ([lta()]) over `n_profiles` and `n_group_classes`, every other `lta()`
#'   argument in `...` held fixed; `model`, when given, crosses structures.
#'   Profile family only; returns a `latents_transition_enumeration`.
#' @param ... Further arguments of [multilpa()], or of [lta()] with `time`,
#'   used for every candidate: typically `n_starts`, `max_iter`, `missing`,
#'   `categorical`, `centering`, `scaling`, `drop_empty` and `boundary`. The
#'   group-class families take only `n_starts`, `max_iter`, `tol` and
#'   `min_variance`.
#' @param cores Maximum number of parallel workers (default 1, serial). With
#'   `seed`, parallel and serial runs give the same result. On Windows, and
#'   on macOS with Apple's Accelerate BLAS, each worker is a separate R
#'   process, so memory grows with `cores`.
#' @param type `"lca"`, `"lpa"` or `"lta"`; `NULL` (default) reads the type
#'   from the other arguments. `"lca"` treats every indicator as categorical,
#'   numeric codes included, as [multilca()] does (indicators in `ordinal` or
#'   `count` keep those distributions, and `categorical` adds to the set).
#'   `"lpa"` needs continuous indicators and no `categorical`, `ordinal` or
#'   `count`. Both need `time = NULL`. `"lta"` needs `time` and reads the
#'   measurement from `categorical`, as [lta()] does, without making numeric
#'   indicators categorical. An explicit type needs `family = "profiles"`.
#' @param verbose Whether to report progress: one message each when the
#'   candidates are prepared, when starts are screened, when they are fitted
#'   and when the grid is done. Default on when `cores > 1`. Not reported
#'   for the group-class families.
#'
#' @return A `multilpa_enumeration` object, or the transition or family
#'   enumeration named under `time` and `family`; the print header names
#'   the resolved type (LPA, LCA, LTA, mixed or family). On every grid:
#'   [get_candidates()] gives one row per candidate (named like
#'   `"profiles3_groups2"`) with its class counts, log likelihood, criteria,
#'   entropy, convergence, boundary and replication, leaving out columns
#'   empty in every row; [get_best()] the best candidate per criterion;
#'   [get_tree()] and [plot_tree()] how classes split; [get_entropy()],
#'   [get_avepp()] and [get_classification()] the classification tables of
#'   every candidate; [plot_criteria()] the criteria; [get_candidate()] one
#'   fitted candidate; [compare_models()] candidates by name (a family grid
#'   is compared in its own candidates table). A profile or
#'   class grid also has the AvePP summaries at each level in its candidates
#'   table, [get_deltas()] and [get_knees()], and a `summary()` that prints
#'   one screen of the evidence; a transition or family grid's `summary()`
#'   is its candidates table. Invalid arguments are refused with the
#'   classed conditions listed on [latents-conditions].
#'
#' @seealso [get_candidate()], [bootstrap_lrt()] to test one class count
#'   against the next, [enumerate_regressions()] for mixture regressions.
#' @examples
#' set.seed(1)
#' d <- data.frame(g = rep(1:10, each = 10), y = rnorm(100))
#' candidates <- enumerate_classes(d, "y", "g", n_profiles = 1:2,
#'                               n_group_classes = 1, n_starts = 2, seed = 1)
#' as.data.frame(candidates)
#' summary(candidates)
#'
#' # Numeric response codes, explicitly categorical:
#' codes <- data.frame(g = rep(seq_len(10), each = 10), item = rep(0:1, 50))
#' classes <- enumerate_classes(codes, "item", "g", type = "lca",
#'                              n_profiles = 1, n_group_classes = 1,
#'                              n_starts = 1, seed = 1)
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
                              time = NULL, ..., cores = 1L, type = NULL,
                              verbose = cores > 1L) {
  data_expression <- substitute(data)
  # A socket worker receives the candidates' closure serialized; an unforced
  # `data` would travel as a promise naming the caller's variable, which the
  # worker cannot find. The transition and family paths never touch `data`
  # before the workers start, so it is forced here for every path.
  force(data)
  cores <- .latents_check_cores(cores)
  stopifnot("`verbose` must be TRUE or FALSE" =
              is.logical(verbose) && length(verbose) == 1L && !is.na(verbose))
  if (missing(id)) {
    .multilpa_missing_id("enumerate_classes",
                         "`enumerate_lpa()` or `enumerate_lca()`")
  }
  # Every candidate's refusal is recorded as that candidate's failure, so a
  # seed set.seed() would not take must be refused here, once, or it comes
  # back as a grid in which every candidate failed.
  .multilpa_check_seed(seed)
  type <- .latents_check_enumeration_type(type, time, family)
  if (!is.null(time)) {
    # Transition models: profiles x group classes (x structures when `model`
    # names codes); every other lta() argument is held fixed through `...`.
    if (!identical(family, "profiles") || !missing(between_variance)) {
      stop(errorCondition(
        "Transition models (`time`) are enumerated for the profile family only.",
        class = "latents_bad_argument", call = NULL))
    }
    result <- .latents_report_once(.lta_enumerate(data, vars, id, time,
                          if (missing(n_profiles)) 2:4 else n_profiles,
                          if (missing(n_group_classes)) 1L else n_group_classes,
                          if (missing(model)) NULL else model, seed, list(...),
                          match.call(), cores = cores, verbose = verbose))
    result$type <- "lta"
    return(result)
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
    result <- .latents_report_once(.additive_enumerate(data, vars, id, family,
                                  n_group_classes, between_variance, seed, list(...),
                                  match.call(), cores = cores))
    result$type <- "family"
    return(result)
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
  .latents_check_data_frame(data, substitute(data))
  .latents_require(
            "`vars` must be a character vector of column names" = is.character(vars),
            "`id` must be a single column name, or NULL for single-level models" =
              single_level || (is.character(id) && length(id) == 1L),
            "`n_profiles` must be one or more positive whole numbers" =
              is.numeric(n_profiles) && length(n_profiles) > 0L &&
              all(is.finite(n_profiles)) && all(n_profiles >= 1) &&
              all(n_profiles == as.integer(n_profiles)),
            "`n_group_classes` must be one or more positive whole numbers" =
              is.numeric(n_group_classes) && length(n_group_classes) > 0L &&
              all(is.finite(n_group_classes)) && all(n_group_classes >= 1) &&
              all(n_group_classes == as.integer(n_group_classes)),
            "`model` must be \"basic\", \"all\", covariance model codes, or NULL" =
              is.null(model) || (is.character(model) &&
                                   length(model) > 0L &&
                                   !anyNA(model)))
  # Every name this guard used to reject is now a formal, so R refuses the call
  # with "matched by multiple actual arguments" before `...` is assembled.
  extra <- list(...)
  if (identical(type, "lca")) {
    extra$categorical <- .latents_lca_categorical(vars, extra)
  }
  if (identical(type, "lpa")) {
    conflicts <- c("categorical", "ordinal", "count")
    conflicts <- conflicts[vapply(extra[conflicts], length, integer(1)) > 0L]
    if (length(conflicts)) {
      stop(errorCondition(sprintf(paste(
        "`type = \"lpa\"` requires continuous indicators and conflicts with %s.",
        "Remove those measurement arguments, or use `type = NULL` for a mixed model."),
        paste(sprintf("`%s`", conflicts), collapse = ", ")),
        class = "latents_bad_argument", call = NULL))
    }
    # Empty declarations carry no distribution under an all-continuous type.
    extra[intersect(names(extra), c("categorical", "ordinal", "count"))] <- NULL
  }
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
  indices <- seq_len(nrow(grid))
  requested <- lapply(grid$model, .multilpa_enumeration_arguments,
                       vars = vars, extra = extra)
  arguments <- lapply(indices, function(i) c(
    list(data = data, vars = vars, id = id, n_profiles = grid$n_profiles[i],
         n_group_classes = grid$n_group_classes[i], seed = seed),
    requested[[i]], extra))
  calls <- lapply(indices, function(i) .latents_candidate_call(
    "multilpa", data_expression, arguments[[i]][setdiff(names(arguments[[i]]), "data")]))
  costs <- grid$n_group_classes * grid$n_profiles *
    .multilpa_pool_width(data, vars, extra)
  fitted <- .multilpa_enumerate_fits(arguments, calls, costs, cores, verbose)
  .latents_refuse_shared_data_error(lapply(fitted, `[[`, "value"))
  # The candidates' conditions come back deferred and are signalled here,
  # candidate by candidate.
  runs <- .latents_report_once(lapply(indices, function(i) {
    conditions <- Filter(function(condition) !inherits(condition, "latents_single_level"),
                          fitted[[i]]$conditions)
    .latents_signal_conditions(conditions)
    warnings <- vapply(Filter(function(condition) inherits(condition, "warning"),
                               conditions), conditionMessage, character(1))
    failed <- inherits(fitted[[i]]$value, "error")
    fit <- if (failed) NULL else fitted[[i]]$value
    error_text <- if (failed) conditionMessage(fitted[[i]]$value) else NA_character_
    row <- cbind(
      data.frame(n_profiles = grid$n_profiles[i], n_group_classes = grid$n_group_classes[i],
        # A latent class model has no covariance structure; naming one (the
        # default the fit records for its absent Gaussian block) misleads.
        model = if (is.null(fit)) grid$model[i] else
          if (length(.multilpa_continuous_names(fit)) == 0L) NA_character_ else
            .latents_comparison_structure(fit),
        log_likelihood = if (is.null(fit)) NA_real_ else fit$log_likelihood,
        n_parameters = if (is.null(fit)) NA_integer_ else fit$n_parameters,
        stringsAsFactors = FALSE),
      .multilpa_enumeration_indices(fit),
      data.frame(
        profile_entropy = if (is.null(fit)) NA_real_ else
          .multilpa_level_relative_entropy(fit, "individuals"),
        group_entropy = if (is.null(fit)) NA_real_ else
          .multilpa_level_relative_entropy(fit, "groups"),
        converged = !is.null(fit) && fit$converged,
        boundary = if (is.null(fit)) NA else .latents_boundary_active(fit),
        min_standardized_eigenvalue = .latents_gaussian_minimum(
          fit, "min_standardized_eigenvalue"),
        min_gaussian_mass = .latents_gaussian_minimum(fit, "gaussian_mass"),
        n_best_replicated = if (is.null(fit)) NA_integer_ else fit$n_best_replicated %||% NA_integer_,
        warnings = paste(unique(warnings), collapse = "; "), error = error_text))
    list(fit = fit, row = row)
  }))
  if (single_level) {
    .multilpa_single_level_notice(
      latent_class = all(vars %in% extra$categorical))
  }
  result <- list(table = do.call(rbind, lapply(runs, `[[`, "row")),
       fits = lapply(runs, `[[`, "fit"), call = match.call())
  class(result) <- "multilpa_enumeration"
  result$type <- type %||% .latents_enumeration_measurement_type(vars, extra)
  names(result$fits) <- .latents_candidate_names(result)
  result
}

#' A candidate's call as the caller would have written it
#'
#' A fit made through `do.call(multilpa, arguments)` records the function's
#' whole body and every argument's value, the data frame included, as its
#' `call`, so a grid held one copy of the data in every candidate. A fit's
#' call is only displayed, never evaluated; this one names the verb and the
#' caller's expression for the data, with the other arguments written out.
#'
#' @param verb The fitting function's name.
#' @param data_expression The caller's expression for `data`.
#' @param arguments The other arguments the candidate was fitted with.
#' @return A call.
#' @noRd
.latents_candidate_call <- function(verb, data_expression, arguments) {
  as.call(c(as.name(verb), list(data = data_expression), arguments))
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
  # Covariances for chosen pairs set the covariance themselves; the choice
  # left is the diagonal they extend, one shared set (EEI) or one per profile
  # (VVI), which `.multilpa_enumeration_arguments()` passes as
  # `variance_model`.
  if (.multilpa_enumeration_patterned(vars, extra)) {
    # A `variance_model` the caller gave is the choice, made once; a `model`
    # beside it contradicts it, as it does without pairs.
    if ("variance_model" %in% names(extra)) {
      if (model_given && !is.null(model)) {
        stop(errorCondition(paste(
          "`model` and `variance_model` both set the covariance structure;",
          "name the structures with `model` alone."),
          class = "latents_bad_argument", call = NULL))
      }
      return(NA_character_)
    }
    # `model = NULL` crosses no structure, with pairs as without.
    if (model_given && is.null(model)) return(NA_character_)
    if (!model_given) return(c("EEI", "VVI"))
    codes <- unique(unlist(lapply(model, \(name) switch(name, basic = , all = c("EEI", "VVI"),
                                                         name)), use.names = FALSE))
    if (!all(codes %in% c("EEI", "VVI"))) {
      stop(errorCondition(paste(
        "With `local_dependence` between continuous indicators the covariance",
        "is the pattern of chosen pairs; `model` can choose only the diagonal",
        "it extends, \"EEI\" or \"VVI\"."),
        class = "latents_bad_argument", call = NULL))
    }
    return(codes)
  }
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

#' Retrieve a fitted model from an enumeration or comparison
#'
#' @description
#' `get_candidate()` retrieves a specified fitted model from an enumeration
#' grid or model comparison for detailed inspection and subsequent analysis.
#' The returned object retains its original model class and supports the
#' summary, diagnostic, plotting, and inference methods available for that
#' model.
#'
#' A candidate can be identified by its name or, for supported enumeration
#' grids, by its class counts and model specification. When several
#' candidates share the same class counts, the covariance structure or
#' candidate name distinguishes them. The function retrieves the stored fit
#' without repeating estimation.
#'
#' @details
#' A candidate is selected by its name, as the grid's `candidate` column shows
#' it (`get_candidate(grid, "profiles3_groups4")`; `classes` instead of
#' `profiles` when every indicator is categorical); by its class counts,
#' `get_candidate(grid, 3, 4)` or `get_candidate(grid, c(3, 4))`; or, on a grid
#' that crosses covariance structures, by its counts and `model`. From a
#' comparison, give the row label: `get_candidate(comparison,
#' "profiles3_groups4")`.
#'
#' @param x An enumeration grid from [enumerate_classes()] (a profile, latent
#'   class, transition or group-class family grid), or a comparison from
#'   [compare_models()]. On a family grid, give the name, or `n_group_classes`
#'   with `model`. A grid of mixture regressions from
#'   [enumerate_regressions()] is not supported (`latents_bad_argument`):
#'   `get_results(grid, "model", n_classes = )` returns one of its models.
#' @param n_profiles Number of individual profiles of the candidate, the
#'   candidate's name (`"profiles3_groups4"`), or its two class counts
#'   (`c(3, 4)`). For a comparison, the row label, as its `model` column shows
#'   it.
#' @param n_group_classes Number of group classes of the candidate. Default 1.
#' @param model The covariance structure of the candidate, needed when
#'   [enumerate_classes()] crossed several structures and the class counts
#'   match several candidates; on a group-class family grid, the family model
#'   as the grid's `model` column shows it.
#' @return The fitted object for the selected candidate, with the class and
#'   methods of its fitting function (for example, `multilpa` for a profile
#'   model).
#' @section Conditions:
#'   `latents_unknown_candidate` for a name or class counts not in the grid,
#'   for counts that match several candidates (the message lists their
#'   names), and for a label that is not a row of a comparison;
#'   `latents_failed_candidate` when the candidate's fit raised an error, so
#'   there is no model to return; `latents_bad_argument` for a name or pair of
#'   counts given together with `n_group_classes` or `model`, for a value that
#'   names no candidate, and for a comparison given no label;
#'   `latents_incomplete_fit` for a comparison that no longer carries its
#'   models. Counts that are not single positive integers, and a `model`
#'   that is not a single string, also raise `latents_bad_argument`.
#' @examples
#' set.seed(1)
#' d <- data.frame(g = rep(1:10, each = 10), y = rnorm(100))
#' candidates <- enumerate_classes(d, "y", "g", n_profiles = 1:2,
#'                               n_group_classes = 1, n_starts = 2, seed = 1)
#' get_candidate(candidates, n_profiles = 2, n_group_classes = 1, model = "VVI")
#' # This grid crosses two structures, so the names carry them:
#' candidates
#' get_candidate(candidates, "profiles2_groups1_VVI")
#' compared <- compare_models(candidates, "profiles1_groups1_VVI",
#'                            "profiles2_groups1_VVI")
#' get_candidate(compared, "profiles2_groups1_VVI")
#' @export
get_candidate <- function(x, n_profiles, n_group_classes = 1L,
                          model = NULL) {
  if (inherits(x, "latents_comparison")) {
    if (missing(n_profiles)) return(.latents_compared_model(x))
    return(.latents_compared_model(x, n_profiles))
  }
  # A grid of mixture regressions counts classes, not profiles and group
  # classes, and has its own accessor; it reached the profile-grid checks
  # below and failed one of them unclassed (review 2026-10-27).
  if (inherits(x, "latents_regression_enumeration")) {
    stop(errorCondition(paste(
      "get_candidate() takes a profile, class, transition or family grid;",
      "take one model from a grid of mixture regressions with",
      "get_results(x, \"model\", n_classes = )."),
      class = "latents_bad_argument", call = NULL))
  }
  .latents_require_fit(x, substitute(x), .latents_enumeration_classes,
                       "an enumeration result (from enumerate_classes(), enumerate_lpa() or enumerate_lca())")
  # Checked once for every grid kind, classed as the other argument checks
  # are: a vector `model` built a message per element on the family grid.
  .latents_require("`model` must be a single model code, or NULL" =
                     is.null(model) || (is.character(model) &&
                                          length(model) == 1L && !is.na(model)))
  # A name, or a pair of counts, names the candidate on its own. (A
  # group-class family grid names candidates without `n_profiles`.)
  if (!missing(n_profiles) && (is.character(n_profiles) || length(n_profiles) == 2L)) {
    if (!missing(n_group_classes) || !is.null(model)) {
      stop(errorCondition(paste(
        "A candidate's name, or its pair of class counts, names it on its",
        "own; give that or the separate counts, not both."),
        class = "latents_bad_argument", call = NULL))
    }
    return(.latents_candidate_at(x, .latents_candidate_row(x, n_profiles)))
  }
  if (inherits(x, "latents_transition_enumeration")) {
    # The same contract as the profile grid below: a fractional count would
    # otherwise be printed truncated, as if a present candidate were absent.
    .latents_check_candidate_counts(n_profiles, n_group_classes)
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
    # A fractional count was reported truncated, as if a present candidate
    # with that many group classes were absent.
    .latents_check_candidate_counts(n_group_classes = n_group_classes)
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
  .latents_check_candidate_counts(n_profiles, n_group_classes)
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
      "%d candidates have %s: %s. Name one, or give its `model`."),
      length(at), described,
      paste(sprintf("\"%s\"", .latents_candidate_names(x)[at]), collapse = ", ")),
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

#' Check the class counts that name a candidate
#' @param n_profiles,n_group_classes The counts [get_candidate()] received
#'   (`n_profiles` 1 for a family grid, whose candidates have none).
#' @return `NULL`, invisibly; `latents_bad_argument` for anything but a
#'   single positive whole number each.
#' @noRd
.latents_check_candidate_counts <- function(n_profiles = 1L, n_group_classes) {
  .latents_require(
    "`n_profiles` must be a single positive integer" =
      is.numeric(n_profiles) && length(n_profiles) == 1L &&
      is.finite(n_profiles) && n_profiles >= 1 &&
      n_profiles == floor(n_profiles),
    "`n_group_classes` must be a single positive integer" =
      is.numeric(n_group_classes) && length(n_group_classes) == 1L &&
      is.finite(n_group_classes) && n_group_classes >= 1 &&
      n_group_classes == floor(n_group_classes))
  invisible(NULL)
}

#' @rdname get_candidate
#' @details `candidate_fit()` is the former name of `get_candidate()`, kept
#'   for compatibility. It takes the same arguments, forwards to
#'   `get_candidate()` and gives no warning.
#' @export
candidate_fit <- function(x, n_profiles, n_group_classes = 1L, model = NULL) {
  # Forwarded as a call, so an argument the caller left out stays missing;
  # the function itself is inserted, so a `get_candidate` defined by the
  # caller cannot stand in for it.
  # The head is a name bound to the function, and the arguments symbols bound
  # to the caller's values, so the call never deparses a body or the data.
  env <- new.env(parent = parent.frame())
  env$get_candidate <- get_candidate
  call <- .latents_bind_values(match.call(), env, parent.frame())
  call[[1L]] <- quote(get_candidate)
  eval(call, env)
}

#' @noRd
#' @export
summary.multilpa_enumeration <- function(object, ...) {
  stopifnot("`object` must be an `multilpa_enumeration` result" =
              inherits(object, "multilpa_enumeration"))
  # 0.9.21 ignored any argument without a word; it is still summarised, and
  # the warning names what was not used (review 2026-10-27).
  .latents_warn_unused_arguments(list(...), "summary()", .multilpa_summary_advice())
  grid <- .multilpa_enumeration_grid(object)
  result <- list(best = .multilpa_enumeration_minima(grid),
                 candidates = grid,
                 n_candidates = nrow(grid),
                 n_converged = sum(grid$converged),
                 n_failed = sum(is.na(grid$log_likelihood)),
                 n_boundary = sum(grid$boundary %in% TRUE),
                 call = object$call)
  result$criteria <- result$best
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
      candidate = if (is.na(best) || is.null(grid$candidate)) NA_character_ else
        grid$candidate[best],
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

#' The candidates table of an enumeration, with its AvePP columns
#'
#' The average posterior probability (AvePP) summaries at each level are
#' computed here from the stored fits, each time the table is read, rather
#' than when the grid is fitted. Grids fitted before these columns existed
#' carry the same fits, so they report the same numbers without refitting. The
#' columns sit after `group_entropy`, beside the entropies they are read with.
#'
#' @param x A `multilpa_enumeration` result.
#' @param cutoff The AvePP below which a class is counted in `*_avepp_below`.
#' @return `x$table` with `avepp_min`, `_max`, `_mean`, `_below` and
#'   the same four for `group_`; `NA` for a candidate whose fit failed.
#' @noRd
.multilpa_enumeration_grid <- function(x, cutoff = 0.8) {
  stopifnot(
    "`cutoff` must be a single number between 0 and 1" =
      is.numeric(cutoff) && length(cutoff) == 1L && is.finite(cutoff) &&
      cutoff > 0 && cutoff <= 1,
    "the enumeration must hold one fit slot per candidate" =
      length(x$fits) == nrow(x$table))
  statistics <- c("min", "max", "mean", "below")
  added <- c(paste0("avepp_", statistics),
             paste0("group_avepp_", statistics))
  # Names an earlier development build stored in the table.
  stale <- c("profile_min_average_posterior", "group_min_average_posterior")
  grid <- x$table[setdiff(names(x$table), c(added, stale))]
  level_summary <- function(level, prefix) {
    rows <- vapply(unname(x$fits), function(fit) {
      if (is.null(fit)) rep(NA_real_, 4L) else
        .multilpa_level_avepp(fit, level, cutoff)
    }, numeric(4L))
    columns <- as.data.frame(t(rows))
    names(columns) <- paste0(prefix, statistics)
    columns[[paste0(prefix, "below")]] <- as.integer(rows[4L, ])
    columns
  }
  columns <- cbind(
    level_summary("individuals", "avepp_"),
    level_summary("groups", "group_avepp_"))
  at <- match("group_entropy", names(grid))
  .latents_named_table(x, cbind(grid[seq_len(at)], columns,
                                grid[setdiff(seq_along(grid), seq_len(at))]))
}

#' Smallest, largest and mean AvePP at one level, and the classes below a cutoff
#'
#' AvePP is the diagonal of the average posterior matrix: for the units
#' modally assigned to a class, their mean posterior of belonging to it.
#' A class with no assigned unit has no AvePP; the smallest,
#' largest and mean are taken over the classes that have one, and such a class
#' is counted as below the cutoff, since nothing is classified into it. A level
#' with no posteriors or a single class classifies nothing and gives `NA`s,
#' the same rule as relative entropy.
#'
#' @param probabilities Posterior probability matrix, or `NULL`.
#' @param cutoff The AvePP below which a class is counted.
#' @param weights `NULL`, or one sampling weight per row: the AvePP is then
#'   the weighted mean, as in the `"classification"` table.
#' @return A numeric vector: `min`, `max`, `mean`, `below`.
#' @noRd
.multilpa_avepp_summary <- function(probabilities, cutoff, weights = NULL) {
  if (is.null(probabilities) || ncol(probabilities) < 2L) {
    return(c(min = NA_real_, max = NA_real_, mean = NA_real_, below = NA_real_))
  }
  avepp <- diag(.multilpa_average_posterior_matrix(probabilities, weights))
  assigned <- avepp[!is.na(avepp)]
  c(min = min(assigned), max = max(assigned), mean = mean(assigned),
    below = sum(is.na(avepp) | avepp < cutoff))
}

#' One fit's classification table with its relative entropy, one-class levels
#' left out
#'
#' @param fit A fitted model of this package.
#' @return The fit's `"classification"` rows at every level with two or more
#'   classes, with that level's `relative_entropy`.
#' @noRd
.multilpa_informative_classification <- function(fit) {
  quality <- .multilpa_classification_table(
    fit, level = .multilpa_default_level(fit, NULL))
  entropy <- .multilpa_entropy_table(fit)
  # A level with one class classifies nothing: its AvePP is 1 by
  # construction and would only pad the table.
  informative <- entropy$level[entropy$n_classes > 1L]
  quality <- quality[quality$level %in% informative, , drop = FALSE]
  quality$relative_entropy <-
    entropy$relative_entropy[match(quality$level, entropy$level)]
  quality
}

#' Per-class classification quality of every candidate in a grid
#'
#' The `"classification"` table of each fitted candidate, keyed by the
#' candidate's class counts and model and carrying that level's relative
#' entropy, so AvePP and entropy are compared across the grid in one table.
#' Candidates whose fit failed have no posteriors and contribute no rows, and
#' neither does a level with a single class.
#'
#' @param x A `multilpa_enumeration` result.
#' @return A base `data.frame`, one row per candidate, level and class.
#' @noRd
.multilpa_enumeration_classification <- function(x) {
  fitted <- which(!vapply(x$fits, is.null, logical(1L)))
  if (length(fitted) == 0L) {
    return(data.frame(candidate = character(), n_profiles = integer(),
                      n_group_classes = integer(), model = character(),
                      level = character(), class = integer(), n_modal = integer(),
                      proportion_modal = numeric(), estimated_n = numeric(),
                      estimated_proportion = numeric(), average_posterior = numeric(),
                      odds_correct_classification = numeric(),
                      relative_entropy = numeric()))
  }
  rows <- lapply(fitted, function(i) {
    quality <- .multilpa_informative_classification(x$fits[[i]])
    keys <- x$table[rep(i, nrow(quality)),
                    c("n_profiles", "n_group_classes", "model"), drop = FALSE]
    cbind(candidate = rep(.latents_candidate_names(x)[i], nrow(quality)),
          keys, quality)
  })
  result <- do.call(rbind, rows)
  row.names(result) <- NULL
  result
}

#' @noRd
#' @export
print.summary_multilpa_enumeration <- function(x, digits = 4L, rows = 10L, ...) {
  stopifnot("`x` must be a `summary_multilpa_enumeration` object" =
              inherits(x, "summary_multilpa_enumeration"))
  .multilpa_check_print_arguments(digits, rows)
  cat(sprintf("Class enumeration: %d candidates, %d converged, %d failed to fit\n",
              x$n_candidates, x$n_converged, x$n_failed))
  identified <- c("n_profiles", "n_group_classes", "model")
  best <- x$best %||% x$criteria
  present <- best[!is.na(best$n_profiles), identified, drop = FALSE]
  disagreement <- unique(present)
  cat(sprintf("%d distinct candidate(s) are minimal under some criterion.\n",
              nrow(disagreement)))
  cat("No candidate is selected automatically. Choose one convention and keep it.\n")
  grid <- x$candidates
  # What a reader must know before reading any minimum: whether the minima lie
  # inside the grid, whether the fits behind them were found more than once,
  # and whether their classes are separated. The tables themselves have their
  # own verbs; the print is the reading of them (review 2026-11-02: the print
  # was seven truncated tables and pointed to a 5,000-line list).
  .multilpa_print_enumeration_warnings(grid, best)
  shown <- best[, intersect(c("criterion", "convention", "candidate", "value"),
                            names(best)), drop = FALSE]
  cat("\n-- lowest value of each criterion ", strrep("-", 32L), "\n", sep = "")
  # Every criterion is shown by default (the table is the reason to print);
  # a `rows` the caller gives limits it as it limits every other table, and
  # `rows = 0` asks for the headings alone.
  limit <- if (missing(rows)) nrow(shown) else rows
  if (limit > 0) print(utils::head(shown, limit), digits = digits, row.names = FALSE)
  if (nrow(shown) > limit) cat(sprintf("   ... %d more rows.\n", nrow(shown) - limit))
  knees <- x$tables$knees
  if (is.data.frame(knees) && nrow(knees) > 0L) {
    cat("\n-- knees of ", knees$criterion[1L], " ", strrep("-", 40L), "\n", sep = "")
    if (rows > 0) {
      print(utils::head(knees[, intersect(c("direction", "line", "knee_candidate",
                                             "knee_without_largest", "stable"),
                                           names(knees)), drop = FALSE], rows),
            row.names = FALSE)
    }
    if (nrow(knees) > rows) cat(sprintf("   ... %d more rows.\n", nrow(knees) - rows))
  }
  cat(paste0("\nTables: get_candidates(x), get_best(x), get_deltas(x), get_knees(x),\n",
             "        get_tree(x), get_classification(x), get_entropy(x), get_checks(x).\n",
             "Plots:  plot_criteria(x), plot_deltas(x), plot_avepp(x), plot_tree(x), plot_all(x).\n"))
  invisible(x)
}

#' The warnings a candidate grid needs before its minima are read
#'
#' * Every criterion lowest at an edge of the grid: no criterion has a minimum
#'   inside it, so the grid does not bracket one.
#' * Candidates whose best log likelihood was reached by one start only.
#' * Candidates with a class whose average posterior is below 0.8.
#' @param grid The candidates table.
#' @param best The table of the lowest value of each criterion.
#' @return `NULL`, invisibly. Called for the printing.
#' @noRd
.multilpa_print_enumeration_warnings <- function(grid, best) {
  lines <- character()
  converged <- grid[grid$converged %in% TRUE, , drop = FALSE]
  chosen <- unique(best$candidate[!is.na(best$candidate)])
  if (length(chosen) > 0L) {
    # The edges are those of the grid that was asked for: a largest count that
    # failed to fit does not make the next one the edge, and a minimum at the
    # largest count is at the edge even when the others failed.
    edge <- function(column) {
      values <- grid[[column]]
      if (is.null(values) || length(unique(values)) < 2L) return(rep(FALSE, length(chosen)))
      at <- grid[[column]][match(chosen, grid$candidate)]
      at == max(values) | (at == min(values) & min(values) > 1L)
    }
    at_edge <- edge("n_profiles") | edge("n_group_classes")
    if (all(at_edge %in% TRUE)) {
      lines <- c(lines, sprintf(paste(
        "Every criterion is lowest at the edge of the grid (%s): the grid does",
        "not bracket a minimum. Extend it, or choose on the knees and the",
        "classification quality."), paste(chosen, collapse = ", ")))
    }
  }
  once <- converged$candidate[converged$n_best_replicated %in% 1L]
  if (length(once) > 0L) {
    named <- if (length(once) > 6L) c(utils::head(once, 6L), "...") else once
    lines <- c(lines, sprintf(paste(
      "%d of %d candidates reached their best log likelihood from one start",
      "only (%s); refit those you consider with more starts."),
      length(once), nrow(converged), paste(named, collapse = ", ")))
  }
  below <- (converged$avepp_below %||% 0) > 0 | (converged$group_avepp_below %||% 0) > 0
  if (any(below %in% TRUE)) {
    lines <- c(lines, sprintf(
      "%d candidates have a class with average posterior below 0.8 (%s).",
      sum(below %in% TRUE), paste(converged$candidate[below %in% TRUE], collapse = ", ")))
  }
  if (length(lines) == 0L) return(invisible(NULL))
  cat("\n")
  cat(paste0("! ", strwrap(lines, width = 76L, exdent = 2L, simplify = FALSE) |>
               vapply(paste, character(1), collapse = "\n  ")), sep = "\n")
  invisible(NULL)
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
  keep <- c("candidate", "n_profiles",
            if (varies("n_group_classes")) "n_group_classes",
            if (varies("model")) "model",
            "log_likelihood", "n_parameters", "aic", bic, "icl_individual",
            "profile_entropy",
            if (!all(is.na(grid$group_entropy))) "group_entropy",
            paste0("avepp_", c("min", "mean", "below")),
            if (!all(is.na(grid$group_avepp_min)))
              paste0("group_avepp_", c("min", "mean", "below")),
            "converged", "boundary")
  compact <- grid[, intersect(keep, names(grid)), drop = FALSE]
  if (same_levels("bic")) names(compact)[names(compact) == "bic_groups"] <- "bic"
  compact
}

#' @noRd
#' @export
as.data.frame.summary_multilpa_enumeration <- function(x, row.names = NULL, optional = FALSE, ...) {
  stopifnot("`x` must be an object of class `summary_multilpa_enumeration`" = inherits(x, "summary_multilpa_enumeration"))
  .multilpa_coerce(x, row.names, list(...))
}

#' Does an enumeration fit covariances for chosen continuous pairs?
#' @param vars The indicators.
#' @param extra The other arguments of the enumeration.
#' @return A single logical.
#' @noRd
.multilpa_enumeration_patterned <- function(vars, extra) {
  pairs <- extra$local_dependence %||% character()
  if (length(pairs) == 0L) return(FALSE)
  items <- .latents_parse_pairs(pairs)$item_1
  any(!items %in% c(extra$categorical, extra$ordinal, extra$count))
}

#' The `multilpa()` arguments of one enumerated structure
#' @param code A structure code, or `NA`.
#' @param vars,extra As for `.multilpa_enumeration_patterned()`.
#' @return A named list.
#' @noRd
.multilpa_enumeration_arguments <- function(code, vars, extra) {
  if (.multilpa_enumeration_patterned(vars, extra) && !is.na(code)) {
    return(list(variance_model = if (identical(code, "EEI")) "equal" else "varying"))
  }
  .multilpa_structure_arguments(code)
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
  # The model conditions on the covariates of its direct effects, so a draw
  # keeps them as observed.
  covariates <- object$extra_data$direct$covariates
  if (!is.null(covariates)) result[colnames(covariates)] <- as.data.frame(covariates)
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
    # A direct effect moves each row's mean by its covariates.
    direct <- object$extra_data$direct
    if (!is.null(direct)) {
      values <- values + tcrossprod(direct$covariates, .latents_direct_matrix(
        direct, object$direct_effects, dimension))
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
    result <- cbind(result, as.data.frame(drawn, stringsAsFactors = FALSE,
                                         check.names = FALSE))
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
#' A covariate both models name must also enter both with the same values:
#' otherwise the two were fitted to different data, and the pair is refused
#' as incomparable before the likelihood check would refuse it as data that
#' do not reproduce a fit (naming `data` the caller may not have supplied).
#'
#' @param null_model,alternative_model Fitted `multilpa_covariates` models.
#' @return `NULL`, invisibly; raises `latents_bad_nesting`, or
#'   `latents_incomparable_models` for a shared covariate with other values.
#' @noRd
.multilpa_check_covariate_nesting <- function(null_model, alternative_model) {
  refuse <- function(message, class = "latents_bad_nesting") {
    stop(errorCondition(message, class = class, call = NULL))
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
  # The values, judged as compare_models(test = "lrt") judges them (column
  # spaces, .latents_lrt_spans()). Group class 1's profile design holds the
  # profile covariates under either slope rule and has one row per
  # observation whatever the class counts, so it compares a profile step and
  # a group-class step alike; the stacked designs differ in rows when the
  # group classes do.
  different_values <- function(inner, outer) {
    !is.null(inner) && !is.null(outer) && !.latents_lrt_spans(inner, outer)
  }
  if (different_values(null_model$profile_design[[1L]],
                       alternative_model$profile_design[[1L]])) {
    refuse(paste(
      "The models name the same profile covariates with different values, so",
      "they were not fitted to the same data. Fit both to one data frame."),
      class = "latents_incomparable_models")
  }
  if (length(null_model$group_covariates) > 0L &&
      different_values(null_model$group_design, alternative_model$group_design)) {
    refuse(paste(
      "The models name the same group covariates with different values, so",
      "they were not fitted to the same data. Fit both to one data frame."),
      class = "latents_incomparable_models")
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

#' Parametric bootstrap likelihood-ratio test for nested mixture models
#'
#' @description
#' `bootstrap_lrt()` tests whether a more flexible model improves on a
#' simpler one nested in it, typically one more profile or group class. The
#' statistic is twice the difference in maximised log likelihoods; its null
#' distribution is simulated by generating data from the fitted null model
#' and refitting both models to each data set, since the chi-squared
#' approximation does not hold when the number of mixture components
#' changes. A small p-value is
#' evidence against the null in favour of this alternative; it does not say
#' the alternative has the right number of classes.
#'
#' @param null_model A converged fit under the null: a profile or class
#'   model with or without membership covariates, an [lta()] fit, or a
#'   group-class family fit from [multilpa()]. See *Supported comparisons*.
#' @param alternative_model A converged, more flexible model fitted to the
#'   same data, in which `null_model` is nested.
#' @param data The data frame both models were fitted to. Required for
#'   transition models, whose simulations keep the observed groups,
#'   occasions and covariates; otherwise taken from the fits when omitted,
#'   and checked against them when supplied.
#' @param iter Number of simulated data sets (default 199; at least 2). More
#'   replicates reduce the Monte Carlo error of the p-value but do not cure
#'   estimation problems.
#' @param n_starts,max_iter,tol Controls for each simulated refit; `NULL`
#'   (default) keeps each model's own controls and search strategy. A
#'   supplied value applies to both models.
#' @param seed Optional integer seed. The caller's random-number state is
#'   restored on exit.
#'
#' @return A `multilpa_bootstrap_lrt` object: the observed statistic, the
#'   bootstrap p-value and its Monte Carlo standard error, `n_valid` and
#'   `n_invalid` replicates, `fixed` (the measurement blocks held in every
#'   refit) and per-replicate diagnostics, including how many starts of
#'   each refit replicated its best likelihood. `as.data.frame()` returns the one-row test
#'   summary, `get_results(result, "replicates")` one row per replicate,
#'   `summary()` both, and `plot()` the simulated distribution with the
#'   observed statistic.
#'
#' @section Supported comparisons:
#' Both models must be of one kind: profile or class models, transition
#' models, or group-class family fits. Profile and class models, with or
#' without membership covariates, must differ by exactly one profile
#' (class) or exactly one group class, and share their observations, group
#' layout, covariance structure, centring and scaling, missing-data
#' handling and `drop_empty`. `centering = "grand"` is repeated in each
#' refit; person-centred or person-scaled fits are refused, because the
#' model gives the group means removed by those transformations no
#' distribution. Covariate models must share their covariates and slope
#' specification (a one-group-class null may be compared with an
#' alternative that adds group covariates or group-class-specific slopes);
#' covariates are held at their observed values when simulating, and both
#' originals must be fully converged, logits included.
#'
#' Transition models may differ in class counts or in restrictions on
#' measurement, transitions, Markov order or covariate effects, with
#' compatible indicators, grouping, occasion layout and measurement
#' distributions. Group-class family fits may differ in the number of group
#' classes or in restrictions on the means and variance components, the
#' alternative being the less restricted one. A profile or class fit with
#' `fixed` (including a [fit_staged()] result) keeps the same held blocks in
#' every replicate, so the p-value is conditional on that measurement
#' solution; a covariate fit's held blocks are re-estimated in each
#' replicate. Fits with a prior, sampling weights or a noise component are
#' refused.
#'
#' @section Calculation:
#' The p-value is `(1 + b) / (V + 1)`, with `b` the number of valid
#' simulated statistics at least as large as the observed one and `V` the
#' number of valid replicates. A replicate is invalid when a refit fails or
#' does not converge, or when the alternative's likelihood falls below the
#' null's beyond numerical tolerance; such replicates are left out with a
#' warning, and the p-value is `NA` when fewer than half of the requested
#' replicates are valid. Leaving them out conditions the reference
#' distribution on successful estimation. Report the invalid share, and look at whether the
#' refits replicated their likelihoods across starts. For covariate models
#' a refit counts as converged when its likelihood has converged even if a
#' membership logit is still moving, as happens with an empty or separated
#' class in the alternative; `logits_settled` records whether the logits
#' converged too.
#'
#' Fits with `missing = "fiml"` are supported when both models use it; each
#' simulated data set keeps the original missing-value pattern. That is
#' calibrated only when missingness is independent of the indicators and
#' the classes given any fixed covariates: a fixed pattern does not
#' reproduce an MAR mechanism that depends on observed indicators, nor an
#' MNAR one, so valid FIML estimation does not by itself make this test
#' valid.
#'
#' @seealso [enumerate_classes()] for fitting the candidates;
#'   [compare_models()] for information-criterion comparisons.
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
                                 iter = 199L, n_starts = NULL, max_iter = NULL,
                                 tol = NULL, seed = NULL) {
  .latents_require_fit(null_model, substitute(null_model), .latents_model_classes)
  .latents_require_fit(alternative_model, substitute(alternative_model),
                       .latents_model_classes)
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
    return(.lta_bootstrap_lrt(null_model, alternative_model, data, iter,
                              n_starts, max_iter, tol, seed, match.call()))
  }
  if (inherits(null_model, "multilpa_additive") ||
      inherits(alternative_model, "multilpa_additive")) {
    return(.additive_bootstrap_lrt(null_model, alternative_model, data, iter,
                                   n_starts, max_iter, tol, seed, match.call()))
  }
  covariate_family <- inherits(null_model, "multilpa_covariates")
  .latents_require(
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
  .latents_check_data_frame(data)
  .latents_require(
    "`iter` must be a single whole number of at least two" =
      is.numeric(iter) && length(iter) == 1L && is.finite(iter) && iter >= 2L &&
      iter == as.integer(iter),
    "`n_starts` must be NULL or a single positive whole number" =
      is.null(n_starts) || (is.numeric(n_starts) && length(n_starts) == 1L &&
        is.finite(n_starts) && n_starts >= 1 && n_starts == as.integer(n_starts)),
    "`max_iter` must be NULL or a single positive whole number" =
      is.null(max_iter) || (is.numeric(max_iter) && length(max_iter) == 1L &&
        is.finite(max_iter) && max_iter >= 1 && max_iter == as.integer(max_iter)),
    "`tol` must be NULL or a single positive number" =
      is.null(tol) || (is.numeric(tol) && length(tol) == 1L && is.finite(tol) &&
        tol > 0))
  # Each model is refitted as it was fitted (its search, starts, stopping
  # rule, acceleration, selection and boundary policy), so the reference
  # distribution comes from the procedure that produced the statistic: a
  # bootstrap that searched less hard than the observed fits rejected a true
  # null at 0.15 instead of 0.05 (SEARCH-PLAN.md).
  controls <- lapply(list(null = null_model, alternative = alternative_model),
                     .latents_refit_controls, n_starts = n_starts,
                     max_iter = max_iter, tol = tol)
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
  scaling_for <- function(model) model$scaling %||% "none"
  drop_empty_for <- function(model) model$drop_empty %||% TRUE
  same_fields <- all(vapply(fields, function(field)
    identical(null_model[[field]], alternative_model[[field]]), logical(1))) &&
    identical(null_model$covariance_pattern, alternative_model$covariance_pattern) &&
    identical(.latents_direct_key(null_model), .latents_direct_key(alternative_model))
  if (!same_fields || !identical(extra_spec(null_model), extra_spec(alternative_model)) ||
      !identical(structure_for(null_model), structure_for(alternative_model)) ||
      !identical(centering_for(null_model), centering_for(alternative_model)) ||
      !identical(scaling_for(null_model), scaling_for(alternative_model)) ||
      !identical(drop_empty_for(null_model), drop_empty_for(alternative_model)) ||
      !identical(missing_for(null_model), missing_for(alternative_model))) {
    stop(errorCondition(paste(
      "Models must use the same observations, group layout, covariance",
      "structure, local dependence, direct effects, centering and scaling",
      "mode, missing-data handling and `drop_empty`."),
        class = "latents_incomparable_models", call = NULL))
  }
  if (covariate_family) .multilpa_check_covariate_nesting(null_model, alternative_model)
  if (identical(centering_for(null_model), "person") ||
      identical(scaling_for(null_model), "person")) {
    stop(errorCondition(paste(
      "A person-centred or person-scaled fit removes each group's baseline or",
      "spread, but this bootstrap does not model those for simulation. Compare",
      "fits without person centering or scaling."),
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
    # Groups are matched by value, as parameter_inference() matches them, so
    # the same identifiers stored as double rather than integer still agree.
    if (!all(c(model$vars, model$id) %in% names(data)) ||
        nrow(data) != model$n_observations ||
        !identical(match(data[[model$id]], model$group_values),
                   model$group_index)) {
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
  reversal_window <- 100 * max(controls$null$tol, controls$alternative$tol) *
    (1 + abs(null_model$log_likelihood))
  observed <- 2 * (alternative_model$log_likelihood - null_model$log_likelihood)
  if (observed < -reversal_window) {
    stop(errorCondition(sprintf(
      paste("Alternative has lower likelihood by %.3g, beyond the %.3g that",
            "`tol = %.3g` can explain; improve its optimization first."),
      -observed, reversal_window, max(controls$null$tol, controls$alternative$tol)),
      class = "latents_reversed_likelihood", call = NULL))
  }
  observed <- max(0, observed)
  replicates <- .latents_blrt_replicates(iter, function(i) {
    simulated <- if (covariate_family) .multilpa_cov_simulate(null_model, data) else
      .multilpa_simulate(null_model)
    if (fiml) {
      simulated <- .multilpa_carry_missingness(simulated, data, null_model$vars)
    }
    models <- Map(function(model, control) {
      ## The covariance and the measurement extensions come from the shared
      ## refit rule: refitting from `variance_model` and `covariance_model`
      ## alone returned a VEI replicate as VVI, two parameters wider, and
      ## dropped direct effects and covariance patterns, so the reference
      ## distribution belonged to a different pair of models than the
      ## statistic compared against it.
      family <- .multilpa_refit_measurement(model)
      do.call(multilpa, c(
        list(data = simulated, vars = model$vars, id = model$id,
             n_profiles = model$n_profiles,
             n_group_classes = model$n_group_classes),
        family,
        list(n_starts = control$n_starts, max_iter = control$max_iter,
             tol = control$tol, search = control$search,
             n_screen = control$n_screen, acceleration = control$acceleration,
             select_start = control$select_start, boundary = control$boundary,
             min_variance = model$min_variance,
             centering = model$centering %||% "none",
             scaling = model$scaling %||% "none",
             drop_empty = model$drop_empty %||% TRUE,
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
    }, list(null_model, alternative_model), controls)
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
       iter = iter, n_valid = tally$n_valid, n_invalid = tally$n_invalid,
       replicates = replicates,
       fixed = constraint$fixed,
       null_profiles = null_model$n_profiles,
       null_group_classes = null_model$n_group_classes,
       alternative_profiles = alternative_model$n_profiles,
       alternative_group_classes = alternative_model$n_group_classes,
       call = match.call(),
       refit_controls = .latents_refit_controls_frame(controls))
  class(result) <- "multilpa_bootstrap_lrt"
  result
}

#' @noRd
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
  n_invalid <- x$n_invalid %||% (x$iter - x$n_valid)
  cat(sprintf("p-value: %s (Monte Carlo SE %s) from %d of %d valid replicates (%d invalid)\n",
              format(x$p_value, digits = 4L), format(x$monte_carlo_se, digits = 3L),
              x$n_valid, x$iter, n_invalid))
  if (length(x$fixed %||% character()) > 0L) {
    cat(sprintf(paste("Every fit held %s fixed at one measurement solution,",
                      "so the p-value is conditional on it.\n"),
                paste(x$fixed, collapse = ", ")))
  }
  if (is.na(x$p_value)) {
    cat(paste("The p-value is withheld because fewer than half of the replicates",
              "were valid; get_results(x, \"replicates\") lists them.\n"))
  } else if (n_invalid > 0L) {
    cat(paste("Replicates with failed or unconverged refits, or an alternative",
              "likelihood below its null beyond tolerance, are left out of",
              "the p-value; get_results(x, \"replicates\") lists them.\n"))
  }
  invisible(x)
}

#' @noRd
#' @export
summary.multilpa_bootstrap_lrt <- function(object, ...) {
  stopifnot("`object` must be an `multilpa_bootstrap_lrt` result" =
              inherits(object, "multilpa_bootstrap_lrt"))
  # 0.9.21 ignored any argument without a word; it is still summarised, and
  # the warning names what was not used (review 2026-10-27).
  .latents_warn_unused_arguments(list(...), "summary()", .multilpa_summary_advice())
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
             n_invalid = x$n_invalid %||% (x$iter - x$n_valid),
             # Named rather than flagged, because which blocks were held is what
             # the conditional p-value is conditional on.
             fixed = if (length(x$fixed %||% character()) == 0L) NA_character_
               else paste(x$fixed, collapse = ", "),
             null_family = x$null_family %||% "profiles",
             alternative_family = x$alternative_family %||% "profiles",
             row.names = NULL, stringsAsFactors = FALSE)
}

#' @noRd
#' @export
print.summary_multilpa_bootstrap_lrt <- function(x, digits = 4L, rows = 10L,
                                                 ...) {
  stopifnot("`x` must be a `summary_multilpa_bootstrap_lrt` object" =
              inherits(x, "summary_multilpa_bootstrap_lrt"))
  .multilpa_check_print_arguments(digits, rows)
  cat("Parametric bootstrap likelihood-ratio comparison\n")
  # `n_errors` counts every invalid replicate: one whose refit raised an
  # error, did not converge, or reversed the likelihood.
  cat(sprintf(paste("%d replicate(s) reached a parameter boundary; %d were",
                    "invalid (an error, no convergence or a reversed",
                    "likelihood).\n"),
              x$n_boundary, x$n_errors))
  cat("The reference distribution is simulated, not chi-square.\n")
  .multilpa_print_tables(x$tables, rows = rows, digits = digits)
  .multilpa_print_table_footer(x$tables)
  invisible(x)
}

#' @noRd
#' @export
as.data.frame.multilpa_bootstrap_lrt <- function(x, row.names = NULL, optional = FALSE, ...) {
  stopifnot("`x` must be an object of class `multilpa_bootstrap_lrt`" = inherits(x, "multilpa_bootstrap_lrt"))
  .multilpa_coerce(x, row.names, list(...))
}

#' @noRd
#' @export
as.data.frame.summary_multilpa_bootstrap_lrt <- function(x, row.names = NULL, optional = FALSE, ...) {
  stopifnot("`x` must be an object of class `summary_multilpa_bootstrap_lrt`" = inherits(x, "summary_multilpa_bootstrap_lrt"))
  .multilpa_coerce(x, row.names, list(...))
}

#' @noRd
#' @export
plot.multilpa_bootstrap_lrt <-function(x, main = NULL, subtitle = NULL,
                                        palette = NULL, style = NULL, ...) {
  stopifnot("`x` must be an `multilpa_bootstrap_lrt` result" =
              inherits(x, "multilpa_bootstrap_lrt"),
            "`palette` must be NULL or a vector of colours" =
              is.null(palette) || (is.character(palette) &&
                                     length(palette) >= 1L && !anyNA(palette)))
  # `style` and named constants through `...` tuned the base-graphics drawing
  # of 0.9.21; the plot is ggplot now, so they have nothing to set.
  .latents_warn_unused_arguments(
    c(if (is.null(style)) list() else list(style = style), list(...)),
    "plot()",
    "Style the returned plot with ggplot2, for example `+ ggplot2::theme()`.")
  colours <- if (is.null(palette)) .gg_okabe_ito[c(1L, 2L)] else
    rep_len(palette, 2L)
  values <- x$replicates$statistic
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    stop(errorCondition(
      "No replicate produced a finite statistic; there is nothing to plot.",
      class = "latents_nothing_to_plot", call = NULL))
  }
  .gg_require()
  breaks <- pretty(range(c(values, x$statistic, 0)), n = 12L)
  # The label sits on the side of the line with more room.
  right_half <- x$statistic > mean(range(breaks))
  drawn <- ggplot2::ggplot(data.frame(statistic = values), ggplot2::aes(statistic)) +
    ggplot2::geom_histogram(breaks = breaks, fill = colours[1L],
                            colour = "white", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = x$statistic, colour = colours[2L],
                        linetype = "dashed", linewidth = 0.9) +
    ggplot2::annotate("text", x = x$statistic, y = Inf,
                      label = if (right_half) "observed " else " observed",
                      hjust = if (right_half) 1 else 0, vjust = 1.5,
                      colour = colours[2L], fontface = "bold", size = 3.4) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
    .gg_theme() +
    ggplot2::labs(x = "Simulated likelihood-ratio statistic", y = "Replicates") +
    .gg_titles(main, subtitle, "Simulated null distribution", sprintf(
      "%d of %d replicates valid; observed %.3f, p = %s", x$n_valid, x$iter,
      x$statistic, format(x$p_value, digits = 3L)))
  .gg_draw_on_call(drawn)
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

#' Raise a data error that failed every candidate of an enumeration
#'
#' A defect in the data (a column that does not exist, an indicator that is
#' constant or entirely missing) fails every candidate for one reason that no
#' choice of model can change. `multilpa()` refuses such data with a classed
#' error; an enumeration used to return a grid of failed cells instead, with
#' no condition, so a misspelled column cost every fit and came back as a
#' table of `NA`s. An invalid argument that every candidate shares (weights,
#' `min_variance`) is refused the same way. Failures that differ between
#' candidates, or that concern the model rather than the data, stay in the
#' grid as before.
#'
#' @param values Each candidate's result: a fit or an error condition.
#' @return `NULL`, invisibly, unless every candidate failed with the same data
#'   error, which is then signalled.
#' @noRd
.latents_refuse_shared_data_error <- function(values) {
  failures <- Filter(function(value) inherits(value, "error"), values)
  if (length(failures) == 0L || length(failures) < length(values)) {
    return(invisible(NULL))
  }
  first <- failures[[1L]]
  # An invalid argument (a `weights` column that is not there, a negative
  # `min_variance`) is shared by every candidate in the same way (API-2).
  data_errors <- c("latents_bad_data", "latents_bad_categorical",
                   "latents_missing_data", "latents_bad_weights",
                   "latents_bad_argument")
  shared <- inherits(first, data_errors) && all(vapply(failures, function(value) {
    identical(class(value), class(first)) &&
      identical(conditionMessage(value), conditionMessage(first))
  }, logical(1)))
  if (shared) stop(.latents_bounded_condition(first))
  invisible(NULL)
}
