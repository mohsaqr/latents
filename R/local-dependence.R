# Local dependence between two categorical indicators: the pair is modelled
# jointly within each class. The two items become one categorical indicator
# whose categories are the full grid of their categories, so each class has its own
# joint distribution of the pair, with no independence assumed between them.
# For two nominal items this is the class-specific two-variable log-linear
# association Latent GOLD fits with `item1 <-> item2 | Cluster`, which has the
# same (C1 C2 - 1) free probabilities per class. Two continuous indicators
# instead get a covariance of their own within each class, every other
# covariance staying zero: the covariance Latent GOLD fits with the same term.

#' Read `local_dependence` into pairs of indicators
#'
#' @param local_dependence A character vector of pairs written `"a:b"`.
#' @return A data frame: `item_1`, `item_2` and `name` (`"a:b"`), one row
#'   per pair; `latents_bad_argument` for a malformed or overlapping pair.
#' @noRd
.latents_parse_pairs <- function(local_dependence) {
  refuse <- function(message) {
    stop(errorCondition(message, class = "latents_bad_argument", call = NULL))
  }
  if (!(is.character(local_dependence) && !anyNA(local_dependence))) {
    refuse("`local_dependence` names pairs of indicators, such as \"item1:item2\".")
  }
  parts <- strsplit(local_dependence, ":", fixed = TRUE)
  if (!all(lengths(parts) == 2L) || any(!nzchar(unlist(parts)))) {
    refuse(sprintf(paste(
      "Each pair in `local_dependence` is two indicator names joined by a colon,",
      "such as \"item1:item2\"; %s is not."),
      paste(sprintf("\"%s\"", local_dependence[lengths(parts) != 2L |
                                                 vapply(parts, \(p) any(!nzchar(p)), logical(1))]),
            collapse = ", ")))
  }
  pairs <- data.frame(item_1 = vapply(parts, `[[`, "", 1L),
                      item_2 = vapply(parts, `[[`, "", 2L),
                      stringsAsFactors = FALSE)
  pairs$name <- paste(pairs$item_1, pairs$item_2, sep = ":")
  if (any(pairs$item_1 == pairs$item_2)) refuse("A pair in `local_dependence` names one indicator twice.")
  items <- c(pairs$item_1, pairs$item_2)
  if (anyDuplicated(items) > 0L) {
    refuse(sprintf(paste(
      "An indicator can be in one pair of `local_dependence` only; %s is in more.",
      "A dependence among three items needs their joint table."),
      paste(sprintf("`%s`", unique(items[duplicated(items)])), collapse = ", ")))
  }
  pairs
}

#' Mark each pair as categorical or continuous, refusing a mixed pair
#' @param pairs From `.latents_parse_pairs()`.
#' @param categorical The categorical indicators.
#' @return `pairs` with a `kind` column.
#' @noRd
.latents_classify_pairs <- function(pairs, categorical) {
  first <- pairs$item_1 %in% categorical
  second <- pairs$item_2 %in% categorical
  if (any(first != second)) {
    stop(errorCondition(sprintf(paste(
      "A pair in `local_dependence` joins items of the same kind: two",
      "categorical items, modelled jointly, or two continuous ones, given a",
      "covariance. %s mixes the two."),
      paste(sprintf("\"%s\"", pairs$name[first != second]), collapse = ", ")),
      class = "latents_bad_argument", call = NULL))
  }
  pairs$kind <- ifelse(first, "categorical", "continuous")
  pairs
}

#' The free covariance entries of continuous local dependence
#'
#' @param pairs From `.latents_parse_pairs()`, continuous pairs.
#' @param vars,categorical,ordinal,count The fit's indicators by kind.
#' @param conflicts A named logical vector: settings that set the covariance
#'   another way, refused when any is `TRUE`.
#' @return A logical continuous-by-continuous matrix: `TRUE` on the diagonal
#'   and for each pair.
#' @noRd
.latents_covariance_pattern <- function(pairs, vars, categorical, ordinal, count,
                                        conflicts) {
  continuous <- setdiff(vars, c(categorical, ordinal, count))
  items <- c(pairs$item_1, pairs$item_2)
  absent <- setdiff(items, continuous)
  if (length(absent) > 0L) {
    stop(errorCondition(sprintf(paste(
      "`local_dependence` names %s, which %s not a continuous indicator among",
      "`vars`; pair two categorical items, or two continuous ones."),
      paste(sprintf("`%s`", absent), collapse = ", "),
      if (length(absent) > 1L) "are" else "is"), class = "latents_bad_argument", call = NULL))
  }
  if (any(conflicts)) {
    stop(errorCondition(sprintf(paste(
      "A covariance for chosen pairs of continuous indicators sets the",
      "covariance itself; it is not combined with %s."),
      paste(sprintf("`%s`", names(conflicts)[conflicts]), collapse = ", ")),
      class = "latents_bad_argument", call = NULL))
  }
  pattern <- diag(length(continuous)) > 0
  dimnames(pattern) <- list(continuous, continuous)
  pattern[cbind(pairs$item_1, pairs$item_2)] <- TRUE
  pattern[cbind(pairs$item_2, pairs$item_1)] <- TRUE
  pattern
}

#' The pairs a covariance pattern frees, as `"a:b"` names
#' @param pattern From `.latents_covariance_pattern()`.
#' @return A character vector, one name per pair.
#' @noRd
.latents_pattern_pairs <- function(pattern) {
  pairs <- which(pattern & upper.tri(pattern), arr.ind = TRUE)
  paste(rownames(pattern)[pairs[, "row"]], colnames(pattern)[pairs[, "col"]], sep = ":")
}

#' The covariance label of a fit with covariances for chosen pairs
#'
#' The diagonal structure the pattern extends, EEI with one shared covariance
#' matrix and VVI with one per profile, and the pairs: `"VVI + y1:y2"`. The
#' fit is estimated as a full covariance, so its structure field reads VVV or
#' EEE, which would describe a model with every covariance free.
#' @param fit A fit with `covariance_pattern`.
#' @return A single string.
#' @noRd
.latents_pattern_label <- function(fit) {
  sprintf("%s + %s", if (identical(fit$variance_model, "equal")) "EEI" else "VVI",
          paste(.latents_pattern_pairs(fit$covariance_pattern), collapse = ", "))
}

#' Hold a start's covariances to a pattern
#' @param parameters A start, with `covariances` under a full model.
#' @param pattern From `.latents_covariance_pattern()`, or `NULL`.
#' @return The start, its covariances zero outside the pattern.
#' @noRd
.latents_mask_covariances <- function(parameters, pattern) {
  if (is.null(pattern) || is.null(parameters$covariances)) return(parameters)
  masked <- parameters$covariances
  masked[!array(pattern, dim(masked))] <- 0
  parameters$covariances <- masked
  parameters
}

#' The number of covariances a pattern holds at zero, per covariance matrix
#' @param pattern From `.latents_covariance_pattern()`.
#' @return An integer.
#' @noRd
.latents_pattern_zeros <- function(pattern) {
  sum(lower.tri(pattern) & !pattern)
}

#' Replace the paired indicators by their joint indicator
#'
#' The joint indicator takes the place of the first item of its pair, so the
#' indicators keep their order.
#' @param vars Indicator names.
#' @param pairs From `.latents_parse_pairs()`.
#' @return The indicator names after the join.
#' @noRd
.latents_joined_names <- function(vars, pairs) {
  vars[match(pairs$item_1, vars)] <- pairs$name
  setdiff(vars, pairs$item_2)
}

#' Add the joint indicator of each pair to a data frame
#'
#' Each joint value is the two items' values joined by a colon, missing when
#' either item is. A frame that already holds a joint indicator (the frame a
#' fit stores, or one this function has joined) is left as it is, so joining
#' is idempotent. A row where one item of a pair is observed and the other is
#' not keeps the observed item, labelled by `.latents_partial_label()`.
#' @param data A data frame.
#' @param pairs From `.latents_parse_pairs()`.
#' @param min_categories Minimum observed categories of each original item;
#'   fitting requires two, prediction subsets may have none.
#' @param rebuild_existing Regenerate an existing joint column from its source
#'   items when a fitted model has already established their meaning. Initial
#'   fitting leaves this false so existing caller columns remain protected.
#' @return The data frame with one more column per pair.
#' @noRd
.latents_join_pairs <- function(data, pairs, min_categories = 0L,
                                 rebuild_existing = FALSE) {
  if (!is.null(pairs$kind)) pairs <- pairs[pairs$kind == "categorical", , drop = FALSE]
  joined <- Reduce(function(frame, index) {
    pair <- pairs[index, ]
    has_items <- all(c(pair$item_1, pair$item_2) %in% names(frame))
    prepared <- pair$name %in% names(frame) &&
      identical(attr(frame[[pair$name]], "latents_pair"), c(pair$item_1, pair$item_2))
    if (pair$name %in% names(frame) && !has_items &&
        (!any(c(pair$item_1, pair$item_2) %in% names(frame)) || prepared)) {
      # A fitted pair may be supplied as its reconstructed joint column.
      # Restore provenance before interpreting its partial-row markers.
      attr(frame[[pair$name]], "latents_pair") <- c(pair$item_1, pair$item_2)
      return(frame)
    }
    if (!has_items) {
      stop(errorCondition(sprintf(
        "`data` needs `%s` and `%s` for the pair \"%s\" of `local_dependence`.",
        pair$item_1, pair$item_2, pair$name), class = "latents_bad_data", call = NULL))
    }
    first <- frame[[pair$item_1]]
    second <- frame[[pair$item_2]]
    first_levels <- .multilpa_categorical_levels(first, pair$item_1, min_categories)
    second_levels <- .multilpa_categorical_levels(second, pair$item_2, min_categories)
    # A joint category is written "first:second", so a category that already
    # holds a colon makes two responses one label: ("x:y", "z") and
    # ("x", "y:z") would be counted as the same response.
    colons <- Filter(function(item) {
      any(grepl(":", unique(as.character(frame[[item]])), fixed = TRUE))
    }, c(pair$item_1, pair$item_2))
    if (length(colons) > 0L) {
      stop(errorCondition(sprintf(paste(
        "The categories of %s contain \":\", which joins the two items of the",
        "pair \"%s\"; recode them without a colon."),
        paste(sprintf("`%s`", colons), collapse = " and "), pair$name),
        class = "latents_bad_data", call = NULL))
    }
    observed <- !is.na(first) & !is.na(second)
    labels <- ifelse(observed, paste(as.character(first), as.character(second), sep = ":"),
                     NA_character_)
    # Base rbind() drops custom factor attributes but keeps the first frame's
    # metadata. Trust that fallback only when the column still contains the
    # joint responses of its source items; stale metadata cannot authorize
    # replacing a newly assigned, unrelated column.
    frame_pairs <- attr(frame, "latents_pairs")
    prepared <- prepared || (pair$name %in% names(frame) &&
      identical(frame_pairs[[pair$name]], c(pair$item_1, pair$item_2)) &&
      is.factor(frame[[pair$name]]) &&
      identical(as.character(frame[[pair$name]]), labels))
    # Only columns produced by this preparation may be regenerated beside
    # their source items. A caller column with the same name is never replaced.
    if (pair$name %in% names(frame) && !prepared && !rebuild_existing) {
      stop(errorCondition(sprintf(paste(
        "The generated joint indicator `%s` conflicts with an existing column;",
        "rename that column before using `local_dependence`."), pair$name),
        class = "latents_bad_data", call = NULL))
    }
    # A row observing one item of the pair keeps that item: it is labelled
    # "a::" (first item only) or ":::b" (second item only), and the
    # categorical engine sums the missing item out of the pair's joint
    # distribution, which is the observed-data likelihood under MAR (review
    # 2026-10-26 round 2, C2-2). Two or more colons never occur in a joint
    # category, because the items' own categories hold none.
    first_only <- !is.na(first) & is.na(second)
    second_only <- is.na(first) & !is.na(second)
    labels[first_only] <- .latents_partial_label(first[first_only], "first")
    labels[second_only] <- .latents_partial_label(second[second_only], "second")
    # Categories ordered by the first item's categories, then the second's;
    # the labels of half-observed rows follow, first item's then second's.
    grid <- expand.grid(second = as.character(second_levels),
                        first = as.character(first_levels), stringsAsFactors = FALSE)
    combinations <- paste(grid$first, grid$second, sep = ":")
    halves <- c(.latents_partial_label(first_levels, "first"),
                .latents_partial_label(second_levels, "second"))
    frame[[pair$name]] <- factor(labels, levels = c(combinations, halves))
    attr(frame[[pair$name]], "latents_pair") <- c(pair$item_1, pair$item_2)
    frame_pairs[[pair$name]] <- c(pair$item_1, pair$item_2)
    attr(frame, "latents_pairs") <- frame_pairs
    frame
  }, seq_len(nrow(pairs)), init = data)
  joined
}

#' The label of a half-observed pair row
#'
#' `"a::"` when only the first item is observed (its category `a`), `":::b"`
#' when only the second is (category `b`). A joint category holds exactly one
#' colon, so these labels cannot be mistaken for one, nor for each other even
#' when a category is the empty string.
#' @param value Categories of the observed item.
#' @param item `"first"` or `"second"`.
#' @return A character vector.
#' @noRd
.latents_partial_label <- function(value, item = c("first", "second")) {
  item <- match.arg(item)
  value <- as.character(value)
  if (identical(item, "first")) paste0(value, "::") else paste0(":::", value)
}

#' Which category labels of a joint item mark half-observed rows
#'
#' A label is a half-observed marker only when the labels around it have the
#' joint item's shape: every other label holds exactly one colon (two items'
#' categories), and the marker's category is one of the matching item's
#' categories in those labels. A column that is not a joined pair never has
#' that shape, so its labels are all ordinary categories.
#' @param labels Category labels (factor levels).
#' @return A logical vector, `TRUE` for a half-observed marker.
#' @noRd
.latents_partial_levels <- function(labels) {
  labels <- as.character(labels)
  distinct <- unique(labels[!is.na(labels)])
  .latents_partial_distinct(distinct)[match(labels, distinct)] %in% TRUE
}

#' `.latents_partial_levels()` on distinct, non-missing labels
#' @noRd
.latents_partial_distinct <- function(labels) {
  marked <- grepl("::", labels, fixed = TRUE)
  if (!any(marked)) return(marked)
  joint <- labels[!marked]
  colons <- lengths(regmatches(joint, gregexpr(":", joint, fixed = TRUE)))
  if (length(joint) == 0L || any(colons != 1L)) return(rep(FALSE, length(labels)))
  first_halves <- sub(":.*$", "", joint)
  second_halves <- sub("^[^:]*:", "", joint)
  second_only <- marked & startsWith(labels, ":::")
  category <- ifelse(second_only, substring(labels, 4L), sub("::$", "", labels))
  valid <- !grepl(":", category, fixed = TRUE) &
    ifelse(second_only, category %in% second_halves, category %in% first_halves)
  if (!all(valid[marked])) return(rep(FALSE, length(labels)))
  marked
}

#' Half-observed rows of a joint item, as sets of joint categories
#'
#' @param value The joint column (labels, `NA` where the pair is unobserved).
#' @param grid The joint categories, in code order.
#' @return `NULL` when no row is half observed; otherwise a list with `index`
#'   (one integer per row: the row's set, `NA` for every other row), `labels`
#'   (the marker of each set, `"a::"` or `":::b"`) and `sets` (for each set
#'   the codes of the joint categories that agree with the observed item).
#' @noRd
.latents_partial_map <- function(value, grid) {
  value <- as.character(value)
  levels_seen <- unique(value[!is.na(value)])
  marked <- .latents_partial_levels(c(grid, levels_seen))[-seq_along(grid)]
  labels <- levels_seen[marked]
  if (length(labels) == 0L) return(NULL)
  # Deterministic set order: first-item markers, then second-item ones, each
  # in the order of the grid's own categories.
  first_halves <- sub(":.*$", "", grid)
  second_halves <- sub("^[^:]*:", "", grid)
  second_only <- startsWith(labels, ":::")
  category <- ifelse(second_only, substring(labels, 4L), sub("::$", "", labels))
  rank <- ifelse(second_only, length(grid) + match(category, second_halves),
                 match(category, first_halves))
  ordering <- order(rank, labels)
  labels <- labels[ordering]
  second_only <- second_only[ordering]
  category <- category[ordering]
  sets <- lapply(seq_along(labels), function(set) {
    which(if (second_only[[set]]) second_halves == category[[set]] else
      first_halves == category[[set]])
  })
  list(index = match(value, labels), labels = labels, sets = sets)
}

#' Put the markers of half-observed pair rows back into a rebuilt column
#'
#' A fit's code matrix holds `NA` for a row that observes one item of a
#' joined pair, with the row's marker in its `"latents_partial"` attribute.
#' Data rebuilt from the codes write the marker back, so a refit sees the
#' observed item as the fit did.
#' @param value The column rebuilt from the codes.
#' @param codes The fit's code matrix.
#' @param name The column's name.
#' @return `value`, with the markers in the half-observed rows.
#' @noRd
.latents_restore_partial <- function(value, codes, name) {
  map <- attr(codes, "latents_partial", exact = TRUE)[[name]]
  if (is.null(map)) return(value)
  half <- !is.na(map$index)
  if (!is.factor(value)) {
    value <- as.character(value)
    value[half] <- map$labels[map$index[half]]
    return(value)
  }
  restored <- factor(as.character(value),
                     levels = union(levels(value), map$labels))
  restored[half] <- map$labels[map$index[half]]
  restored
}

#' Fit a model whose `local_dependence` pairs are joined first
#'
#' Called by [multilpa()] when `local_dependence` is given: the pairs are
#' joined in `data`, `vars` and `categorical`, and the model is fitted on the
#' joint indicators, with the caller's own arguments otherwise unchanged.
#' @param call The caller's `match.call()`.
#' @param env The calling function's environment, holding the matched values.
#' @param local_dependence The caller's pairs.
#' @return The fit, carrying `local_dependence` (the pairs) and the caller's
#'   `call`.
#' @noRd
.latents_fit_local_dependence <- function(call, env, local_dependence) {
  pairs <- .latents_parse_pairs(local_dependence)
  arguments <- mget(setdiff(names(call)[-1L], "local_dependence"), envir = env)
  if (!is.null(arguments$family) && !identical(arguments$family, "profiles")) {
    stop(errorCondition(paste(
      "`local_dependence` pairs the indicators of a profile or class model;",
      "the group-class families have no categorical items to pair."),
      class = "latents_bad_argument", call = NULL))
  }
  items <- c(pairs$item_1, pairs$item_2)
  categorical <- arguments$categorical %||% character()
  # Joining removes original item names and deduplicates the new name list;
  # preserve the ordinary identity checks before either transformation.
  if (anyNA(arguments$vars) || anyDuplicated(arguments$vars) ||
      !all(arguments$vars %in% names(arguments$data)) ||
      any(arguments$id %in% arguments$vars)) {
    stop(errorCondition(
      "Supply unique existing indicators and a distinct group column.",
      class = "latents_bad_data", call = NULL))
  }
  if (anyDuplicated(categorical) || !all(categorical %in% arguments$vars)) {
    stop(errorCondition("`categorical` must name distinct indicators listed in `vars`.",
                        class = "latents_bad_categorical", call = NULL))
  }
  absent <- setdiff(items, arguments$vars)
  if (length(absent) > 0L) {
    stop(errorCondition(sprintf(
      "`local_dependence` names %s, which %s not among `vars`.",
      paste(sprintf("`%s`", absent), collapse = ", "),
      if (length(absent) > 1L) "are" else "is"), class = "latents_bad_argument", call = NULL))
  }
  pairs <- .latents_classify_pairs(pairs, categorical)
  # Once joined, a paired item is no longer among `vars`, so the checks that
  # keep an indicator from also being a covariate cannot see it; they are made
  # here, with the classes those checks raise.
  covariates <- c(arguments$profile_covariates, arguments$group_covariates,
                  if (length(arguments$direct_effects) > 0L)
                    .latents_parse_direct_effects(arguments$direct_effects)$covariate)
  reused <- intersect(covariates, items)
  if (length(reused) > 0L) {
    stop(errorCondition(sprintf(
      "%s %s an indicator of a `local_dependence` pair, so it cannot also be a covariate.",
      paste(sprintf("`%s`", reused), collapse = ", "),
      if (length(reused) > 1L) "are" else "is"),
      class = c("multilpa_indicator_reused", "latents_bad_covariate"), call = NULL))
  }
  joined_pairs <- pairs[pairs$kind == "categorical", , drop = FALSE]
  continuous_pairs <- pairs[pairs$kind == "continuous", , drop = FALSE]
  arguments$local_dependence <- if (nrow(continuous_pairs) > 0L) continuous_pairs$name
  # Apply the same native-value and variation requirements as unpaired items.
  # The join itself permits constant/single-row subsets used by prediction.
  arguments$data <- .latents_join_pairs(arguments$data, joined_pairs, min_categories = 2L)
  # Each joint category's two levels, read from the data rather than parsed
  # back out of the label.
  categories <- lapply(stats::setNames(seq_len(nrow(joined_pairs)), joined_pairs$name),
                       function(index) {
    first <- as.character(arguments$data[[joined_pairs$item_1[index]]])
    second <- as.character(arguments$data[[joined_pairs$item_2[index]]])
    label <- as.character(arguments$data[[joined_pairs$name[index]]])
    unique(data.frame(label = label, level_1 = first, level_2 = second,
                      stringsAsFactors = FALSE)[!is.na(label), , drop = FALSE])
  })
  original <- list(vars = arguments$vars, categorical = categorical)
  arguments$vars <- .latents_joined_names(arguments$vars, joined_pairs)
  arguments$categorical <- .latents_joined_names(categorical, joined_pairs)
  fit <- do.call(multilpa, arguments)
  fit$local_dependence <- pairs
  fit$local_dependence_categories <- categories
  # The caller's own indicators, so the fit compares with one of the same data
  # whose items were not paired.
  fit$original_vars <- original$vars
  fit$original_categorical <- original$categorical
  fit$call <- call
  fit
}

#' A fit's categorical responses as labels, one column per original item
#'
#' A joint indicator is split back into its two items through the categories
#' the fit stored, so a model with paired items and one without describe the
#' same data alike.
#' @param fit A fitted model.
#' @return A character matrix, one row per observation and one column per
#'   categorical item of the caller, or `NULL` without categorical items.
#' @noRd
.latents_categorical_labels <- function(fit) {
  # A general lta() fit keeps its categorical codes as `codes`.
  codes <- fit$categorical_data %||% fit$codes
  if (is.null(codes) || ncol(codes) == 0L) return(NULL)
  labels <- lapply(colnames(codes), function(name) {
    # A half-observed pair row comes back as its marker, which the stored
    # categories split into the observed item and a missing one.
    .latents_restore_partial(fit$categorical_levels[[name]][codes[, name]], codes, name)
  })
  names(labels) <- colnames(codes)
  pairs <- fit$local_dependence
  if (!is.null(pairs)) {
    labels <- unlist(lapply(names(labels), function(name) {
      if (!name %in% pairs$name) return(stats::setNames(list(labels[[name]]), name))
      categories <- fit$local_dependence_categories[[name]]
      at <- match(labels[[name]], categories$label)
      pair <- pairs[pairs$name == name, ]
      stats::setNames(list(categories$level_1[at], categories$level_2[at]),
                      c(pair$item_1, pair$item_2))
    }), recursive = FALSE)
    labels <- labels[fit$original_categorical]
  }
  do.call(cbind, labels)
}

#' The joint and independent probabilities of each locally dependent pair
#'
#' For every pair and profile, every combination of the two items'
#' categories: its probability under the fitted joint distribution, the
#' product of the two marginal probabilities it would have if the items were
#' independent within the profile, and their ratio. A ratio far from one in a
#' profile is the dependence the joint indicator captures.
#' @param x A fitted model with `local_dependence`.
#' @return A data frame, one row per pair, profile and combination: `pair`,
#'   `profile`, `item_1`, `category_1`, `item_2`, `category_2`,
#'   `probability`, `independent` and `ratio`.
#' @noRd
.latents_local_dependence_table <- function(x) {
  pairs <- x$local_dependence
  if (is.null(pairs) || nrow(pairs) == 0L) {
    stop(errorCondition(paste(
      "This fit has no locally dependent pairs; fit with",
      "`local_dependence = \"item1:item2\"` to model one."),
      class = "latents_no_local_dependence", call = NULL))
  }
  kinds <- pairs$kind %||% rep("categorical", nrow(pairs))
  rows <- lapply(seq_len(nrow(pairs)), function(index) {
    pair <- pairs[index, ]
    if (identical(kinds[index], "continuous")) {
      # The pair's covariance in each profile, and the correlation it implies.
      continuous <- .multilpa_continuous_names(x)
      first <- match(pair$item_1, continuous)
      second <- match(pair$item_2, continuous)
      return(do.call(rbind, lapply(seq_len(x$n_profiles), function(profile) {
        covariance <- matrix(x$covariances[, , profile], length(continuous))
        data.frame(pair = pair$name, kind = "continuous", profile = profile,
                   item_1 = pair$item_1, category_1 = NA_character_,
                   item_2 = pair$item_2, category_2 = NA_character_,
                   probability = NA_real_, independent = NA_real_, ratio = NA_real_,
                   covariance = covariance[first, second],
                   correlation = covariance[first, second] /
                     sqrt(covariance[first, first] * covariance[second, second]),
                   stringsAsFactors = FALSE)
      })))
    }
    probabilities <- x$response_probabilities[[pair$name]]
    labels <- x$categorical_levels[[pair$name]]
    categories <- x$local_dependence_categories[[pair$name]]
    at <- match(labels, categories$label)
    first_items <- categories$level_1[at]
    second_items <- categories$level_2[at]
    # A combination nobody gave is still a category of the joint item, at the
    # probability floor, but it is absent from the categories read from the
    # data. Its label is "first:second", and a category holding ":" is refused
    # when the pair is joined, so the label's two halves are its two levels.
    unseen <- is.na(at)
    if (any(unseen)) {
      first_items[unseen] <- sub(":.*$", "", labels[unseen])
      second_items[unseen] <- sub("^[^:]*:", "", labels[unseen])
    }
    do.call(rbind, lapply(seq_len(nrow(probabilities)), function(profile) {
      joint <- probabilities[profile, ]
      marginal_1 <- tapply(joint, first_items, sum)[first_items]
      marginal_2 <- tapply(joint, second_items, sum)[second_items]
      independent <- unname(marginal_1 * marginal_2)
      data.frame(pair = pair$name, kind = "categorical", profile = profile,
                 item_1 = pair$item_1, category_1 = first_items,
                 item_2 = pair$item_2, category_2 = second_items,
                 probability = unname(joint), independent = independent,
                 ratio = unname(joint) / independent, covariance = NA_real_,
                 correlation = NA_real_, stringsAsFactors = FALSE)
    }))
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}
