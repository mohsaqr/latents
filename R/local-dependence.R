# Local dependence between two categorical indicators: the pair is modelled
# jointly within each class. The two items become one categorical indicator
# whose categories are their observed combinations, so each class has its own
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
#' is idempotent. Rows where one item of a pair is observed and the other is
#' not lose the pair, which is reported with `latents_partial_pair`.
#' @param data A data frame.
#' @param pairs From `.latents_parse_pairs()`.
#' @return The data frame with one more column per pair.
#' @noRd
.latents_join_pairs <- function(data, pairs) {
  if (!is.null(pairs$kind)) pairs <- pairs[pairs$kind == "categorical", , drop = FALSE]
  joined <- Reduce(function(frame, index) {
    pair <- pairs[index, ]
    has_items <- all(c(pair$item_1, pair$item_2) %in% names(frame))
    if (pair$name %in% names(frame) && !has_items) return(frame)
    if (!has_items) {
      stop(errorCondition(sprintf(
        "`data` needs `%s` and `%s` for the pair \"%s\" of `local_dependence`.",
        pair$item_1, pair$item_2, pair$name), class = "latents_bad_data", call = NULL))
    }
    first <- frame[[pair$item_1]]
    second <- frame[[pair$item_2]]
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
    partial <- xor(is.na(first), is.na(second))
    if (any(partial)) {
      warning(warningCondition(sprintf(paste(
        "%d rows observe only one item of the pair \"%s\"; the pair is treated",
        "as missing in those rows."), sum(partial), pair$name),
        class = "latents_partial_pair", call = NULL))
    }
    labels <- ifelse(observed, paste(as.character(first), as.character(second), sep = ":"),
                     NA_character_)
    # Categories ordered by the first item's categories, then the second's.
    order_of <- function(value) {
      if (is.factor(value)) levels(droplevels(value)) else sort(unique(value[!is.na(value)]))
    }
    grid <- expand.grid(second = as.character(order_of(second)),
                        first = as.character(order_of(first)), stringsAsFactors = FALSE)
    combinations <- paste(grid$first, grid$second, sep = ":")
    frame[[pair$name]] <- factor(labels, levels = intersect(combinations, labels))
    frame
  }, seq_len(nrow(pairs)), init = data)
  joined
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
  arguments$data <- .latents_join_pairs(arguments$data, joined_pairs)
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
    fit$categorical_levels[[name]][codes[, name]]
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
