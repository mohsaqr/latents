# Local dependence between two categorical indicators: the pair is modelled
# jointly within each class. The two items become one categorical indicator
# whose categories are their observed combinations, so each class has its own
# joint distribution of the pair, with no independence assumed between them.
# For two nominal items this is the class-specific two-variable log-linear
# association Latent GOLD fits with `item1 <-> item2 | Cluster`, which has the
# same (C1 C2 - 1) free probabilities per class.

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
  continuous <- setdiff(items, categorical)
  if (length(continuous) > 0L) {
    stop(errorCondition(sprintf(paste(
      "`local_dependence` pairs categorical indicators; %s %s not declared in",
      "`categorical`."), paste(sprintf("`%s`", continuous), collapse = ", "),
      if (length(continuous) > 1L) "are" else "is"),
      class = "latents_bad_argument", call = NULL))
  }
  arguments$data <- .latents_join_pairs(arguments$data, pairs)
  # Each joint category's two levels, read from the data rather than from the
  # label, so a level that itself holds a colon is not split.
  categories <- lapply(stats::setNames(seq_len(nrow(pairs)), pairs$name), function(index) {
    first <- as.character(arguments$data[[pairs$item_1[index]]])
    second <- as.character(arguments$data[[pairs$item_2[index]]])
    label <- as.character(arguments$data[[pairs$name[index]]])
    unique(data.frame(label = label, level_1 = first, level_2 = second,
                      stringsAsFactors = FALSE)[!is.na(label), , drop = FALSE])
  })
  original <- list(vars = arguments$vars, categorical = categorical)
  arguments$vars <- .latents_joined_names(arguments$vars, pairs)
  arguments$categorical <- .latents_joined_names(categorical, pairs)
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
  codes <- fit$categorical_data
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
  rows <- lapply(seq_len(nrow(pairs)), function(index) {
    pair <- pairs[index, ]
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
      data.frame(pair = pair$name, profile = profile, item_1 = pair$item_1,
                 category_1 = first_items, item_2 = pair$item_2,
                 category_2 = second_items, probability = unname(joint),
                 independent = independent, ratio = unname(joint) / independent,
                 stringsAsFactors = FALSE)
    }))
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}
