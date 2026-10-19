# How classes split across the solutions of an enumeration grid.

#
# A tree is a sequence of solutions (rows). Each row divides the units into its
# classes by posterior share, and the band between two rows carries the
# posterior mass the two classes share: the sum over units of p(a) * p(b).
# At the individual level the rows are the numbers of classes at one number of
# group classes; at the group level they are the numbers of group classes at
# one number of classes; or the rows are candidates named by the caller. Every
# tree of one call classifies the same units, so the classes of one tree can
# be matched to those of another and drawn in the same colour.

utils::globalVariables(c("row_facet", "col_facet", "row", "xmin", "xmax",
                         "fill", "label", "text_colour", "group", "x", "y"))

#' The trees a call asks for
#'
#' @param x A `multilpa_enumeration`.
#' @param level `"individuals"` or `"groups"`.
#' @param n_profiles,n_group_classes,model Optional filters: the one that
#'   indexes trees selects trees, the other selects rows.
#' @param candidates `NULL`, or two or more candidate names drawn as the rows
#'   of one tree, in that order.
#' @return A list of trees, each a list of `name`, `label`, `row_facet`,
#'   `col_facet` and `rows` (a data frame: `at`, the grid row; `candidate`;
#'   `key`, the count matched across trees; `label`, the row's axis label).
#' @noRd
.latents_tree_specs <- function(x, level, n_profiles, n_group_classes, model,
                                candidates) {
  table <- .subset2(x, "table")
  names <- .latents_candidate_names(x)
  if (!is.null(candidates)) {
    if (!is.null(n_profiles) || !is.null(n_group_classes) || !is.null(model)) {
      stop(errorCondition(paste(
        "`candidates` names the rows of the tree itself; it is not combined",
        "with `n_profiles`, `n_group_classes` or `model`."),
        class = "latents_bad_argument", call = NULL))
    }
    if (!(is.character(candidates) && length(candidates) >= 2L &&
          !anyNA(candidates) && !anyDuplicated(candidates))) {
      stop(errorCondition(paste(
        "`candidates` names two or more different candidates, such as",
        "c(\"classes3_groups3\", \"classes4_groups3\"); `from` and `to` take",
        "class counts, such as from = c(3, 3), to = c(4, 3)."),
        class = "latents_bad_argument", call = NULL))
    }
    at <- vapply(candidates, \(name) as.integer(.latents_candidate_row(x, name)),
                 integer(1), USE.NAMES = FALSE)
    invisible(lapply(at, \(row) .latents_candidate_at(x, row)))
    unconverged <- !table$converged[at] %in% TRUE
    if (any(unconverged)) {
      stop(errorCondition(sprintf(
        "%s did not converge, so its classification is not a solution to compare.",
        paste(sprintf("\"%s\"", names[at][unconverged]), collapse = ", ")),
        class = "latents_unconverged", call = NULL))
    }
    if (identical(level, "groups") && all(table$n_group_classes[at] %in% 1L)) {
      stop(errorCondition(paste(
        "These candidates have one group class, so their groups have no tree."),
        class = c("latents_no_tree", "latents_nothing_to_plot"), call = NULL))
    }
    # Each row is labelled by its counts (and structure, on a grid crossed
    # with structures), as the axis of every other tree.
    word <- .latents_class_word(x)
    counts_label <- \(count, noun) {
      sprintf("%d %s", count, if (count == 1L) .latents_singular(noun) else noun)
    }
    crossed <- length(unique(table$model[!is.na(table$model)])) > 1L
    labels <- vapply(at, \(row) paste(c(
      counts_label(table$n_profiles[row], word),
      if (!all(table$n_group_classes %in% 1L)) {
        counts_label(table$n_group_classes[row], "group classes")
      },
      if (crossed) table$model[row]), collapse = ", "), character(1))
    return(list(list(
      name = "candidates", label = "", row_facet = "", col_facet = "",
      rows = data.frame(at = at, candidate = names[at], key = NA_integer_,
                        label = labels, stringsAsFactors = FALSE))))
  }
  structure <- ifelse(is.na(table$model), "", table$model)
  fitted <- !vapply(.subset2(x, "fits"), is.null, logical(1))
  usable <- table$converged %in% TRUE & fitted
  .latents_check_tree_filter(model, structure[usable], "model")
  .latents_check_tree_filter(n_profiles, table$n_profiles[usable], "n_profiles")
  .latents_check_tree_filter(n_group_classes, table$n_group_classes[usable],
                             "n_group_classes")
  keep <- usable &
    (if (is.null(model)) TRUE else structure %in% model) &
    (if (is.null(n_profiles)) TRUE else table$n_profiles %in% n_profiles) &
    (if (is.null(n_group_classes)) TRUE else
      table$n_group_classes %in% n_group_classes)
  groups_level <- identical(level, "groups")
  # The count that indexes trees, and the count that runs down each tree.
  tree_count <- if (groups_level) table$n_profiles else table$n_group_classes
  row_count <- if (groups_level) table$n_group_classes else table$n_profiles
  # One count asked for down the rows is one level of every tree: the rows
  # then run across the trees, comparing that level's solutions.
  fixed <- unique(if (groups_level) n_group_classes else n_profiles)
  across <- length(fixed) == 1L
  if (across) {
    row_count <- tree_count
    tree_count <- rep(0L, nrow(table))
  }
  single_level <- all(table$n_group_classes %in% 1L)
  if (groups_level && single_level) {
    stop(errorCondition(paste(
      "This grid has one group class, so its groups have no tree.",
      "A group-level tree needs two or more numbers of group classes."),
      class = c("latents_no_tree", "latents_nothing_to_plot"), call = NULL))
  }
  word <- .latents_class_word(x)
  # What the rows count, for the axis and the messages, and the level held
  # fixed when the rows run across the trees.
  rows_noun <- if (groups_level != across) "group classes" else word
  fixed_noun <- if (groups_level) "group classes" else word
  fixed_label <- if (across) sprintf("%d %s", fixed, if (fixed == 1L)
    .latents_singular(fixed_noun) else fixed_noun)
  keys <- unique(data.frame(structure = structure[keep],
                            count = tree_count[keep], stringsAsFactors = FALSE))
  keys <- keys[order(match(keys$structure, unique(structure)), keys$count), ,
               drop = FALSE]
  # A structure is named only when the trees span more than one.
  several_structures <- length(unique(keys$structure)) > 1L
  trees <- lapply(seq_len(nrow(keys)), \(index) {
    rows_in <- which(keep & structure == keys$structure[index] &
                       tree_count == keys$count[index])
    rows_in <- rows_in[order(row_count[rows_in])]
    count <- keys$count[index]
    count_label <- if (across) fixed_label else if (groups_level) {
      sprintf("%d %s", count, if (count == 1L) .latents_singular(word) else word)
    } else if (single_level) "" else {
      sprintf("%d group class%s", count, if (count == 1L) "" else "es")
    }
    tree_structure <- if (several_structures) keys$structure[index] else ""
    list(
      name = paste(c(if (nzchar(tree_structure)) tree_structure,
                     if (across) {
                       sprintf("%s%d", if (groups_level) "groups" else word, fixed)
                     } else if (groups_level) sprintf("%s%d", word, count) else
                       if (!single_level) sprintf("groups%d", count)),
                   collapse = "_"),
      label = paste(c(if (nzchar(tree_structure)) tree_structure,
                      if (!across && nzchar(count_label)) count_label),
                    collapse = ", "),
      row_facet = if (across) "" else count_label, col_facet = tree_structure,
      rows = data.frame(at = rows_in, candidate = names[rows_in],
                        key = row_count[rows_in],
                        label = as.character(row_count[rows_in]),
                        stringsAsFactors = FALSE))
  })
  trees <- Filter(\(tree) nrow(tree$rows) >= 2L, trees)
  if (length(trees) == 0L) {
    stop(errorCondition(sprintf(
      "A tree needs at least two numbers of %s for one model.", rows_noun),
      class = c("latents_no_tree", "latents_nothing_to_plot"), call = NULL))
  }
  # A tree with no structure and no count label is the grid's only tree.
  if (length(trees) == 1L && !nzchar(trees[[1L]]$name)) trees[[1L]]$name <- "tree"
  attr(trees, "rows_noun") <- rows_noun
  attr(trees, "fixed") <- if (across) fixed_label
  attr(trees, "fixed_solution") <- if (across) {
    sprintf("%d-%s", fixed, gsub(" ", "-", .latents_singular(fixed_noun)))
  }
  trees
}

#' The candidates `from` and `to` name, or the caller's `candidates`
#'
#' `from` and `to` each name one candidate as a grid names it everywhere: by
#' its name, by its class counts `c(n_profiles, n_group_classes)`, or by one
#' count for one group class.
#' @param x The grid.
#' @param candidates The caller's `candidates`.
#' @param from,to The caller's `from` and `to`.
#' @return A character vector of candidate names, or `NULL`.
#' @noRd
.latents_tree_candidates <- function(x, candidates, from, to, model = NULL) {
  if (is.null(from) && is.null(to)) return(candidates)
  if (is.null(from) || is.null(to) || !is.null(candidates)) {
    stop(errorCondition(paste(
      "`from` and `to` name the two solutions to compare, both of them and",
      "without `candidates`."), class = "latents_bad_argument", call = NULL))
  }
  names <- .latents_candidate_names(x)
  table <- .subset2(x, "table")
  resolve <- \(value) {
    # On a grid crossed with structures, `model` picks the structure of a
    # candidate given by its class counts.
    if (is.null(model) || is.character(value)) {
      return(names[.latents_candidate_row(x, value)])
    }
    counts <- c(value, 1)[seq_len(2L)]
    at <- which(table$n_profiles == counts[1L] & table$n_group_classes == counts[2L] &
                  table$model %in% model)
    if (length(at) != 1L) {
      stop(errorCondition(sprintf(
        "No single candidate has %s classes, %s group classes and `model` %s.",
        counts[1L], counts[2L], paste(model, collapse = ", ")),
        class = "latents_unknown_candidate", call = NULL))
    }
    names[at]
  }
  c(resolve(from), resolve(to))
}

#' Check the arguments of a tree before anything is drawn
#'
#' @param level The caller's `level`.
#' @param n_profiles,n_group_classes,model The caller's filters.
#' @param combine The caller's `combine`.
#' @return The matched `level`.
#' @noRd
.latents_check_tree_arguments <- function(level, n_profiles, n_group_classes,
                                          model, combine = FALSE) {
  refuse <- \(message) {
    stop(errorCondition(message, class = "latents_bad_argument", call = NULL))
  }
  levels <- c("individuals", "groups")
  if (identical(level, levels)) level <- levels[1L]
  if (!(is.character(level) && length(level) == 1L && level %in% levels)) {
    refuse("`level` must be \"individuals\" or \"groups\".")
  }
  counts_ok <- \(value) is.null(value) || (is.numeric(value) && length(value) >= 1L &&
    all(is.finite(value)) && all(value >= 1) && all(value == floor(value)))
  if (!counts_ok(n_profiles)) refuse("`n_profiles` must be whole numbers of 1 or more.")
  if (!counts_ok(n_group_classes)) {
    refuse("`n_group_classes` must be whole numbers of 1 or more.")
  }
  if (!(is.null(model) || (is.character(model) && length(model) >= 1L && !anyNA(model)))) {
    refuse("`model` must name one or more covariance structures.")
  }
  if (!(isTRUE(combine) || isFALSE(combine))) refuse("`combine` must be TRUE or FALSE.")
  level
}

#' The singular of a class word
#' @param word `"profiles"`, `"classes"` or `"group classes"`.
#' @return The singular.
#' @noRd
.latents_singular <- function(word) {
  c(profiles = "profile", classes = "class", "group classes" = "group class")[[word]]
}

#' Refuse a tree filter that names nothing in the grid
#' @param value The caller's filter, or `NULL`.
#' @param available The values present among converged candidates.
#' @param argument The argument's name.
#' @return `NULL`, invisibly.
#' @noRd
.latents_check_tree_filter <- function(value, available, argument) {
  if (is.null(value)) return(invisible(NULL))
  available <- sort(unique(available))
  missing_values <- setdiff(value, available)
  if (length(missing_values) > 0L) {
    stop(errorCondition(sprintf(
      "`%s` names %s, which no converged candidate has. The grid has %s.",
      argument, paste(missing_values, collapse = ", "),
      paste(available, collapse = ", ")),
      class = "latents_unknown_candidate", call = NULL))
  }
  invisible(NULL)
}

#' The posteriors a tree divides, with the weight of each unit
#' @param fit A fitted candidate.
#' @param level `"individuals"` or `"groups"`.
#' @return A list of `posteriors` (units x classes) and `weights`.
#' @noRd
.latents_tree_posteriors <- function(fit, level) {
  sampling <- fit$sampling_weights
  if (identical(level, "groups")) {
    # The groups the fit counts, weighted by their sampling weights, as in
    # the fit's effective group counts.
    posteriors <- fit$group_posteriors
    counted <- .multilpa_informative_groups(fit) %||% rep(TRUE, nrow(posteriors))
    weights <- if (is.null(sampling)) rep(1, nrow(posteriors)) else
      if (is.null(names(sampling)) || is.null(rownames(posteriors))) unname(sampling) else
        unname(sampling)[match(rownames(posteriors), names(sampling))]
    return(list(posteriors = posteriors[counted, , drop = FALSE],
                weights = weights[counted]))
  }
  if (isTRUE(fit$noise)) {
    stop(errorCondition(paste(
      "A fit with a noise component has no tree of its classes yet: the noise",
      "component's share would be left out of every row."),
      class = "latents_unsupported_noise", call = NULL))
  }
  posteriors <- fit$subject_posteriors
  counted <- .multilpa_informative_rows(fit)
  weights <- if (is.null(sampling)) rep(1, nrow(posteriors)) else
    unname(sampling)[fit$group_index]
  list(posteriors = posteriors[counted, , drop = FALSE], weights = weights[counted])
}

#' Every ordering of `seq_len(n)`, one per row
#' @param n A small positive integer.
#' @return An integer matrix with `factorial(n)` rows.
#' @noRd
.latents_permutations <- function(n) {
  if (n == 1L) return(matrix(1L))
  smaller <- .latents_permutations(n - 1L)
  do.call(rbind, lapply(seq_len(n), \(first) {
    cbind(first, matrix(setdiff(seq_len(n), first)[smaller], ncol = n - 1L))
  }))
}

#' Match the classes of one solution to those of another
#'
#' The one-to-one assignment with the largest total shared posterior mass:
#' exact for square matrices of up to eight classes, and greedy otherwise
#' (the largest overlap first, ties broken by class number). A class whose
#' assigned reference class shares no mass with it is left unmatched.
#' @param overlap A classes-by-reference-classes matrix of shared mass.
#' @return For each class (row), the matched reference class, or `0`.
#' @noRd
.latents_match_classes <- function(overlap) {
  k <- nrow(overlap)
  matched <- if (k == ncol(overlap) && k <= 8L) {
    orderings <- .latents_permutations(k)
    totals <- rowSums(matrix(overlap[cbind(rep(seq_len(k), each = nrow(orderings)),
                                           as.vector(orderings))],
                             nrow = nrow(orderings)))
    # The first of equal totals, so ties resolve the same way every time.
    orderings[which.max(totals), ]
  } else .latents_greedy_match(overlap)
  shared <- overlap[cbind(seq_len(k), pmax(matched, 1L))]
  as.integer(ifelse(matched > 0L & shared > 0, matched, 0L))
}

#' Greedy one-to-one matching, the largest overlap first
#' @param overlap A classes-by-reference-classes matrix of shared mass.
#' @return For each class (row), the matched reference class, or `0`.
#' @noRd
.latents_greedy_match <- function(overlap) {
  cells <- which(!is.na(overlap), arr.ind = TRUE)
  cells <- cells[order(-overlap[cells], cells[, 1L], cells[, 2L]), , drop = FALSE]
  Reduce(\(matched, index) {
    class <- cells[index, 1L]
    reference <- cells[index, 2L]
    if (matched[class] == 0L && !reference %in% matched) matched[class] <- reference
    matched
  }, seq_len(nrow(cells)), init = integer(nrow(overlap)))
}

#' Make the colours of one row distinct
#'
#' The first class with a colour keeps it; a later class with the same colour
#' takes the first palette colour the row does not use yet. With more classes
#' than palette colours, the extra classes repeat palette colours.
#' @param colours The row's colours, indexed by class.
#' @return The colours, distinct while the palette lasts.
#' @noRd
.latents_distinct_colours <- function(colours) {
  repeated <- which(duplicated(colours))
  if (length(repeated) == 0L) return(colours)
  spare <- setdiff(.gg_okabe_ito, colours)
  replacement <- rep_len(c(spare, .gg_okabe_ito), length(repeated))
  colours[repeated] <- replacement
  colours
}

#' Shares, flows, positions and lineage colours of one tree
#'
#' Each row's classes are ordered by the flow-weighted centre of their parents,
#' so bands cross only where units really are reshuffled. A class that is its
#' dominant parent's main continuation keeps the parent's colour; any other
#' class takes a new one (Zappia and Oshlack, 2018).
#' @param x The grid.
#' @param tree One element of `.latents_tree_specs()`.
#' @param level `"individuals"` or `"groups"`.
#' @return The tree with `posteriors` and `weights` (one per row), `shares`,
#'   `flows` (between consecutive rows), `segments` (per row: `class`,
#'   `xmin`, `xmax`, `share`) and `colours` (per row, indexed by class).
#' @noRd
.latents_tree_layout <- function(x, tree, level) {
  gap <- 0.03
  drawn <- lapply(tree$rows$at, \(at) {
    .latents_tree_posteriors(.latents_candidate_at(x, at), level)
  })
  posteriors <- lapply(drawn, `[[`, "posteriors")
  weights <- drawn[[1L]]$weights
  total <- sum(weights)
  shares <- lapply(posteriors, \(p) colSums(p * weights) / total)
  flows <- lapply(seq_along(posteriors)[-1L], \(row) {
    crossprod(posteriors[[row - 1L]] * weights, posteriors[[row]]) / total
  })
  segments_of <- \(share, ordering) {
    widths <- share[ordering]
    k <- length(widths)
    right <- cumsum(widths) + gap * (seq_len(k) - 1) - gap * (k - 1) / 2
    data.frame(class = ordering, xmin = right - widths, xmax = right,
               share = widths)
  }
  first_order <- order(-shares[[1L]], seq_along(shares[[1L]]))
  segments <- Reduce(\(previous, row) {
    parent <- previous[[length(previous)]]
    flow <- flows[[row - 1L]]
    centres <- (parent$xmin + parent$xmax)[order(parent$class)] / 2
    child_centre <- colSums(flow * centres) / colSums(flow)
    c(previous, list(segments_of(shares[[row]],
                                 order(child_centre, seq_along(child_centre)))))
  }, seq_along(posteriors)[-1L], init = list(segments_of(shares[[1L]], first_order)))
  first_colours <- character(length(shares[[1L]]))
  first_colours[first_order] <- rep_len(.gg_okabe_ito, length(first_order))
  colours <- Reduce(\(state, row) {
    flow <- flows[[row - 1L]]
    dominant_parent <- apply(flow, 2L, which.max)
    main_child <- apply(flow, 1L, which.max)
    continues <- main_child[dominant_parent] == seq_len(ncol(flow))
    fresh <- rep_len(.gg_okabe_ito,
                     state$used + sum(!continues))[-seq_len(state$used)]
    row_colours <- character(ncol(flow))
    row_colours[continues] <- state$colours[[row - 1L]][dominant_parent[continues]]
    row_colours[!continues] <- fresh
    list(colours = c(state$colours, list(row_colours)),
         used = state$used + sum(!continues))
  }, seq_along(posteriors)[-1L],
  init = list(colours = list(first_colours), used = length(first_colours)))$colours
  # Rows with one number of classes are one level of different solutions:
  # each class takes the colour of the first row's class it matches.
  if (length(unique(lengths(shares))) == 1L) {
    colours <- c(colours[1L], lapply(seq_along(posteriors)[-1L], \(row) {
      matched <- .latents_match_classes(
        crossprod(posteriors[[row]] * weights, posteriors[[1L]]))
      .latents_distinct_colours(
        ifelse(matched > 0L, colours[[1L]][pmax(matched, 1L)], colours[[row]]))
    }))
  }
  colours <- lapply(colours, .latents_distinct_colours)
  c(tree, list(posteriors = posteriors, weights = weights, shares = shares,
               flows = flows, segments = segments, colours = colours))
}

#' Give every tree the colours of the matching classes of the first tree
#'
#' A row of a later tree whose key (number of classes, or of group classes)
#' also occurs in the first tree has its classes matched to that row's
#' classes by shared posterior mass, and takes their colours. Other rows keep
#' their own lineage colours.
#' @param layouts Trees from `.latents_tree_layout()`.
#' @return The trees, recoloured.
#' @noRd
.latents_tree_anchor <- function(layouts) {
  if (length(layouts) < 2L) return(layouts)
  reference <- layouts[[1L]]
  c(list(reference), lapply(layouts[-1L], \(tree) {
    tree$colours <- lapply(seq_along(tree$colours), \(row) {
      at <- match(tree$rows$key[row], reference$rows$key)
      own <- tree$colours[[row]]
      if (is.na(at) || length(reference$weights) != length(tree$weights)) return(own)
      overlap <- crossprod(tree$posteriors[[row]] * tree$weights,
                           reference$posteriors[[at]])
      matched <- .latents_match_classes(overlap)
      .latents_distinct_colours(
        ifelse(matched > 0L, reference$colours[[at]][pmax(matched, 1L)], own))
    })
    tree
  }))
}

#' The flows of every tree, one row per pair of classes in consecutive rows
#'
#' @inheritParams .latents_tree_specs
#' @return A data frame: `tree`, `level`, `from_candidate`, `from_class`,
#'   `from_share`, `to_candidate`, `to_class`, `to_share`, `mass` (the share
#'   of all units the two classes share), `share_of_from` and `share_of_to`.
#' @noRd
.latents_tree_table <- function(x, n_group_classes = NULL, n_profiles = NULL,
                                model = NULL, candidates = NULL, from = NULL,
                                to = NULL, level = c("individuals", "groups")) {
  level <- .latents_check_tree_arguments(level, n_profiles, n_group_classes, model)
  candidates <- .latents_tree_candidates(x, candidates, from, to, model)
  if (!is.null(from)) model <- NULL
  trees <- lapply(.latents_tree_specs(x, level, n_profiles, n_group_classes,
                                      model, candidates),
                  \(tree) .latents_tree_layout(x, tree, level))
  frame <- do.call(rbind, lapply(trees, \(tree) {
    do.call(rbind, lapply(seq_along(tree$flows), \(index) {
      flow <- tree$flows[[index]]
      cells <- expand.grid(from_class = seq_len(nrow(flow)),
                           to_class = seq_len(ncol(flow)))
      from_share <- tree$shares[[index]][cells$from_class]
      to_share <- tree$shares[[index + 1L]][cells$to_class]
      mass <- flow[cbind(cells$from_class, cells$to_class)]
      data.frame(tree = tree$name, level = level,
                 from_candidate = tree$rows$candidate[index],
                 from_class = cells$from_class, from_share = from_share,
                 to_candidate = tree$rows$candidate[index + 1L],
                 to_class = cells$to_class, to_share = to_share, mass = mass,
                 share_of_from = mass / from_share, share_of_to = mass / to_share,
                 stringsAsFactors = FALSE)
    }))
  }))
  rownames(frame) <- NULL
  frame
}

#' The data of drawn trees: bars and bands, with each tree's facets
#' @param trees Anchored layouts.
#' @return A list of `bars` and `bands` data frames.
#' @noRd
.gg_tree_frames <- function(trees) {
  half <- 0.16
  parts <- lapply(trees, \(tree) {
    facets <- data.frame(tree = tree$label, row_facet = tree$row_facet,
                         col_facet = tree$col_facet, stringsAsFactors = FALSE)
    position <- tree$rows$position
    bars <- do.call(rbind, lapply(seq_along(tree$segments), \(row) {
      segment <- tree$segments[[row]]
      cbind(segment, facets, row = position[row],
            fill = tree$colours[[row]][segment$class],
            stringsAsFactors = FALSE)
    }))
    bands <- do.call(rbind, lapply(seq_along(tree$flows), \(index) {
      parent <- tree$segments[[index]]
      child <- tree$segments[[index + 1L]]
      flow <- tree$flows[[index]]
      links <- expand.grid(from = parent$class, to = child$class)
      links$mass <- flow[cbind(links$from, links$to)]
      links <- links[links$mass > 0, , drop = FALSE]
      # Outflows stack inside a parent in child order and inflows inside a
      # child in parent order, so bands tile both bars without overlap.
      links$from_rank <- match(links$from, parent$class)
      links$to_rank <- match(links$to, child$class)
      links <- links[order(links$from_rank, links$to_rank), , drop = FALSE]
      links$top_right <- parent$xmin[links$from_rank] +
        stats::ave(links$mass, links$from, FUN = cumsum)
      links$top_left <- links$top_right - links$mass
      links <- links[order(links$to_rank, links$from_rank), , drop = FALSE]
      links$bottom_right <- child$xmin[links$to_rank] +
        stats::ave(links$mass, links$to, FUN = cumsum)
      links$bottom_left <- links$bottom_right - links$mass
      do.call(rbind, lapply(seq_len(nrow(links)), \(i) {
        band <- .gg_flow_band(links$top_left[i], links$top_right[i],
                              links$bottom_left[i], links$bottom_right[i],
                              y_top = position[index] + half,
                              y_bottom = position[index + 1L] - half)
        cbind(band, facets,
              group = sprintf("%s-%d-%d-%d", tree$name, index, links$from[i],
                              links$to[i]),
              fill = tree$colours[[index]][links$from[i]],
              stringsAsFactors = FALSE)
      }))
    }))
    list(bars = bars, bands = bands)
  })
  list(bars = do.call(rbind, lapply(parts, `[[`, "bars")),
       bands = do.call(rbind, lapply(parts, `[[`, "bands")))
}

#' Draw one tree, or several as panels of one plot
#'
#' @param trees Anchored layouts, all drawn in this plot.
#' @param titles A list of `main`, `subtitle`, `caption` and `y`.
#' @return A ggplot object.
#' @noRd
.gg_tree_plot <- function(trees, titles) {
  half <- 0.16
  # Rows sit at their count, so trees that start at different counts align
  # when drawn together; named candidates sit in the order given.
  keys <- sort(unique(unlist(lapply(trees, \(tree) tree$rows$key))))
  trees <- lapply(trees, \(tree) {
    tree$rows$position <- if (anyNA(tree$rows$key)) seq_len(nrow(tree$rows)) else
      match(tree$rows$key, keys)
    tree
  })
  row_labels <- if (length(keys) == 0L) trees[[1L]]$rows$label else
    as.character(keys)
  frames <- .gg_tree_frames(trees)
  bars <- frames$bars
  bands <- frames$bands
  columns <- length(unique(bars$col_facet))
  # A share is written only where its segment is wide enough to hold it; the
  # more trees side by side, the narrower every segment.
  min_share <- 0.025 + 0.04 * (columns - 1L)
  bars$label <- ifelse(bars$share < min_share, "", .gg_percent(bars$share))
  bars$text_colour <- .gg_ink(bars$fill)
  row_levels <- unique(unlist(lapply(trees, \(tree) tree$row_facet)))
  col_levels <- unique(unlist(lapply(trees, \(tree) tree$col_facet)))
  bars$row_facet <- factor(bars$row_facet, levels = row_levels)
  bands$row_facet <- factor(bands$row_facet, levels = row_levels)
  bars$col_facet <- factor(bars$col_facet, levels = col_levels)
  bands$col_facet <- factor(bands$col_facet, levels = col_levels)
  several_rows <- length(row_levels) > 1L
  several_cols <- length(col_levels) > 1L
  facets <- if (several_rows && several_cols) {
    ggplot2::facet_grid(rows = ggplot2::vars(row_facet),
                        cols = ggplot2::vars(col_facet))
  } else if (several_rows) {
    ggplot2::facet_grid(rows = ggplot2::vars(row_facet))
  } else if (several_cols) {
    ggplot2::facet_wrap(ggplot2::vars(col_facet), nrow = 1L)
  }
  ggplot2::ggplot() +
    ggplot2::geom_polygon(data = bands,
                          ggplot2::aes(x = x, y = y, group = group, fill = fill),
                          alpha = 0.4, colour = NA) +
    ggplot2::geom_rect(data = bars,
                       ggplot2::aes(xmin = xmin, xmax = xmax,
                                    ymin = row - half, ymax = row + half,
                                    fill = fill),
                       colour = "white", linewidth = 0.9) +
    ggplot2::geom_text(data = bars,
                       ggplot2::aes(x = (xmin + xmax) / 2, y = row,
                                    label = label, colour = text_colour),
                       size = 3, fontface = "bold") +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_colour_identity() +
    facets +
    ggplot2::scale_y_reverse(breaks = seq_along(row_labels), labels = row_labels,
                             expand = ggplot2::expansion(add = 0.25)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(0)) +
    .gg_titles(titles$main, titles$subtitle, titles$default_main,
               titles$default_subtitle, caption = titles$caption) +
    ggplot2::labs(x = NULL, y = titles$y) +
    .gg_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.spacing = ggplot2::unit(1.5, "lines"))
}

#' The trees of an enumeration grid, drawn
#'
#' @inheritParams .latents_tree_specs
#' @param combine `TRUE` draws every tree as a panel of one plot.
#' @param main,subtitle Titles; `NULL` uses the view's own.
#' @return A ggplot object for one tree or `combine = TRUE`; otherwise a
#'   `latents_plots` list named by tree.
#' @noRd
.gg_view_tree <- function(x, level = "individuals", n_profiles = NULL,
                          n_group_classes = NULL, model = NULL,
                          candidates = NULL, combine = FALSE, main = NULL,
                          subtitle = NULL) {
  specs <- .latents_tree_specs(x, level, n_profiles, n_group_classes, model,
                               candidates)
  trees <- .latents_tree_anchor(lapply(specs, \(tree) {
    .latents_tree_layout(x, tree, level)
  }))
  word <- .latents_class_word(x)
  groups_level <- identical(level, "groups")
  noun <- if (groups_level) "group classes" else word
  units <- if (groups_level) "groups" else "cases"
  fixed <- attr(specs, "fixed")
  rows_noun <- attr(specs, "rows_noun") %||% noun
  default_main <- if (!is.null(candidates)) {
    sprintf("How the %s of these solutions correspond", noun)
  } else if (!is.null(fixed)) {
    sprintf("The %s compared across numbers of %s", fixed, rows_noun)
  } else sprintf("How %s split as more are added", noun)
  default_subtitle <- paste0(
    if (!is.null(fixed)) {
      sprintf("Each row is the %s solution at one number of %s, divided by share of %s\n",
              attr(specs, "fixed_solution"), rows_noun, units)
    } else {
      sprintf("Each row is one solution, divided into its %s by share of %s\n",
              noun, units)
    },
    if (!is.null(fixed) || all(lengths(trees[[1L]]$shares) ==
                                 length(trees[[1L]]$shares[[1L]]))) {
      sprintf("A %s takes the colour of the %s it matches in the first row",
              .latents_singular(noun), .latents_singular(noun))
    } else if (length(trees) > 1L) {
      sprintf("A colour marks the same %s in every tree; a new colour marks one that split off",
              .latents_singular(noun))
    } else {
      sprintf("Colour follows a %s across solutions; a new colour marks one that split off",
              .latents_singular(noun))
    })
  titles <- list(
    main = main, subtitle = subtitle, default_subtitle = default_subtitle,
    # Wrapped by hand: a ggplot caption does not wrap, and runs off the plot.
    caption = if (!is.null(candidates) || !is.null(fixed)) {
      paste0("A band's width is the share of units two classes have in common.\n",
             "Bands that cross come from units classified differently, or from ",
             "uncertain classification.")
    } else {
      paste0("A band that forks: a class split in two (nested, stable solutions).\n",
             "Bands that merge or cross: units shared between classes, by ",
             "reassignment or by uncertain classification.")
    },
    y = if (!is.null(candidates)) NULL else
      sprintf("Number of %s", rows_noun))
  if (length(trees) == 1L || isTRUE(combine)) {
    titles$default_main <- if (length(trees) == 1L && nzchar(trees[[1L]]$label)) {
      sprintf("%s: %s", default_main, trees[[1L]]$label)
    } else default_main
    return(.gg_tree_plot(trees, titles))
  }
  plots <- lapply(trees, \(tree) {
    titles$default_main <- sprintf("%s: %s", default_main, tree$label)
    .gg_tree_plot(list(tree), titles)
  })
  .gg_plots(stats::setNames(plots, vapply(trees, `[[`, "", "name")))
}
