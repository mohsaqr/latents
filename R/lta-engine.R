# General latent transition engine. Every transition from occasion t - 1 to t
# of group j has its own multinomial logit, driven by a design row Z[j, t, ]:
# an intercept alone gives the homogeneous model; occasion indicators give a
# transition matrix per occasion; covariate columns (fixed or changing over
# occasions) give covariate-dependent transitions. The initial profile
# distribution is a multinomial logit on a group-level design. The reference
# outcome of a transition logit is staying in the origin profile, so a
# coefficient is the log odds of moving to that profile rather than staying.
# Measurement is the multilpa() model, invariant across occasions or with its
# own means and variances at each occasion.

#' Covariate columns as a numeric matrix, factors expanded to indicators
#' @return A matrix with one row per data row and named columns, or `NULL`.
#' @noRd
.lta_covariate_matrix <- function(data, covariates) {
  if (length(covariates) == 0L) return(NULL)
  missing_columns <- setdiff(covariates, names(data))
  if (length(missing_columns) > 0L) {
    stop(errorCondition(sprintf("Covariates %s are not columns of `data`.",
                                paste(sprintf("`%s`", missing_columns), collapse = ", ")),
                        class = "latents_bad_data", call = NULL))
  }
  # A factor level no row takes (common after subset()) would give an
  # all-zero indicator, refused as collinear; glm() drops it, and so does this
  # (TRANS-7). Every row is read, so a refit or a check on the same data
  # builds the same columns.
  frame <- droplevels(data[, covariates, drop = FALSE])
  if (anyNA(frame)) {
    stop(errorCondition("Transition and initial covariates must be complete.",
                        class = "latents_bad_data", call = NULL))
  }
  design <- stats::model.matrix(stats::reformulate(covariates), frame)
  design <- design[, colnames(design) != "(Intercept)", drop = FALSE]
  if (!all(is.finite(design))) {
    stop(errorCondition("Covariates must be finite.", class = "latents_bad_data",
                        call = NULL))
  }
  design
}

#' Initial and transition designs on the sequence layout
#'
#' @param layout Sequence layout.
#' @param covariate_rows Row-level covariate matrix for transitions, or `NULL`.
#' @param initial_rows Row-level covariate matrix for the initial
#'   distribution, or `NULL`; read at each group's first occasion.
#' @param transitions `"homogeneous"` (one intercept) or `"occasion"` (one
#'   intercept per occasion from the second on).
#' @return A list with `initial` (groups x p0), `transition` (a list over
#'   occasions 2..T of groups x p matrices; transition covariates are read at
#'   the destination occasion, zero past a group's last occasion, where no
#'   move is made; `.lta_check_covariate_gaps()` refuses skipped positions
#'   inside a span) and the column names.
#' @param order Transition order. With order 2 the first-order logits govern
#'   only the move into occasion 2 and the second-order logits the later moves,
#'   each on its own subset of columns, so each block is checked for rank on
#'   its own rows and columns.
#' @param used `NULL` (every row), or the rows the likelihood reads, from
#'   `.lta_used_rows()`: `initial` (one logical per group) and `moves`
#'   (groups by occasions). The rank is judged on those rows only, so a
#'   covariate that varies only in rows the likelihood never reads (a unit of
#'   weight zero, a move after a group's last observed occasion) is refused
#'   as the constant it is where it counts (R2-LTA-5).
#' @details A rank-deficient initial or transition design is refused with
#'   `latents_bad_covariate`, as `r3step()` refuses one: the logits of an
#'   aliased column are not identified (any split between the aliased columns
#'   fits equally well), yet they would be reported as estimates and counted
#'   as free parameters in AIC and BIC.
#' @noRd
.lta_designs <- function(layout, covariate_rows, initial_rows, transitions,
                         order = 1L, used = NULL) {
  stopifnot("`order` must be 1 or 2" = length(order) == 1L && order %in% c(1, 2))
  n_groups <- nrow(layout$slot)
  first_row <- layout$slot[cbind(seq_len(n_groups), 1L)]
  initial <- cbind(matrix(1, n_groups, 1L, dimnames = list(NULL, "(Intercept)")),
                   if (!is.null(initial_rows)) initial_rows[first_row, , drop = FALSE])
  occasion_names <- paste0("occasion_", seq_len(layout$n_occasions)[-1L])
  transition <- lapply(seq_len(layout$n_occasions)[-1L], function(t) {
    intercepts <- if (identical(transitions, "occasion")) {
      matrix(as.numeric(seq_len(layout$n_occasions)[-1L] == t), n_groups,
             length(occasion_names), byrow = TRUE,
             dimnames = list(NULL, occasion_names))
    } else matrix(1, n_groups, 1L, dimnames = list(NULL, "(Intercept)"))
    if (is.null(covariate_rows)) return(intercepts)
    rows <- layout$slot[, t]
    values <- matrix(0, n_groups, ncol(covariate_rows),
                     dimnames = list(NULL, colnames(covariate_rows)))
    observed <- !is.na(rows)
    values[observed, ] <- covariate_rows[rows[observed], , drop = FALSE]
    cbind(intercepts, values)
  })
  used <- used %||% list(initial = rep(TRUE, n_groups), moves = layout$within)
  if (any(used$initial)) {
    .lta_check_design_rank(initial[used$initial, , drop = FALSE],
                           "initial-state design (`initial_covariates`)", 1L)
  }
  .lta_check_transition_rank(layout, transition, order, used$moves)
  list(initial = initial, transition = transition)
}

#' Number of leading intercept columns of a transition design
#' @param terms Transition design column names.
#' @param n_occasions Occasions of the layout.
#' @return 1 for `"(Intercept)"`, otherwise one occasion intercept per move.
#' @noRd
.lta_n_intercepts <- function(terms, n_occasions) {
  if (identical(terms[1L], "(Intercept)")) 1L else as.integer(n_occasions - 1L)
}

#' The centring and scaling of each logit design's covariate columns
#'
#' The logit models are fitted on covariates centred at their mean and divided
#' by their standard deviation over the rows the likelihood reads (each group's
#' first occasion for the initial design, the moves made for the transition
#' designs). The leading intercept columns span the constant on every row (one
#' intercept, or occasion indicators of which exactly one is 1), so the
#' centred design spans the same models and the map between the two
#' coefficient charts is exact and linear:
#' raw slope = chart slope / scale, raw intercept = chart intercept -
#' sum(centre / scale * chart slope). Without it a covariate held at 1e8 is
#' almost collinear with the intercept for the quasi-Newton finish, which then
#' stops at a lower likelihood and reports convergence.
#'
#' @param designs Raw designs from `.lta_designs()`.
#' @param layout Sequence layout.
#' @param used Active initial and move masks from `.lta_used_rows()`, or NULL.
#' @return A list with `initial` and `transition`, each holding
#'   `n_intercepts`, `centre` and `scale` (one per covariate column).
#' @noRd
.lta_design_chart <- function(designs, layout, used = NULL) {
  used <- used %||% list(initial = rep(TRUE, nrow(designs$initial)),
                         moves = layout$within)
  part <- function(design, n_intercepts) {
    covariates <- design[, seq_len(ncol(design)) > n_intercepts, drop = FALSE]
    centre <- colMeans(covariates)
    spread <- sqrt(colMeans(sweep(covariates, 2L, centre, "-")^2))
    # A column with no spread was refused by the rank check; a unit scale
    # keeps the map defined should one reach here by another route.
    list(n_intercepts = n_intercepts, centre = centre,
         scale = ifelse(is.finite(spread) & spread > 0, spread, 1))
  }
  moves <- do.call(rbind, lapply(seq_along(designs$transition), function(index) {
    designs$transition[[index]][used$moves[, index + 1L], , drop = FALSE]
  }))
  list(initial = part(designs$initial[used$initial, , drop = FALSE], 1L),
       transition = part(moves, .lta_n_intercepts(colnames(moves), layout$n_occasions)))
}

#' Designs in the fitting chart: covariate columns centred and scaled
#' @param designs Raw designs.
#' @param chart From `.lta_design_chart()`.
#' @return Designs of the same shape and column names.
#' @noRd
.lta_chart_designs <- function(designs, chart) {
  standardise <- function(design, part) {
    covariates <- seq_len(ncol(design)) > part$n_intercepts
    design[, covariates] <- sweep(sweep(design[, covariates, drop = FALSE], 2L,
                                        part$centre, "-"), 2L, part$scale, "/")
    design
  }
  list(initial = standardise(designs$initial, chart$initial),
       transition = lapply(designs$transition, standardise, part = chart$transition))
}

#' The linear map between coefficient charts of one design
#' @param part One element of `.lta_design_chart()`.
#' @param direction `"to_raw"` (chart to supplied units) or `"to_chart"`.
#' @return A square matrix M with coefficients_out = M coefficients_in, over
#'   all columns of the design.
#' @noRd
.lta_chart_matrix <- function(part, direction = c("to_raw", "to_chart")) {
  direction <- match.arg(direction)
  n_covariates <- length(part$centre)
  n_columns <- part$n_intercepts + n_covariates
  map <- diag(n_columns)
  if (n_covariates == 0L) return(map)
  covariates <- part$n_intercepts + seq_len(n_covariates)
  intercepts <- seq_len(part$n_intercepts)
  if (identical(direction, "to_raw")) {
    map[cbind(covariates, covariates)] <- 1 / part$scale
    map[intercepts, covariates] <- matrix(-part$centre / part$scale,
                                          part$n_intercepts, n_covariates, byrow = TRUE)
  } else {
    map[cbind(covariates, covariates)] <- part$scale
    map[intercepts, covariates] <- matrix(part$centre, part$n_intercepts,
                                          n_covariates, byrow = TRUE)
  }
  map
}

#' Map the initial and transition coefficients between charts
#'
#' The same matrix applies to every class, profile and origin, over the first
#' (design-column) dimension. A column a block does not use (an occasion
#' intercept outside a second-order block) is moved with the rest, and moved
#' back exactly, without entering any likelihood.
#'
#' @param parameters A parameter list with `initial`, `transition` and
#'   optionally `transition2` coefficient arrays.
#' @param chart From `.lta_design_chart()`.
#' @param direction `"to_raw"` or `"to_chart"`.
#' @return `parameters` with those arrays mapped; dimnames and attributes kept.
#' @noRd
.lta_map_coefficients <- function(parameters, chart, direction = c("to_raw", "to_chart")) {
  direction <- match.arg(direction)
  apply_map <- function(coefficients, part) {
    if (is.null(coefficients)) return(NULL)
    map <- .lta_chart_matrix(part, direction)
    stopifnot("coefficient rows must match the design columns" =
                dim(coefficients)[1L] == nrow(map))
    mapped <- map %*% matrix(coefficients, nrow(map))
    attributes(mapped) <- attributes(coefficients)
    mapped
  }
  parameters$initial <- apply_map(parameters$initial, chart$initial)
  parameters$transition <- apply_map(parameters$transition, chart$transition)
  if (!is.null(parameters$transition2)) {
    parameters$transition2 <- apply_map(parameters$transition2, chart$transition)
  }
  parameters
}

#' Jacobian of the estimation-scale vector from the fitting chart to the raw one
#'
#' theta_raw = J theta_chart, with J the identity on every coordinate but the
#' initial and transition logits, whose blocks are the free-column submatrix
#' of `.lta_chart_matrix()` (every covariate column is free, so a free
#' intercept depends only on itself and the covariates).
#'
#' @param fit A fit or its specification (raw designs, `layout`, `order`,
#'   `n_profiles`, `n_group_classes`, `stayer`).
#' @param chart From `.lta_design_chart()`.
#' @param n_theta Length of the estimation-scale vector.
#' @param direction `"to_raw"` gives J, `"to_chart"` its inverse.
#' @return An `n_theta` x `n_theta` matrix.
#' @noRd
.lta_chart_jacobian <- function(fit, chart, n_theta, direction = c("to_raw", "to_chart")) {
  direction <- match.arg(direction)
  n_profiles <- fit$n_profiles
  n_classes <- fit$n_group_classes
  n_movers <- sum(!(fit$stayer %||% rep(FALSE, n_classes)))
  terms <- colnames(fit$designs$transition[[1L]])
  second <- (fit$order %||% 1L) >= 2L
  n_intercepts <- fit$layout$n_occasions - 1L
  first <- .lta_transition_columns(terms, second, 1L, n_intercepts)
  later <- .lta_transition_columns(terms, second, 2L, n_intercepts)
  transition_map <- .lta_chart_matrix(chart$transition, direction)
  blocks <- c(
    rep(list(.lta_chart_matrix(chart$initial, direction)), n_classes * (n_profiles - 1L)),
    rep(list(transition_map[first, first, drop = FALSE]),
        n_movers * n_profiles * (n_profiles - 1L)),
    if (second) rep(list(transition_map[later, later, drop = FALSE]),
                    n_movers * n_profiles^2 * (n_profiles - 1L)))
  sizes <- vapply(blocks, nrow, integer(1))
  offset <- n_theta - sum(sizes)
  stopifnot("the logit blocks must fit in the estimation vector" = offset >= 0L)
  jacobian <- diag(n_theta)
  starts <- offset + cumsum(c(0L, sizes[-length(sizes)]))
  # Fills disjoint diagonal blocks of one matrix in place.
  invisible(lapply(seq_along(blocks), function(index) {
    rows <- starts[index] + seq_len(sizes[index])
    jacobian[rows, rows] <<- blocks[[index]]
  }))
  jacobian
}

#' Refuse rank-deficient transition designs
#'
#' Only the moves the likelihood reads are checked: past a group's last
#' occasion its design rows are zero placeholders, and a move after a group's
#' last observed occasion, or of a unit of weight zero, does not enter the
#' likelihood. An occasion intercept none of those moves uses is flat and not
#' counted (`.lta_identified_columns()`), so it is left out rather than read
#' as an aliased column.
#'
#' @param layout Sequence layout.
#' @param transition The list of per-occasion transition designs.
#' @param order Transition order.
#' @param moves_used Groups by occasions: which moves enter the likelihood.
#' @return `NULL`, invisibly; raises `latents_bad_covariate`.
#' @noRd
.lta_check_transition_rank <- function(layout, transition, order,
                                       moves_used = layout$within) {
  terms <- colnames(transition[[1L]])
  moves <- lapply(seq_along(transition), function(index) {
    transition[[index]][moves_used[, index + 1L], , drop = FALSE]
  })
  check <- function(design, what, n_intercepts) {
    if (nrow(design) == 0L) return(invisible(NULL))
    live <- seq_len(ncol(design)) > n_intercepts | colSums(design != 0) > 0
    .lta_check_design_rank(design[, live, drop = FALSE], what,
                           sum(live[seq_len(n_intercepts)]))
  }
  # The intercept columns lead: one "(Intercept)", or one indicator per
  # occasion, which sum to one on every move row.
  n_intercepts <- if (identical(terms[1L], "(Intercept)")) 1L else layout$n_occasions - 1L
  if (order < 2L) {
    check(do.call(rbind, moves),
          "transition design (`transitions`, `transition_covariates`)", n_intercepts)
    return(invisible(NULL))
  }
  # With occasion intercepts, the first block keeps only occasion_2 and the
  # second block the other n_occasions - 2.
  occasion <- n_intercepts > 1L
  first <- .lta_transition_columns(terms, TRUE, 1L, layout$n_occasions - 1L)
  check(moves[[1L]][, first, drop = FALSE],
        "first-order transition design (the move into occasion 2)", 1L)
  second <- .lta_transition_columns(terms, TRUE, 2L, layout$n_occasions - 1L)
  check(do.call(rbind, moves[-1L])[, second, drop = FALSE],
        "second-order transition design (the moves after occasion 2)",
        if (occasion) n_intercepts - 1L else 1L)
  invisible(NULL)
}

#' Refuse a rank-deficient logit design, naming it and its aliased columns
#'
#' The rule of `.multilpa_r3step_check_rank()`: columns are put on a common
#' scale first, so a covariate held in large or small units is judged by its
#' pattern, not by its magnitude, and the QR rank (default tolerance) is
#' compared with the number of columns. The covariate columns are centred
#' before that: the leading intercept columns span the constant vector (one
#' intercept, or occasion indicators that sum to one on every row), so
#' centring leaves the column space and the rank unchanged, and it keeps a
#' covariate held at a large origin (`w + 1e8`) from looking like the
#' intercept to the QR tolerance.
#'
#' @param design A design matrix with named columns, intercept columns first.
#' @param what The design's name, for the message.
#' @param n_intercepts Number of leading intercept columns.
#' @return `design`, invisibly; raises `latents_bad_covariate`.
#' @noRd
.lta_check_design_rank <- function(design, what, n_intercepts) {
  stopifnot("`design` must be a numeric matrix" = is.matrix(design) && is.numeric(design),
            "`what` must be a single string" = is.character(what) && length(what) == 1L,
            "`n_intercepts` must lie between 1 and the column count" =
              length(n_intercepts) == 1L && n_intercepts >= 1L &&
              n_intercepts <= ncol(design))
  covariates <- seq_len(ncol(design)) > n_intercepts
  centred <- design
  centred[, covariates] <- sweep(design[, covariates, drop = FALSE], 2L,
                                 colMeans(design[, covariates, drop = FALSE]), "-")
  decomposition <- qr(sweep(centred, 2L, .multilpa_design_scale(centred), "/"))
  if (decomposition$rank < ncol(design)) {
    aliased <- colnames(design)[decomposition$pivot[-seq_len(decomposition$rank)]]
    # Occasion intercepts are generated, so a covariate aliased with them is
    # one that takes a single value at each occasion (the occasion number).
    hint <- if (any(grepl("^occasion_[0-9]+$", colnames(design)))) paste(
      " With transitions = \"occasion\", a covariate that is constant within",
      "each occasion is collinear with the occasion intercepts.") else ""
    stop(errorCondition(sprintf(paste0(
      "The %s is rank deficient (rank %d of %d columns): %s cannot be told ",
      "apart from the other columns, so its logits are not identified. Remove ",
      "constant or collinear covariates.%s"),
      what, decomposition$rank, ncol(design),
      paste(sprintf("`%s`", aliased), collapse = ", "), hint),
      class = "latents_bad_covariate", call = NULL))
  }
  invisible(design)
}

#' Refuse transition covariates at grid positions a group skipped
#'
#' Under `occasions = "grid"` a move into a skipped position has no data row,
#' so its covariates have no value. Any fill-in (zero, the last value) would be
#' an assumption the likelihood then depends on: with zero, shifting a
#' covariate by a constant changed the fit. The move is refused instead.
#'
#' @param layout Sequence layout.
#' @param groups Prepared groups (`ids`, `index`).
#' @param time_values Each row's position.
#' @param transition_covariates Transition covariate names.
#' @return `NULL`, invisibly; raises `latents_bad_data`.
#' @noRd
.lta_check_covariate_gaps <- function(layout, groups, time_values,
                                      transition_covariates) {
  if (length(transition_covariates) == 0L || layout$n_occasions < 2L) {
    return(invisible(NULL))
  }
  # Inside a group's span, a slot without a row is a skipped position. The
  # first slot is always observed, so only destinations can be skipped.
  skipped <- which(layout$within & is.na(layout$slot), arr.ind = TRUE)
  if (nrow(skipped) == 0L) return(invisible(NULL))
  positions <- sort(unique(time_values))
  first <- as.integer(tapply(match(time_values, positions), groups$index, min))
  cells <- unique(data.frame(
    group = groups$ids[skipped[, "row"]],
    position = positions[first[skipped[, "row"]] + skipped[, "col"] - 1L]))
  shown <- utils::head(cells, 5L)
  stop(errorCondition(sprintf(paste(
    "Transition covariates have no value where a group skips a grid position:",
    "%s%s. Use occasions = \"observed\", or add those rows (indicators may be",
    "missing under missing = \"fiml\") with their covariate values."),
    paste(sprintf("group %s at %s", shown$group, shown$position), collapse = ", "),
    if (nrow(cells) > nrow(shown)) sprintf(" and %d more", nrow(cells) - nrow(shown)) else ""),
    class = "latents_bad_data", call = NULL))
}

#' Log transition probabilities of one class at one occasion
#'
#' @param design Groups x p design at that occasion.
#' @param coefficients p x (K - 1) x K array for this class: for origin k, the
#'   logits of the K - 1 destinations other than k, in profile order.
#' @return A groups x K x K array of log probabilities, origin by destination.
#' @noRd
.lta_log_transition <- function(design, coefficients) {
  n_profiles <- dim(coefficients)[3L]
  n_groups <- nrow(design)
  out <- array(0, c(n_groups, n_profiles, n_profiles))
  invisible(lapply(seq_len(n_profiles), function(k) {
    others <- setdiff(seq_len(n_profiles), k)
    # Staying is the reference category: .multilpa_log_softmax() puts it last.
    log_p <- .multilpa_log_softmax(design, matrix(coefficients[, , k], ncol(design)))
    out[, k, others] <<- log_p[, seq_along(others), drop = FALSE]
    out[, k, k] <<- log_p[, length(others) + 1L]
  }))
  out
}

#' Log initial probabilities of one class: groups x K (reference: last profile)
#' @noRd
.lta_log_initial <- function(design, coefficients) {
  .multilpa_log_softmax(design, coefficients)
}

#' Posterior occupancy and pair probabilities of one class, per group and
#' occasion, weighted by the group's class posterior
#' @return A list with `posterior` (rows x K), `initial` (groups x K) and
#'   `pairs` (list over occasions 2..T of groups x K x K arrays).
#' @noRd
.lta_moments <- function(pass, emission, layout, log_transition, weights,
                         n_observations) {
  n_profiles <- ncol(pass$alpha[[1L]])
  n_groups <- nrow(layout$slot)
  occupancy <- lapply(seq_len(layout$n_occasions), function(t) {
    exp(pass$alpha[[t]] + pass$beta[[t]] - pass$log_scaled) * weights
  })
  posterior <- matrix(0, n_observations, n_profiles)
  invisible(lapply(seq_len(layout$n_occasions), function(t) {
    rows <- layout$slot[, t]
    observed <- !is.na(rows)
    if (any(observed)) posterior[rows[observed], ] <<- occupancy[[t]][observed, , drop = FALSE]
  }))
  pairs <- lapply(seq_len(layout$n_occasions)[-1L], function(t) {
    active <- weights * layout$within[, t]
    from <- pass$alpha[[t - 1L]] - pass$log_scaled
    to <- pass$beta[[t]] + emission[[t]]
    # The pair log probability is formed whole before exponentiating, so the
    # scaled forward and backward factors cannot overflow apart.
    joint <- array(from, c(n_groups, n_profiles, n_profiles)) +
      aperm(array(to, c(n_groups, n_profiles, n_profiles)), c(1L, 3L, 2L)) +
      log_transition[[t - 1L]]
    exp(joint) * active
  })
  list(posterior = posterior, initial = occupancy[[1L]], pairs = pairs)
}

#' Pair-state moments back to profile moments
#'
#' Occupancy and initial counts sum over the previous profile; the first move
#' becomes first-order counts (origin j, destination l); every later move is
#' kept by origin pair (i, j) and destination l.
#' @noRd
.lta_reduce_moments <- function(moments, n_states) {
  n_pairs <- n_states * n_states
  current <- rep(seq_len(n_states), times = n_states)
  collapse <- function(values) {
    # Sum the columns of each current profile.
    vapply(seq_len(n_states), function(j) rowSums(values[, current == j, drop = FALSE]),
           numeric(nrow(values)))
  }
  pair_to <- function(i, j) (j - 1L) * n_states + seq_len(n_states)
  first <- moments$pairs[[1L]]
  first_order <- array(0, c(dim(first)[1L], n_states, n_states))
  invisible(lapply(seq_len(n_pairs), function(s) {
    i <- (s - 1L) %/% n_states + 1L
    j <- (s - 1L) %% n_states + 1L
    first_order[, j, ] <<- first_order[, j, ] + first[, s, pair_to(i, j)]
  }))
  later <- lapply(moments$pairs[-1L], function(pair) {
    out <- array(0, c(dim(pair)[1L], n_pairs, n_states))
    invisible(lapply(seq_len(n_pairs), function(s) {
      i <- (s - 1L) %/% n_states + 1L
      j <- (s - 1L) %% n_states + 1L
      out[, s, ] <<- pair[, s, pair_to(i, j)]
    }))
    out
  })
  list(posterior = matrix(collapse(moments$posterior), nrow(moments$posterior)),
       initial = matrix(collapse(moments$initial), nrow(moments$initial)),
       pairs = list(first_order), pairs2 = later)
}

#' Measurement log densities, rows x K, invariant or per occasion
#' @noRd
.lta_log_density <- function(x, codes, parameters, occasion_of_row, extra = NULL) {
  n_profiles <- nrow(parameters$measurement[[1L]]$means)
  out <- matrix(0, nrow(x), n_profiles)
  moments <- vector("list", length(parameters$measurement))
  invisible(lapply(seq_along(parameters$measurement), function(t) {
    rows <- if (length(parameters$measurement) == 1L) seq_len(nrow(x)) else
      which(occasion_of_row == t)
    if (length(rows) == 0L) return(NULL)
    # Missing indicators are integrated out and full covariances used by the
    # cross-sectional moments, which the M-step reuses.
    measurement <- .latents_measurement_log_density(
      x[rows, , drop = FALSE], parameters$measurement[[t]],
      if (is.null(codes)) NULL else codes[rows, , drop = FALSE],
      if (is.null(extra)) NULL else .latents_extra_rows(extra, rows),
      offset = FALSE, n_profiles = n_profiles)
    out[rows, ] <<- measurement$log_density
    moments[t] <<- list(measurement$moments)
  }))
  structure(out, moments = moments)
}

#' Second-order transitions as a first-order chain on profile pairs
#'
#' State (i, j) at occasion t holds the profile at t - 1 (i) and at t (j);
#' index (i - 1) K + j. Occasion 1 occupies only the diagonal pairs (j, j)
#' with the initial probabilities; the move to occasion 2 uses the
#' first-order logits; every later move (i, j) -> (j, l) uses the
#' second-order logits of origin pair (i, j). Emissions read the current
#' profile j.
#' @return A list with the augmented `log_initial` (groups x K^2) and
#'   `log_transition` (list of groups x K^2 x K^2 arrays, -Inf where a move is
#'   impossible).
#' @noRd
.lta_augment <- function(log_initial, first_order, second_order) {
  n_groups <- nrow(log_initial)
  n_states <- ncol(log_initial)
  n_pairs <- n_states * n_states
  pair <- function(i, j) (i - 1L) * n_states + j
  initial <- matrix(-Inf, n_groups, n_pairs)
  initial[, pair(seq_len(n_states), seq_len(n_states))] <- log_initial
  expand <- function(per_origin) {
    # per_origin(i, j) gives the groups x K log probabilities of the next
    # profile l from pair (i, j).
    out <- array(-Inf, c(n_groups, n_pairs, n_pairs))
    cells <- expand.grid(i = seq_len(n_states), j = seq_len(n_states))
    invisible(lapply(seq_len(nrow(cells)), function(r) {
      i <- cells$i[r]
      j <- cells$j[r]
      out[, pair(i, j), pair(j, seq_len(n_states))] <<- per_origin(i, j)
    }))
    out
  }
  transition <- c(
    list(expand(function(i, j) matrix(first_order[[1L]][, j, ], n_groups))),
    lapply(second_order, function(log_p) {
      expand(function(i, j) matrix(log_p[, pair(i, j), ], n_groups))
    }))
  list(log_initial = initial, log_transition = transition)
}

#' Log transitions of a stayer class: never leaving the current profile
#' @noRd
.lta_identity_transition <- function(n_groups, n_profiles) {
  out <- array(-Inf, c(n_groups, n_profiles, n_profiles))
  invisible(lapply(seq_len(n_profiles), function(k) out[, k, k] <<- 0))
  out
}

#' Second-order version: from pair (i, j) the next profile is j
#' @noRd
.lta_identity_transition2 <- function(n_groups, n_profiles) {
  out <- array(-Inf, c(n_groups, n_profiles^2, n_profiles))
  invisible(lapply(seq_len(n_profiles^2), function(pair) {
    out[, pair, (pair - 1L) %% n_profiles + 1L] <<- 0
  }))
  out
}

#' Log second-order transition probabilities of one class at one occasion:
#' groups x K^2 origin pairs x K destinations (reference: staying in j)
#' @noRd
.lta_log_transition2 <- function(design, coefficients) {
  n_states <- dim(coefficients)[2L] + 1L
  n_groups <- nrow(design)
  out <- array(0, c(n_groups, n_states * n_states, n_states))
  invisible(lapply(seq_len(n_states * n_states), function(s) {
    j <- (s - 1L) %% n_states + 1L
    others <- setdiff(seq_len(n_states), j)
    log_p <- .multilpa_log_softmax(design, matrix(coefficients[, , s], ncol(design)))
    out[, s, others] <<- log_p[, seq_along(others), drop = FALSE]
    out[, s, j] <<- log_p[, length(others) + 1L]
  }))
  out
}

#' E-step of the general transition model
#' @noRd
.lta_expectation <- function(x, codes, layout, designs, parameters, occasion_of_row,
                             sampling_weights = NULL, extra = NULL) {
  n_types <- length(parameters$group_probabilities)
  n_groups <- nrow(layout$slot)
  log_density <- .lta_log_density(x, codes, parameters, occasion_of_row, extra)
  gaussian_moments <- attr(log_density, "moments")
  attr(log_density, "moments") <- NULL
  emission <- .multilpa_sequence_emission(log_density, layout, validate = FALSE)
  second <- !is.null(parameters$transition2)
  n_states <- nrow(parameters$measurement[[1L]]$means)
  # Emissions on the augmented pair states read the current profile.
  augmented_emission <- if (second) lapply(emission$log_density, function(block) {
    block[, rep(seq_len(n_states), times = n_states), drop = FALSE]
  }) else NULL
  stayer <- parameters$stayer %||% rep(FALSE, n_types)
  per_class <- lapply(seq_len(n_types), function(h) {
    log_initial <- .lta_log_initial(designs$initial,
                                    matrix(parameters$initial[, , h], ncol(designs$initial)))
    log_transition <- lapply(designs$transition, function(design) {
      if (stayer[h]) return(.lta_identity_transition(nrow(design), n_states))
      .lta_log_transition(design, array(parameters$transition[, , , h],
                                        dim(parameters$transition)[1:3]))
    })
    if (!second) {
      return(list(log_initial = log_initial, log_transition = log_transition,
                  pass = .latents_forward_backward(emission$log_density, layout,
                                               log_initial, log_transition)))
    }
    log_transition2 <- lapply(designs$transition[-1L], function(design) {
      if (stayer[h]) return(.lta_identity_transition2(nrow(design), n_states))
      .lta_log_transition2(design, array(parameters$transition2[, , , h],
                                         dim(parameters$transition2)[1:3]))
    })
    augmented <- .lta_augment(log_initial, log_transition, log_transition2)
    list(log_initial = log_initial, log_transition = log_transition,
         log_transition2 = log_transition2, augmented = augmented,
         pass = .latents_forward_backward(augmented_emission, layout,
                                      augmented$log_initial, augmented$log_transition))
  })
  scores <- matrix(vapply(per_class, function(c) c$pass$log_scaled, numeric(n_groups)),
                   n_groups, n_types)
  weighted <- sweep(scores, 2L, log(parameters$group_probabilities), "+")
  group_log_likelihood <- .multilpa_log_sum_exp(weighted) + emission$offset
  group_posteriors <- exp(weighted - .multilpa_row_max(weighted))
  group_posteriors <- group_posteriors / rowSums(group_posteriors)
  # Sampling weights scale each person's posteriors, and so every moment the
  # M-step and the scores read, into weighted counts; pseudo maximum
  # likelihood then needs nothing else. Integer weights equal duplication.
  if (!is.null(sampling_weights)) {
    group_posteriors <- group_posteriors * sampling_weights
  }
  moments <- lapply(seq_len(n_types), function(h) {
    if (!second) {
      return(.lta_moments(per_class[[h]]$pass, emission$log_density, layout,
                          per_class[[h]]$log_transition, group_posteriors[, h],
                          nrow(x)))
    }
    .lta_reduce_moments(.lta_moments(per_class[[h]]$pass, augmented_emission, layout,
                                     per_class[[h]]$augmented$log_transition,
                                     group_posteriors[, h], nrow(x)), n_states)
  })
  joint <- lapply(moments, `[[`, "posterior")
  subject_posteriors <- Reduce(`+`, joint)
  log_likelihood <- sum((sampling_weights %||% 1) * group_log_likelihood)
  if (!is.finite(log_likelihood) || any(!is.finite(subject_posteriors))) {
    stop("Non-finite likelihood or posterior probabilities.")
  }
  list(log_likelihood = log_likelihood, group_log_likelihood = group_log_likelihood,
       group_posteriors = group_posteriors, subject_posteriors = subject_posteriors,
       joint = joint, moments = moments, per_class = per_class,
       gaussian_moments = gaussian_moments, weighted = !is.null(sampling_weights))
}

#' M-step of the general transition model
#' @noRd
.lta_maximization <- function(x, codes, layout, designs, expectation, parameters,
                              occasion_of_row, variance_model, min_variance,
                              n_categories, min_probability,
                              covariance_model = "diagonal", structure = NULL,
                              extra = NULL) {
  n_types <- length(parameters$group_probabilities)
  n_profiles <- ncol(expectation$subject_posteriors)
  # Measurement: the cross-sectional M-step, on all rows or occasion by occasion.
  measurement <- lapply(seq_along(parameters$measurement), function(t) {
    rows <- if (length(parameters$measurement) == 1L) seq_len(nrow(x)) else
      which(occasion_of_row == t)
    sub <- list(subject_posteriors = expectation$subject_posteriors[rows, , drop = FALSE],
                group_posteriors = expectation$group_posteriors,
                joint = lapply(expectation$joint, function(j) j[rows, , drop = FALSE]),
                gaussian_moments = expectation$gaussian_moments[[t]])
    # An occasion reached only by units of weight zero carries no weighted
    # membership: its block is absent from the likelihood (flat, uncounted
    # by .lta_identified_blocks()), so it keeps its values, as a trailing
    # occasion with nothing observed does in effect (TRANS-3).
    if (isTRUE(expectation$weighted) && !(sum(sub$subject_posteriors) > 0)) {
      return(parameters$measurement[[t]])
    }
    # Weighted posteriors sum to the weighted row count, which the shared
    # variance blocks divide by; unweighted, that count is the row count.
    if (isTRUE(expectation$weighted)) sub$n_total <- sum(sub$subject_posteriors)
    updated <- .multilpa_maximization(x[rows, , drop = FALSE], sub, variance_model,
                                      min_variance, covariance_model,
                                      if (is.null(codes)) NULL else codes[rows, , drop = FALSE],
                                      n_categories, min_probability,
                                      structure = structure,
                                      previous = parameters$measurement[[t]],
                                      extra = .latents_extra_rows(extra, rows))
    updated[intersect(c("means", "variances", "covariances", "response_probabilities",
                        "ordinal_intercepts", "ordinal_locations", "count_means",
                        "count_dispersion"),
                      names(updated))]
  })
  initial <- array(vapply(seq_len(n_types), function(h) {
    counts <- expectation$moments[[h]]$initial
    keep <- rowSums(counts) > 0
    .multilpa_weighted_logits(designs$initial[keep, , drop = FALSE],
                              counts[keep, , drop = FALSE],
                              matrix(parameters$initial[, , h],
                                     ncol(designs$initial)))$coefficients
  }, numeric(ncol(designs$initial) * (n_profiles - 1L))),
  c(ncol(designs$initial), n_profiles - 1L, n_types))
  p <- ncol(designs$transition[[1L]])
  second <- !is.null(parameters$transition2)
  # With second-order transitions the first-order logits govern only the move
  # into occasion 2.
  stacked_design <- if (second) designs$transition[[1L]] else
    do.call(rbind, designs$transition)
  stayer <- parameters$stayer %||% rep(FALSE, n_types)
  transition <- array(vapply(seq_len(n_types), function(h) {
    if (stayer[h]) return(as.vector(parameters$transition[, , , h]))
    vapply(seq_len(n_profiles), function(k) {
      others <- setdiff(seq_len(n_profiles), k)
      # Counts of moves out of k, destination columns with staying last.
      counts <- do.call(rbind, lapply(expectation$moments[[h]]$pairs, function(pair) {
        matrix(pair[, k, c(others, k)], nrow(pair))
      }))
      keep <- rowSums(counts) > 0
      if (!any(keep)) return(as.vector(parameters$transition[, , k, h]))
      as.vector(.multilpa_weighted_logits(
        stacked_design[keep, , drop = FALSE], counts[keep, , drop = FALSE],
        matrix(parameters$transition[, , k, h], p))$coefficients)
    }, numeric(p * (n_profiles - 1L)))
  }, numeric(p * (n_profiles - 1L) * n_profiles)),
  c(p, n_profiles - 1L, n_profiles, n_types))
  transition2 <- if (!second) NULL else {
    later_design <- do.call(rbind, designs$transition[-1L])
    array(vapply(seq_len(n_types), function(h) {
      if (stayer[h]) return(as.vector(parameters$transition2[, , , h]))
      vapply(seq_len(n_profiles^2), function(pair) {
        j <- (pair - 1L) %% n_profiles + 1L
        others <- setdiff(seq_len(n_profiles), j)
        counts <- do.call(rbind, lapply(expectation$moments[[h]]$pairs2, function(moves) {
          matrix(moves[, pair, c(others, j)], nrow(moves))
        }))
        keep <- rowSums(counts) > 0
        if (!any(keep)) return(as.vector(parameters$transition2[, , pair, h]))
        as.vector(.multilpa_weighted_logits(
          later_design[keep, , drop = FALSE], counts[keep, , drop = FALSE],
          matrix(parameters$transition2[, , pair, h], p))$coefficients)
      }, numeric(p * (n_profiles - 1L)))
    }, numeric(p * (n_profiles - 1L) * n_profiles^2)),
    c(p, n_profiles - 1L, n_profiles^2, n_types))
  }
  dimnames(transition) <- list(colnames(designs$transition[[1L]]), NULL, NULL, NULL)
  attr(transition, "n_occasion_intercepts") <- layout$n_occasions - 1L
  c(list(measurement = measurement, initial = initial, transition = transition),
    if (second) list(transition2 = transition2),
    list(group_probabilities = if (isTRUE(expectation$weighted))
      colSums(expectation$group_posteriors) / sum(expectation$group_posteriors) else
        colMeans(expectation$group_posteriors),
         stayer = parameters$stayer))
}

#' Each row's sequence weight, from the occasion layout
#' @param layout The sequence layout (`slot`: sequences by occasions of row
#'   indices).
#' @param sampling_weights `NULL`, or one weight per sequence.
#' @param n_rows Number of rows.
#' @return `NULL` without weights; otherwise one weight per row, `0` for a row
#'   the layout does not place.
#' @noRd
.lta_row_weights <- function(layout, sampling_weights, n_rows) {
  if (is.null(sampling_weights)) return(NULL)
  placed <- !is.na(layout$slot)
  weights <- numeric(n_rows)
  weights[layout$slot[placed]] <- sampling_weights[row(layout$slot)[placed]]
  weights
}

#' EM for the general transition model from one start
#'
#' `pause` and `resume` are as for `.latents_em()`; with either, the result
#' also carries the `checkpoint`.
#' @noRd
.lta_em <- function(x, codes, layout, designs, parameters, occasion_of_row,
                    variance_model, min_variance, n_categories, min_probability,
                    max_iter, tol, covariance_model = "diagonal", structure = NULL,
                    sampling_weights = NULL, extra = NULL, pause = Inf,
                    resume = NULL, acceleration = "squarem") {
  # The inner logit Newton steps of the transition M-step leave more rounding
  # in the likelihood than a closed-form step, hence the looser guard.
  fit <- .latents_em(
    parameters,
    evaluate = function(point) {
      .lta_expectation(x, codes, layout, designs, point, occasion_of_row,
                       sampling_weights, extra)
    },
    maximize = function(point, expectation) {
      .lta_maximization(x, codes, layout, designs, expectation, point,
                        occasion_of_row, variance_model, min_variance,
                        n_categories, min_probability, covariance_model,
                        structure, extra)
    },
    max_iter = max_iter, tol = tol,
    settings = .latents_em_settings(
      decrease_tolerance = 1e-8, accelerate = .lta_accelerator(acceleration, general = TRUE),
      # A unit-free stop (review 2026-10-26 round 2, K2-6): the observed values
      # of the rows whose sequence has positive weight.
      stop_scale = .latents_observed_values(
        x, codes, .latents_extra_matrix(extra),
        row_weights = .lta_row_weights(layout, sampling_weights, nrow(x)))),
    pause = pause, resume = resume)
  result <- list(parameters = fit$state, expectation = fit$expectation,
                 converged = fit$converged, iterations = fit$iterations,
                 history = fit$history)
  if (!is.null(fit$checkpoint)) result$checkpoint <- fit$checkpoint
  result
}

#' Starting values: the homogeneous initializer, converted to logit
#' coefficients with zero covariate effects
#' @param row_weights Optional sampling weights in observation order, as for
#'   `.multilpa_initialize()`.
#' @noRd
.lta_initialize <- function(x, group_index, layout, designs, n_profiles, n_types,
                            variance_model, min_variance, start_index, codes,
                            n_categories, min_probability, n_measurement,
                            order = 1L, mover_stayer = 0L,
                            covariance_model = "diagonal", structure = NULL,
                            extra = NULL, kind = "kmeans", row_weights = NULL) {
  start <- .multilpa_transition_initialize(x, group_index, n_profiles, n_types,
                                           variance_model, min_variance, start_index,
                                           covariance_model, codes, n_categories,
                                           min_probability, extra, kind = kind,
                                           row_weights = row_weights)
  if (!is.null(structure)) {
    start <- .multilpa_project_start(start, structure, nrow(x), min_variance)
  }
  block <- start[intersect(c("means", "variances", "covariances",
                             "response_probabilities", "ordinal_intercepts",
                             "ordinal_locations", "count_means", "count_dispersion"),
                           names(start))]
  p0 <- ncol(designs$initial)
  p <- ncol(designs$transition[[1L]])
  initial <- array(0, c(p0, n_profiles - 1L, n_types))
  transition <- array(0, c(p, n_profiles - 1L, n_profiles, n_types),
                      dimnames = list(colnames(designs$transition[[1L]]), NULL, NULL, NULL))
  attr(transition, "n_occasion_intercepts") <- layout$n_occasions - 1L
  # Intercept columns carry the starting logits; covariate columns start at 0.
  intercepts_initial <- which(colnames(designs$initial) == "(Intercept)")
  intercepts_transition <- if (identical(colnames(designs$transition[[1L]])[1L],
                                        "(Intercept)")) 1L else
    seq_len(layout$n_occasions - 1L)
  invisible(lapply(seq_len(n_types), function(h) {
    shares <- start$initial_probabilities[h, ]
    initial[intercepts_initial, , h] <<- log(shares[-n_profiles] / shares[n_profiles])
    tp <- matrix(start$transition_probabilities[, , h], n_profiles)
    lapply(seq_len(n_profiles), function(k) {
      others <- setdiff(seq_len(n_profiles), k)
      logits <- log(tp[k, others] / tp[k, k])
      transition[intercepts_transition, , k, h] <<-
        matrix(logits, length(intercepts_transition), length(others), byrow = TRUE)
    })
  }))
  # Second order starts as the first-order chain: each origin pair (i, j)
  # takes the first-order logits of its current profile j.
  transition2 <- if (order < 2L) NULL else {
    array(vapply(seq_len(n_types), function(h) {
      vapply(seq_len(n_profiles^2), function(pair) {
        as.vector(transition[, , (pair - 1L) %% n_profiles + 1L, h])
      }, numeric(p * (n_profiles - 1L)))
    }, numeric(p * (n_profiles - 1L) * n_profiles^2)),
    c(p, n_profiles - 1L, n_profiles^2, n_types))
  }
  c(list(measurement = rep(list(block), n_measurement), initial = initial,
         transition = transition),
    if (order >= 2L) list(transition2 = transition2),
    list(group_probabilities = start$group_probabilities,
         stayer = c(rep(FALSE, n_types - mover_stayer), rep(TRUE, mover_stayer))))
}

#' Fit through the general engine from its explicit internal arguments
#' @inheritParams .lta_setup_general
#' @return A fitted multilpa_lta object.
#' @noRd
.lta_fit_general <- function(data, vars, id, time, n_profiles, n_group_classes,
                             variance_model, n_starts, max_iter, tol, min_variance,
                             seed, categorical, min_probability, occasions,
                             transitions, transition_covariates, initial_covariates,
                             measurement_model, select_start, call, order = 1L,
                             mover_stayer = FALSE, missing = "error",
                             covariance_model = "diagonal", structure = NULL,
                             weights = NULL, ordinal = character(),
                             count = character(), count_model = "poisson",
                             count_dispersion = "varying", drop_empty = TRUE,
                             search = "kmeans", n_screen = NULL,
                             acceleration = "squarem", cores = 1L) {
  plan <- do.call(.lta_setup_general, as.list(environment()), quote = TRUE)
  .lta_run_plan(plan, cores)
}

#' Prepare the general transition model
#'
#' Called by `lta()` when transitions vary by occasion or covariate, the
#' initial distribution depends on covariates, or measurement varies by
#' occasion.
#' @return A prepared start plan for the general transition engine.
#' @noRd
.lta_setup_general <- function(data, vars, id, time, n_profiles, n_group_classes,
                             variance_model, n_starts, max_iter, tol, min_variance,
                             seed, categorical, min_probability, occasions,
                             transitions, transition_covariates, initial_covariates,
                             measurement_model, select_start, call, order = 1L,
                             mover_stayer = FALSE, missing = "error",
                             covariance_model = "diagonal", structure = NULL,
                             weights = NULL, ordinal = character(),
                             count = character(), count_model = "poisson",
                             count_dispersion = "varying", drop_empty = TRUE,
                             search = "kmeans", n_screen = NULL,
                             acceleration = "squarem", cores = 1L) {
  n_movers <- as.integer(n_group_classes)
  # A stayer class is one more group class, whose transitions are fixed.
  n_group_classes <- n_movers + as.integer(mover_stayer)
  # The indicator checks (constancy, an observed value per indicator) judge
  # only the rows of positive sampling weight, as in multilpa() (ML-2).
  measurement <- .multilpa_prepare_indicators(data, vars, categorical, missing,
                                              min_probability,
                                              other = c(ordinal, count),
                                              weights = weights)
  # Ordinal categories and the count checks are judged on the rows of
  # positive weight, as in multilpa() (review 2026-10-27, CORE-2).
  extra <- .latents_set_count_model(
    .latents_prepare_extra(data, ordinal, count, missing,
                           contributing = .multilpa_contributing_rows(data, weights)),
    count_model, count_dispersion)
  x <- measurement$x
  codes <- measurement$codes
  n_categories <- measurement$n_categories
  time_values <- .multilpa_time_values(data, time, id, vars)
  groups <- .multilpa_prepare_groups(data[[id]])
  layout <- .multilpa_sequence_layout(groups$index, time_values, groups$n, occasions)
  if (layout$n_occasions < 2L) {
    stop(errorCondition(
      "No group is observed at two occasions, so no transition can be estimated.",
      class = "latents_bad_transition", call = NULL))
  }
  if (order >= 2L && layout$n_occasions < 3L) {
    stop(errorCondition("Second-order transitions require at least three occasions.",
                        class = "latents_bad_transition", call = NULL))
  }
  sampling_weights <- .latents_sampling_weights(data, weights, groups$index, groups$n)
  # Identification is judged on the groups that enter the likelihood: a group
  # of weight zero supplies neither a distinct row nor a group for a class to
  # hold (ML-2, as multilpa()'s `carrying`).
  carrying <- if (is.null(sampling_weights)) rep(TRUE, groups$n) else sampling_weights > 0
  # `NULL`: a millionth of the smallest indicator variance over the rows of
  # positive weight (review 2026-10-26 round 2, K2-2).
  min_variance <- .latents_resolve_min_variance(
    min_variance, x,
    if (is.null(sampling_weights)) NULL else sampling_weights[groups$index])
  distinct_rows <- nrow(unique(cbind(x, codes, extra$ordinal, extra$count)[
    carrying[groups$index], , drop = FALSE]))
  if (n_profiles > distinct_rows || n_group_classes > sum(carrying)) {
    stop(errorCondition(paste(
      "Profiles cannot exceed the distinct observed indicator rows, and group",
      "classes (including stayers) cannot exceed the number of groups, both",
      "counted over the groups of positive weight."),
      class = "latents_unidentified", call = NULL))
  }
  .lta_check_covariate_gaps(layout, groups, time_values, transition_covariates)
  occupied <- .lta_occupied_occasions(layout, x, codes, extra,
                                      if (is.null(sampling_weights)) NULL else
                                        sampling_weights[groups$index])
  if (!any(occupied[-1L])) {
    stop(errorCondition(paste(
      "No occasion after the first has an observed indicator, so no transition",
      "can be estimated."), class = "latents_bad_transition", call = NULL))
  }
  # An empty occasion before the last observed one is not flat (unlike a
  # trailing one): with occasion-specific transitions its moves in and out
  # enter only through their product, so their logits are not separately
  # identified. The fit is still well defined, so this warns, not refuses;
  # the logits of the move into it are not counted (.lta_identified_columns()).
  inner_empty <- which(!occupied & seq_along(occupied) < max(which(occupied)))
  if (identical(transitions, "occasion") && length(inner_empty) > 0L) {
    warning(warningCondition(sprintf(paste(
      "No row observes an indicator at occasion %s, which later occasions",
      "follow. With transitions = \"occasion\", the moves into and out of an",
      "empty occasion enter the likelihood only through their product, so only",
      "the product is identified: the logits of the move into it are not",
      "counted in n_parameters, AIC or BIC and have no standard error, and",
      "those of the move out of it are reported given them."),
      paste(inner_empty, collapse = ", ")),
      class = "latents_unidentified_occasion", call = NULL))
  }
  occasion_of_row <- integer(nrow(x))
  occasion_of_row[layout$slot[!is.na(layout$slot)]] <- col(layout$slot)[!is.na(layout$slot)]
  used <- .lta_used_rows(layout, x, codes, extra, sampling_weights)
  raw_designs <- .lta_designs(layout, .lta_covariate_matrix(data, transition_covariates),
                              .lta_covariate_matrix(data, initial_covariates), transitions,
                              order = order,
                              used = used)
  # EM and the quasi-Newton finish run on centred and scaled covariates, so the
  # search depends neither on a covariate's units nor on its origin (LTA-F5);
  # the coefficients are mapped back to the supplied units at the end.
  chart <- .lta_design_chart(raw_designs, layout, used)
  designs <- .lta_chart_designs(raw_designs, chart)
  n_measurement <- if (identical(measurement_model, "occasion")) layout$n_occasions else 1L
  .latents_local_seed(seed)
  centers <- .multilpa_initialization_center(x, sampling_weights[groups$index])
  x <- sweep(x, 2L, centers, "-")
  categorical_names <- categorical
  if (!is.null(n_categories)) names(n_categories) <- categorical_names
  # What the parameter maps and the scores read, before a fit exists.
  spec <- list(n_profiles = as.integer(n_profiles),
               n_group_classes = as.integer(n_group_classes),
               continuous = measurement$continuous, categorical = categorical_names,
               n_categories = n_categories, variance_model = variance_model,
               measurement = vector("list", n_measurement), designs = designs,
               layout = layout, occasion_of_row = occasion_of_row, x = x,
               codes = codes, group_index = groups$index, n_groups = groups$n,
               group_ids = groups$ids, order = as.integer(order),
               covariance_model = covariance_model, structure = structure,
               stayer = c(rep(FALSE, n_movers), rep(TRUE, as.integer(mover_stayer))),
               sampling_weights = sampling_weights, extra_data = extra)
  # Each start runs to completion or is recorded as failed; the best is
  # chosen by the shared rule (see .latents_run_starts()).
  # Every start is prepared before any is run; EM draws no random numbers, so
  # the draws are those of preparing and running each start in turn.
  screen <- .latents_search_screen(search, n_starts, n_screen, max_iter)
  prepare <- function(start_index) {
    .lta_initialize(x, groups$index, layout, designs, n_profiles,
                    n_group_classes, variance_model, min_variance,
                    start_index, codes, n_categories, min_probability,
                    n_measurement, order, as.integer(mover_stayer),
                    covariance_model, structure, extra,
                    kind = .latents_start_kind(start_index, n_starts),
                    row_weights = sampling_weights[groups$index])
  }
  initial <- .latents_prepare_starts(
    if (is.null(screen)) n_starts else screen$n_pool, prepare)
  context <- as.list(environment())
  context[c("prepare", "initial")] <- NULL
  .lta_plan(context, initial, .lta_run_general_start, .lta_finalise_general,
             c("x", "codes", "layout", "designs", "occasion_of_row",
               "variance_model", "min_variance", "n_categories", "min_probability",
               "max_iter", "tol", "covariance_model", "structure",
               "sampling_weights", "extra", "acceleration", "spec"))
}

#' Run one general transition start, including its numerical finish
#' @param start_index Original start index.
#' @param start Starting values, or NULL when resuming.
#' @param pause,resume Screening budget and checkpoint.
#' @param context Forced EM and numerical-finish inputs from setup.
#' @return One fitted or screened attempt.
#' @noRd
.lta_run_general_start <- function(start_index, start, pause = Inf, resume = NULL, context) {
  with(context, {
    # EM to a loose tolerance, then a quasi-Newton finish on the exact
    # likelihood: EM alone crawls when profiles are weakly separated.
    em <- .lta_em(x, codes, layout, designs, start %||% resume$state, occasion_of_row,
                  variance_model, min_variance, n_categories, min_probability,
                  max_iter, max(tol, 1e-6), covariance_model, structure,
                  sampling_weights, extra, pause = pause, resume = resume,
                  acceleration = acceleration)
    # A screened start is ranked on its EM likelihood; the finish waits for
    # the final stage.
    if (is.finite(pause)) return(em)
    # The quasi-Newton finish works on diagonal measurement; a covariance
    # structure is maximized by EM alone, to the full tolerance.
    if (!identical(covariance_model, "diagonal") || !is.null(structure)) {
      if (tol < 1e-6 && max_iter > 0L) {
        em <- .lta_em(x, codes, layout, designs, em$parameters, occasion_of_row,
                      variance_model, min_variance, n_categories, min_probability,
                      max_iter, tol, covariance_model, structure, sampling_weights,
                      extra, acceleration = acceleration)
      }
      return(em)
    }
    if (max_iter == 0L) return(em)
    finish <- .lta_quasi_newton(spec, em$parameters, min_variance, tol)
    expectation <- .lta_expectation(x, codes, layout, designs, finish$parameters,
                                    occasion_of_row, sampling_weights, extra)
    list(parameters = finish$parameters, expectation = expectation,
         converged = isTRUE(finish$converged), iterations = em$iterations,
         history = c(em$history, expectation$log_likelihood),
         max_gradient = finish$gradient)
  })
}

#' Finalise a general transition fit from its selected starts
#' @param starts The shared start-selection result.
#' @param context Forced setup inputs.
#' @return A fitted multilpa_lta object.
#' @noRd
.lta_finalise_general <- function(starts, context) {
  with(context, {
    found <- starts
    attempts <- starts$attempts
    valid <- starts$valid
    scores <- starts$scores
    converged <- starts$converged
    best_start <- starts$best_start
    best <- attempts[[best_start]]
    # An origin no expected move leaves has flat logits that EM left at their
    # start values: they are fixed at zero (equal probabilities), flagged and
    # not counted (TRANS-1). The likelihood does not move, but the expectation
    # is recomputed so every table reads the parameters reported.
    empty_origins <- .lta_empty_origins(best$expectation, best$parameters$stayer,
                                        n_profiles)
    if (any(empty_origins$first) || any(empty_origins$second)) {
      best$parameters <- .lta_zero_empty_origins(best$parameters, empty_origins)
      best$expectation <- .lta_expectation(x, codes, layout, designs, best$parameters,
                                           occasion_of_row, sampling_weights, extra)
    }
    # The weighted posteriors are counts; each person is classified by their
    # own posterior, while the likelihood stays the weighted one.
    weighted_expectation <- best$expectation
    if (!is.null(sampling_weights)) {
      best$expectation <- .lta_expectation(x, codes, layout, designs, best$parameters,
                                           occasion_of_row, extra = extra)
      best$expectation$log_likelihood <- weighted_expectation$log_likelihood
    }
    # From the fitting chart to the units the covariates arrived in: the same
    # model, so the expectation and the likelihood carry over unchanged.
    best$parameters <- .lta_map_coefficients(best$parameters, chart, "to_raw")
    designs <- raw_designs
    result <- .lta_result(best, attempts, scores, valid, converged, best_start, centers,
                          data, vars, id, time, groups, layout, designs, x, codes,
                          measurement, occasion_of_row, n_profiles, n_group_classes,
                          variance_model, transitions, transition_covariates,
                          initial_covariates, measurement_model, min_variance,
                          min_probability, call, order, mover_stayer, missing,
                          covariance_model, structure,
                          final_stage = .latents_final_stage(found),
                          sampling_weights = sampling_weights, occupied = occupied, used = used,
                          empty_origins = empty_origins)
    result$starts <- .latents_search_columns(result$starts, found, n_starts)
    record <- .latents_search_record(search, n_starts, screen, select_start, "allow")
    result[names(record)] <- record
    result$max_iter <- as.integer(max_iter)
    result$tol <- tol
    result$arguments$search <- search
    result$arguments$n_screen <- n_screen
    result$acceleration <- acceleration
    result$arguments$acceleration <- acceleration
    if (!is.null(sampling_weights)) {
      result$weights <- weights
      result$sampling_weights <- stats::setNames(sampling_weights, groups$ids)
      result$arguments$weights <- weights
    }
    result$ordinal <- ordinal
    result$count <- count
    result$arguments$ordinal <- ordinal
    result$arguments$count <- count
    result$arguments$count_model <- count_model
    result$arguments$count_dispersion <- count_dispersion
    if (!is.null(extra)) {
      result$extra_data <- extra
      profile_names <- paste0("profile_", seq_len(n_profiles))
      result$measurement <- lapply(result$measurement, .latents_label_extra,
                                   extra = extra, profile_names = profile_names)
      result$n_parameters <- result$n_parameters +
        .lta_identified_blocks(n_measurement, occupied) *
        .latents_extra_n_parameters(extra, n_profiles)
      result$aic <- -2 * result$log_likelihood + 2 * result$n_parameters
      # Ordinal and count columns can be the only ones a row observes; a unit
      # of weight zero still observes nothing that counts (ML-1).
      positive <- if (is.null(sampling_weights)) TRUE else
        sampling_weights[groups$index] > 0
      result$informative <- (result$informative |
        rowSums(!is.na(.latents_extra_matrix(extra))) > 0L) & positive
      result$n_informative <- sum(result$informative)
      result$n_informative_groups <- length(unique(groups$index[result$informative]))
      result$bic <- -2 * result$log_likelihood +
        log(result$n_informative_groups) * result$n_parameters
      result$bic_individual <- -2 * result$log_likelihood +
        log(result$n_informative) * result$n_parameters
      .latents_warn_poisson_limit(do.call(cbind, lapply(result$measurement,
                                                        `[[`, "count_dispersion")))
    }
    result$drop_empty <- drop_empty
    result$arguments$drop_empty <- drop_empty
    if (isFALSE(drop_empty)) {
      # The sample sizes of 0.9.37 and earlier: every row and every group,
      # except those of a unit of weight zero, which is in no count under
      # either convention (ML-1; the units `nobs()` and `BIC()` count).
      result$n_informative <- sum(.multilpa_carrying_rows(result, result$n_observations))
      result$bic <- -2 * result$log_likelihood +
        log(.multilpa_group_sample_size(result)) * result$n_parameters
      result$bic_individual <- -2 * result$log_likelihood +
        log(result$n_informative) * result$n_parameters
    }
    # This diagnostic concerns Gaussian measurement, so rows without any
    # Gaussian observation cannot support its spread, regardless of drop_empty.
    # Use the weighted E-step memberships that actually estimated each density;
    # the reported posteriors have already been converted to probabilities.
    gaussian_rows <- rowSums(!is.na(x)) > 0L
    profile_parameters <- .latents_profile_parameters(
      length(measurement$continuous), variance_model, covariance_model, structure)
    invisible(lapply(seq_len(n_measurement), function(index) {
      # The block of an occasion nothing of positive weight observes is flat
      # and not counted (.lta_identified_blocks()), so it is not thin either.
      if (n_measurement > 1L && !occupied[index]) return(NULL)
      rows <- gaussian_rows & (n_measurement == 1L | occasion_of_row == index)
      .latents_warn_thin_profiles(
        colSums(weighted_expectation$subject_posteriors[rows, , drop = FALSE]),
        profile_parameters,
        context = if (n_measurement > 1L)
          sprintf("at occasion %s", index) else NULL)
    }))
    result$gaussian_diagnostics <- do.call(rbind, lapply(seq_len(n_measurement), function(index) {
      rows <- n_measurement == 1L | occasion_of_row == index
      .latents_gaussian_diagnostics(
        result$measurement[[index]], x[rows, , drop = FALSE],
        result$subject_posteriors[rows, , drop = FALSE],
        weights = if (is.null(sampling_weights)) NULL else sampling_weights[groups$index[rows]],
        min_variance = min_variance,
        context = if (n_measurement == 1L) "profiles" else sprintf("occasion %s", index))
    }))
    result
  })
}

#' The design rows the likelihood reads
#'
#' A group's initial design row enters the likelihood when the group has
#' positive sampling weight and observes something; its move into occasion t
#' when, in addition, the move lies inside its span and the group observes
#' something at t or later. Past a group's last observed occasion every
#' transition row sums to one in the forward recursion, so those moves are
#' flat, as are all rows of a unit of weight zero (ML-2, R2-LTA-5).
#'
#' @param layout Sequence layout.
#' @param x,codes Continuous and categorical indicators, one row per
#'   observation.
#' @param extra Prepared ordinal and count indicators, or `NULL`.
#' @param unit_weights Sampling weight of each group, or `NULL`.
#' @return A list with `initial` (one logical per group) and `moves` (a
#'   logical groups-by-occasions matrix).
#' @noRd
.lta_used_rows <- function(layout, x, codes, extra, unit_weights = NULL) {
  observed <- rowSums(!is.na(cbind(x, codes, .latents_extra_matrix(extra)))) > 0L
  slots <- layout$slot
  present <- !is.na(slots)
  seen <- matrix(FALSE, nrow(slots), ncol(slots))
  seen[present] <- observed[slots[present]]
  last <- max.col(seen * rep(seq_len(ncol(slots)), each = nrow(slots)), "last")
  last[rowSums(seen) == 0L] <- 0L
  positive <- if (is.null(unit_weights)) rep(TRUE, nrow(slots)) else unit_weights > 0
  list(initial = positive & last > 0L,
       moves = layout$within & outer(last, seq_len(ncol(slots)), ">=") & positive)
}

#' Which occasions have at least one observed indicator
#'
#' An occasion counts when some group's row at that position of its sequence
#' observes an indicator (continuous, categorical, ordinal or count) and
#' carries positive sampling weight; this is the row rule of `n_informative`.
#'
#' @param layout Sequence layout.
#' @param x,codes Continuous and categorical indicators, one row per
#'   observation.
#' @param extra Prepared ordinal and count indicators, or `NULL`.
#' @param row_weights Sampling weight of each row, or `NULL`.
#' @return A logical vector, one per occasion of the layout.
#' @noRd
.lta_occupied_occasions <- function(layout, x, codes, extra, row_weights = NULL) {
  observed <- rowSums(!is.na(cbind(x, codes, .latents_extra_matrix(extra)))) > 0L
  if (!is.null(row_weights)) observed <- observed & row_weights > 0
  vapply(seq_len(layout$n_occasions), function(t) {
    rows <- layout$slot[, t]
    any(observed[rows[!is.na(rows)]])
  }, logical(1))
}

#' Transition design columns whose logits are identified
#'
#' AIC, BIC and `n_parameters` count only these columns, and Wald inference
#' is computed over them alone (`.lta_identified_coordinates()`).
#'
#' With `transitions = "occasion"`, the intercept `occasion_t` enters only the
#' move into occasion t, and the move into an occasion at which no row
#' observes an indicator is not counted:
#'
#' * past the last observed occasion the forward recursion sums each
#'   transition row to one, so the move into t never enters any group's
#'   likelihood: its intercepts are flat and keep their start values;
#' * an empty occasion followed by an observed one enters only through the
#'   product of the moves into and out of it (a run of r empty occasions
#'   joins r + 1 transition matrices of which only the product, K (K - 1)
#'   logits, is identified), so the moves into the empty occasions are not
#'   counted and the move out of the run carries the product;
#' * when occasion 1 itself is empty, the moves up to the first observed
#'   occasion f join the initial distribution, of which only the distribution
#'   at f is identified, so they are not counted.
#'
#' With `order = 2` the second-order logits govern the moves into occasions 3
#' and later; when nothing is observed from occasion 3 on, every one of them
#' is flat and the whole block (covariate columns included) is not counted.
#' Covariate columns are shared across moves and stay counted otherwise.
#'
#' The two interior rules (an empty occasion followed by an observed one, and
#' an empty first occasion) are applied to first-order chains only. Under
#' second order the state two steps back still enters the move out of the
#' run, so more directions are identified than the product rule says (a
#' numerical Hessian had rank 20 where the rule counted 19); there every
#' logit of those moves stays counted and the count is an upper bound, which
#' the `latents_unidentified_occasion` warning states. Moves past the last
#' observed occasion are flat under either order.
#'
#' @param terms Transition design column names.
#' @param n_occasions Occasions of the layout.
#' @param transitions `"homogeneous"` or `"occasion"`.
#' @param order Transition order.
#' @param occupied Logical, one per occasion, from `.lta_occupied_occasions()`.
#' @return A list with `first` and `second`, the counted column indices of the
#'   first- and second-order blocks.
#' @noRd
.lta_identified_columns <- function(terms, n_occasions, transitions, order, occupied) {
  stopifnot("`occupied` must have one value per occasion" =
              is.logical(occupied) && length(occupied) == n_occasions && any(occupied))
  n_intercepts <- n_occasions - 1L
  # Intercept column j is the move into occasion j + 1.
  destination <- seq_len(n_intercepts) + 1L
  last_observed <- max(which(occupied))
  uncounted <- if (identical(transitions, "occasion") && order < 2L) {
    which(!occupied[destination] | destination <= min(which(occupied)))
  } else if (identical(transitions, "occasion")) {
    which(destination > last_observed)
  } else integer(0)
  second <- if (order >= 2L && max(which(occupied)) < 3L) integer(0) else
    setdiff(.lta_transition_columns(terms, order >= 2L, 2L, n_intercepts), uncounted)
  list(first = setdiff(.lta_transition_columns(terms, order >= 2L, 1L, n_intercepts),
                       uncounted),
       second = second)
}

#' Transition origins that no expected move leaves
#'
#' A profile that no unit of positive weight occupies before a later occasion
#' (in a group class) has no expected move out of it, so the likelihood does
#' not depend on its outgoing logits: they are flat, as the homogeneous
#' engine's uniform rows are (`latents_empty_transition_row`). Under second
#' order the same holds for an origin pair of profiles. Read from the final
#' weighted expectation, so a unit of weight zero informs no origin.
#'
#' @param expectation The weighted expectation of the fit's parameters.
#' @param stayer One logical per group class; a stayer class has no free
#'   transition logits, so none of its origins is reported.
#' @param n_profiles Number of profiles.
#' @return A list with `first` (profiles by group classes) and `second`
#'   (profile pairs by group classes, `NULL` without second order), logical.
#' @noRd
.lta_empty_origins <- function(expectation, stayer, n_profiles) {
  empty <- function(moves, n_origins) {
    vapply(seq_along(expectation$moments), function(h) {
      if (isTRUE(stayer[h])) return(rep(FALSE, n_origins))
      # Arrays are units x origins x destinations, one per move.
      totals <- Reduce(`+`, lapply(expectation$moments[[h]][[moves]], function(pair) {
        apply(pair, 2L, sum)
      }), numeric(n_origins))
      !(totals > 0)
    }, logical(n_origins))
  }
  second <- !is.null(expectation$moments[[1L]]$pairs2)
  list(first = matrix(empty("pairs", n_profiles), n_profiles),
       second = if (second) matrix(empty("pairs2", n_profiles^2), n_profiles^2) else NULL)
}

#' Fix the logits of empty transition origins at zero
#'
#' Their likelihood is flat, so EM leaves them at their start values, which
#' would then be reported as estimates. Zero logits give the destinations of
#' an empty origin equal probabilities, the homogeneous engine's uniform row,
#' whatever the start.
#'
#' @param parameters Parameters on the engine's chart.
#' @param empty The list from `.lta_empty_origins()`.
#' @return The parameters with those logits set to zero.
#' @noRd
.lta_zero_empty_origins <- function(parameters, empty) {
  zero <- function(coefficients, origins) {
    if (is.null(origins) || !any(origins)) return(coefficients)
    # `origins` is origin x group class, the last two dimensions of the
    # coefficients; repeating each cell over the first two marks every logit
    # of an empty origin.
    stopifnot("empty origins must match the transition coefficients" =
                identical(dim(origins), dim(coefficients)[3:4]))
    coefficients[rep(origins, each = prod(dim(coefficients)[1:2]))] <- 0
    coefficients
  }
  parameters$transition <- zero(parameters$transition, empty$first)
  if (!is.null(parameters$transition2)) {
    parameters$transition2 <- zero(parameters$transition2, empty$second)
  }
  parameters
}

#' Measurement blocks the likelihood depends on
#'
#' Under occasion-specific measurement, every row of an occasion with no
#' observed indicator has log density 0 whatever that occasion's means,
#' variances or response probabilities, so its block is flat and not counted.
#'
#' @param n_measurement Measurement blocks (1, or one per occasion).
#' @param occupied Logical, one per occasion.
#' @return The number of counted measurement blocks.
#' @noRd
.lta_identified_blocks <- function(n_measurement, occupied) {
  if (n_measurement == 1L) 1L else sum(occupied)
}

#' Free parameters of the general transition model
#' @noRd
.lta_count_parameters <- function(n_profiles, n_types, d, n_categories,
                                  variance_model, n_measurement, p0, p, order = 1L,
                                  n_stayers = 0L, covariance_model = "diagonal",
                                  structure = NULL, p_second = p) {
  n_movers <- n_types - n_stayers
  # One measurement block as the cross-sectional counter counts it, without
  # its group-class and profile-share terms (one class: K - 1 shares).
  per_measurement <- .multilpa_count_parameters(
    n_profiles, 1L, d, n_categories, variance_model, covariance_model, structure) -
    (n_profiles - 1L)
  as.integer(n_measurement * per_measurement + (n_types - 1L) +
               n_types * p0 * (n_profiles - 1L) +
               n_movers * p * (n_profiles - 1L) * n_profiles +
               if (order >= 2L) n_movers * p_second * (n_profiles - 1L) * n_profiles^2 else 0L)
}

#' Assemble a `multilpa_lta` fit
#' @noRd
.lta_result <- function(best, attempts, scores, valid, converged, best_start,
                        centers, data, vars, id, time, groups, layout, designs, x,
                        codes, measurement, occasion_of_row, n_profiles,
                        n_group_classes, variance_model, transitions,
                        transition_covariates, initial_covariates,
                        measurement_model, min_variance, min_probability, call,
                        order = 1L, mover_stayer = FALSE, missing = "error",
                        covariance_model = "diagonal", structure = NULL,
                        final_stage = rep(TRUE, length(attempts)),
                        sampling_weights = NULL, occupied = NULL, used = NULL,
                        empty_origins = NULL) {
  occupied <- occupied %||% rep(TRUE, layout$n_occasions)
  profile_names <- paste0("profile_", seq_len(n_profiles))
  class_names <- c(paste0("group_class_", seq_len(n_group_classes - mover_stayer)),
                   if (mover_stayer) "stayers")
  parameters <- best$parameters
  parameters$measurement <- lapply(parameters$measurement, function(block) {
    if (!is.null(block$means)) {
      block$means <- sweep(block$means, 2L, centers, "+")
      dimnames(block$means) <- dimnames(block$variances) <-
        list(profile_names, measurement$continuous)
    }
    block
  })
  expectation <- best$expectation
  free <- .lta_identified_columns(colnames(designs$transition[[1L]]), layout$n_occasions,
                                  transitions, order, occupied)
  n_parameters <- .lta_count_parameters(n_profiles, n_group_classes, ncol(x),
                                        measurement$n_categories, variance_model,
                                        .lta_identified_blocks(length(parameters$measurement),
                                                               occupied),
                                        ncol(designs$initial), length(free$first), order,
                                        as.integer(mover_stayer), covariance_model,
                                        structure, length(free$second))
  # The flat logits of an empty origin (.lta_empty_origins()) are not counted:
  # K - 1 destinations on each counted column, per empty origin and class.
  n_parameters <- as.integer(n_parameters -
    (n_profiles - 1L) * (sum(empty_origins$first) * length(free$first) +
                           sum(empty_origins$second) * length(free$second)))
  log_likelihood <- expectation$log_likelihood
  subject_posteriors <- expectation$subject_posteriors
  group_posteriors <- expectation$group_posteriors
  # A row with no observed indicator adds nothing to the likelihood, so it does
  # not enlarge the individual sample size or count as classified.
  observed_per_row <- rowSums(!is.na(x)) +
    if (is.null(codes)) 0L else rowSums(!is.na(codes))
  # A unit of sampling weight zero contributes nothing either, so its rows are
  # not classified and do not count in any sample size (ML-1).
  if (!is.null(sampling_weights)) {
    observed_per_row[sampling_weights[groups$index] <= 0] <- 0L
  }
  n_informative <- sum(observed_per_row > 0L)
  n_informative_groups <- length(unique(groups$index[observed_per_row > 0L]))
  dimnames(subject_posteriors) <- list(NULL, profile_names)
  dimnames(group_posteriors) <- list(groups$ids, class_names)
  names(parameters$group_probabilities) <- class_names
  dimnames(parameters$initial) <- list(colnames(designs$initial),
                                       profile_names[-n_profiles], class_names)
  # Dimension 2 indexes the destinations other than the origin, in profile
  # order; the transition table names them.
  dimnames(parameters$transition) <- list(colnames(designs$transition[[1L]]),
                                          paste0("destination_", seq_len(n_profiles - 1L)),
                                          profile_names, class_names)
  if (!is.null(parameters$transition2)) {
    pairs <- expand.grid(current = profile_names, previous = profile_names,
                         stringsAsFactors = FALSE)
    dimnames(parameters$transition2) <- list(
      colnames(designs$transition[[1L]]),
      paste0("destination_", seq_len(n_profiles - 1L)),
      paste0(pairs$previous, "->", pairs$current), class_names)
  }
  probability_boundary <- .lta_probability_boundary(expectation, layout,
                                                    parameters$stayer, used = used)
  boundary <- probability_boundary || any(vapply(parameters$measurement, function(block) {
    .multilpa_covariance_boundary(block, min_variance)
  }, logical(1)))
  result <- list(
    call = call, vars = vars, continuous = measurement$continuous,
    categorical = setdiff(vars, c(measurement$continuous, measurement$other)),
    categorical_levels = measurement$encoded$levels,
    id = id, time = time, group_ids = groups$ids, group_index = groups$index,
    # The identifiers as supplied (typed), for the assignments' data columns.
    group_values = groups$values,
    time_values = .multilpa_time_values(data, time, id, vars),
    n_profiles = as.integer(n_profiles), n_group_classes = as.integer(n_group_classes),
    n_groups = groups$n, n_observations = nrow(x), n_occasions = layout$n_occasions,
    transitions = transitions, transition_covariates = transition_covariates,
    initial_covariates = initial_covariates, measurement_model = measurement_model,
    variance_model = variance_model,
    measurement = parameters$measurement, initial_coefficients = parameters$initial,
    transition_coefficients = parameters$transition,
    second_order_coefficients = parameters$transition2, order = as.integer(order),
    missing = missing, covariance_model = covariance_model,
    covariance_structure = structure,
    stayer = stats::setNames(parameters$stayer, class_names),
    mover_stayer = isTRUE(mover_stayer),
    group_probabilities = parameters$group_probabilities,
    designs = designs, layout = layout, occasion_of_row = occasion_of_row,
    centers = centers, x = x, codes = codes,
    n_categories = if (is.null(measurement$n_categories)) NULL else
      stats::setNames(measurement$n_categories, setdiff(vars, c(measurement$continuous, measurement$other))),
    subject_posteriors = subject_posteriors, group_posteriors = group_posteriors,
    expectation = expectation,
    log_likelihood = log_likelihood, n_parameters = n_parameters,
    aic = -2 * log_likelihood + 2 * n_parameters,
    bic = -2 * log_likelihood + log(n_informative_groups) * n_parameters,
    bic_individual = -2 * log_likelihood + log(n_informative) * n_parameters,
    n_informative = n_informative, informative = observed_per_row > 0L,
    n_informative_groups = n_informative_groups,
    converged = isTRUE(best$converged), iterations = best$iterations,
    log_likelihood_history = best$history, boundary = boundary,
    probability_boundary = probability_boundary,
    # Which occasions observe anything: they decide which logits and
    # measurement blocks are identified (.lta_identified_coordinates()).
    occupied = occupied,
    # Origins (profiles; pairs under second order) by group class that no
    # expected move leaves: their logits are flat (.lta_empty_origins()).
    empty_transition_rows = matrix(empty_origins$first %||% FALSE, n_profiles,
                                   n_group_classes,
                                   dimnames = list(profile_names, class_names)),
    empty_second_order_rows = if (is.null(empty_origins$second)) NULL else
      matrix(empty_origins$second, n_profiles^2, n_group_classes,
             dimnames = list(NULL, class_names)),
    starts = data.frame(start = seq_along(attempts),
                        log_likelihood = ifelse(valid, scores, NA_real_),
                        converged = converged,
                        iterations = vapply(seq_along(attempts), function(i) {
                          if (valid[i]) as.integer(attempts[[i]]$iterations %||% NA_integer_)
                          else NA_integer_
                        }, integer(1)),
                        error = vapply(attempts, function(a) a$error %||% NA_character_,
                                       character(1)),
                        selected = seq_along(attempts) == best_start),
    n_best_replicated = sum(valid & final_stage & abs(scores - log_likelihood) <=
                              1e-6 * (1 + abs(log_likelihood))),
    min_variance = min_variance, min_probability = min_probability,
    # The specification, so the fit can be refitted (bootstrap) or nested.
    arguments = list(vars = vars, id = id, time = time, n_profiles = as.integer(n_profiles),
                     n_group_classes = as.integer(n_group_classes - mover_stayer),
                     variance_model = variance_model,
                     categorical = setdiff(vars, c(measurement$continuous, measurement$other)),
                     occasions = layout$occasions, transitions = transitions,
                     transition_covariates = transition_covariates,
                     initial_covariates = initial_covariates,
                     measurement = measurement_model, order = as.integer(order),
                     mover_stayer = isTRUE(mover_stayer), missing = missing,
                     covariance_model = covariance_model,
                     model = structure, min_variance = min_variance,
                     min_probability = min_probability))
  class(result) <- "multilpa_lta"
  if (!result$converged) {
    warning(warningCondition("The best start did not converge; increase max_iter.",
                             class = "latents_unconverged", call = NULL))
  }
  if (boundary) {
    warning(warningCondition(if (probability_boundary) paste(
      "A transition or initial probability is estimated at zero (a move or",
      "start that never occurs); this is a boundary fit and its standard",
      "errors are withheld.") else
      paste("A variance or covariance eigenvalue reached min_variance;",
            "this is a bound-active fit."),
      class = "latents_boundary", call = NULL))
  }
  if (any(result$empty_transition_rows) || any(result$empty_second_order_rows)) {
    warning(warningCondition(paste(
      "A profile is never occupied before a later occasion, so no move leaves",
      "it and its outgoing logits are not identified. They are fixed at zero",
      "(equal probabilities), are not counted in n_parameters, AIC or BIC, and",
      "are flagged estimated = FALSE in get_results(fit, \"transitions\") and",
      "parameter_inference(fit)."),
      class = "latents_empty_transition_row", call = NULL))
  }
  # As in the homogeneous engine: a class that has emptied (a stayer class no
  # one belongs to, say) is reported, not left to look like an estimate.
  if (any(colSums(subject_posteriors) < 1) || any(colSums(group_posteriors) < 1)) {
    warning(warningCondition(
      "A profile or group class has effective membership below one.",
      class = "latents_small_classes", call = NULL))
  }
  result
}
