#' Refuse transition specifications with a known short-panel likelihood ridge
#'
#' These are structural results for models without membership covariates.
#' Two occasions identify a latent joint distribution, not a mixture of
#' transition chains or a separate stayer population. Three occasions do not
#' identify multiple occasion-specific chains, or a mixture of second-order
#' chains. Covariate-dependent identification is not decided by this check.
#'
#' @param n_occasions Number of occasions in the fitted layout.
#' @param n_movers Number of freely moving group classes.
#' @param transitions Homogeneous or occasion-specific transitions.
#' @param order First- or second-order dependence.
#' @param mover_stayer Whether a separate stayer class is fitted.
#' @param has_covariates Whether initial or transition covariates are fitted.
#' @param n_profiles Optional state count for the necessary latent-path
#'   dimension check.
#' @return NULL invisibly, or a latents_unidentified error.
#' @noRd
.lta_check_short_panel <- function(n_occasions, n_movers, transitions, order = 1L,
                                    mover_stayer = FALSE, has_covariates = FALSE,
                                    n_profiles = NULL) {
  stopifnot(length(n_occasions) == 1L, n_occasions >= 2L,
            length(n_movers) == 1L, n_movers >= 1L,
            transitions %in% c("homogeneous", "occasion"), order %in% c(1L, 2L),
            is.logical(mover_stayer), length(mover_stayer) == 1L,
            is.logical(has_covariates), length(has_covariates) == 1L)
  if (has_covariates) return(invisible(NULL))
  classes <- n_movers + as.integer(mover_stayer)
  dimension_problem <- FALSE
  if (!is.null(n_profiles)) {
    stopifnot(length(n_profiles) == 1L, n_profiles >= 2L)
    first_columns <- if (transitions == "homogeneous" || order >= 2L) 1L else n_occasions - 1L
    second_columns <- if (transitions == "homogeneous") 1L else n_occasions - 2L
    chain_parameters <- classes - 1L + classes * (n_profiles - 1L) +
      n_movers * n_profiles * (n_profiles - 1L) * first_columns +
      if (order >= 2L) n_movers * n_profiles^2 * (n_profiles - 1L) * second_columns else 0L
    # The observed likelihood depends on these coefficients only through the
    # K^T latent path masses. Logarithms avoid overflowing that upper bound.
    dimension_problem <- log1p(chain_parameters) > n_occasions * log(n_profiles) + 1e-12
  }
  problem <- if (n_occasions == 2L && classes > 1L) {
    "Two occasions cannot identify multiple group classes or a separate stayer class."
  } else if (n_occasions == 3L && classes > 1L && order >= 2L) {
    "Three occasions cannot identify multiple group classes in a second-order chain."
  } else if (n_occasions == 3L && n_movers > 1L && transitions == "occasion") {
    "Three occasions cannot identify multiple freely moving group classes with occasion-specific transitions."
  } else if (dimension_problem) {
    "The initial, transition and group-class coordinates exceed the dimension of the latent path distribution."
  } else NULL
  if (!is.null(problem)) {
    stop(errorCondition(paste(problem,
      "Use one freely moving group class, a supported simpler transition specification, or a longer panel."),
      class = "latents_unidentified", call = NULL))
  }
  invisible(NULL)
}

#' Weighted expected outgoing moves by origin and group class
#' @param expectation The final weighted E-step result.
#' @param stayer Logical indicator of structural stayer classes.
#' @param n_profiles Number of profiles.
#' @param sampling_weights Optional group weights when the supplied moments
#'   are unweighted, as in an older saved fit's reporting expectation.
#' @return Counts aggregated over applicable moves and by occasion, separately
#'   for first-order origins and second-order origin pairs. Stayers have Inf
#'   internally so that their fixed rows are never flagged as weak estimates.
#' @noRd
.lta_origin_support <- function(expectation, stayer, n_profiles, sampling_weights = NULL) {
  stopifnot(is.list(expectation$moments), length(n_profiles) == 1L, n_profiles >= 2L)
  n_classes <- length(expectation$moments)
  stayer <- stayer %||% rep(FALSE, n_classes)
  stopifnot(is.logical(stayer), length(stayer) == n_classes,
            is.null(sampling_weights) || is.numeric(sampling_weights))
  second <- !is.null(expectation$moments[[1L]]$pairs2)
  by_move <- function(block, n_origins) {
    count <- length(expectation$moments[[1L]][[block]])
    if (identical(block, "pairs") && second) count <- min(count, 1L)
    lapply(seq_len(count), function(t) {
      matrix(vapply(seq_len(n_classes), function(h) {
        if (stayer[h]) return(rep(Inf, n_origins))
        pairs <- expectation$moments[[h]][[block]][[t]]
        if (!is.null(sampling_weights)) {
          stopifnot(length(sampling_weights) == dim(pairs)[1L])
          pairs <- pairs * sampling_weights
        }
        apply(pairs, 2L, sum)
      }, numeric(n_origins)), n_origins, n_classes)
    })
  }
  first <- by_move("pairs", n_profiles)
  second_moves <- if (second) by_move("pairs2", n_profiles^2) else NULL
  list(first = Reduce(`+`, first, matrix(0, n_profiles, n_classes)),
       second = if (second) Reduce(`+`, second_moves, matrix(0, n_profiles^2, n_classes)),
       first_by_occasion = first, second_by_occasion = second_moves)
}

#' Origins with fewer than one expected outgoing transition
#'
#' This is a weak-support diagnostic, not a declaration that a nonzero row is
#' structurally unidentified. Its formally free coordinates remain counted.
#' @param support Counts from .lta_origin_support().
#' @param transitions Homogeneous or occasion-specific transitions.
#' @return Logical matrices for first- and second-order origins.
#' @noRd
.lta_weak_origins <- function(support, transitions) {
  stopifnot(is.list(support), transitions %in% c("homogeneous", "occasion"))
  weak <- function(counts) counts > 0 & counts < 1
  select <- function(total, by_occasion) {
    if (is.null(total)) return(NULL)
    if (transitions == "occasion") {
      # A row can be supported in other occasions yet have no outgoing mass
      # at this one. Entirely empty origins use the separate exact-zero rule.
      scarce <- Reduce(`|`, lapply(by_occasion, function(counts) counts < 1),
                        matrix(FALSE, nrow(total), ncol(total)))
      return(scarce & total > 0 & is.finite(total))
    }
    weak(total)
  }
  list(first = select(support$first, support$first_by_occasion),
       second = select(support$second, support$second_by_occasion))
}

#' Warn when reported transition origins have negligible expected support
#' @param first,second Logical weak-support origin matrices.
#' @return NULL invisibly, with a classed warning when any origin is weak.
#' @noRd
.lta_warn_weak_origins <- function(first, second = NULL) {
  stopifnot(is.logical(first), is.null(second) || is.logical(second))
  if (any(first) || any(second)) {
    warning(warningCondition(paste(
      "Some transition origins have fewer than one expected outgoing move.",
      "Their rows have weak support; inspect origin_expected_count and supported",
      "in the transition tables. Nonzero rows remain in the likelihood and its",
      "formal parameter count, but Wald inference is withheld."),
      class = "latents_weak_transition_support", call = NULL))
  }
  invisible(NULL)
}

#' Preserve an entire logit design term, including trailing dots
#' @param names Character coefficient names in block.owner.outcome.term form.
#' @return Character design terms.
#' @noRd
.lta_coefficient_term <- function(names) {
  stopifnot(is.character(names))
  sub("^[^.]*\\.[^.]*\\.[^.]*\\.", "", names)
}

#' Split coefficient names without discarding a trailing empty component
#' @param names Character coefficient names.
#' @return A list of character name components.
#' @noRd
.lta_name_pieces <- function(names) {
  stopifnot(is.character(names))
  lapply(strsplit(paste0(names, "\r"), ".", fixed = TRUE), function(parts) {
    parts[length(parts)] <- sub("\r$", "", parts[length(parts)])
    parts
  })
}

#' Check transition logits in a unit-free design chart
#' @param coefficients Logit coefficients, design columns by destinations.
#' @param design The contributing initial or transition design rows.
#' @param limit Maximum absolute unit-free coefficient.
#' @return Whether a coefficient is numerically extreme.
#' @noRd
.lta_extreme_logits <- function(coefficients, design, limit = 20) {
  stopifnot(is.matrix(coefficients), is.matrix(design), nrow(coefficients) == ncol(design))
  if (nrow(design) == 0L) return(FALSE)
  basis <- .multilpa_design_basis(design)
  scaled <- if (is.null(basis)) coefficients * .multilpa_design_scale(design) else
    basis %*% coefficients
  any(!is.finite(scaled)) || any(abs(scaled) > limit)
}

#' Weak-origin diagnostics on a fitted general transition model
#' @param fit A general transition fit.
#' @return A list of weak first- and second-order origin matrices, recomputed
#'   for older saved fits that do not carry the diagnostic.
#' @noRd
.lta_fit_weak_origins <- function(fit) {
  stopifnot(inherits(fit, "multilpa_lta"))
  if (!is.null(fit$weak_transition_rows)) {
    return(list(first = fit$weak_transition_rows, second = fit$weak_second_order_rows))
  }
  support <- .lta_origin_support(fit$expectation, fit$stayer, fit$n_profiles,
                                 unname(fit$sampling_weights))
  .lta_weak_origins(support, fit$transitions %||% "homogeneous")
}

#' Support counts displayed beside general transition probabilities
#' @param fit A general transition fit.
#' @param kind First- or second-order transition block.
#' @param class_index Group-class index.
#' @param occasion Move index for an occasion-specific first-order table.
#' @return A list with expected origin counts and supported-row flags.
#' @noRd
.lta_table_support <- function(fit, kind, class_index, occasion = NULL) {
  stopifnot(inherits(fit, "multilpa_lta"), kind %in% c("first", "second"))
  support <- fit$origin_support %||%
    .lta_origin_support(fit$expectation, fit$stayer, fit$n_profiles,
                        unname(fit$sampling_weights))
  counts <- if (is.null(occasion)) support[[kind]][, class_index] else
    support[[paste0(kind, "_by_occasion")]][[occasion]][, class_index]
  weak <- if (is.null(occasion)) .lta_fit_weak_origins(fit)[[kind]][, class_index] else
    counts > 0 & counts < 1
  if (isTRUE(fit$stayer[class_index])) {
    return(list(count = rep(NA_real_, length(counts)), supported = rep(TRUE, length(counts))))
  }
  list(count = counts, supported = counts >= 1 & !weak)
}
