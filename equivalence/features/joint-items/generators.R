# The generating design of the joint-items kit, shared by make-kit.R (which
# writes the kit) and the reference study in the tests (which draws many more
# cheap datasets the same way). Every dataset comes from simulab; see make-kit.R.

# -- generating parameters -------------------------------------------------------

#' P(category 1) of binary items: one row per class, one column per item
binary_profile <- function(k, m) {
  rows <- list(c(0.85, 0.80, 0.75, 0.80, 0.85, 0.75),
               c(0.20, 0.25, 0.30, 0.15, 0.20, 0.25),
               c(0.80, 0.75, 0.20, 0.25, 0.80, 0.20),
               c(0.25, 0.20, 0.80, 0.85, 0.25, 0.80))
  do.call(rbind, lapply(rows[seq_len(k)], \(row) row[seq_len(m)]))
}

#' A class-specific conditional table P(second | first): rho on a diagonal
#' association plus (1 - rho) times a class-specific base distribution
association <- function(rho, base, c1) {
  k <- length(rho)
  c2 <- ncol(base)
  cond <- array(0, dim = c(k, c1, c2))
  invisible(lapply(seq_len(k), \(j) lapply(seq_len(c1), \(l) {
    onehot <- as.numeric(seq_len(c2) == ((l - 1L) %% c2) + 1L)
    cond[j, l, ] <<- rho[j] * onehot + (1 - rho[j]) * base[j, ]
  })))
  cond
}

#' Draw one case
#'
#' @param seed The simulab seed; the redraw of the second items uses seed + 7919.
#' @param n Rows.
#' @param proportions Class proportions.
#' @param drawn Named list of k x C probability matrices, one per item simulab
#'   draws (first items of pairs and the unpaired items).
#' @param pairs List of pairs, each list(first, second, cond = k x C1 x C2).
#' @param order Column order of the items.
draw_case <- function(seed, n, proportions, drawn, pairs, order) {
  k <- length(proportions)
  widest <- max(vapply(drawn, ncol, integer(1)))
  probabilities <- array(0, dim = c(k, length(drawn), widest))
  invisible(lapply(seq_along(drawn), \(i) {
    probabilities[, i, seq_len(ncol(drawn[[i]]))] <<- drawn[[i]]
  }))
  simulated <- as.data.frame(simulab::simulate_lca(
    n = n, probabilities = probabilities, proportions = proportions, seed = seed))
  membership <- as.integer(sub("Class ", "", simulated$latent_class))
  items <- stats::setNames(lapply(seq_along(drawn), \(i) simulated[[paste0("item_", i)]]),
                           names(drawn))
  set.seed(seed + 7919L)
  joint <- lapply(pairs, function(pair) {
    first <- items[[pair$first]]
    probs <- t(vapply(seq_len(n), \(i) pair$cond[membership[i], first[i], ],
                      numeric(dim(pair$cond)[3L])))
    u <- stats::runif(n)
    items[[pair$second]] <<- 1L + as.integer(rowSums(u > t(apply(probs, 1L, cumsum))))
    marginal <- drawn[[pair$first]]
    table <- array(0, dim = dim(pair$cond))
    invisible(lapply(seq_len(k), \(j) table[j, , ] <<- marginal[j, ] * pair$cond[j, , ]))
    table
  })
  names(joint) <- vapply(pairs, \(pair) paste(pair$first, pair$second, sep = ":"), "")
  list(items = as.data.frame(items[order]), membership = membership,
       truth = list(proportions = proportions, drawn = drawn,
                    conditional = stats::setNames(lapply(pairs, `[[`, "cond"), names(joint)),
                    joint = joint))
}

two_by_two <- function(rho = c(0.6, 0.5)) {
  association(rho, rbind(c(0.7, 0.3), c(0.3, 0.7), c(0.5, 0.5), c(0.6, 0.4))[seq_along(rho), ,
                                                                            drop = FALSE], 2L)
}
binary <- function(p) cbind(p, 1 - p)

generate <- list(
  bin_small = function(seed, n = 80) {
    drawn <- c(list(p1 = binary(c(0.75, 0.25))),
               stats::setNames(lapply(1:3, \(i) binary(binary_profile(2, 3)[, i])),
                               paste0("u", 1:3)))
    draw_case(seed, n, c(0.55, 0.45), drawn,
              list(list(first = "p1", second = "p2", cond = two_by_two())),
              c("p1", "p2", "u1", "u2", "u3"))
  },
  bin_mod = function(seed) {
    generate$bin_small(seed, n = 800)
  },
  cat3_mod = function(seed, n = 800, others = 4L) {
    k <- 3L
    drawn <- c(list(p1 = rbind(c(0.6, 0.3, 0.1), c(0.15, 0.35, 0.5), c(0.3, 0.4, 0.3))),
               stats::setNames(lapply(seq_len(others), \(i) binary(binary_profile(k, others)[, i])),
                               paste0("u", seq_len(others))))
    cond <- association(c(0.6, 0.4, 0.5),
                        rbind(c(0.5, 0.3, 0.2), c(0.2, 0.3, 0.5), c(0.33, 0.34, 0.33)), 3L)
    draw_case(seed, n, c(0.4, 0.35, 0.25), drawn,
              list(list(first = "p1", second = "p2", cond = cond)),
              c("p1", "p2", paste0("u", seq_len(others))))
  },
  bin4_mod = function(seed, n = 800) {
    drawn <- c(list(p1 = binary(c(0.7, 0.3))),
               stats::setNames(lapply(1:4, \(i) binary(binary_profile(2, 4)[, i])),
                               paste0("u", 1:4)))
    cond <- array(0, dim = c(2, 2, 4))
    cond[1, 1, ] <- c(0.5, 0.3, 0.15, 0.05); cond[1, 2, ] <- c(0.05, 0.15, 0.3, 0.5)
    cond[2, 1, ] <- c(0.4, 0.4, 0.1, 0.1); cond[2, 2, ] <- c(0.1, 0.2, 0.3, 0.4)
    draw_case(seed, n, c(0.6, 0.4), drawn,
              list(list(first = "p1", second = "p2", cond = cond)),
              c("p1", "p2", paste0("u", 1:4)))
  },
  twopair_mod = function(seed, n = 800) {
    k <- 3L
    drawn <- list(p1 = binary(c(0.8, 0.2, 0.5)), q1 = binary(c(0.3, 0.75, 0.8)),
                  u1 = binary(binary_profile(k, 2)[, 1]), u2 = binary(binary_profile(k, 2)[, 2]))
    draw_case(seed, n, c(0.4, 0.35, 0.25), drawn,
              list(list(first = "p1", second = "p2", cond = two_by_two(c(0.6, 0.4, 0.5))),
                   list(first = "q1", second = "q2", cond = two_by_two(c(0.3, 0.6, 0.5)))),
              c("p1", "p2", "q1", "q2", "u1", "u2"))
  },
  k4_large = function(seed) {
    k <- 4L
    drawn <- c(list(p1 = binary(c(0.8, 0.2, 0.7, 0.3))),
               stats::setNames(lapply(1:5, \(i) binary(binary_profile(k, 5)[, i])),
                               paste0("u", 1:5)))
    draw_case(seed, 5000, c(0.3, 0.3, 0.2, 0.2), drawn,
              list(list(first = "p1", second = "p2", cond = two_by_two(c(0.6, 0.5, 0.4, 0.3)))),
              c("p1", "p2", paste0("u", 1:5)))
  },
  huge = function(seed) {
    generate$cat3_mod(seed, n = 20000, others = 6L)
  },
  sparse = function(seed, n = 600) {
    drawn <- c(list(p1 = rbind(c(0.5, 0.3, 0.2), c(0.2, 0.3, 0.5))),
               stats::setNames(lapply(1:3, \(i) binary(binary_profile(2, 3)[, i])),
                               paste0("u", 1:3)))
    # Structural zeros: (1, 3) and (3, 1) never occur in either class.
    cond <- array(0, dim = c(2, 3, 3))
    cond[1, 1, ] <- c(0.7, 0.3, 0); cond[1, 2, ] <- c(0.3, 0.5, 0.2); cond[1, 3, ] <- c(0, 0.4, 0.6)
    cond[2, 1, ] <- c(0.4, 0.6, 0); cond[2, 2, ] <- c(0.2, 0.3, 0.5); cond[2, 3, ] <- c(0, 0.2, 0.8)
    draw_case(seed, n, c(0.5, 0.5), drawn,
              list(list(first = "p1", second = "p2", cond = cond)),
              c("p1", "p2", "u1", "u2", "u3"))
  },
  perfect = function(seed) {
    drawn <- c(list(p1 = rbind(c(0.6, 0.3, 0.1), c(0.1, 0.3, 0.6))),
               stats::setNames(lapply(1:3, \(i) binary(binary_profile(2, 3)[, i])),
                               paste0("u", 1:3)))
    cond <- array(0, dim = c(2, 3, 3))
    invisible(lapply(1:2, \(j) cond[j, , ] <<- diag(3)))
    draw_case(seed, 600, c(0.5, 0.5), drawn,
              list(list(first = "p1", second = "p2", cond = cond)),
              c("p1", "p2", "u1", "u2", "u3"))
  },
  boundary = function(seed) {
    drawn <- c(list(p1 = binary(c(0.75, 0.3))),
               stats::setNames(lapply(1:3, \(i) binary(binary_profile(2, 3)[, i])),
                               paste0("u", 1:3)))
    # The second item's category 2 never occurs in class 1.
    cond <- array(0, dim = c(2, 2, 3))
    cond[1, 1, ] <- c(0.8, 0, 0.2); cond[1, 2, ] <- c(0.3, 0, 0.7)
    cond[2, 1, ] <- c(0.3, 0.4, 0.3); cond[2, 2, ] <- c(0.1, 0.3, 0.6)
    draw_case(seed, 700, c(0.5, 0.5), drawn,
              list(list(first = "p1", second = "p2", cond = cond)),
              c("p1", "p2", "u1", "u2", "u3"))
  },
  partial = function(seed, n = 800) {
    drawn <- c(list(p1 = binary(c(0.75, 0.3))),
               stats::setNames(lapply(1:3, \(i) binary(binary_profile(2, 3)[, i])),
                               paste0("u", 1:3)))
    cond <- array(0, dim = c(2, 2, 3))
    cond[1, 1, ] <- c(0.6, 0.3, 0.1); cond[1, 2, ] <- c(0.1, 0.2, 0.7)
    cond[2, 1, ] <- c(0.5, 0.3, 0.2); cond[2, 2, ] <- c(0.2, 0.2, 0.6)
    case <- draw_case(seed, n, c(0.5, 0.5), drawn,
                      list(list(first = "p1", second = "p2", cond = cond)),
                      c("p1", "p2", "u1", "u2", "u3"))
    set.seed(seed + 31L)
    n <- nrow(case$items)
    holes <- sample(n, round(0.18 * n))
    case$items$p2[holes[seq_len(round(0.08 * n))]] <- NA
    case$items$p1[holes[round(0.08 * n) + seq_len(round(0.05 * n))]] <- NA
    case$items$u1[holes[round(0.13 * n) + seq_len(round(0.05 * n))]] <- NA
    case
  },
  indep = function(seed, n = 1000) {
    drawn <- c(list(p1 = binary(c(0.7, 0.3))),
               stats::setNames(lapply(1:3, \(i) binary(binary_profile(2, 3)[, i])),
                               paste0("u", 1:3)))
    draw_case(seed, n, c(0.5, 0.5), drawn,
              list(list(first = "p1", second = "p2", cond = two_by_two(c(0, 0)))),
              c("p1", "p2", "u1", "u2", "u3"))
  },
  strings = function(seed) {
    drawn <- c(list(p1 = rbind(c(0.6, 0.3, 0.1), c(0.1, 0.3, 0.6))),
               stats::setNames(lapply(1:3, \(i) binary(binary_profile(2, 3)[, i])),
                               paste0("u", 1:3)))
    cond <- association(c(0.5, 0.4), rbind(c(0.7, 0.3), c(0.3, 0.7)), 3L)
    case <- draw_case(seed, 600, c(0.5, 0.5), drawn,
                      list(list(first = "p1", second = "p2", cond = cond)),
                      c("p1", "p2", "u1", "u2", "u3"))
    case$labels <- list(p1 = c("low", "mid", "high"), p2 = c("no", "yes"),
                        u1 = c("agree", "disagree"), u2 = c("agree", "disagree"),
                        u3 = c("agree", "disagree"))
    case
  },
  numeric = function(seed) {
    drawn <- c(list(p1 = binary(c(0.7, 0.25))),
               stats::setNames(lapply(1:3, \(i) binary(binary_profile(2, 3)[, i])),
                               paste0("u", 1:3)))
    cond <- array(0, dim = c(2, 2, 3))
    cond[1, 1, ] <- c(0.6, 0.3, 0.1); cond[1, 2, ] <- c(0.1, 0.3, 0.6)
    cond[2, 1, ] <- c(0.3, 0.3, 0.4); cond[2, 2, ] <- c(0.5, 0.3, 0.2)
    case <- draw_case(seed, 600, c(0.5, 0.5), drawn,
                      list(list(first = "p1", second = "p2", cond = cond)),
                      c("p1", "p2", "u1", "u2", "u3"))
    # Codes whose concatenation without a separator would collide:
    # 1 & 10 vs 11 & 0, and 1 & 1 vs 11 & "".
    case$labels <- list(p1 = c(1, 11), p2 = c(0, 10, 1))
    case
  },
  weights = function(seed) {
    case <- generate$bin_small(seed, 400)
    set.seed(seed + 53L)
    case$weights <- sample(1:3, nrow(case$items), replace = TRUE)
    case
  },
  multilevel = function(seed) {
    groups <- 60L
    size <- 15L
    group_class <- rep(c(1L, 2L), c(30L, 30L))
    drawn <- c(list(p1 = binary(c(0.75, 0.25))),
               stats::setNames(lapply(1:3, \(i) binary(binary_profile(2, 3)[, i])),
                               paste0("u", 1:3)))
    pairs <- list(list(first = "p1", second = "p2", cond = two_by_two()))
    order <- c("p1", "p2", "u1", "u2", "u3")
    parts <- lapply(1:2, function(h) {
      draw_case(seed + h, sum(group_class == h) * size, list(c(0.8, 0.2), c(0.25, 0.75))[[h]],
                drawn, pairs, order)
    })
    case <- list(items = rbind(parts[[1]]$items, parts[[2]]$items),
                 membership = c(parts[[1]]$membership, parts[[2]]$membership),
                 truth = c(parts[[1]]$truth[c("drawn", "conditional", "joint")],
                           list(group_proportions = rbind(c(0.8, 0.2), c(0.25, 0.75)))))
    case$group <- rep(seq_len(groups), each = size)
    case
  }
)

specification <- list(
  bin_small = list(k = 2L, replicates = 5L, mplus = TRUE),
  bin_mod = list(k = 2L, replicates = 5L, mplus = TRUE),
  cat3_mod = list(k = 3L, replicates = 5L, mplus = TRUE),
  bin4_mod = list(k = 2L, replicates = 5L, mplus = TRUE),
  twopair_mod = list(k = 3L, replicates = 5L, mplus = TRUE, pairs = c("p1:p2", "q1:q2")),
  k4_large = list(k = 4L, replicates = 5L, mplus = TRUE),
  huge = list(k = 3L, replicates = 1L, mplus = FALSE),
  sparse = list(k = 2L, replicates = 4L, mplus = TRUE),
  perfect = list(k = 2L, replicates = 3L, mplus = TRUE),
  boundary = list(k = 2L, replicates = 2L, mplus = TRUE),
  partial = list(k = 2L, replicates = 3L, mplus = TRUE, missing = TRUE),
  indep = list(k = 2L, replicates = 3L, mplus = TRUE),
  strings = list(k = 2L, replicates = 2L, mplus = TRUE),
  numeric = list(k = 2L, replicates = 2L, mplus = TRUE),
  weights = list(k = 2L, replicates = 2L, mplus = TRUE),
  multilevel = list(k = 2L, h = 2L, replicates = 2L, mplus = FALSE)
)

