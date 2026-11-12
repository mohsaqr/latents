# Datasets of the argument sweep: bundled data, classic public data and
# simulated data with a known structure. Every simulation has its own seed.

#' Run code under a seed without disturbing the caller's stream
sweep_seeded <- function(seed, code) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE, after = FALSE)
  set.seed(seed, kind = "Mersenne-Twister")
  code()
}

#' Continuous datasets
#' @return Named list of data.frames holding only numeric indicators.
sweep_continuous <- function() {
  blobs <- sweep_seeded(101L, function() {
    centre <- rbind(c(0, 0, 0), c(4, 1, -3), c(-2, 5, 3))
    profile <- rep(1:3, c(100L, 80L, 60L))
    out <- as.data.frame(centre[profile, ] + matrix(stats::rnorm(240L * 3L), 240L) %*%
                           chol(matrix(c(1, 0.4, 0.2, 0.4, 1.5, 0.3, 0.2, 0.3, 0.8), 3L)))
    stats::setNames(out, c("a", "b", "c"))
  })
  iris <- stats::setNames(datasets::iris[1:4], c("sepal_l", "sepal_w", "petal_l", "petal_w"))
  list(faithful = datasets::faithful, iris = iris, srl = latents::srl, blobs = blobs)
}

#' Punch MCAR holes into a data.frame
#' @param data data.frame.
#' @param share Share of cells set to `NA`.
#' @param seed Seed.
#' @param columns Columns that may lose values.
#' @return The data with missing cells.
sweep_punch <- function(data, share = 0.1, seed = 7L, columns = names(data)) {
  sweep_seeded(seed, function() {
    data[columns] <- lapply(data[columns], function(column) {
      column[stats::runif(length(column)) < share] <- NA
      column
    })
    data
  })
}

#' Categorical datasets
#' @return Named list of data.frames of factors.
sweep_categorical <- function() {
  environment <- new.env()
  utils::data(list = c("carcinoma", "values"), package = "poLCA", envir = environment)
  as_factors <- function(data) as.data.frame(lapply(data, factor))
  nominal <- sweep_seeded(202L, function() {
    class <- sample(1:3, 500L, TRUE, c(0.5, 0.3, 0.2))
    probabilities <- list(c(0.7, 0.2, 0.1), c(0.1, 0.2, 0.7), c(0.2, 0.6, 0.2))
    items <- lapply(1:5, function(j) {
      shift <- (j - 1L) %% 3L
      vapply(class, function(k) {
        sample(c("low", "mid", "high"), 1L, prob = probabilities[[(k + shift - 1L) %% 3L + 1L]])
      }, character(1))
    })
    as_factors(stats::setNames(as.data.frame(items), paste0("item", 1:5)))
  })
  esm <- latents::student_esm
  binary <- names(esm)[vapply(esm, is.factor, logical(1))]
  list(carcinoma = as_factors(environment$carcinoma), values = as_factors(environment$values),
       nominal = nominal, esm = esm[binary[seq_len(min(5L, length(binary)))]])
}

#' A two-level continuous dataset with two group classes and two profiles
#' @return data.frame: id, y1, y2, y3.
sweep_two_level <- function() {
  sweep_seeded(303L, function() {
    groups <- 60L
    size <- sample(6:12, groups, TRUE)
    class <- rep(1:2, each = groups / 2L)
    id <- rep(seq_len(groups), size)
    share <- c(0.8, 0.25)[class[id]]
    profile <- 1L + (stats::runif(length(id)) > share)
    centre <- rbind(c(0, 0, 0), c(3, -2, 2))
    out <- as.data.frame(centre[profile, ] + matrix(stats::rnorm(length(id) * 3L), length(id)))
    cbind(id = id, stats::setNames(out, c("y1", "y2", "y3")))
  })
}

#' A two-level categorical dataset
#' @return data.frame: id and four binary factors.
sweep_two_level_categorical <- function() {
  sweep_seeded(404L, function() {
    groups <- 80L
    size <- sample(8:14, groups, TRUE)
    class <- rep(1:2, each = groups / 2L)
    id <- rep(seq_len(groups), size)
    profile <- 1L + (stats::runif(length(id)) > c(0.85, 0.2)[class[id]])
    yes <- rbind(c(0.9, 0.8, 0.85, 0.75), c(0.15, 0.2, 0.1, 0.3))
    items <- lapply(1:4, function(j) {
      factor(ifelse(stats::runif(length(id)) < yes[profile, j], "yes", "no"), c("no", "yes"))
    })
    cbind(id = id, stats::setNames(as.data.frame(items), paste0("q", 1:4)))
  })
}

#' Count data with two Poisson classes and an overdispersed column
#' @return data.frame: two Poisson counts, one negative-binomial count, one continuous.
sweep_counts <- function() {
  sweep_seeded(505L, function() {
    class <- rep(1:2, c(180L, 120L))
    data.frame(visits = stats::rpois(300L, c(1.5, 6)[class]),
               posts = stats::rpois(300L, c(4, 0.8)[class]),
               clicks = stats::rnbinom(300L, mu = c(2, 9)[class], size = 2),
               score = stats::rnorm(300L, c(0, 2)[class]))
  })
}

#' Ordinal data with two classes
#' @return data.frame of four ordered factors with four categories.
sweep_ordinal <- function() {
  sweep_seeded(606L, function() {
    class <- rep(1:2, c(220L, 180L))
    probabilities <- rbind(c(0.45, 0.3, 0.15, 0.1), c(0.08, 0.17, 0.3, 0.45))
    items <- lapply(1:4, function(j) {
      factor(vapply(class, function(k) sample(1:4, 1L, prob = probabilities[k, ]), integer(1)),
             levels = 1:4, ordered = TRUE)
    })
    stats::setNames(as.data.frame(items), paste0("o", 1:4))
  })
}

#' Draw a first-order chain of profiles per unit
#' @param units Number of units.
#' @param occasions Number of occasions.
#' @param initial Initial probabilities.
#' @param transition Transition matrix.
#' @return units x occasions integer matrix of states.
sweep_chain <- function(units, occasions, initial, transition) {
  first <- sample(seq_along(initial), units, TRUE, initial)
  Reduce(function(states, occasion) {
    previous <- states[, ncol(states)]
    cbind(states, vapply(previous, function(s) {
      sample(seq_along(initial), 1L, prob = transition[s, ])
    }, integer(1)))
  }, seq_len(occasions - 1L), matrix(first, ncol = 1L))
}

#' Longitudinal datasets for transition models
#' @return Named list of sets: data, vars, id, time, and the indicator roles.
sweep_longitudinal <- function() {
  sweep_seeded(707L, function() {
    transition <- rbind(c(0.8, 0.2), c(0.3, 0.7))
    states <- as.vector(t(sweep_chain(150L, 4L, c(0.45, 0.55), transition)))
    gaussian <- data.frame(id = rep(seq_len(150L), each = 4L), wave = rep(1:4, 150L),
                           y1 = stats::rnorm(600L, c(0, 3)[states]),
                           y2 = stats::rnorm(600L, c(0, -2)[states]))
    three <- as.vector(t(sweep_chain(200L, 5L, c(0.5, 0.3, 0.2),
                                     rbind(c(0.7, 0.2, 0.1), c(0.15, 0.7, 0.15), c(0.1, 0.2, 0.7)))))
    three_state <- data.frame(id = rep(seq_len(200L), each = 5L), wave = rep(1:5, 200L),
                              y1 = stats::rnorm(1000L, c(0, 3, 6)[three]),
                              y2 = stats::rnorm(1000L, c(0, -3, 2)[three]))
    yes <- rbind(c(0.9, 0.85, 0.8, 0.9), c(0.15, 0.1, 0.2, 0.25))
    item_states <- as.vector(t(sweep_chain(200L, 4L, c(0.6, 0.4), transition)))
    items <- as.data.frame(lapply(1:4, function(j) {
      factor(ifelse(stats::runif(800L) < yes[item_states, j], "yes", "no"), c("no", "yes"))
    }))
    binary <- cbind(id = rep(seq_len(200L), each = 4L), wave = rep(1:4, 200L),
                    stats::setNames(items, paste0("q", 1:4)))
    mixed <- cbind(gaussian, visits = stats::rpois(600L, c(1.5, 5)[states]),
                   q = factor(ifelse(stats::runif(600L) < c(0.8, 0.2)[states], "yes", "no")),
                   o = factor(vapply(states, function(s) {
                     sample(1:3, 1L, prob = rbind(c(0.6, 0.3, 0.1), c(0.1, 0.3, 0.6))[s, ])
                   }, integer(1)), levels = 1:3, ordered = TRUE))
    ragged <- gaussian[stats::runif(nrow(gaussian)) > 0.15 | gaussian$wave == 1L, ]
    hours <- latents::study_hours
    list(
      gaussian = list(data = gaussian, vars = c("y1", "y2"), id = "id", time = "wave"),
      three_state = list(data = three_state, vars = c("y1", "y2"), id = "id", time = "wave"),
      ragged = list(data = ragged, vars = c("y1", "y2"), id = "id", time = "wave"),
      binary = list(data = binary, vars = paste0("q", 1:4), id = "id", time = "wave",
                    categorical = paste0("q", 1:4)),
      mixed = list(data = mixed, vars = c("y1", "visits", "q", "o"), id = "id", time = "wave",
                   categorical = "q", count = "visits", ordinal = "o"),
      hours = list(data = hours, vars = c("hours", "sleep"), id = "student", time = "week"))
  })
}

#' Regression data with two latent classes, nested in units
#' @return data.frame: unit, x (numeric), arm (factor), z (unit covariate), and a
#'   Gaussian (y), binary (b), count (n) and ordered (grade) outcome.
sweep_regression <- function() {
  sweep_seeded(808L, function() {
    units <- 60L
    per_unit <- 10L
    z <- stats::rnorm(units)
    unit_class <- 1L + (stats::runif(units) < stats::plogis(-0.2 + 1.2 * z))
    unit <- rep(seq_len(units), each = per_unit)
    class <- unit_class[unit]
    x <- stats::rnorm(units * per_unit)
    arm <- factor(sample(c("control", "treated"), units * per_unit, TRUE))
    treated <- as.integer(arm == "treated")
    eta <- cbind(1 + 2 * x + 0.5 * treated, 6 - 1.5 * x - 1 * treated)[cbind(seq_along(class), class)]
    latent <- cbind(-0.5 + 1.2 * x, 1.5 - 1 * x)[cbind(seq_along(class), class)]
    grade <- cut(latent + stats::rlogis(length(class)), c(-Inf, -0.5, 1, Inf),
                 labels = c("low", "mid", "high"), ordered_result = TRUE)
    data.frame(unit = unit, x = x, arm = arm, z = z[unit],
               y = eta + stats::rnorm(length(class), 0, c(0.8, 1.4)[class]),
               b = stats::rbinom(length(class), 1L, stats::plogis(latent)),
               n = stats::rpois(length(class), exp(c(0.3, 1.6)[class] + c(0.5, -0.3)[class] * x)),
               grade = grade)
  })
}
