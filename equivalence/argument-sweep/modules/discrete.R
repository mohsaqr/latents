# The count and ordinal arguments of multilpa() (`count`, `count_model`,
# `count_dispersion`, `ordinal`) alone and mixed with continuous and
# categorical indicators. References: closed forms and MASS::glm.nb at one
# profile, flexmix's Poisson EM continued from the package's posteriors, the
# adjacent-category structure checked on the reported probabilities, and the
# likelihood rebuilt from the result tables with its local-maximum check.

sweep_module_discrete <- function() {
  # flexmix's logLik() is an S4 method found only with the package attached.
  suppressPackageStartupMessages(library(flexmix))
  counts <- sweep_counts()
  ordinal <- sweep_ordinal()
  fit <- function(data, k, vars, ..., n_starts = 10L) {
    lpa(data, vars, n_profiles = k, seed = 41L, n_starts = n_starts, tol = 1e-10,
        max_iter = 20000L, ...)
  }
  ll_of <- function(x) sweep_table(x, "information_criteria")$log_likelihood
  par_of <- function(x) sweep_table(x, "information_criteria")$n_parameters

  poisson_cases <- lapply(sweep_grid(k = 1:3, vars = c("visits", "visits+posts",
                                                         "visits+posts+score")),
                          function(option) {
    sweep_define(sprintf("poisson[k=%d,vars=%s]", option$k, option$vars), function() {
      vars <- strsplit(option$vars, "+", fixed = TRUE)[[1L]]
      count <- setdiff(vars, "score")
      model <- fit(counts, option$k, vars, count = count)
      table <- sweep_table(model, "count_means")
      rows <- rbind(
        sweep_likelihood_checks(model, counts),
        sweep_check("parameter count", par_of(model),
                    (option$k - 1L) + option$k * length(count) +
                      2L * option$k * ("score" %in% vars), 0))
      if (option$k == 1L) {
        return(rbind(rows, sweep_check("one profile: Poisson means are the sample means",
                                       table$mean[match(count, table$indicator)],
                                       vapply(counts[count], mean, numeric(1)), 1e-8)))
      }
      if ("score" %in% vars) return(rows)
      components <- lapply(count, function(v) {
        flexmix::FLXMRglm(stats::as.formula(sprintf("%s ~ 1", v)), family = "poisson")
      })
      continued <- suppressWarnings(flexmix::flexmix(
        ~ 1, data = counts, k = option$k, model = components, cluster = sweep_posteriors(model),
        control = list(tolerance = 1e-12, iter.max = 20000L, minprior = 0)))
      rbind(rows, sweep_check("flexmix's EM from the package's posteriors: log likelihood",
                              ll_of(model), as.numeric(logLik(continued)), 1e-6))
    })
  })

  negative_binomial_cases <- lapply(
    sweep_grid(k = 1:2, dispersion = c("varying", "equal"), vars = c("clicks", "clicks+visits")),
    function(option) {
      sweep_define(sprintf("negative-binomial[k=%d,dispersion=%s,vars=%s]", option$k,
                           option$dispersion, option$vars), function() {
        vars <- strsplit(option$vars, "+", fixed = TRUE)[[1L]]
        model <- fit(counts, option$k, vars, count = vars, count_model = "negative_binomial",
                     count_dispersion = option$dispersion)
        table <- sweep_table(model, "count_means")
        dispersions <- if (identical(option$dispersion, "equal")) 1L else option$k
        rows <- rbind(
          sweep_likelihood_checks(model, counts),
          sweep_check("parameter count", par_of(model),
                      (option$k - 1L) + length(vars) * (option$k + dispersions), 0),
          if (identical(option$dispersion, "equal")) {
            sweep_check("equal dispersion is one value per indicator",
                        as.numeric(tapply(table$dispersion, table$indicator,
                                          function(d) diff(range(d)))),
                        rep(0, length(vars)), 1e-10)
          })
        if (option$k > 1L) return(rows)
        reference <- lapply(vars, function(v) {
          MASS::glm.nb(stats::as.formula(sprintf("%s ~ 1", v)), data = counts,
                       control = stats::glm.control(epsilon = 1e-12, maxit = 200L))
        })
        rbind(rows,
              sweep_check("one profile: glm.nb's log likelihood", ll_of(model),
                          sum(vapply(reference, function(m) as.numeric(stats::logLik(m)),
                                     numeric(1))), 1e-7),
              sweep_check("one profile: glm.nb's means",
                          table$mean[match(vars, table$indicator)],
                          vapply(reference, function(m) exp(stats::coef(m)[[1L]]), numeric(1)),
                          1e-6),
              sweep_check("one profile: glm.nb's dispersion (1 / theta)",
                          table$dispersion[match(vars, table$indicator)],
                          vapply(reference, function(m) 1 / m$theta, numeric(1)), 1e-3))
      })
    })

  ordinal_cases <- lapply(sweep_grid(k = 1:3, items = 2:4), function(option) {
    sweep_define(sprintf("ordinal[k=%d,items=%d]", option$k, option$items), function() {
      vars <- names(ordinal)[seq_len(option$items)]
      model <- fit(ordinal, option$k, vars, ordinal = vars)
      table <- sweep_table(model, "ordinal")
      categories <- nlevels(ordinal[[1L]])
      # Adjacent-category logit: against the first profile, the log odds of
      # each next category move by the difference of locations, whatever the category.
      structure <- do.call(rbind, lapply(vars, function(v) {
        rows <- table[table$indicator == v, , drop = FALSE]
        log_p <- sweep_wide(transform(rows, log_p = log(probability)), "profile", "category",
                            "log_p")
        steps <- t(apply(log_p, 1L, diff))
        location <- rows$location[match(rownames(log_p), rows$profile)]
        observed <- sweep(steps, 2L, steps[1L, ])
        expected <- matrix(location - location[[1L]], nrow(steps), ncol(steps))
        data.frame(observed = as.numeric(observed), expected = as.numeric(expected))
      }))
      rows <- rbind(
        sweep_likelihood_checks(model, ordinal),
        sweep_check("parameter count", par_of(model),
                    (option$k - 1L) + length(vars) * ((option$k - 1L) + (categories - 1L)), 0),
        sweep_check("adjacent-category log odds differ by the locations", structure$observed,
                    structure$expected, 1e-6))
      if (option$k > 1L) return(rows)
      rbind(rows, sweep_check("one profile: probabilities are the observed proportions",
                              table$probability,
                              unlist(lapply(ordinal[vars], function(v) prop.table(table(v)))),
                              1e-7))
    })
  })

  equivalence_cases <- list(
    sweep_define("ordinal-with-two-categories-is-categorical", function() {
      binary <- as.data.frame(lapply(ordinal, function(v) {
        factor(as.integer(v) > 2L, c(FALSE, TRUE), c("low", "high"), ordered = TRUE)
      }))
      as_ordinal <- fit(binary, 2L, names(binary), ordinal = names(binary))
      as_categorical <- fit(binary, 2L, names(binary), categorical = names(binary))
      rbind(sweep_check("log likelihood", ll_of(as_ordinal), ll_of(as_categorical), 1e-7),
            sweep_check("parameter count", par_of(as_ordinal), par_of(as_categorical), 0))
    }),
    sweep_define("ordinal-input-types", function() {
      vars <- names(ordinal)
      base <- fit(ordinal, 2L, vars, ordinal = vars)
      integers <- as.data.frame(lapply(ordinal, as.integer))
      sweep_check("integer codes against ordered factors: log likelihood", ll_of(base),
                  ll_of(fit(integers, 2L, vars, ordinal = vars)), 1e-10)
    }),
    sweep_define("count-input-types", function() {
      base <- fit(counts, 2L, c("visits", "posts"), count = c("visits", "posts"))
      doubles <- transform(counts, visits = as.numeric(visits), posts = as.numeric(posts))
      sweep_check("double against integer counts: log likelihood", ll_of(base),
                  ll_of(fit(doubles, 2L, c("visits", "posts"), count = c("visits", "posts"))), 0)
    }),
    sweep_define("poisson-is-the-negative-binomial-limit", function() {
      # The Poisson fit is a negative-binomial fit at dispersion zero, so the
      # negative binomial's maximum cannot be lower.
      poisson <- fit(counts, 2L, c("visits", "posts"), count = c("visits", "posts"))
      wider <- fit(counts, 2L, c("visits", "posts"), count = c("visits", "posts"),
                   count_model = "negative_binomial")
      sweep_at_least("negative binomial at least the Poisson likelihood", ll_of(wider),
                     ll_of(poisson), 1e-7)
    }))

  mixed <- cbind(counts, ordinal[seq_len(nrow(counts)), 1:2],
                 item = sweep_categorical()$nominal$item1[seq_len(nrow(counts))])
  mixed_cases <- lapply(sweep_grid(k = 1:3, missing = c("complete", "fiml"),
                                   weights = c("none", "weighted")), function(option) {
    sweep_define(sprintf("mixed-indicators[k=%d,%s,%s]", option$k, option$missing, option$weights),
                 function() {
      vars <- c("score", "visits", "o1", "o2", "item")
      data <- if (identical(option$missing, "fiml")) sweep_punch(mixed, 0.1, columns = vars) else mixed
      w <- sweep_seeded(17L, function() stats::runif(nrow(data), 0.5, 2))
      arguments <- c(
        list(cbind(data, w = w), option$k, vars, count = "visits", ordinal = c("o1", "o2"),
             categorical = "item"),
        if (identical(option$missing, "fiml")) list(missing = "fiml"),
        if (identical(option$weights, "weighted")) list(weights = "w"))
      model <- do.call(fit, arguments)
      observed <- rowSums(!is.na(data[vars])) > 0L
      scaled <- if (identical(option$weights, "weighted")) w * nrow(data) / sum(w) else NULL
      sweep_likelihood_checks(model, data, weights = scaled, n_bic = sum(observed),
                              empty = !observed)
    })
  })

  refusal_cases <- list(sweep_define("refusals", function() {
    rbind(
      sweep_refuses("a negative count",
                    function() lpa(transform(counts, visits = visits - 3L), "visits", 2L,
                                   count = "visits")),
      sweep_refuses("a fractional count",
                    function() lpa(transform(counts, visits = visits + 0.5), "visits", 2L,
                                   count = "visits")),
      sweep_refuses("a count named in categorical too",
                    function() lpa(counts, "visits", 2L, count = "visits", categorical = "visits")),
      sweep_refuses("a count outside vars",
                    function() lpa(counts, "visits", 2L, count = "posts")),
      sweep_refuses("an unknown count model",
                    function() lpa(counts, "visits", 2L, count = "visits", count_model = "zip")),
      sweep_warns("count_dispersion with a Poisson model warns",
                  function() lpa(counts, "visits", 2L, count = "visits", n_starts = 2L, seed = 1L,
                                 count_dispersion = "equal"), "latents_unused_argument"),
      sweep_warns("count_model without a count indicator warns",
                  function() lpa(counts, "score", 2L, n_starts = 2L, seed = 1L,
                                 count_model = "negative_binomial"), "latents_unused_argument"),
      sweep_refuses("an ordinal indicator with one category",
                    function() lpa(transform(ordinal, o1 = factor(1L)), names(ordinal), 2L,
                                   ordinal = names(ordinal))),
      sweep_refuses("a start with ordinal indicators",
                    function() lpa(ordinal, names(ordinal), 2L, ordinal = names(ordinal),
                                   start = list(means = matrix(0, 2L, 0L)))))
  }))

  c(poisson_cases, negative_binomial_cases, ordinal_cases, equivalence_cases, mixed_cases,
    refusal_cases)
}
