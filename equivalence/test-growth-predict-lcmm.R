# predict() for growth mixture models against lcmm 2.2.1 (Proust-Lima,
# Philipps and Liquet 2017). Our fit is given hlme's own estimates, so the
# predictions are compared at identical parameters and neither optimizer
# matters:
#   "class_response"    predictY(marg = TRUE)                     pred_m<k>
#   "response"          hlme's fitted pred_m (prior-weighted over classes)
#   "posterior"         predictClass()                            prob<k>
#   "class_conditional" predictY(marg = FALSE)                    pred_ss<k>
#   "conditional"       predictY(marg = FALSE)                    pred_ss
#   "random_effects"    predictRE() is their posterior-weighted mean

skip_if_not_installed("lcmm")

predict_data <- local({
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  set.seed(5)
  n <- 150L
  times <- 0:5
  out <- data.frame(id = rep(seq_len(n), each = length(times)),
                    time = rep(times, n), age = rep(stats::rnorm(n), each = length(times)))
  class <- 1L + (stats::runif(n) < stats::plogis(-0.2 + 0.8 * out$age[!duplicated(out$id)]))
  intercept <- stats::rnorm(n)
  slope <- stats::rnorm(n, 0, 0.3)
  person_class <- class[out$id]
  out$y <- c(2, 5)[person_class] + c(0.6, -0.3)[person_class] * out$time +
    intercept[out$id] + slope[out$id] * out$time + stats::rnorm(nrow(out), 0, 0.7)
  out <- out[stats::runif(nrow(out)) > 0.15, ]
  if (!is.null(old)) assign(".Random.seed", old, globalenv())
  out
})
new_persons <- subset(predict_data, id <= 12)

hlme_fit <- function(random, ng, nwg) {
  one <- lcmm::hlme(y ~ time, random = random, subject = "id", ng = 1,
                    data = predict_data, verbose = FALSE)
  hlme <- lcmm::hlme
  model <- bquote(hlme(y ~ time, mixture = ~ time, random = .(random),
                       subject = "id", ng = .(ng), nwg = .(nwg), classmb = ~ age,
                       data = predict_data, verbose = FALSE))
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (!is.null(old)) assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(2)
  eval(bquote(lcmm::gridsearch(.(model), rep = 10, maxiter = 30, minit = one)))
}

configurations <- list(
  list(label = "random intercept and slope, equal covariance", random = ~ 1 + time,
       hlme_random = ~ time, ng = 2L, nwg = FALSE, covariance = "equal"),
  list(label = "proportional covariance (nwg)", random = ~ 1 + time,
       hlme_random = ~ time, ng = 2L, nwg = TRUE, covariance = "proportional"),
  list(label = "random intercept, three classes", random = ~ 1,
       hlme_random = ~ 1, ng = 3L, nwg = FALSE, covariance = "equal"))

invisible(lapply(configurations, function(cfg) {
  test_that(sprintf("growth predict() matches lcmm: %s", cfg$label), {
    theirs <- hlme_fit(cfg$hlme_random, cfg$ng, cfg$nwg)
    expect_equal(theirs$conv, 1)
    ours <- mixture_regression(y ~ time, predict_data, n_classes = cfg$ng, id = "id",
                               class_level = "group", membership = ~ age,
                               random = cfg$random, random_covariance = cfg$covariance,
                               variance = "equal", n_starts = 2, seed = 1,
                               vcov_type = "none")
    ours$params <- params_from_hlme(ours, theirs)
    classes <- seq_len(cfg$ng)

    marginal <- lcmm::predictY(theirs, newdata = new_persons, var.time = "time")$pred
    class_response <- predict(ours, new_persons, type = "class_response")
    expect_equal(class_response$fitted, as.vector(marginal), tolerance = 1e-8)

    # hlme's fitted marginal mean weighs the classes by each person's prior.
    response <- predict(ours, predict_data, type = "response")
    expect_equal(response$fitted, theirs$pred$pred_m, tolerance = 1e-8)

    posterior <- predict(ours, new_persons, type = "posterior")
    their_posterior <- lcmm::predictClass(theirs, newdata = new_persons, subject = "id")
    expect_equal(as.matrix(posterior[sprintf("posterior_class_%d", classes)]),
                 as.matrix(their_posterior[sprintf("prob%d", classes)]),
                 tolerance = 1e-6, ignore_attr = TRUE)

    subject <- lcmm::predictY(theirs, newdata = new_persons, var.time = "time",
                              marg = FALSE, subject = "id")$pred
    class_conditional <- predict(ours, new_persons, type = "class_conditional")
    expect_equal(class_conditional$fitted,
                 unlist(subject[sprintf("pred_ss%d", classes)], use.names = FALSE),
                 tolerance = 1e-6)
    expect_equal(predict(ours, new_persons, type = "conditional")$fitted,
                 subject$pred_ss, tolerance = 1e-6)

    effects <- predict(ours, new_persons, type = "random_effects")
    terms <- setdiff(names(effects), c("id", "class", "posterior"))
    weighted <- vapply(terms, function(term) {
      as.vector(tapply(effects$posterior * effects[[term]], effects$id, sum))
    }, numeric(length(unique(effects$id))))
    their_effects <- lcmm::predictRE(theirs, newdata = new_persons, subject = "id")
    expect_equal(unname(weighted), unname(as.matrix(their_effects[-1L])),
                 tolerance = 1e-6, ignore_attr = TRUE)
  })
}))
