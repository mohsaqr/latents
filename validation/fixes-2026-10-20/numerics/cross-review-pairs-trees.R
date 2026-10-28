pkgload::load_all(quiet = TRUE)
# Read-only cross-review: reproduce data-frame provenance across common base
# transforms, and inspect posterior/weight units after complete-case omission.
record <- function(label, expression) {
  cat("\n", label, "\n", sep = "")
  tryCatch(print(force(expression)), error = function(condition) {
    cat(paste(class(condition), collapse = "/"), ": ", conditionMessage(condition), "\n", sep = "")
  })
}
d <- expand.grid(a = 0:1, b = 0:1, c = 0:1)
pairs <- .latents_parse_pairs("a:b")
joined <- .latents_join_pairs(d, pairs)
record("row-bound prepared data", identical(
  as.character(.latents_join_pairs(rbind(joined[1:3, ], joined[4:8, ]), pairs)[["a:b"]]),
  as.character(joined[["a:b"]])))
record("serialized prepared subset", identical(
  as.character(.latents_join_pairs(unserialize(serialize(joined[1:4, ], NULL)), pairs)[["a:b"]]),
  as.character(joined[["a:b"]][1:4])))
record("raw joint-only fitting", lca(joined[c("a:b", "c")], c("a", "b", "c"),
  1, local_dependence = "a:b", n_starts = 1, seed = 1))
record("raw joint-only prediction", {
  fit <- lca(d, c("a", "b", "c"), 1, local_dependence = "a:b", n_starts = 1, seed = 1)
  all.equal(predict(fit, joined[c("a:b", "c")], type = "posterior"),
            predict(fit, d, type = "posterior"))
})
set.seed(43)
ids <- rep(sample(101:140), rep(c(2, 3, 4, 5), 10))
truth <- rep(rep(c(-3, 3), each = 20), rep(c(2, 3, 4, 5), 10))
frame <- data.frame(id = ids, y = stats::rnorm(length(ids), truth, .7),
                    x = sequence(rep(c(2, 3, 4, 5), 10)))
frame$w <- rep(seq_len(40), rep(c(2, 3, 4, 5), 10))
frame$y[c(1, 2, 30)] <- NA_real_
frame$x[45] <- NA_real_
invisible(lapply(c("observation", "group"), function(level) {
  fit <- suppressWarnings(mixture_regression(y ~ x, frame, 2, id = "id",
    class_level = level, weights = "w", missing = "omit", n_starts = 2, seed = 1))
  draw <- .latents_tree_posteriors(fit, "individuals")
  reference <- if (level == "group") fit$spec$sampling_weights else fit$spec$row_weights
  record(paste("omitted rows", level, "weight alignment"),
         list(posterior_rows = nrow(draw$posteriors), weight_length = length(draw$weights),
              weights_agree = identical(draw$weights, reference),
              shares = colSums(draw$posteriors * draw$weights) / sum(draw$weights)))
  stopifnot(nrow(draw$posteriors) == length(draw$weights), identical(draw$weights, reference))
}))
invisible(lapply(c("observation", "group"), function(level) {
  grid <- suppressWarnings(enumerate_regressions(y ~ x, frame, n_classes = 1:2,
    id = "id", class_level = level, weights = "w", missing = "omit",
    n_starts = 2, seed = 1))
  flows <- get_results(grid, "tree")
  candidate <- .latents_tree_fit(grid, 2L)
  draw <- .latents_tree_posteriors(candidate, "individuals")
  expected <- colSums(draw$posteriors * draw$weights) / sum(draw$weights)
  actual <- flows$to_share[order(flows$to_class)]
  record(paste("public omitted-row weighted tree", level),
         list(shares = actual, reference = expected, max_error = max(abs(actual - expected))))
  stopifnot(max(abs(actual - expected)) < 1e-12)
}))
fit <- lca(d, c("a", "b", "c"), 1, local_dependence = "a:b", n_starts = 1, seed = 1)
transforms <- list(subset = subset(joined, a == 0),
                   merge = merge(joined, data.frame(c = 0:1, z = 1:2)),
                   rowbound_then_columns = rbind(joined, joined)[c("a", "b", "a:b", "c")])
invisible(lapply(names(transforms), function(name) {
  record(paste("prepared public prediction after", name), {
    posterior <- predict(fit, transforms[[name]], type = "posterior")
    list(rows = nrow(posterior), expected_rows = nrow(transforms[[name]]))
  })
}))
