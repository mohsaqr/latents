pkgload::load_all('.', quiet = TRUE)
set.seed(91)
# Unit mapping, zero weights and complete-case deletion at the tree boundary.
d <- data.frame(id = rep(1:18, rep(2:4, 6)),
                x = sequence(rep(2:4, 6)))
d$y <- ifelse(d$id <= 8, -4, 4) + .3 * d$x + rnorm(nrow(d), sd = .4)
d$w <- rep(c(0, 1, 8), length.out = 18)[d$id]
d$y[1] <- NA
models <- list(
  observation = mixture_regression(y ~ x, d, 2, id = 'id', weights = 'w', missing = 'omit',
                                    n_starts = 2, seed = 2),
  group = mixture_regression(y ~ x, d, 2, id = 'id', class_level = 'group',
                             weights = 'w', missing = 'omit', n_starts = 2, seed = 2),
  two_level = mixture_regression(y ~ x, d, 2, id = 'id', n_group_classes = 2,
                                 weights = 'w', missing = 'omit', n_starts = 2, seed = 2))
invisible(lapply(models, function(fit) {
  levels <- if (identical(fit$spec$nesting, 'two-level')) c('individuals', 'groups') else
    'individuals'
  invisible(lapply(levels, function(level) {
    draw <- .latents_tree_posteriors(fit, level)
    reference <- if (level == 'groups' || identical(fit$spec$nesting, 'group'))
      fit$spec$sampling_weights else fit$spec$row_weights
    stopifnot(nrow(draw$posteriors) == length(reference),
              identical(draw$weights, reference), all(is.finite(draw$posteriors)))
    print(list(nesting = fit$spec$nesting, level = level, units = length(reference),
               zero_weights = sum(reference == 0),
               shares = colSums(draw$posteriors * reference) / sum(reference)))
  }))
}))
# Weighted multilevel growth is explicitly unsupported. Unweighted clustered
# growth trees must count persons at one level, clusters at the other.
g <- d[!is.na(d$y), ]
g$cluster <- ceiling(g$id / 3)
growth <- suppressWarnings(mixture_regression(y ~ x, g, 2, id='id',
  class_level='group', random=~1, cluster='cluster', n_group_classes=2,
  n_starts=2, seed=2))
invisible(lapply(c('individuals', 'groups'), function(level) {
  draw <- .latents_tree_posteriors(growth, level)
  expected <- if (level == 'groups') growth$spec$n_clusters else growth$spec$n_groups
  stopifnot(nrow(draw$posteriors) == expected, identical(draw$weights, rep(1, expected)))
  print(list(growth_level=level, units=expected))
}))
# Ambiguous categories, source aliases, missing source columns and prediction subsets.
types <- list(numeric=c(2, 10), integer=1:2, logical=c(FALSE, TRUE),
              character=c('a','b'), factor=factor(c('a','b'), levels=c('b','a')),
              ordered=ordered(c('a','b'), levels=c('b','a')))
invisible(lapply(types, function(value) {
  frame <- data.frame(a=value[rep(1:2, each=2)], b=rep(0:1,2), c=rep(0:1,2))
  joined <- .latents_join_pairs(frame, .latents_parse_pairs('a:b'))
  stopifnot(nlevels(joined[['a:b']])==4L,
            identical(.latents_join_pairs(joined, .latents_parse_pairs('a:b')),joined))
  fit <- lca(frame,c('a','b','c'),1,local_dependence='a:b',n_starts=1,seed=1)
  stopifnot(all(predict(fit, frame[1,,drop=FALSE], type='posterior')==1))
}))
# Audit duplication/identity checks before the pair is removed from vars.
a <- expand.grid(a=0:1,b=0:1,c=0:1)
probe <- tryCatch(lca(a,c('a','a','b','c'),1,local_dependence='a:b',n_starts=1,seed=1),
                   error=function(e)e)
cat('Duplicate original vars result:', class(probe)[1], '\n')
cat('ADVERSARIAL REVIEW COMPLETE\n')
