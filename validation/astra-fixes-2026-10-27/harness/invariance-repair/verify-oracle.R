pkgload::load_all('.', quiet = TRUE)
source('equivalence/harness/registry.R')
source('equivalence/harness/suites/invariance.R')
saved <- readRDS('validation/astra-fixes-2026-10-27/harness/invariance-repair/oracle-fixture.rds')
fit <- saved$original
expected <- .expect_excluded_gaussian_groups(saved$tables, fit, saved$extra,
                                             'school', c('y1','y2'), 150L)
# An independent exhaustive oracle: enumerate every row-profile path and
# group class for each of the four groups (2^10 x 2 joint possibilities).
by_group <- split(seq_len(nrow(saved$extra)), saved$extra$school)
checks <- lapply(names(by_group), function(group) {
  at <- by_group[[group]]
  paths <- as.matrix(expand.grid(rep(list(seq_len(nrow(fit$means))), length(at))))
  joint <- vapply(seq_along(fit$group_probabilities), function(h) {
    apply(paths, 1, function(path) {
      log(fit$group_probabilities[h]) + sum(log(fit$profile_probabilities[h,path])) +
        sum(stats::dnorm(saved$extra$y1[at], fit$means[path,1], sqrt(fit$variances[path,1]), log=TRUE)) +
        sum(stats::dnorm(saved$extra$y2[at], fit$means[path,2], sqrt(fit$variances[path,2]), log=TRUE))
    })
  }, numeric(nrow(paths)))
  mass <- exp(joint-max(joint)); total<-sum(mass); mass<-mass/total
  likelihood <- max(joint)+log(total)
  group_posterior <- colSums(mass)
  posterior <- vapply(seq_len(nrow(fit$means)),function(k) {
    vapply(seq_along(at),function(i)sum(mass[paths[,i]==k,,drop=FALSE]),numeric(1))
  },numeric(length(at)))
  want_group <- expected$group_posteriors[as.character(expected$group_posteriors$group)==group,]
  want_rows <- expected$posteriors[expected$posteriors$row %in% (150L+at),]
  c(group=max(abs(want_group$posterior-group_posterior)),
    likelihood=max(abs(want_group$log_likelihood-likelihood)),
    posterior=max(abs(want_rows$posterior-as.vector(posterior))))
})
errors <- do.call(rbind, checks)
print(errors)
stopifnot(all(errors < 1e-12))
old <- saved$tables
old$model$n_observations <- nrow(saved$padded)
old$model$n_groups <- length(unique(saved$padded$school))
expected$model <- old$model
before <- .invariance_compare('oracle check',saved$actual,old)
after <- .invariance_compare('oracle check',saved$actual,expected)
print(names(after))
print(table(before$ledger$status))
print(table(after$ledger$status))
stopifnot(sum(before$ledger$status=='partial')==11L,
          all(!after$ledger$status %in% c('partial','not covered')),
          all(after$comparisons$agrees))
