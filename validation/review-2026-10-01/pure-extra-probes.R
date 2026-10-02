pkgload::load_all('.', quiet = TRUE)
set.seed(990)
n <- 240L
z <- rep(1:2, each = n / 2)
d <- data.frame(g = rep(seq_len(n / 3), each = 3), x = rnorm(n),
                o = ordered(vapply(z, function(k) sample(1:4, 1L, prob=if(k==1L) c(.7,.2,.08,.02) else c(.02,.08,.2,.7)), integer(1))),
                k = rnbinom(n, mu = c(2, 12)[z], size = 3),
                a = factor(ifelse(runif(n) < c(.15,.8)[z], 'a', 'b')))
specs <- list(ordinal = list(vars='o', ordinal='o'),
              poisson = list(vars='k', count='k'),
              nb = list(vars='k', count='k', count_model='negative_binomial'),
              categorical = list(vars='a', categorical='a'),
              mixed = list(vars=c('o','k','a'), ordinal='o', count='k', categorical='a'))
results <- lapply(names(specs), function(name) {
  lapply(c(FALSE, TRUE), function(covariates) {
    args <- c(list(data=d, id='g', n_profiles=2, n_group_classes=1,
                   n_starts=1, seed=5), specs[[name]],
              if (covariates) list(profile_covariates='x'))
    f <- do.call(multilpa, args)
    tables <- get_results(f, 'all')
    stopifnot(all(vapply(tables, is.data.frame, logical(1))), nrow(tables$profiles)==0L,
              length(coef(f,scale='unconstrained'))==f$n_parameters)
    if (!covariates) stopifnot(all(is.finite(predict(f,d,type='density')$log_density)))
    invisible(summary(f))
    sim <- if (covariates) .multilpa_cov_simulate(f,d) else .multilpa_simulate(f)
    stopifnot(nrow(sim)==n, all(f$vars %in% names(sim)))
    refit_args <- if (covariates) args[setdiff(names(args), c('data','n_starts','seed'))] else .multilpa_refit_arguments(f)
    r <- do.call(multilpa, c(list(data=d,n_starts=1,seed=5), refit_args))
    stopifnot(r$n_parameters==f$n_parameters, isTRUE(all.equal(r$log_likelihood,f$log_likelihood)))
    cat(name, covariates, 'OK', length(coef(f,scale='unconstrained')), '\n')
    list(spec=name, covariates=covariates, n_parameters=f$n_parameters)
  })
})
saveRDS(results,'validation/review-2026-10-01/pure-extra-probes.rds')
