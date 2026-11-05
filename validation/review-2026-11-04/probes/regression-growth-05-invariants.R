suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(5)
# Duplicates vs integer weights, row order, unit rescale, NA-omit alignment, unused factor levels.
N <- 60
d <- do.call(rbind, lapply(1:N, function(i) { k <- rbinom(1, 1, .5); t <- 0:3; b <- rnorm(1, 0, .6); data.frame(id = i, t = t, y = c(0, 2)[k + 1] + c(.8, -.4)[k + 1] * t + b + rnorm(4, sd = .5), v = rep(sample(c(1, 2), 1), 4), cnt = rpois(4, exp(c(.2, 1.4)[k + 1] + .2 * t))) }))
dup <- do.call(rbind, lapply(1:N, function(i) { a <- d[d$id == i, ]; do.call(rbind, lapply(seq_len(a$v[1]), function(r) { a$id <- i * 10 + r; a })) }))
fw <- mixture_regression(y ~ t, d, 2, id = "id", class_level = "group", random = "intercept", weights = "v", seed = 1, n_starts = 2, tol = 1e-10)
fd <- mixture_regression(y ~ t, dup, 2, id = "id", class_level = "group", random = "intercept", seed = 1, n_starts = 2, tol = 1e-10)
cat(sprintf("growth: weights LL %.8f dup LL %.8f; max coef diff %.2e\n", fw$log_likelihood, fd$log_likelihood, max(abs(coef(fw) - coef(fd)))))
print(as.data.frame(get_results(fw, "fit"))[c("n_persons","n_observations","bic","bic_rows","sabic")]); print(as.data.frame(get_results(fd, "fit"))[c("n_persons","n_observations","bic","bic_rows","sabic")])
# Row order and label invariance
o <- sample(nrow(d)); f0 <- mixture_regression(y ~ t, d, 2, id = "id", class_level = "group", random = "intercept", seed = 1, n_starts = 2, tol = 1e-10)
f1 <- mixture_regression(y ~ t, d[o, ], 2, id = "id", class_level = "group", random = "intercept", seed = 1, n_starts = 2, tol = 1e-10)
cat("row order coef diff", max(abs(coef(f0) - coef(f1))), " se", max(abs(sqrt(diag(vcov(f0))) - sqrt(diag(vcov(f1))))), "\n")
rr <- get_results(f1, "random_effects"); r0 <- get_results(f0, "random_effects"); cat("RE diff", max(abs(rr[["(Intercept)"]] - r0[["(Intercept)"]])), "\n")
ind <- get_results(f1, "individual"); hand <- d[o, ]; cat("individual aligned", all(ind$id == hand$id & ind$t == hand$t & abs(ind$observed - hand$y) < 1e-12), "\n")
# Rescale y*1000 and t*12
d3 <- transform(d, y = y * 1000 + 5e4, t = t * 12 + 2000)
f3 <- mixture_regression(y ~ t, d3, 2, id = "id", class_level = "group", random = "intercept", seed = 1, n_starts = 2, tol = 1e-10)
cat("rescale LL shift", f3$log_likelihood - f0$log_likelihood, " expected", -nrow(d) * log(1000), "\n")
cat("slope ratio", coef(f3)[["coefficient.class_1.t"]] / coef(f0)[["coefficient.class_1.t"]], " (1000/12 =", 1000/12, ")  se ratio", sqrt(vcov(f3)["coefficient.class_1.t","coefficient.class_1.t"] / vcov(f0)["coefficient.class_1.t","coefficient.class_1.t"]), "\n")
print(get_results(f0, "random")[c("class","parameter","estimate","std_error")]); print(get_results(f3, "random")[c("class","parameter","estimate","std_error")])
# NA omit: row column points to original rows; recovery with truth
dn <- d; dn$truth <- rep(sample(1:2, N, TRUE), each = 4); dn$y[c(3, 50, 100)] <- NA
fn <- suppressWarnings(mixture_regression(y ~ t, dn, 2, id = "id", class_level = "group", random = "intercept", missing = "omit", seed = 1, n_starts = 2))
pn <- predict(fn); cat("predict rows skip NA rows:", !any(c(3, 50, 100) %in% pn$row), " nrow", nrow(pn), "\n")
print(get_results(fn, "recovery", data = dn, truth = "truth"))
# Unused factor level in membership / random / formula after subsetting
dd <- d; dd$grp <- factor(ifelse(dd$v == 1, "a", "b"), levels = c("a", "b", "zz"))
fm <- withCallingHandlers(mixture_regression(y ~ t, dd, 2, id = "id", class_level = "group", membership = ~ grp, seed = 1, n_starts = 2),
  warning = function(w) { cat("  warning class", class(w)[1], ":", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
print(colnames(fm$spec$w)); cat("k", fm$n_parameters, "length coef", length(coef(fm)), "\n"); print(get_results(fm, "membership"))
fg <- mixture_regression(y ~ t + grp, dd, 2, id = "id", class_level = "group", seed = 1, n_starts = 2)
cat("formula factor unused dropped: columns", colnames(fg$spec$x), " k", fg$n_parameters, "\n")
dt <- dd; dt$grp2 <- factor(dt$grp); fmd <- mixture_regression(y ~ t, dt, 2, id = "id", class_level = "group", membership = ~ grp2, seed = 1, n_starts = 2)
cat("LL with unused level", fm$log_likelihood, " without", fmd$log_likelihood, "  BIC", get_results(fm, "fit")$bic, get_results(fmd, "fit")$bic, " k", fm$n_parameters, fmd$n_parameters, "\n")
dr <- dd
fr <- tryCatch(withCallingHandlers(mixture_regression(y ~ t, dr, 2, id = "id", class_level = "group", random = ~ 0 + grp, random_covariance = "equal", seed = 1, n_starts = 2),
  warning = function(w) { cat("  warning class", class(w)[1], ":", conditionMessage(w), "\n"); invokeRestart("muffleWarning") }), error = function(e) { cat("  error", class(e)[1], conditionMessage(e), "\n"); NULL })
if (!is.null(fr)) { print(colnames(fr$spec$random_design)); cat("random unused k", fr$n_parameters, "\n"); print(get_results(fr, "random")) }
