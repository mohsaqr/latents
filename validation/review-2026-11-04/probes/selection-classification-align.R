suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(6)
ng <- 60; nper <- 5
g <- sample(1:2, ng, TRUE); id <- rep(seq_len(ng), each = nper)
k <- vapply(g[id], function(h) sample(1:3, 1, prob = if (h == 1) c(.6,.3,.1) else c(.1,.3,.6)), 1L)
d <- data.frame(id = id, y1 = rnorm(length(id), c(-2,0,2)[k]), y2 = rnorm(length(id), c(0,2,-2)[k], c(1,.7,1.2)[k]), z = rnorm(length(id)), gz = rnorm(ng)[id])
d$a <- ifelse(runif(nrow(d)) < c(.2,.5,.85)[k], "x", "y")
d$o <- factor(pmin(4, pmax(1, round(2.5 + c(-1,0,1)[k] + rnorm(nrow(d), 0, .7)))), ordered = TRUE)
d$cnt <- rnbinom(nrow(d), size = 3, mu = c(1,4,9)[k])
loglik_check <- function(x) {
  fams <- inherits(x, "multilpa_covariates")
  data <- .multilpa_resolve_data(x, NULL)
  if (fams) { .multilpa_cov_check_likelihood(x, if (is.null(x$stored_data)) d else data); return("ok (cov check)") }
  continuous <- .multilpa_continuous_names(x)
  X <- if (length(continuous) == 0L) matrix(numeric(0), nrow(d), 0L) else .multilpa_center_like(x, as.matrix(d[continuous]))
  codes <- if (length(x$categorical %||% character()) == 0L) NULL else .multilpa_encode_categorical(d[, x$categorical, drop = FALSE])$codes
  extra <- if (is.null(x$extra_data)) NULL else .latents_prepare_extra_like(x, d)
  e <- .multilpa_expectation(X, x$group_index, x, codes, extra = extra)
  list(ll = e$log_likelihood, post = e$subject_posteriors %||% e$posteriors)
}
perm_check <- function(label, f, porder = c(3L,1L,2L), gorder = c(2L,1L)) {
  cov <- inherits(f, "multilpa_covariates")
  p <- if (cov) .multilpa_cov_permute(f, porder, gorder) else .multilpa_permute_group_classes(.multilpa_permute_profiles(f, porder), gorder)
  out <- if (cov) { tryCatch({.multilpa_cov_check_likelihood(p, d); "LL reproduced"}, error = function(e) paste("LL NOT reproduced:", conditionMessage(e))) } else {
    ll <- loglik_check(p); sprintf("LL diff %.2e", ll$ll - f$log_likelihood) }
  # assignment consistency
  ok1 <- identical(as.integer(p$subject_profiles), max.col(p$subject_posteriors, "first"))
  ok2 <- identical(as.integer(p$group_classes), max.col(p$group_posteriors, "first"))
  ic <- all.equal(unlist(get_results(p,"information_criteria",format="wide")), unlist(get_results(f,"information_criteria",format="wide")), tolerance = 1e-10)
  # classification table follows
  a <- get_results(f, "classification"); b <- get_results(p, "classification")
  cls <- all.equal(a$average_posterior[a$level=="individuals"][porder], b$average_posterior[b$level=="individuals"], tolerance = 1e-10)
  cat(sprintf("%-34s %s | modal ok %s/%s | IC same %s | classification follows %s\n", label, out, ok1, ok2, isTRUE(ic), isTRUE(cls)))
  invisible(p)
}
base <- function(...) multilpa(d, id = "id", n_profiles = 3, n_group_classes = 2, n_starts = 4, seed = 1, ...)
f1 <- base(vars = c("y1","y2"), model = "VVI"); perm_check("VVI", f1)
f2 <- base(vars = c("y1","y2"), model = "VVV"); perm_check("VVV", f2)
f3 <- base(vars = c("y1","y2","a"), categorical = "a", model = "VVI"); perm_check("mixed categorical", f3)
f4 <- base(vars = c("y1","o"), ordinal = "o", model = "VVI"); perm_check("ordinal", f4)
f5 <- base(vars = c("y1","cnt"), count = "cnt", count_model = "negative_binomial", model = "VVI"); perm_check("NB count", f5)
f6 <- base(vars = c("y1","y2"), model = "VVI", direct_effects = "y1 ~ z"); perm_check("direct effect", f6)
f7 <- base(vars = c("y1","y2"), variance_model = "varying", local_dependence = "y1:y2"); perm_check("pair covariance", f7)
c1 <- base(vars = c("y1","y2"), model = "VVI", profile_covariates = "z"); perm_check("profile covariate", c1)
c2 <- base(vars = c("y1","y2"), model = "VVI", profile_covariates = "z", profile_slopes = "group_class"); perm_check("slopes by group class", c2)
c3 <- base(vars = c("y1","y2"), model = "VVI", group_covariates = "gz", profile_covariates = "z"); perm_check("group covariate", c3)
# Field-level audit: which fields have a profile or group-class dimension and did NOT change
audit <- function(f, p, K, G) {
  n <- names(f)
  cand <- n[vapply(n, function(nm) { v <- f[[nm]]; (is.array(v) || is.numeric(v) || is.list(v)) && !identical(v, p[[nm]]) == FALSE }, TRUE)]
  dims <- vapply(cand, function(nm) { v <- f[[nm]]; d0 <- if (is.array(v)) dim(v) else length(v); any(d0 %in% c(K, G)) }, TRUE)
  cand[dims]
}
pp <- .multilpa_permute_group_classes(.multilpa_permute_profiles(f4, c(3L,1L,2L)), c(2L,1L))
cat("Unchanged fields with a dimension of 3 or 2 (plain fit, ordinal):\n"); print(audit(f4, pp, 3L, 2L))
qq <- .multilpa_cov_permute(c2, c(3L,1L,2L), c(2L,1L))
cat("Unchanged fields, covariate:\n"); print(audit(c2, qq, 3L, 2L))
# End to end: a refit on shuffled rows aligned to the reference
ord <- sample(ng); d2 <- do.call(rbind, split(d, d$id)[ord])
fs <- multilpa(d2, c("y1","y2"), "id", n_profiles = 3, n_group_classes = 2, model = "VVI", n_starts = 4, seed = 7)
al <- align_labels(fs, f1)
cat("LL same", all.equal(al$log_likelihood, fs$log_likelihood), "\n")
print(round(rbind(f1 = as.vector(f1$means), aligned = as.vector(al$means)), 3))
print(round(rbind(f1 = as.vector(f1$profile_probabilities), aligned = as.vector(al$profile_probabilities)), 3))
# Noise
fn <- lpa(d, c("y1","y2"), n_profiles = 3, model = "VVI", noise = TRUE, n_starts = 3, seed = 1)
pn <- .multilpa_permute_profiles(fn, c(3L,1L,2L))
cat("noise: modal", all(pn$subject_profiles[fn$subject_profiles == 0] == 0), " noise class stays; entropy same", isTRUE(all.equal(get_results(fn,"entropy"), get_results(pn,"entropy"))), "\n")
# Label-permutation invariance of diagnostics numbers and of the residuals table
r0 <- get_results(f3, "residuals"); r1 <- get_results(perm_check("mixed again", f3), "residuals")
s <- function(r) sort(round(r$statistic[is.finite(r$statistic)], 8)); cat("residual stats invariant:", isTRUE(all.equal(s(r0), s(r1))), "\n")
# BCH invariance
b0 <- get_results(f1, "bch_weights", level = "individuals"); pb <- .multilpa_permute_profiles(f1, c(3L,1L,2L)); b1 <- get_results(pb, "bch_weights", level = "individuals")
cat("BCH invariant under relabel:", isTRUE(all.equal(sort(b0$weight), sort(b1$weight))), "\n")
# sensitivity agreement by hand
sens <- sensitivity(f1, seeds = 1:3); print(as.data.frame(sens))
