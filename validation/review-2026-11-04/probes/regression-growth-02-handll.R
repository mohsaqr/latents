suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(2)
lse <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }
softmax <- function(s) exp(s) / sum(exp(s))
# Data: groups with row-level classes shifted by a group class
J <- 40; nj <- 12
d <- data.frame(g = rep(sample(1000 + seq_len(J)), each = nj)); n <- nrow(d)
d$x <- rnorm(n); d$w <- rnorm(n); d$gv <- rep(rnorm(J), each = nj)
h <- rep(rbinom(J, 1, plogis(d$gv[seq(1, n, nj)])), each = nj)
cl <- rbinom(n, 1, ifelse(h == 1, .8, .2))
d$y <- ifelse(cl == 1, 2 + d$x, -1 - d$x) + rnorm(n, sd = .8)
d$cnt <- rpois(n, exp(ifelse(cl == 1, 1.5 + .4 * d$x, .2 - .3 * d$x)))
d <- d[sample(n), ]
# 1. two-level with membership and group covariates (hand likelihood)
f <- mixture_regression(y ~ x, d, 2, id = "g", n_group_classes = 2, membership = ~ w, group_membership = ~ gv, seed = 1, n_starts = 4, tol = 1e-10)
p <- f$params; print(p$class_logits); print(p$delta)
ids <- sort(unique(d$g))
ll <- sum(vapply(ids, function(j) {
  at <- d[d$g == j, ]
  by_h <- vapply(1:2, function(hh) {
    sum(vapply(seq_len(nrow(at)), function(i) {
      pr <- softmax(c(c(hh == 1, hh == 2, at$w[i]) %*% p$class_logits))
      dens <- vapply(1:2, function(k) dnorm(at$y[i], sum(c(1, at$x[i]) * p$beta[, k]), sqrt(p$sigma2[k]), log = TRUE), 0)
      lse(log(pr) + dens)
    }, 0))
  }, 0)
  lse(log(softmax(c(c(1, at$gv[1]) %*% p$delta))) + by_h)
}, 0))
cat(sprintf("two-level LL reported %.10f hand %.10f diff %.2e\n", f$log_likelihood, ll, f$log_likelihood - ll))
ft <- get_results(f, "fit"); k_hand <- 4 + 2 + 1*2 + (2 + 1)*1
cat(sprintf("k %d hand %d ; bic %.6f hand(groups) %.6f; bic_rows %.6f hand %.6f; sabic %.6f hand %.6f\n", ft$n_parameters, k_hand, ft$bic, -2*ll + k_hand*log(J), ft$bic_rows, -2*ll + k_hand*log(n), ft$sabic, -2*ll+k_hand*log((J+2)/24)))
cat("nobs", nobs(f), " BIC()", BIC(f), "\n")
# hand score vs numDeriv SE for two-level
th <- coef(f); nll <- function(t) -.mixture_expectation(f$spec, .mixture_unpack(f$spec, t, f$params))$log_likelihood
H <- numDeriv::hessian(nll, th); se_n <- sqrt(diag(solve(H))); se <- sqrt(diag(vcov(f)))
cat("two-level max rel SE diff vs numDeriv", max(abs(se / se_n - 1)), "\n")
# Table labels
print(get_results(f, "classes")); print(get_results(f, "group_classes"))
cat("rows", n, "sum n_assigned", sum(get_results(f, "classes")$n_assigned), "\n")
# ICL hand: H(rho) + sum_h rho_h sum_i H(k_i|h,y)
ex <- f$expectation; r <- ex$rho; H0 <- -sum(ifelse(r > 0, r * log(r), 0))
w_ent <- sum(vapply(1:2, function(hh) {
  jt <- ex$tau_by_group_class[[hh]]; wt <- r[f$spec$group_index, hh]; c1 <- jt / wt
  sum(wt * rowSums(-c1 * log(c1)))
}, 0))
cat(sprintf("ICL %.6f hand %.6f\n", ft$icl, ft$bic + 2 * (H0 + w_ent)))

# 2. group nesting with membership covariate, poisson: hand LL and logit reference
d2 <- d; d2$gv2 <- d2$gv
cnt_g <- rep(rbinom(J, 1, plogis(d$gv[match(ids, d$g)])), each = nj)
dg <- data.frame(g = rep(ids, each = nj), x = rnorm(n)); dg$gv <- rep(rnorm(J), each = nj)
cc <- rep(rbinom(J, 1, plogis(1.2 * dg$gv[seq(1, n, nj)])), each = nj)
dg$cnt <- rpois(n, exp(ifelse(cc == 1, 1.5 + .4 * dg$x, .2 - .3 * dg$x)))
f2 <- mixture_regression(cnt ~ x, dg, 2, family = "poisson", id = "g", class_level = "group", membership = ~ gv, seed = 1, n_starts = 3, tol = 1e-10)
p <- f2$params
ll <- sum(vapply(ids, function(j) {
  at <- dg[dg$g == j, ]
  lse(log(softmax(c(c(1, at$gv[1]) %*% p$gamma))) + vapply(1:2, function(k) sum(dpois(at$cnt, exp(cbind(1, at$x) %*% p$beta[, k]), log = TRUE)), 0))
}, 0))
cat(sprintf("group poisson LL %.10f hand %.10f diff %.2e; BIC %.5f hand %.5f\n", f2$log_likelihood, ll, f2$log_likelihood - ll, get_results(f2, "fit")$bic, -2*ll + f2$n_parameters*log(J)))
print(get_results(f2, "membership"))
cat("gamma col1 zero:", all(p$gamma[, 1] == 0), "  sizes", colSums(f2$expectation$group_tau), "\n")
th <- coef(f2); nll <- function(t) -.mixture_expectation(f2$spec, .mixture_unpack(f2$spec, t, f2$params))$log_likelihood
se_n <- sqrt(diag(solve(numDeriv::hessian(nll, th))))
cat("group poisson max rel SE diff", max(abs(sqrt(diag(vcov(f2))) / se_n - 1)), "\n")
# 3. Observation nesting with id, H = 1: BIC unit, robust cluster hand
f3 <- mixture_regression(y ~ x, d, 2, id = "g", seed = 1, n_starts = 3, vcov_type = "robust", tol = 1e-10)
ft3 <- get_results(f3, "fit"); cat("obs+id: nesting", ft3$nesting, " nobs", nobs(f3), " bic unit rows:", isTRUE(all.equal(ft3$bic, ft3$bic_rows)), "\n")
print(get_results(f3, "classes")[c("class", "share", "n_assigned", "sigma", "sigma_std_error")])
