suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(3)
lse <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }
softmax <- function(s) exp(s) / sum(exp(s))
dmvn <- function(y, m, V) { U <- chol(V); r <- backsolve(U, y - m, transpose = TRUE); -0.5 * (length(y) * log(2 * pi) + 2 * sum(log(diag(U))) + sum(r^2)) }
N <- 90
mk <- function(N, slope = c(.8, -.5), int = c(0, 2)) {
  do.call(rbind, lapply(seq_len(N), function(i) {
    m <- sample(4:7, 1); t <- sort(sample(0:8, m)); k <- rbinom(1, 1, .4) + 1
    b <- MASS::mvrnorm(1, c(0, 0), matrix(c(1, .3, .3, .25), 2))
    data.frame(id = i, t = t, k = k, x = rnorm(m), cov = rep(rnorm(1), m), y = int[k] + slope[k] * t + b[1] + b[2] * t + rnorm(m, sd = .7))
  }))
}
d <- mk(N); d$c1 <- d$cov
d <- d[sample(nrow(d)), ]
# One class vs lme4 and nlme (ML)
f1 <- mixture_regression(y ~ t, d, 1, id = "id", class_level = "group", random = ~ 1 + t, seed = 1, tol = 1e-12)
m <- nlme::lme(y ~ t, random = ~ 1 + t | id, d, method = "ML", control = nlme::lmeControl(opt = "optim", msMaxIter = 500, tolerance = 1e-10))
m4 <- lme4::lmer(y ~ t + (1 + t | id), d, REML = FALSE, control = lme4::lmerControl(optCtrl = list(ftol_abs = 1e-12, xtol_abs = 1e-12)))
tb <- as.data.frame(f1)
cat(sprintf("LL latents %.8f nlme %.8f lme4 %.8f\n", f1$log_likelihood, logLik(m), logLik(m4)))
cat("fixed diff nlme", max(abs(tb$estimate - nlme::fixef(m))), " se ratio nlme", tb$std_error / sqrt(diag(vcov(m))), " se ratio lme4", tb$std_error / sqrt(diag(as.matrix(vcov(m4)))), "\n")
print(get_results(f1, "random")[c("parameter", "term", "estimate", "std_error", "conf_low", "conf_high")]); print(nlme::VarCorr(m))
cat("intervals nlme (natural sd, log-scale wald):\n"); print(nlme::intervals(m, which = "var-cov")$reStruct$id)
ft <- get_results(f1, "fit"); cat(sprintf("k %d lme4 %d; bic %.5f lme4 BIC(rows) %.5f, bic_rows %.5f\n", ft$n_parameters, attr(logLik(m4), "df"), ft$bic, BIC(m4), ft$bic_rows))
# Empirical Bayes vs lme4/nlme
re <- predict(f1, type = "random_effects"); rr <- nlme::ranef(m); 
cat("EB diff vs nlme", max(abs(re[["(Intercept)"]] - rr[as.character(re$id), "(Intercept)"])), max(abs(re$t - rr[as.character(re$id), "t"])), "\n")
cond <- predict(f1, type = "conditional"); cat("conditional vs nlme fitted", max(abs(cond$fitted - fitted(m, level = 1)[cond$row])), "\n")
# Two-class hand likelihood: each covariance option, membership covariate
configs <- list(
  varying = list(random_covariance = "varying", variance = "varying", random_diagonal = FALSE),
  equal = list(random_covariance = "equal", variance = "equal", random_diagonal = FALSE),
  proportional = list(random_covariance = "proportional", variance = "varying", random_diagonal = FALSE),
  diagonal = list(random_covariance = "varying", variance = "varying", random_diagonal = TRUE))
for_each <- lapply(names(configs), function(nm) {
  cfg <- configs[[nm]]
  f <- suppressWarnings(do.call(mixture_regression, c(list(y ~ t, d, 2, id = "id", class_level = "group", random = ~ 1 + t, membership = ~ cov, seed = 1, n_starts = 3, tol = 1e-10), cfg)))
  p <- f$params; G <- .growth_covariances(f$spec, p)
  ll <- sum(vapply(sort(unique(d$id)), function(i) {
    a <- d[d$id == i, ]; Z <- cbind(1, a$t)
    lse(log(softmax(c(c(1, a$cov[1]) %*% p$gamma))) + vapply(1:2, function(k) dmvn(a$y, drop(Z %*% p$beta[, k]), Z %*% G[[k]] %*% t(Z) + diag(p$sigma2[k], nrow(a))), 0))
  }, 0))
  th <- coef(f)
  nll <- function(t) -.growth_expectation(f$spec, f$stats, .growth_unpack(f$spec, t, p))$log_likelihood
  se_n <- sqrt(diag(solve(numDeriv::hessian(nll, th, method.args = list(d = 1e-3, r = 4)))))
  rnd <- get_results(f, "random")
  # delta for sd and correlation: numeric hand on the packed vcov
  V <- vcov(f); ks <- if (identical(cfg$random_covariance, "equal")) 1L else 1:2
  delta_err <- max(unlist(lapply(ks, function(k) {
    ff <- function(t) { g <- .growth_covariances(f$spec, .growth_unpack(f$spec, t, p))[[k]]; c(diag(g), sqrt(diag(g)), if (!cfg$random_diagonal) g[1, 2] / sqrt(g[1, 1] * g[2, 2])) }
    J <- numDeriv::jacobian(ff, th); se <- sqrt(diag(J %*% V %*% t(J)))
    lab <- if (identical(cfg$random_covariance, "equal")) "all" else paste0("class_", k)
    abs(rnd$std_error[rnd$class == lab] / se - 1)
  })))
  cat(sprintf("%-12s LL %.8f hand %.8f diff %.1e; k %d length(theta) %d; SE vs numDeriv max rel %.2e; random-table delta max rel %.2e\n", nm, f$log_likelihood, ll, f$log_likelihood - ll, f$n_parameters, length(th), max(abs(sqrt(diag(V)) / se_n - 1)), delta_err))
  invisible(f)
})
f <- for_each[[1]]
fitt <- get_results(f, "fit"); cat(sprintf("bic persons %.5f hand %.5f; bic_rows %.5f hand %.5f\n", fitt$bic, -2*fitt$log_likelihood + fitt$n_parameters*log(N), fitt$bic_rows, -2*fitt$log_likelihood + fitt$n_parameters*log(nrow(d))))
print(get_results(f, "membership")); print(get_results(f, "classes"))
# Predictions marginalise (hand)
p <- f$params; pr <- exp(.mixture_log_softmax(f$spec$w, p$gamma))
pred <- predict(f)
hand <- vapply(seq_len(nrow(d)), function(r) { i <- match(d$id[r], f$spec$group_levels); pi <- softmax(c(c(1, d$cov[r]) %*% p$gamma)); sum(pi * (cbind(1, d$t[r]) %*% p$beta)) }, 0)
cat("predict response vs hand (row alignment via kept_rows)", max(abs(pred$fitted - hand[pred$row])), "\n")
# conditional hand
cond <- predict(f, type = "conditional")
G <- .growth_covariances(f$spec, p)
hand2 <- numeric(nrow(d))
invisible(lapply(sort(unique(d$id)), function(i) {
  at <- which(d$id == i); a <- d[at, ]; Z <- cbind(1, a$t)
  ld <- vapply(1:2, function(k) dmvn(a$y, drop(Z %*% p$beta[, k]), Z %*% G[[k]] %*% t(Z) + diag(p$sigma2[k], nrow(a))), 0)
  lp <- log(softmax(c(c(1, a$cov[1]) %*% p$gamma))) + ld; po <- exp(lp - lse(lp))
  ow <- vapply(1:2, function(k) { V <- Z %*% G[[k]] %*% t(Z) + diag(p$sigma2[k], nrow(a)); drop(Z %*% p$beta[, k]) + drop(Z %*% G[[k]] %*% t(Z) %*% solve(V, a$y - Z %*% p$beta[, k])) }, numeric(nrow(a)))
  hand2[at] <<- drop(ow %*% po)
}))
cat("predict conditional vs hand", max(abs(cond$fitted - hand2[cond$row])), "\n")
tr <- get_results(f, "trajectories"); g <- tr$t[1:3]
cat("trajectory hand", max(abs(tr$estimate[tr$class == "class_2"][1:3] - (p$beta[1, 2] + p$beta[2, 2] * g))), "\n")
Vb <- vcov(f)[c("coefficient.class_2.(Intercept)", "coefficient.class_2.t"), c("coefficient.class_2.(Intercept)", "coefficient.class_2.t")]
cat("traj se hand", abs(tr$std_error[tr$class == "class_2"][2] - sqrt(drop(t(c(1, g[2])) %*% Vb %*% c(1, g[2])))), "\n")
# parameter inference BH set
pi <- parameter_inference(f, adjust = "BH")
print(as.data.frame(pi)[c("block", "class", "term", "estimate", "standard_error", "p_value", "p_adjusted")])
tst <- pi$block %in% c("coefficient", "membership")
cat("BH over tested rows incl intercepts equals table:", isTRUE(all.equal(p.adjust(pi$p_value[tst], "BH"), pi$p_adjusted[tst])), "\n")
