suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(4)
lse <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }
softmax <- function(s) exp(s) / sum(exp(s))
dmvn <- function(y, m, V) { U <- chol(V); r <- backsolve(U, y - m, transpose = TRUE); -0.5 * (length(y) * log(2 * pi) + 2 * sum(log(diag(U))) + sum(r^2)) }
S <- 24; per <- 8
mk <- function() do.call(rbind, lapply(seq_len(S), function(s) {
  h <- rbinom(1, 1, .5); sv <- rnorm(1)
  do.call(rbind, lapply(seq_len(per), function(p) {
    k <- rbinom(1, 1, if (h == 1) .8 else .2) + 1; t <- 0:4; b <- rnorm(1, 0, 1); pc <- rep(rnorm(1), 5)
    data.frame(school = s, id = s * 100 + p, t = t, sv = sv, pc = pc, k = k, y = c(0, 2)[k] + c(.8, -.4)[k] * t + b + rnorm(5, sd = .7))
  }))
}))
d <- mk(); d <- d[sample(nrow(d)), ]
f <- suppressWarnings(mixture_regression(y ~ t, d, 2, id = "id", class_level = "group", random = "intercept", cluster = "school", n_group_classes = 2,
   membership = ~ pc, group_membership = ~ sv, random_covariance = "equal", seed = 1, n_starts = 3, tol = 1e-10))
p <- f$params; print(p$delta); print(p$class_logits); G <- .growth_covariances(f$spec, p)
ll <- sum(vapply(sort(unique(d$school)), function(s) {
  a <- d[d$school == s, ]
  byh <- vapply(1:2, function(h) sum(vapply(unique(a$id), function(i) {
    b <- a[a$id == i, ]; Z <- matrix(1, nrow(b))
    lse(log(softmax(c(c(h == 1, h == 2, b$pc[1]) %*% p$class_logits))) + vapply(1:2, function(k) dmvn(b$y, drop(cbind(1, b$t) %*% p$beta[, k]), Z %*% G[[k]] %*% t(Z) + diag(p$sigma2[k], nrow(b))), 0))
  }, 0)), 0)
  lse(log(softmax(c(c(1, a$sv[1]) %*% p$delta))) + byh)
}, 0))
cat(sprintf("cluster LL %.9f hand %.9f diff %.2e\n", f$log_likelihood, ll, f$log_likelihood - ll))
ft <- get_results(f, "fit"); print(as.data.frame(ft)[c("n_parameters","bic","bic_persons","bic_rows","sabic","icl")])
cat("hand k: beta4+sig1(equal? no)...", length(coef(f)), " nobs", nobs(f), " bic hand", -2*ll + f$n_parameters*log(S), " persons hand", -2*ll + f$n_parameters*log(S*per), "\n")
th <- coef(f); nll <- function(t) -.growth_expectation(f$spec, f$stats, .growth_unpack(f$spec, t, p))$log_likelihood
# independent: gradient of LL vs reported total score in data coordinates; SE via Jacobian of the numeric gradient of LL
score <- colSums(.growth_scores(f$spec, f$stats, p))
cat("score vs numDeriv grad max abs", max(abs(score + numDeriv::grad(nll, th))), "\n")
se_n <- sqrt(diag(solve(numDeriv::hessian(nll, th, method.args = list(d = 1e-3, r = 4)))))
cat("cluster SE vs numDeriv Hessian max rel", max(abs(sqrt(diag(vcov(f))) / se_n - 1)), "\n")
print(get_results(f, "group_classes"))
print(get_results(f, "membership"))
# predict(response) marginal over group class and classes vs hand
pr <- predict(f); hh <- vapply(seq_len(nrow(d)), function(r) { s <- d[r, ]; eta <- softmax(c(c(1, s$sv) %*% p$delta)); pk <- Reduce(`+`, lapply(1:2, function(h) eta[h] * softmax(c(c(h == 1, h == 2, s$pc) %*% p$class_logits)))); sum(pk * c(p$beta[1, 1] + p$beta[2, 1] * s$t, p$beta[1, 2] + p$beta[2, 2] * s$t)) }, 0)
cat("cluster predict vs hand", max(abs(pr$fitted - hh[pr$row])), "\n")
# Correlated: AR and Brownian hand LL + SE
d2 <- do.call(rbind, lapply(1:80, function(i) { t <- sort(sample(1:8, 5)); k <- rbinom(1, 1, .5) + 1; b <- rnorm(1, 0, .8)
  W <- 1.0^2 * exp(-.6 * abs(outer(t, t, "-"))); e <- drop(rnorm(5) %*% chol(W + diag(.2, 5)))
  data.frame(id = i, t = t, y = c(0, 2)[k] + c(.6, -.3)[k] * t + b + e) }))
for (rc in c("ar", "brownian")) {
  fc <- withCallingHandlers(mixture_regression(y ~ t, d2, 2, id = "id", class_level = "group", random = "intercept", residual_correlation = rc, random_covariance = "equal", variance = "equal", seed = 1, n_starts = 2),
     warning = function(w) { cat("  warning:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
  p <- fc$params; G <- .growth_covariances(fc$spec, p)
  ll <- sum(vapply(unique(d2$id), function(i) {
    b <- d2[d2$id == i, ]; Z <- matrix(1, nrow(b)); W <- if (rc == "ar") p$correlation[["sd"]]^2 * exp(-p$correlation[["rate"]] * abs(outer(b$t, b$t, "-"))) else p$correlation[["sd"]]^2 * outer(b$t, b$t, pmin)
    lse(log(softmax(c(p$gamma))) + vapply(1:2, function(k) dmvn(b$y, drop(cbind(1, b$t) %*% p$beta[, k]), Z %*% G[[k]] %*% t(Z) + W + diag(p$sigma2[k], nrow(b))), 0))
  }, 0))
  th <- coef(fc); nll <- function(t) -.growth_expectation(fc$spec, fc$stats, .growth_unpack(fc$spec, t, p))$log_likelihood
  cat(sprintf("%s LL %.8f hand %.8f k %d length %d conv %s boundary {%s}\n", rc, fc$log_likelihood, ll, fc$n_parameters, length(th), fc$converged, paste(fc$residual_boundary, collapse = ",")))
  se_n <- sqrt(diag(solve(numDeriv::hessian(nll, th, method.args = list(d = 1e-3, r = 4)))))
  cat("  SE vs numDeriv max rel", max(abs(sqrt(diag(vcov(fc))) / se_n - 1), na.rm = TRUE), " score vs grad", max(abs(colSums(.growth_scores(fc$spec, fc$stats, p)) + numDeriv::grad(nll, th))), "\n")
  print(get_results(fc, "random"))
  re <- predict(fc, type = "random_effects"); b1 <- d2[d2$id == re$id[1], ]; Z <- matrix(1, nrow(b1)); W <- if (rc == "ar") p$correlation[["sd"]]^2 * exp(-p$correlation[["rate"]] * abs(outer(b1$t, b1$t, "-"))) else p$correlation[["sd"]]^2 * outer(b1$t, b1$t, pmin)
  V <- Z %*% G[[1]] %*% t(Z) + W + diag(p$sigma2[1], nrow(b1)); eb <- drop(G[[1]] %*% t(Z) %*% solve(V, b1$y - cbind(1, b1$t) %*% p$beta[, 1]))
  cat("  EB correlated class1 person1 reported", re[1, "(Intercept)"], "hand", eb, "\n")
}
