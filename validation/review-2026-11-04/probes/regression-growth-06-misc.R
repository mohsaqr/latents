suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(6)
lse <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }
softmax <- function(s) exp(s) / sum(exp(s))
dmvn <- function(y, m, V) { U <- chol(V); r <- backsolve(U, y - m, transpose = TRUE); -0.5 * (length(y) * log(2 * pi) + 2 * sum(log(diag(U))) + sum(r^2)) }
## A. cluster SE vs an independent difference of the reported LL, step scaled to the SE, with Richardson
S <- 20; per <- 6
d <- do.call(rbind, lapply(seq_len(S), function(s) { h <- rbinom(1, 1, .5)
  do.call(rbind, lapply(seq_len(per), function(p) { k <- rbinom(1, 1, if (h == 1) .8 else .2) + 1; t <- 0:3; b <- rnorm(1, 0, 1)
    data.frame(school = s, id = s * 100 + p, t = t, y = c(0, 2)[k] + c(.8, -.4)[k] * t + b + rnorm(4, sd = .7)) })) }))
f <- suppressWarnings(mixture_regression(y ~ t, d, 2, id = "id", class_level = "group", random = "intercept", cluster = "school", n_group_classes = 2, random_covariance = "equal", seed = 1, n_starts = 3, tol = 1e-10))
th <- coef(f); p <- f$params; se <- sqrt(diag(vcov(f)))
ll <- function(t) .growth_expectation(f$spec, f$stats, .growth_unpack(f$spec, t, p))$log_likelihood
hess <- function(h) { n <- length(th); H <- matrix(0, n, n); f0 <- ll(th)
  for (i in seq_len(n)) for (j in seq_len(i)) { # small fixed-size finite-difference table; vectorising adds nothing
    ei <- replace(numeric(n), i, h * se[i]); ej <- replace(numeric(n), j, h * se[j])
    H[i, j] <- H[j, i] <- (ll(th + ei + ej) - ll(th + ei - ej) - ll(th - ei + ej) + ll(th - ei - ej)) / (4 * h^2 * se[i] * se[j]) }
  H }
H1 <- hess(.2); H2 <- hess(.1); H <- (4 * H2 - H1) / 3
cat("cluster: SE vs Richardson 2nd difference of LL, max rel", max(abs(sqrt(diag(solve(-H))) / se - 1)), "\n")
print(round(sqrt(diag(solve(-H))) / se - 1, 6))
## B. Robust sandwich for a growth fit, by hand from per-person numerical gradients
N <- 50
dg <- do.call(rbind, lapply(1:N, function(i) { k <- rbinom(1, 1, .5); t <- 0:4; b <- rnorm(1, 0, .6); data.frame(id = i, t = t, y = c(0, 2)[k + 1] + c(.8, -.4)[k + 1] * t + b + rnorm(5, sd = .5)) }))
fr <- mixture_regression(y ~ t, dg, 2, id = "id", class_level = "group", random = "intercept", seed = 1, n_starts = 2, vcov_type = "robust", tol = 1e-10)
fo <- mixture_regression(y ~ t, dg, 2, id = "id", class_level = "group", random = "intercept", seed = 1, n_starts = 2, tol = 1e-10)
th <- coef(fr); p <- fr$params
person_ll <- function(t, i) { g <- .growth_unpack(fr$spec, t, p); G <- .growth_covariances(fr$spec, g); a <- dg[dg$id == i, ]; Z <- matrix(1, nrow(a))
  lse(log(softmax(c(g$gamma))) + vapply(1:2, function(k) dmvn(a$y, drop(cbind(1, a$t) %*% g$beta[, k]), Z %*% G[[k]] %*% t(Z) + diag(g$sigma2[k], nrow(a))), 0)) }
S1 <- t(vapply(1:N, function(i) numDeriv::grad(person_ll, th, i = i), numeric(length(th))))
bread <- vcov(fo); hand <- bread %*% (crossprod(S1) * N / (N - 1)) %*% bread
cat("growth robust SE vs hand sandwich max rel", max(abs(sqrt(diag(vcov(fr))) / sqrt(diag(hand)) - 1)), "\n")
print(fr$inference$vcov_type); print(get_results(fr, "fit")[c("bic")])
## C. Two-class binomial counts, NB with covariates: hand LL
dc <- data.frame(x = rnorm(400)); k <- rbinom(400, 1, .4)
dc$s <- rbinom(400, 10, plogis(ifelse(k == 1, 1 + dc$x, -1 - .5 * dc$x))); dc$f <- 10 - dc$s
fb <- mixture_regression(cbind(s, f) ~ x, dc, 2, family = "binomial", seed = 1, n_starts = 3, tol = 1e-10); q <- fb$params
hb <- sum(vapply(1:400, function(i) lse(log(softmax(c(q$gamma))) + dbinom(dc$s[i], 10, plogis(c(1, dc$x[i]) %*% q$beta), log = TRUE)), 0))
cat(sprintf("binomial counts LL %.9f hand %.9f; k %d\n", fb$log_likelihood, hb, fb$n_parameters)); print(as.data.frame(fb)[c("class","term","estimate","exp_estimate","p_value","p_adjusted")])
s <- as.data.frame(fb); ok <- s$term != "(Intercept)"; cat("BH non-intercepts:", isTRUE(all.equal(p.adjust(s$p_value[ok], "BH"), s$p_adjusted[ok])), " intercept NA:", all(is.na(s$p_adjusted[!ok])), "\n")
dc$nb <- rnbinom(400, size = ifelse(k == 1, 5, 1.5), mu = exp(ifelse(k == 1, 2 + .3 * dc$x, .5 - .4 * dc$x)))
fn <- mixture_regression(nb ~ x, dc, 2, family = "negative_binomial", seed = 1, n_starts = 3, tol = 1e-10); q <- fn$params
hn <- sum(vapply(1:400, function(i) lse(log(softmax(c(q$gamma))) + dnbinom(dc$nb[i], size = 1 / q$sigma2, mu = exp(c(1, dc$x[i]) %*% q$beta), log = TRUE)), 0))
cat(sprintf("NB two-class LL %.9f hand %.9f; k %d\n", fn$log_likelihood, hn, fn$n_parameters)); print(get_results(fn, "classes"))
## D. A column called `class` as id or cluster
dcl <- dg; names(dcl)[names(dcl) == "id"] <- "class"; dcl$class <- dcl$class + 500
fc <- mixture_regression(y ~ t, dcl, 2, id = "class", class_level = "group", seed = 1, n_starts = 2)
cat("LCGA, id = 'class', groups table:\n"); print(head(as.data.frame(get_results(fc, "groups")), 3)); cat(" names:", names(get_results(fc, "groups")), "\n"); print(head(as.data.frame(get_results(fc, "assignments")), 3))
fg <- mixture_regression(y ~ t, dcl, 2, id = "class", class_level = "group", random = "intercept", seed = 1, n_starts = 2)
cat("growth, id = 'class', assignments:\n"); a <- get_results(fg, "assignments"); cat(" names:", names(a), "\n"); print(head(as.data.frame(a), 3)); print(head(a, 3))
print(head(as.data.frame(get_results(fg, "random_effects")), 3)); print(head(as.data.frame(get_results(fg, "individual")), 3))
dsc <- d; names(dsc)[names(dsc) == "school"] <- "class"
fgc <- suppressWarnings(mixture_regression(y ~ t, dsc, 2, id = "id", class_level = "group", random = "intercept", cluster = "class", n_group_classes = 2, random_covariance = "equal", seed = 1, n_starts = 2))
a <- get_results(fgc, "assignments"); cat("cluster = 'class' names:", names(a), "\n"); print(head(a, 3)); cat(" a$class[1:3]:", a$class[1:3], "\n"); print(table(a$class)[1:3])
## E. Factor level present only in weight-zero rows in membership (the N-row guard)
dw <- dg; dw$sw <- 1; dw$grp <- factor(ifelse(dw$id %% 2 == 0, "a", "b"), levels = c("a", "b", "c")); dw$grp[dw$id %in% 1:4] <- "c"; dw$sw[dw$id %in% 1:4] <- 0
fw <- tryCatch(withCallingHandlers(mixture_regression(y ~ t, dw, 2, id = "id", class_level = "group", membership = ~ grp, weights = "sw", seed = 1, n_starts = 2),
  warning = function(w) { cat("  warning", class(w)[1], ":", substr(conditionMessage(w), 1, 90), "\n"); invokeRestart("muffleWarning") }), error = function(e) { cat("  error", class(e)[1], ":", substr(conditionMessage(e), 1, 120), "\n"); NULL })
fw2 <- mixture_regression(y ~ t, dw, 2, id = "id", class_level = "group", membership = ~ grp, weights = "sw", seed = 1, n_starts = 2, vcov_type = "none")
print(colnames(fw2$spec$w)); cat("k", fw2$n_parameters, " spec$n", fw2$spec$n, "\n")
dw3 <- subset(dw, sw > 0); dw3$grp <- factor(dw3$grp); fw3 <- mixture_regression(y ~ t, dw3, 2, id = "id", class_level = "group", membership = ~ grp, seed = 1, n_starts = 2, vcov_type = "none")
cat("same model without those rows: k", fw3$n_parameters, " LL", fw3$log_likelihood, " vs", fw2$log_likelihood, " BIC", get_results(fw3, "fit")$bic, get_results(fw2, "fit")$bic, "\n")
## F. Recovery with data in another row order (row level)
dr <- data.frame(x = rnorm(300)); kk <- rbinom(300, 1, .5); dr$truth <- kk; dr$y <- ifelse(kk == 1, 3 + dr$x, -3 - dr$x) + rnorm(300, sd = .3)
fx <- mixture_regression(y ~ x, dr, 2, seed = 1, n_starts = 2)
print(as.data.frame(get_results(fx, "recovery", data = dr, truth = "truth"))); print(as.data.frame(get_results(fx, "recovery", data = dr[order(dr$x), ], truth = "truth")))
## G. Class order by rows in a group nesting: sizes in rows vs share of groups
du <- do.call(rbind, lapply(1:60, function(i) { k <- as.integer(i > 40); m <- if (k == 1) 12 else 2; data.frame(id = i, t = seq_len(m), y = c(0, 6)[k + 1] + rnorm(m, sd = .4)) }))
fu <- mixture_regression(y ~ t, du, 2, id = "id", class_level = "group", seed = 1, n_starts = 2)
print(get_results(fu, "classes")[c("class", "share", "n_assigned")]); print(get_results(fu, "membership"))
