suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL)
q <- function(e) { w <- character(); v <- withCallingHandlers(tryCatch(e, error = function(x) structure(list(msg = conditionMessage(x), cls = class(x)[1]), class = "failed")),
  warning = function(c) { w <<- c(w, class(c)[1]); invokeRestart("muffleWarning") }); attr(v, "w") <- unique(w); v }
one <- function(label, code) { cat("\n=====", label, "\n"); r <- tryCatch(code, error = function(e) cat("SCRIPT ERROR:", conditionMessage(e), "\n")); invisible(r) }

one("C1b person centring, 12 seeds, k2 vs k3 picked; share of zero rows in a spike", {
  r <- lapply(1101:1112, function(s) { set.seed(s)
    J <- 80; sizes <- c(rep(1L, 12), sample(5:9, J - 12, TRUE)); g <- rep(seq_len(J), sizes); n <- length(g)
    z <- rbinom(n, 1, .5); off <- rep(rnorm(J, 0, 1.5), sizes)
    d <- data.frame(g = g, a = off + rnorm(n, c(-1, 1)[z + 1], .7), b = off + rnorm(n, c(1, -1)[z + 1], .7), c = rnorm(n, c(0, .5)[z + 1], 1))
    d2 <- subset(d, g %in% which(sizes > 1)); v <- c("a", "b", "c")
    fit <- function(D, k) q(multilpa(D, v, "g", n_profiles = k, n_group_classes = 1, centering = "person", n_starts = 5, seed = 1))
    a2 <- fit(d, 2); a3 <- fit(d, 3); b2 <- fit(d2, 2); b3 <- fit(d2, 3)
    data.frame(seed = s, with_k3_minus_k2 = round(a3$bic - a2$bic, 1), without_k3_minus_k2 = round(b3$bic - b2$bic, 1), k3_min_var = signif(min(a3$variances), 2), spike = min(a3$variances) < 1e-3, warned = paste(attr(a3, "w"), collapse = ",")) })
  print(do.call(rbind, r), row.names = FALSE)
})

one("R2b unused factor level, what failed", {
  set.seed(801); d <- do.call(rbind, lapply(1:60, function(i) { gr <- sample(c("a", "b"), 1); k <- rbinom(1, 1, ifelse(gr == "a", .7, .3)); data.frame(id = i, grp = gr, t = 0:3, y = c(0, 4)[k + 1] + c(.5, -.5)[k + 1] * (0:3) + rnorm(4, sd = .8)) }))
  d1 <- d; d1$grp <- factor(d1$grp, levels = c("a", "b", "zz")); d2 <- d; d2$grp <- factor(d2$grp)
  m <- function(D) q(mixture_regression(y ~ t, D, 2, id = "id", class_level = "group", membership = ~ grp, n_starts = 3, seed = 1))
  a <- m(d1); b <- m(d2)
  if (inherits(a, "failed")) cat("a FAILED", a$cls, a$msg, "\n")
  if (inherits(b, "failed")) cat("b FAILED", b$cls, b$msg, "\n")
  if (!inherits(a, "failed") && !inherits(b, "failed")) { cat(sprintf("unused: k %d LL %.4f BIC %.3f [%s]\ndropped: k %d LL %.4f BIC %.3f [%s]\n", a$n_parameters, a$log_likelihood, a$bic, paste(attr(a, "w"), collapse = ","), b$n_parameters, b$log_likelihood, b$bic, paste(attr(b, "w"), collapse = ",")))
    cat("names membership:", names(get_results(a, "membership")), "\n"); print(get_results(a, "membership")); print(get_results(b, "membership")) }
})

one("I3 staged bootstrap with a foreign measurement", {
  sim <- function(s, shift) { set.seed(s); G <- 60; m <- 8; h <- rep(rbinom(G, 1, .5), each = m); k <- rbinom(G * m, 1, c(.25, .75)[h + 1]) + 1L
    data.frame(g = rep(seq_len(G), each = m), y1 = rnorm(G * m, c(-1.5, 1.5)[k] + shift), y2 = rnorm(G * m, c(1, -1)[k])) }
  d <- sim(1, 0); v <- c("y1", "y2")
  for (cfg in list(list("other sample shifted by 3", sim(2, 3)), list("other sample, same population", sim(3, 0)), list("same data", d))) {
    meas <- q(multilpa(cfg[[2]], v, "g", n_profiles = 2, n_group_classes = 1, n_starts = 3, seed = 1))
    st <- q(fit_staged(d, v, "g", n_profiles = 2, n_group_classes = 2, measurement = meas, n_starts = 2, seed = 1))
    b <- q(parameter_inference(st, method = "bootstrap", iter = 30, n_starts = 2, seed = 1))
    if (inherits(b, "failed")) { cat(cfg[[1]], "FAILED", b$cls, b$msg, "\n"); next }
    out <- b$estimate < b$conf_low | b$estimate > b$conf_high
    cat(sprintf("%s: rows outside own interval %d of %d; warnings %s\n", cfg[[1]], sum(out, na.rm = TRUE), length(out), paste(attr(b, "w"), collapse = ","))) }
})

one("T1 state reached only at the last wave", {
  set.seed(1301); n <- 300; TT <- 4; v <- paste0("y", 1:5); mu <- rbind(0, c(2.5, 2.5, 0, 0, 0), 5)
  s <- matrix(1L, n, TT); s[, 1] <- sample(1:2, n, TRUE); for (t in 2:(TT - 1)) s[, t] <- ifelse(runif(n) < .2, 3L - s[, t - 1], s[, t - 1]); s[, TT] <- ifelse(runif(n) < .5, 3L, s[, TT - 1])
  d <- do.call(rbind, lapply(seq_len(TT), function(t) { y <- mu[s[, t], ] + matrix(rnorm(n * 5), n, 5); colnames(y) <- v; data.frame(id = seq_len(n), time = t, y) }))
  f <- q(lta(d, v, "id", 3, "time", n_starts = 5, seed = 1)); top <- which.max(rowSums(f$means)); tb <- get_results(f, "transitions")
  print(subset(tb, as.character(from) %in% c(top, paste0("profile_", top)))[, c("from", "to", "probability", "expected_count", "estimated")], row.names = FALSE); cat("n_parameters", f$n_parameters, "warnings:", attr(f, "w"), "\n")
})

one("T2 stayer class with two occasions", {
  set.seed(1401); n <- 300; P <- matrix(c(.8, .2, .3, .7), 2, byrow = TRUE)
  d <- do.call(rbind, lapply(seq_len(n), function(i) { s <- integer(2); s[1] <- sample(2, 1, prob = c(.6, .4)); s[2] <- sample(2, 1, prob = P[s[1], ]); data.frame(id = i, time = 1:2, y1 = rnorm(2, c(0, 3)[s]), y2 = rnorm(2, c(1, -2)[s], .8), y3 = rnorm(2, c(0, 2)[s])) }))
  v <- c("y1", "y2", "y3"); a <- q(lta(d, v, "id", 2, "time", n_starts = 6, seed = 1)); b <- q(lta(d, v, "id", 2, "time", mover_stayer = TRUE, n_starts = 6, seed = 1))
  cat(sprintf("plain LL %.6f k %d | stayer LL %.6f k %d [warnings: %s]\n", a$log_likelihood, a$n_parameters, b$log_likelihood, b$n_parameters, paste(attr(b, "w"), collapse = ",")))
  print(get_results(b, "group_classes")); w <- q(parameter_inference(b)); if (inherits(w, "failed")) cat("Wald refused", w$cls, "\n") else cat("largest Wald SE", signif(max(w$standard_error, na.rm = TRUE), 4), "\n")
})

one("T3 bootstrap with one longer person", {
  set.seed(1501); P <- matrix(c(.8, .2, .3, .7), 2, byrow = TRUE); path <- function(i, TT) { s <- integer(TT); s[1] <- sample(2, 1, prob = c(.6, .4)); for (t in 2:TT) s[t] <- sample(2, 1, prob = P[s[t - 1], ]); data.frame(id = i, time = seq_len(TT), y1 = rnorm(TT, c(0, 3)[s]), y2 = rnorm(TT, c(1, -2)[s], .8)) }
  d <- do.call(rbind, c(lapply(1:60, path, TT = 3), list(path(61, 5)))); o <- q(lta(d, c("y1", "y2"), "id", 2, "time", transitions = "occasion", n_starts = 2, seed = 1))
  b <- q(parameter_inference(o, method = "bootstrap", iter = 30, n_starts = 1, seed = 1)); if (inherits(b, "failed")) cat("FAILED:", b$cls, "|", b$msg, "\n") else cat("completed rows", nrow(b), "\n")
})

one("3S1 r3step at groups level near separation, robust vs observed", {
  r <- lapply(2113:2124, function(s) { set.seed(s); G <- 25; m <- 10; z <- rnorm(G); gc <- 1L + (runif(G) < plogis(4 * z)); g <- rep(seq_len(G), each = m); tr <- 1L + (runif(G * m) < c(.15, .85)[gc[g]])
    d <- data.frame(g = g, z = rep(z, each = m), a = rnorm(G * m, c(-1.5, 1.5)[tr]), b = rnorm(G * m, c(-1.5, 1.5)[tr])); f <- q(multilpa(d, c("a", "b"), "g", n_profiles = 2, n_group_classes = 2, n_starts = 3, seed = 1))
    rb <- q(r3step(f, d, "z", level = "groups")); ob <- q(r3step(f, d, "z", level = "groups", vcov_type = "observed")); if (inherits(rb, "failed") || inherits(ob, "failed")) return(data.frame(seed = s, slope = NA, rob_se = NA, rob_p = NA, obs_se = NA, obs_p = NA, warned = "refused"))
    k <- which(rb$term == "z")[1]; data.frame(seed = s, slope = round(rb$estimate[k], 2), rob_se = round(rb$standard_error[k], 2), rob_p = signif(rb$p_value[k], 2), obs_se = round(ob$standard_error[k], 2), obs_p = signif(ob$p_value[k], 2), warned = paste(attr(rb, "w"), collapse = ",")) })
  print(do.call(rbind, r), row.names = FALSE)
})

one("3S2 step_one estimated on staged vs joint fit", {
  set.seed(2201); G <- 100; m <- 10; n <- G * m; gc <- rep(1:2, length.out = G); g <- rep(seq_len(G), each = m); tr <- 1L + (runif(n) < c(.3, .7)[gc[g]])
  d <- data.frame(g = g, a = rnorm(n, c(-.8, .8)[tr]), b = rnorm(n, c(-.7, .7)[tr]), c = rnorm(n, c(.7, -.7)[tr]), y = rnorm(n, c(0, 1)[tr])); v <- c("a", "b", "c")
  full <- q(multilpa(d, v, "g", n_profiles = 2, n_group_classes = 2, n_starts = 5, seed = 1, tol = 1e-10)); st <- q(fit_staged(d, v, "g", n_profiles = 2, n_group_classes = 2, n_starts = 5, seed = 1))
  for (l in list(list("joint", full), list("staged", st))) { a <- q(three_step(l[[2]], d, "y", method = "bch")); b <- q(three_step(l[[2]], d, "y", method = "bch", step_one = "estimated")); cat(l[[1]], "SE ratio estimated/fixed:", signif(b$standard_error / a$standard_error, 4), " attr step_one:", attr(b, "step_one"), " warnings:", attr(b, "w"), "\n") }
})

one("A1 cross-level criteria on data with no group structure", {
  r <- lapply(3101:3106, function(s) { set.seed(s); J <- 100; n <- 6; k <- sample(1:2, J * n, TRUE); pm <- rbind(c(-1.5, -1), c(1.5, 1))
    D <- data.frame(team = rep(sprintf("t%03d", seq_len(J)), each = n), y1 = rnorm(J * n, pm[k, 1], .6), y2 = rnorm(J * n, pm[k, 2], .6)); v <- c("y1", "y2")
    fx <- function(fam, H) q(multilpa(D, v, "team", n_profiles = 2, n_group_classes = H, family = fam, seed = 1, n_starts = 3, variance_model = "equal", between_variance = "equal"))
    f1 <- fx("full_cross_level", 1); f2 <- fx("full_cross_level", 2); r2 <- fx("restricted_cross_level", 2); p1 <- q(multilpa(D, v, "team", n_profiles = 2, n_group_classes = 1, variance_model = "equal", seed = 1, n_starts = 3)); p2 <- q(multilpa(D, v, "team", n_profiles = 2, n_group_classes = 2, variance_model = "equal", seed = 1, n_starts = 3))
    data.frame(seed = s, full_H2_minus_H1 = round(f2$bic - f1$bic, 1), restr_H2_minus_H1 = round(r2$bic - f1$bic, 1), full_minus_restr_LL = round(f2$log_likelihood - r2$log_likelihood, 1), profile_family_H2_minus_H1 = round(p2$bic - p1$bic, 1)) })
  print(do.call(rbind, r), row.names = FALSE)
})

one("A2 prevalence, full vs restricted, sizes related to class", {
  set.seed(3201); J <- 80; h <- rep(1:2, length.out = J); sz <- c(5, 15)[h]; g <- rep(seq_len(J), sz); n <- length(g); k <- 1L + rbinom(n, 1, c(.2, .8)[h[g]]); pm <- rbind(c(-1.5, -1), c(1.5, 1))
  D <- data.frame(team = sprintf("t%03d", g), y1 = rnorm(n, pm[k, 1], .6), y2 = rnorm(n, pm[k, 2], .6))
  for (fam in c("full_cross_level", "restricted_cross_level")) { f <- q(multilpa(D, c("y1", "y2"), "team", n_profiles = 2, n_group_classes = 2, family = fam, seed = 1, n_starts = 3)); p <- get_results(f, "profiles"); cat(fam, "names:", names(p), "\n"); print(unique(subset(p, !is.na(prevalence))[, intersect(names(p), c("profile", "prevalence"))]), row.names = FALSE); cat("mean posterior over individuals:", round(colMeans(f$subject_posteriors), 3), "\n") }
  cat("true share of profile 2 among individuals:", round(mean(k == 2), 3), " in the average group:", round(mean(tapply(k == 2, g, mean)), 3), "\n")
})
