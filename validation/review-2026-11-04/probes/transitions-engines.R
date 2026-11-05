# Homogeneous and general engines on one specification; scores and vcov against
# finite differences of the plain likelihood; natural Jacobian numerically.
suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
sim <- function(n, T, seed, K = 2) {
  set.seed(seed)
  P <- matrix(c(.8,.2,.3,.7), 2, byrow = TRUE)
  do.call(rbind, lapply(seq_len(n), function(i) {
    s <- integer(T); s[1] <- sample(2, 1, prob = c(.6,.4))
    for (t in 2:T) s[t] <- sample(2, 1, prob = P[s[t-1], ])
    data.frame(id = i, time = seq_len(T), y1 = rnorm(T, c(0, 2)[s]), y2 = rnorm(T, c(1, -1)[s], .8), y3 = rnorm(T, c(0, 1.5)[s]),
               g = factor(rbinom(T, 1, c(.2, .8)[s]), levels = 0:1, labels = c("no", "yes")))
  }))
}
vars <- c("y1","y2","y3")
d <- sim(150, 4, 2)
d$unit <- 1
# 1. same spec: all-one weights send it to the general engine
h1 <- suppressWarnings(lta(d, vars, "id", 2, "time", n_starts = 5, seed = 1, tol = 1e-10))
g1 <- suppressWarnings(lta(d, vars, "id", 2, "time", n_starts = 5, seed = 1, tol = 1e-10, weights = "unit"))
cat(sprintf("H=1 LL hom %.7f gen %.7f; npar %d %d; BIC %.4f %.4f; BICind %.4f %.4f; nobs %d %d; BIC() %.4f %.4f\n",
  h1$log_likelihood, g1$log_likelihood, h1$n_parameters, g1$n_parameters, h1$bic, g1$bic, h1$bic_individual, g1$bic_individual, nobs(h1), nobs(g1), BIC(h1), BIC(g1)))
th <- subset(get_results(h1, "transitions"), TRUE)[, c("from","to","probability")]; tg <- get_results(g1, "transitions")
cat("max |transition diff| ", max(abs(sort(th$probability) - sort(tg$probability))), "\n")
# 2. T=2: occasion transitions = homogeneous; SEs on the logit scale
d2 <- subset(d, time <= 2)
h2 <- suppressWarnings(lta(d2, vars, "id", 2, "time", n_starts = 5, seed = 1))
g2 <- suppressWarnings(lta(d2, vars, "id", 2, "time", transitions = "occasion", n_starts = 5, seed = 1))
cat(sprintf("T=2 LL hom %.7f gen %.7f; npar %d %d; BIC %.4f %.4f\n", h2$log_likelihood, g2$log_likelihood, h2$n_parameters, g2$n_parameters, h2$bic, g2$bic))
if (cor(h2$means[,1], g2$measurement[[1]]$means[,1]) < 0) cat("labels flipped between engines\n")
vh <- vcov(h2, scale = "unconstrained"); vg <- vcov(g2)
print(round(sqrt(diag(vh)), 5)); print(round(sqrt(diag(vg)), 5))
# Natural-scale SE of the transition probabilities from the general logits (delta) vs homogeneous natural
pih <- parameter_inference(h2); print(subset(pih, parameter == "transition_probability")[, c("outcome","term","estimate","standard_error","conf_low","conf_high")])
pig <- parameter_inference(g2); print(subset(pig, block == "transition")[, c("from","to","estimate","standard_error")])
ps <- subset(pig, block == "transition")
pstay <- plogis(ps$estimate)
cat("delta-method SE of move probability from general logit:", round(pstay * (1 - pstay) * ps$standard_error, 6), "\n")
# 3. Homogeneous analytic score vs finite differences of the likelihood, several specs
fd_check <- function(label, f) {
  x <- f$indicator_data; if (is.null(x)) x <- matrix(numeric(0), f$n_observations, 0L)
  x <- matrix(as.numeric(x), nrow(x), ncol(x))
  centers <- if (ncol(x)) colMeans(x, na.rm = TRUE) else numeric(0)
  cx <- if (ncol(x)) sweep(x, 2L, centers, "-") else x
  layout <- .multilpa_sequence_layout(f$group_index, f$time_values, f$n_groups, f$occasions)
  view <- .multilpa_transition_view(f, centers)
  b <- .multilpa_transition_shared_blocks(view, "unconstrained")
  w <- c(measurement = length(b$measurement), group = length(b$group))
  th <- unname(.multilpa_transition_encode(f, "unconstrained", centers))
  ll <- function(v) .multilpa_transition_expectation(cx, layout, .multilpa_transition_decode(v, view, w[["measurement"]], w[["group"]]), f$categorical_data)$log_likelihood
  an <- colSums(.multilpa_transition_group_scores(th, cx, layout, view, w, f$categorical_data))
  fd <- vapply(seq_along(th), function(j) { e <- replace(numeric(length(th)), j, 1e-5); (ll(th + e) - ll(th - e)) / 2e-5 }, 0)
  # numeric Hessian from FD of likelihood alone
  np <- length(th); H <- matrix(0, np, np); h <- 1e-3
  scale <- rep(1, np)
  ll0 <- ll(th)
  for (a in seq_len(np)) for (c in seq_len(a)) { # probe-only double loop
    ea <- replace(numeric(np), a, h); ec <- replace(numeric(np), c, h)
    H[a, c] <- H[c, a] <- (ll(th+ea+ec) - ll(th+ea-ec) - ll(th-ea+ec) + ll(th-ea-ec)) / (4*h*h)
  }
  se_num <- sqrt(diag(solve(-H)))
  se_pkg <- sqrt(diag(vcov(f, scale = "unconstrained")))
  # Natural Jacobian
  J <- .multilpa_transition_jacobian(f, view)
  natural_of <- function(v) {
    p <- .multilpa_transition_decode(v, view, w[["measurement"]], w[["group"]])
    fp <- f; fp$means <- p$means; fp$variances <- p$variances; fp$covariances <- p$covariances
    fp$response_probabilities <- p$response_probabilities; fp$group_probabilities <- p$group_probabilities
    fp$initial_probabilities <- p$initial_probabilities; fp$transition_probabilities <- p$transition_probabilities
    unname(.multilpa_transition_encode(fp, "natural"))
  }
  Jn <- vapply(seq_along(th), function(j) { e <- replace(numeric(np), j, 1e-6); (natural_of(th + e) - natural_of(th - e)) / 2e-6 }, numeric(nrow(J)))
  cat(sprintf("%-28s np %2d max|score-fd| %.2e  max rel |SE pkg - SE fd Hessian| %.2e  max|J-Jfd| %.2e\n", label, np, max(abs(an - fd)), max(abs(se_pkg/se_num - 1)), max(abs(J - Jn))))
}
dsmall <- sim(80, 3, 5)
fits <- list(
  "continuous H=1" = lta(dsmall, vars, "id", 2, "time", n_starts = 2, seed = 1),
  "continuous H=2" = lta(sim(120, 5, 7), vars, "id", 2, "time", n_group_classes = 2, n_starts = 3, seed = 1),
  "equal variance" = lta(dsmall, vars, "id", 2, "time", variance_model = "equal", n_starts = 2, seed = 1),
  "full covariance" = lta(dsmall, c("y1","y2"), "id", 2, "time", covariance_model = "full", n_starts = 2, seed = 1),
  "mixed categorical" = lta(dsmall, c("y1","y2","g"), "id", 2, "time", categorical = "g", n_starts = 2, seed = 1))
ms <- dsmall; ms$y1[seq(1, nrow(ms), 7)] <- NA; ms$y2[seq(2, nrow(ms), 7)] <- NA
fits[["FIML"]] <- lta(ms, vars, "id", 2, "time", missing = "fiml", n_starts = 2, seed = 1)
gp <- dsmall[!(dsmall$id %% 4 == 0 & dsmall$time == 2), ]
fits[["grid gaps"]] <- lta(gp, vars, "id", 2, "time", occasions = "grid", n_starts = 2, seed = 1)
invisible(Map(fd_check, names(fits), fits))
# 4. natural CI bounds of transition probabilities inside [0,1]
p <- parameter_inference(fits[[1]])
print(subset(p, parameter %in% c("transition_probability", "initial_probability", "probability", "group_probability")))
# 5. Weights: integer weight 2 against duplication, general engine
dup <- d; dup2 <- d; dup2$id <- dup2$id + 1000L; both <- rbind(dup, dup2); both$unit <- 1
d$wt <- 2
a <- suppressWarnings(lta(d, vars, "id", 2, "time", weights = "wt", n_starts = 3, seed = 1))
b <- suppressWarnings(lta(both, vars, "id", 2, "time", weights = "unit", n_starts = 3, seed = 1))
cat(sprintf("weights 2 vs duplicates: LL %.6f %.6f; n_par %d %d; BIC %.4f %.4f (log(150) vs log(300) convention); nobs %s %s\n", a$log_likelihood, b$log_likelihood, a$n_parameters, b$n_parameters, a$bic, b$bic, nobs(a), nobs(b)))
