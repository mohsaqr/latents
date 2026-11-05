suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
# 1. Tiny case: direct numerical integration of the group intercept
d <- data.frame(g = c(1, 2, 2, 3, 3, 3, 3), y = c(0.3, -1.2, 0.4, 2.1, 1.7, 3.0, 2.2),
                z = c(1.0, 0.2, 0.7, -0.5, -1.1, 0.1, 0.6))
stats <- .additive_prepare(d, c("y", "z"), "g")
par <- list(means = rbind(c(-0.5, 0.2), c(1.8, -0.3)), between = rbind(c(0.7, 0.2), c(0.3, 1.1)),
            within = rbind(c(0.9, 0.4), c(0.5, 1.3)), weights = c(0.35, 0.65))
direct <- function(values, mu, tt, w) {
  f <- function(b) vapply(b, function(u) prod(dnorm(values, u, sqrt(w))) * dnorm(u, mu, sqrt(tt)), 0)
  integrate(f, mu - 12 * sqrt(tt), mu + 12 * sqrt(tt), rel.tol = 1e-11, subdivisions = 500L)$value
}
groups <- split(d, d$g)
hand <- t(vapply(groups, function(gr) vapply(1:2, function(h) {
  log(direct(gr$y, par$means[h, 1], par$between[h, 1], par$within[h, 1])) +
    log(direct(gr$z, par$means[h, 2], par$between[h, 2], par$within[h, 2]))
}, 0), c(0, 0)))
engine <- .additive_log_density(stats, par)
cat("max |quadrature - closed form| log density:", max(abs(hand - engine)), "\n")
print(unname(engine))
ll_hand <- sum(apply(sweep(hand, 2, log(par$weights), "+"), 1, function(v) log(sum(exp(v)))))
cat("LL hand", ll_hand, "engine", .additive_expectation(stats, par)$log_likelihood, "\n")
# conditional intercept moment by integration (group 3, class 2, indicator y)
gr <- groups[[3]]; mu <- par$means[2, 1]; tt <- par$between[2, 1]; w <- par$within[2, 1]
num <- integrate(function(b) b * vapply(b, function(u) prod(dnorm(gr$y, u, sqrt(w))) * dnorm(u, mu, sqrt(tt)), 0), mu - 12*sqrt(tt), mu + 12*sqrt(tt), rel.tol = 1e-11)$value
den <- direct(gr$y, mu, tt, w)
m2 <- integrate(function(b) b^2 * vapply(b, function(u) prod(dnorm(gr$y, u, sqrt(w))) * dnorm(u, mu, sqrt(tt)), 0), mu - 12*sqrt(tt), mu + 12*sqrt(tt), rel.tol = 1e-11)$value
ex <- .additive_expectation(stats, par)
cat("cond mean quad", num/den, "engine", ex$conditional_mean[3, 1, 2], "\n")
cat("cond var quad", m2/den - (num/den)^2, "engine", ex$conditional_variance[3, 1, 2], "\n")

# 2. A real fit
draw <- function(J = 150, seed = 3, sizes = c(2, 5, 10), w = c(.5, .5), sep = 1) {
  set.seed(seed)
  cl <- sample(seq_along(w), J, TRUE, w)
  n <- sample(sizes, J, TRUE)
  do.call(rbind, lapply(seq_len(J), function(j) {
    mu <- c(-sep, sep)[cl[j]]
    b <- rnorm(2, c(mu, -mu), sqrt(c(.3, .4)))
    data.frame(team = paste0("t", j), type = cl[j], y1 = rnorm(n[j], b[1], 1), y2 = rnorm(n[j], b[2], sqrt(1.5)))
  }))
}
D <- draw()
fit <- multilpa(D, c("y1", "y2"), "team", n_group_classes = 2, family = "additive", seed = 1, n_starts = 3)
print(fit)
print(get_results(fit, "group_classes"))
print(get_results(fit, "parameters"))
# hand likelihood from raw ratings via dense covariance
dense <- function(fit, D) {
  ids <- unique(D$team); ll <- 0
  vapply(ids, function(id) {
    g <- D[D$team == id, ]
    comps <- vapply(1:2, function(h) {
      sum(vapply(1:2, function(r) {
        v <- g[[c("y1","y2")[r]]]; n <- length(v)
        S <- diag(fit$within_variances[h, r], n) + fit$between_variances[h, r]
        ch <- chol(S); e <- backsolve(ch, v - fit$means[h, r], transpose = TRUE)
        -0.5*n*log(2*pi) - sum(log(diag(ch))) - 0.5*sum(e^2)
      }, 0)) + log(fit$group_probabilities[h])
    }, 0)
    m <- max(comps); m + log(sum(exp(comps - m)))
  }, 0)
}
cat("reported LL", fit$log_likelihood, "dense", sum(dense(fit, D)), "\n")
# Effective groups formula and the doc claim
inf <- .additive_inference(fit, "observed", warn_weak = FALSE)
gc <- .additive_effective_groups(fit, inf)
print(gc)
vw <- vcov(fit)[c("weight.group_class_1", "weight.group_class_2"), c("weight.group_class_1", "weight.group_class_2")]
print(vw)
w <- fit$group_probabilities
cat("hand w^2(1-w)/V:", w^2 * (1 - w) / diag(vw), "\n")
cat("perfectly classified groups giving same Var(w) = w(1-w)/V:", w * (1 - w) / diag(vw), "\n")
cat("J w (expected class count):", 150 * w, "\n")
# Logit & log CI by hand
tab <- get_results(fit, "parameters")
wr <- subset(tab, parameter == "weight"); vr <- subset(tab, parameter == "variance")
se <- wr$standard_error; est <- wr$estimate
cat("logit CI hand", plogis(qlogis(est) - qnorm(.975)*se/(est*(1-est))), plogis(qlogis(est) + qnorm(.975)*se/(est*(1-est))), "\nreported low", wr$conf_low, "high", wr$conf_high, "\n")
cat("log CI hand", head(vr$estimate*exp(-qnorm(.975)*vr$standard_error/vr$estimate),3), "reported", head(vr$conf_low,3), "\n")
# independent SE of a variance via profile-free Hessian on log scale: reported log-SE vs vcov unconstrained
vu <- vcov(fit, scale = "unconstrained")
cat("log-scale se vs delta:", sqrt(diag(vu))["log_between_variance.group_class_1.y1"],
    subset(tab, level == "between" & parameter == "variance" & group_class == "group_class_1" & indicator == "y1")$standard_error /
    subset(tab, level == "between" & parameter == "variance" & group_class == "group_class_1" & indicator == "y1")$estimate, "\n")
# conditional intercept table vs formula, row alignment after canonical reorder
ci <- get_results(fit, "conditional_intercepts")
pick <- subset(ci, team == "t7" & indicator == "y2" & group_class == "group_class_2")
g <- D[D$team == "t7", ]; n <- nrow(g)
T2 <- fit$between_variances["group_class_2", "y2"]; W2 <- fit$within_variances["group_class_2", "y2"]; m2 <- fit$means["group_class_2", "y2"]
cat("t7 cond mean table", pick$intercept_mean, "hand", m2 + n*T2/(W2 + n*T2) * (mean(g$y2) - m2), "; var table", pick$intercept_variance, "hand", T2*W2/(W2+n*T2), "\n")

# 3. Weights 2 vs duplicated groups, additive
dupl <- D[D$team %in% paste0("t", 1:50), ]
dupl$team <- paste0(dupl$team, "_d")
D2 <- rbind(D, dupl)
Dw <- D; Dw$w <- ifelse(Dw$team %in% paste0("t", 1:50), 2, 1)
f_dup <- multilpa(D2, c("y1", "y2"), "team", n_group_classes = 2, family = "additive", seed = 1, n_starts = 3, tol = 1e-10)
f_w <- multilpa(Dw, c("y1", "y2"), "team", n_group_classes = 2, family = "additive", seed = 1, n_starts = 3, tol = 1e-10, weights = "w")
cat("additive: means diff", max(abs(f_dup$means - f_w$means)), "between", max(abs(f_dup$between_variances - f_w$between_variances)), "within", max(abs(f_dup$within_variances - f_w$within_variances)), "weights", max(abs(f_dup$group_probabilities - f_w$group_probabilities)), "\n")
cat("LL scale: dup", f_dup$log_likelihood, "weighted", f_w$log_likelihood, "weighted*200/150", f_w$log_likelihood * 200 / 150, "\n")
cat("BIC groups count weighted fit:", (f_w$bic + 2 * f_w$log_likelihood) / f_w$n_parameters, "log(150)", log(150), "\n")

# 4. Label permutation of the likelihood + unit rescale
perm <- list(means = par$means[2:1, ], between = par$between[2:1, ], within = par$within[2:1, ], weights = par$weights[2:1])
cat("perm LL diff", .additive_expectation(stats, par)$log_likelihood - .additive_expectation(stats, perm)$log_likelihood, "\n")
Ds <- D; Ds$y1 <- D$y1 * 10; Ds$y2 <- D$y2 * 10
fs <- multilpa(Ds, c("y1", "y2"), "team", n_group_classes = 2, family = "additive", seed = 1, n_starts = 3)
cat("rescale LL shift", fit$log_likelihood - fs$log_likelihood, "theory", nrow(D) * 2 * log(10), "; eff groups", get_results(fs, "group_classes")$effective_groups, "vs", gc$effective_groups, "\n")

# 5. Comparison with a flat fit: T = 0, one class
one <- multilpa(D, c("y1", "y2"), "team", n_group_classes = 1, family = "additive", seed = 1, n_starts = 1)
flat <- suppressMessages(multilpa(D, c("y1", "y2"), NULL, n_profiles = 1, seed = 1, n_starts = 1))
cat("one-class additive LL", one$log_likelihood, "flat single-profile LL", flat$log_likelihood, "params", one$n_parameters, flat$n_parameters, "bic", one$bic, flat$bic, "bic_individual", one$bic_individual, "\n")

# 6. bootstrap_lrt on a weighted additive
print(tryCatch(bootstrap_lrt(multilpa(Dw, c("y1", "y2"), "team", n_group_classes = 1, family = "additive", seed = 1, n_starts = 1, weights = "w"), f_w, iter = 2, n_starts = 1, seed = 1), error = function(e) class(e)))
sessionInfo()$R.version$version.string
