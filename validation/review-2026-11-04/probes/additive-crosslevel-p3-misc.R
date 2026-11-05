suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
vars <- c("y1", "y2")
# a. One class additive = random-intercept ML per indicator (independent: nlme)
set.seed(11)
J <- 80; n <- sample(2:9, J, TRUE)
D <- do.call(rbind, lapply(seq_len(J), function(j) { b <- rnorm(2, c(1, -2), sqrt(c(.5, .2)))
  data.frame(team = paste0("t", j), y1 = rnorm(n[j], b[1], 1), y2 = rnorm(n[j], b[2], .7)) }))
one <- multilpa(D, vars, "team", n_group_classes = 1, family = "additive", seed = 1, n_starts = 1, tol = 1e-10)
if (requireNamespace("nlme", quietly = TRUE)) {
  fits <- lapply(vars, function(v) nlme::lme(stats::as.formula(paste(v, "~ 1")), random = ~ 1 | team, data = D, method = "ML", control = nlme::lmeControl(opt = "optim", msMaxIter = 500)))
  vc <- lapply(fits, function(f) as.numeric(nlme::VarCorr(f)[, "Variance"]))
  cat("nlme LL sum", sum(vapply(fits, function(f) as.numeric(logLik(f)), 0)), "additive LL", one$log_likelihood, "\n")
  cat("mean nlme", vapply(fits, function(f) nlme::fixef(f)[[1]], 0), "additive", one$means, "\n")
  cat("T nlme", vapply(vc, `[`, 0, 1), "additive", one$between_variances, "\n")
  cat("W nlme", vapply(vc, `[`, 0, 2), "additive", one$within_variances, "\n")
  cat("SE(mean) nlme", vapply(fits, function(f) sqrt(diag(vcov(f))), 0), "additive", subset(get_results(one, "parameters"), parameter == "mean")$standard_error, "\n")
}
# b. No higher point: optim from the estimate, dispersion and additive_dispersion equal (dense-free base R on stats)
draw2 <- function(seed) { set.seed(seed); J <- 120; cl <- sample(1:2, J, TRUE); n <- sample(3:8, J, TRUE)
  do.call(rbind, lapply(seq_len(J), function(j) { b <- rnorm(2, c(0, 0) + c(-1, 1)[cl[j]], sqrt(.3)); s <- c(.7, 1.6)[cl[j]]
    data.frame(team = paste0("t", j), y1 = rnorm(n[j], b[1], s), y2 = rnorm(n[j], b[2], s)) })) }
E <- draw2(21)
indep_ll <- function(m, b, w, p, E) sum(vapply(split(E, E$team), function(g) { n <- nrow(g)
  comp <- vapply(seq_along(p), function(h) log(p[h]) + sum(vapply(1:2, function(r) { v <- g[[vars[r]]]; a <- mean(v); S <- sum((v - a)^2)
    -0.5 * (n * log(2 * pi) + (n - 1) * log(w[h, r]) + log(w[h, r] + n * b[h, r]) + S / w[h, r] + n * (a - m[h, r])^2 / (w[h, r] + n * b[h, r])) }, 0)), 0)
  mx <- max(comp); mx + log(sum(exp(comp - mx))) }, 0))
chk <- function(fam, bv) {
  f <- multilpa(E, vars, "team", n_group_classes = 2, family = fam, between_variance = bv, seed = 1, n_starts = 3, tol = 1e-9)
  obj <- function(th) { p <- .additive_unpack(th, 2L, 2L, f$structure); -indep_ll(p$means, p$between, p$within, p$weights, E) }
  th0 <- coef(f, scale = "unconstrained")
  op <- optim(th0, obj, method = "BFGS", control = list(maxit = 60, reltol = 1e-12))
  cat(fam, bv, "reported LL", f$log_likelihood, "independent at estimate", -obj(th0), "after BFGS from estimate", -op$value, "n_par", f$n_parameters, length(th0), "\n")
  f
}
chk("additive", "equal"); chk("additive", "varying"); dis <- chk("dispersion", "equal"); chk("additive_dispersion", "varying"); chk("additive_dispersion", "equal")

# c. effective_groups: which covariance on which table
set.seed(31); J <- 160; cl <- sample(1:2, J, TRUE, c(.6, .4)); n <- sample(4:9, J, TRUE)
H <- do.call(rbind, lapply(seq_len(J), function(j) { m <- c(-.7, .7)[cl[j]]; b <- rnorm(2, m, .4)
  data.frame(team = paste0("t", j), y1 = b[1] + rt(n[j], 3), y2 = b[2] + rt(n[j], 3)) }))
fh <- multilpa(H, vars, "team", n_group_classes = 2, family = "additive", between_variance = "equal", seed = 1, n_starts = 3)
obs <- .additive_effective_groups(fh, .additive_inference(fh, "observed", warn_weak = FALSE))
rob <- .additive_effective_groups(fh, .additive_inference(fh, "robust", warn_weak = FALSE))
cat("observed eff", obs$effective_groups, "robust eff", rob$effective_groups, "\n")
w <- character()
tab <- withCallingHandlers(get_results(fh, "group_classes"), warning = function(c) { w <<- c(w, class(c)[1]); invokeRestart("muffleWarning") })
cat("group_classes$effective_groups", tab$effective_groups, "warnings", w, "\n")
w <- character()
s <- withCallingHandlers(summary(fh, vcov_type = "robust"), warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
cat("summary(vcov_type = robust)$group_classes eff:", s$group_classes$effective_groups, "\n warn:", w, "\n")

# d. seed = NULL grid: result depends on the RNG state; cores = 1 path uses the caller's stream
set.seed(41); J <- 90; cl <- sample(1:3, J, TRUE); n <- sample(2:5, J, TRUE)
G <- do.call(rbind, lapply(seq_len(J), function(j) { m <- c(-.6, 0, .6)[cl[j]]; b <- rnorm(2, m, .5)
  data.frame(team = paste0("t", j), y1 = rnorm(n[j], b[1], 1), y2 = rnorm(n[j], b[2], 1)) }))
grid <- function(state) { set.seed(state); suppressWarnings(as.data.frame(enumerate_classes(G, vars, "team", family = "additive", between_variance = "varying", n_group_classes = 3:4, n_starts = 2, seed = NULL)))$log_likelihood }
cat("grid LL, set.seed(1):", grid(1), "\n again set.seed(1):", grid(1), "\n set.seed(2):", grid(2), "\n set.seed(3):", grid(3), "\n set.seed(4):", grid(4), "\n")

# e. Worker loss with mocks (no workers): the retry result and on_lost paths
lost_once <- local({ calls <- 0L; function(indices, collected, cores, retry = FALSE) { calls <<- calls + 1L
  lapply(seq_along(indices), function(i) if (!retry && i == 2L) NULL else collected(indices[[i]])) } })
res <- testthat::with_mocked_bindings(
  withCallingHandlers(.latents_map_candidates(1:4, function(i) { set.seed(100 + i); i + runif(1) }, 2L, backend = "fork"),
    latents_worker_lost = function(c) { cat("WARN:", conditionMessage(c), "\n"); invokeRestart("muffleWarning") }),
  .latents_fork_pass = lost_once)
cat("recovered results", unlist(res), "; serial", unlist(lapply(1:4, function(i) { set.seed(100 + i); i + runif(1) })), "\n")
always_lost <- function(indices, collected, cores, retry = FALSE) lapply(seq_along(indices), function(i) if (i == 2L || retry) NULL else collected(indices[[i]]))
gridlost <- testthat::with_mocked_bindings(
  withCallingHandlers(enumerate_classes(G, vars, "team", family = "additive", between_variance = "varying", n_group_classes = 2:4, n_starts = 1, seed = 1, cores = 2, verbose = FALSE),
    warning = function(c) { cat("WARN[", class(c)[1], "]:", substr(conditionMessage(c), 1, 100), "\n"); invokeRestart("muffleWarning") }),
  .latents_fork_pass = always_lost, .latents_socket_pass = always_lost)
print(as.data.frame(gridlost)[, c("n_group_classes", "log_likelihood", "converged", "error")])
print(tryCatch(get_best(gridlost), error = function(e) conditionMessage(e)))
# sensitivity under a lost worker
sfit <- lpa(data.frame(y1 = G$y1, y2 = G$y2), vars, n_profiles = 2, n_starts = 1, seed = 1)
always_lost1 <- function(indices, collected, cores, retry = FALSE) lapply(seq_along(indices), function(i) if (i == 2L) NULL else collected(indices[[i]]))
print(tryCatch(testthat::with_mocked_bindings(
  suppressWarnings(sensitivity(sfit, seeds = 1:3, n_starts = 1, cores = 2)),
  .latents_fork_pass = always_lost1, .latents_socket_pass = always_lost1), error = function(e) paste("ERROR", class(e)[1], conditionMessage(e))))

# f. Selection under ties: deterministic, not order of arrival
ll <- c(-100, -100 + 1e-12, -100 - 1e-5, -99.9999)
print(.multilpa_select_start(ll, c(TRUE, TRUE, TRUE, FALSE), "likelihood", tie = 1e-3))
print(.multilpa_select_start(ll, c(TRUE, TRUE, TRUE, FALSE), "likelihood", tie = 1e-6))
print(.multilpa_select_start(c(-Inf, NaN, -5, -5), c(FALSE, FALSE, FALSE, TRUE), "likelihood", tie = 0))
# g. cross-level: group classes beyond the positive-weight groups
Z <- G; Z$w <- ifelse(Z$team %in% paste0("t", 1:88), 0, 1)
print(tryCatch(suppressWarnings(multilpa(Z, vars, "team", n_profiles = 1, n_group_classes = 3, family = "restricted_cross_level", n_starts = 1, seed = 1, weights = "w"))$group_probabilities, error = function(e) paste(class(e)[1], conditionMessage(e))))
print(tryCatch(suppressWarnings(multilpa(Z, vars, "team", n_group_classes = 3, family = "additive", n_starts = 1, seed = 1, weights = "w"))$group_probabilities, error = function(e) paste(class(e)[1], conditionMessage(e))))
