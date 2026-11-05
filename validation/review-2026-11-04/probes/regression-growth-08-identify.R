suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(8)
warn <- function(expr) withCallingHandlers(expr, warning = function(w) { cat("  [warning ", class(w)[1], "] ", substr(conditionMessage(w), 1, 110), "\n", sep = ""); invokeRestart("muffleWarning") })
# A. Exact likelihood ridge of a random intercept + slope with two occasions
pp <- do.call(rbind, lapply(1:150, function(i) { b <- rnorm(1, 0, 1); data.frame(id = i, t = 0:1, y = b + .5 * (0:1) + rnorm(2, sd = .7)) }))
f <- warn(mixture_regression(y ~ t, pp, 1, id = "id", class_level = "group", random = ~ 1 + t, seed = 1))
p <- f$params; G <- p$random_covariance[[1]]; s2 <- p$sigma2
at <- function(delta) { q <- p; g <- G; g[1, 1] <- g[1, 1] - delta; g[1, 2] <- g[2, 1] <- g[1, 2] + delta; g[2, 2] <- g[2, 2] - 2 * delta; q$random_covariance <- list(g); q$sigma2 <- s2 + delta
  c(delta = delta, LL = .growth_expectation(f$spec, f$stats, q)$log_likelihood, g00 = g[1, 1], g01 = g[1, 2], g11 = g[2, 2], s2 = s2 + delta, min_eig = min(eigen(g)$values)) }
print(round(do.call(rbind, lapply(c(-.1, -.05, 0, .03), at)), 8))
cat("reported: LL", f$log_likelihood, " n_parameters", f$n_parameters, " converged", f$converged, " random_boundary", f$random_boundary, "\n")
I <- f$inference$information; e <- eigen(I / tcrossprod(sqrt(diag(I))), symmetric = TRUE, only.values = TRUE)$values; cat("information correlation-form eigenvalues:", signif(e, 3), "\n")
v0 <- warn(mixture_regression(y ~ t, pp, 1, id = "id", class_level = "group", random = ~ 1 + t, min_variance = 1e-4, seed = 1)); cat("other floor: residual_sd", sqrt(v0$params$sigma2), " vs", sqrt(s2), "; G11", v0$params$random_covariance[[1]][1, 1], "vs", G[1, 1], "\n")
y2 <- pp; y2$y <- y2$y * 1.0000001; v1 <- warn(mixture_regression(y ~ t, y2, 1, id = "id", class_level = "group", random = ~ 1 + t, seed = 1, max_iter = 200, tol = 1e-5)); cat("max_iter = 200, tol = 1e-5: residual_sd", sqrt(v1$params$sigma2), "; G11", v1$params$random_covariance[[1]][1, 1], " converged", v1$converged, "\n")
# And a realistic pre-post with two classes
pc <- do.call(rbind, lapply(1:200, function(i) { k <- rbinom(1, 1, .5); b <- rnorm(1, 0, 1); data.frame(id = i, t = 0:1, y = c(0, 3)[k + 1] + c(.5, -1)[k + 1] * (0:1) + b + rnorm(2, sd = .7)) }))
f2 <- warn(mixture_regression(y ~ t, pc, 2, id = "id", class_level = "group", random = ~ 1 + t, random_covariance = "equal", variance = "equal", seed = 1, n_starts = 3))
print(as.data.frame(get_results(f2, "random"))[c("class", "parameter", "term", "estimate", "std_error", "boundary")])
print(as.data.frame(get_results(f2, "classes"))[c("class", "share", "share_std_error", "residual_sd", "residual_sd_std_error")])
print(as.data.frame(f2)[c("class", "term", "estimate", "std_error")]); print(as.data.frame(get_results(f2, "fit"))[c("n_parameters", "bic", "converged", "random_boundary", "n_best_replicated")])
# B. Group nesting escapes the row-level identification refusals
bb <- data.frame(id = rep(1:150, each = 2)); pk <- rep(rbinom(150, 1, .4), each = 2); bb$b <- rbinom(300, 1, ifelse(pk == 1, .85, .2))
fb <- warn(mixture_regression(b ~ 1, bb, 2, family = "binomial", id = "id", class_level = "group", seed = 1, n_starts = 5))
print(as.data.frame(fb)[c("class", "term", "estimate", "std_error")]); print(as.data.frame(get_results(fb, "classes"))[c("class", "share", "prior", "prior_std_error")]); print(as.data.frame(get_results(fb, "starts"))[c("log_likelihood", "converged", "selected")]); cat("n_best_replicated", fb$n_best_replicated, " k", fb$n_parameters, " df available 2\n")
sg <- data.frame(id = 1:300, x = rnorm(300)); sg$b <- rbinom(300, 1, plogis(.5 * sg$x))
fs <- tryCatch(warn(mixture_regression(b ~ x, sg, 2, family = "binomial", id = "id", class_level = "group", seed = 1, n_starts = 3)), error = function(e) { cat("  refused:", class(e)[1], "\n"); NULL })
if (!is.null(fs)) { print(as.data.frame(fs)[c("class", "term", "estimate", "std_error")]); cat(" k", fs$n_parameters, "\n") }
fo <- tryCatch(warn(mixture_regression(o ~ 1, transform(sg, o = cut(x, 3, labels = FALSE)), 2, family = "ordinal", id = "id", class_level = "group", seed = 1, n_starts = 3)), error = function(e) { cat("  refused:", class(e)[1], "\n"); NULL })
if (!is.null(fo)) { print(as.data.frame(fo)[c("class", "term", "estimate", "std_error")]); cat(" k", fo$n_parameters, "\n") }
# C. Correlated-residual fits: is the reported maximum replicated?
dd <- subset(growth_scores, student <= 150)
runs <- lapply(1:4, function(s) { x <- suppressWarnings(mixture_regression(score ~ wave, dd, 3, id = "student", class_level = "group", random = "intercept", residual_correlation = "ar", random_covariance = "equal", variance = "equal", seed = s, n_starts = 3))
  st <- get_results(x, "starts"); c(seed = s, fit_LL = x$log_likelihood, n_best = x$n_best_replicated, start_LL = max(st$log_likelihood[st$selected]), conv = x$converged) })
print(do.call(rbind, runs))
