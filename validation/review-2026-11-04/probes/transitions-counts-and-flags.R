# n_parameters against the coordinates; FD scores for covariance structures; flags.
suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
sim <- function(n, T, seed) {
  set.seed(seed)
  P <- matrix(c(.8,.2,.3,.7), 2, byrow = TRUE)
  do.call(rbind, lapply(seq_len(n), function(i) {
    s <- integer(T); s[1] <- sample(2, 1, prob = c(.6,.4))
    for (t in 2:T) s[t] <- sample(2, 1, prob = P[s[t-1], ])
    data.frame(id = i, time = seq_len(T), y1 = rnorm(T, c(0, 2)[s]), y2 = rnorm(T, c(1, -1)[s], .8), y3 = rnorm(T, c(0, 1.5)[s]),
               z = rnorm(T), w = rep(rnorm(1), T), cat = factor(sample(c("a","b","c"), T, TRUE, prob = list(c(.6,.3,.1), c(.1,.3,.6))[[1]])),
               ord = sample(1:3, T, TRUE), cnt = rpois(T, c(1, 4)[s]))
  }))
}
d <- sim(120, 4, 2); d$ord <- ordered(d$ord)
vars <- c("y1","y2","y3")
args <- list(
  cov_init = list(vars = vars, initial_covariates = "w", transition_covariates = "z"),
  occasion = list(vars = vars, transitions = "occasion"),
  meas_occ = list(vars = vars, measurement = "occasion"),
  order2 = list(vars = vars, order = 2),
  order2_occ = list(vars = vars, order = 2, transitions = "occasion"),
  stayer = list(vars = vars, mover_stayer = TRUE),
  stayer_order2_occ_meas = list(vars = vars, mover_stayer = TRUE, order = 2, measurement = "occasion"),
  categorical = list(vars = c("y1","cat"), categorical = "cat", transitions = "occasion"),
  ordinal = list(vars = c("y1","ord"), ordinal = "ord"),
  count = list(vars = c("y1","cnt"), count = "cnt"),
  negbin = list(vars = c("y1","cnt"), count = "cnt", count_model = "negative_binomial"),
  negbin_eq_occ = list(vars = c("y1","cnt"), count = "cnt", count_model = "negative_binomial", count_dispersion = "equal", measurement = "occasion"),
  equal = list(vars = vars, variance_model = "equal", transitions = "occasion"),
  full = list(vars = vars, covariance_model = "full", transitions = "occasion"),
  H2_cov = list(vars = vars, n_group_classes = 2, transition_covariates = "z"))
structures <- .multilpa_structures()
args <- c(args, setNames(lapply(structures, function(m) list(vars = vars, model = m, transitions = "occasion")), paste0("model_", structures)))
fits <- lapply(names(args), function(n) {
  a <- args[[n]]
  tryCatch(suppressWarnings(do.call(lta, c(list(data = d, id = "id", n_profiles = 2, time = "time", n_starts = 1, seed = 1, max_iter = 40), a))),
           error = function(e) { cat(n, "ERROR", conditionMessage(e), "\n"); NULL })
})
names(fits) <- names(args)
tab <- do.call(rbind, lapply(names(fits), function(n) {
  f <- fits[[n]]; if (is.null(f)) return(NULL)
  theta <- .lta_theta(f)
  sc <- colSums(.lta_group_scores(theta, f))
  ll <- function(v) .lta_expectation(f$x, f$codes, f$layout, f$designs, .lta_unpack(v, f), f$occasion_of_row, unname(f$sampling_weights), f$extra_data)$log_likelihood
  fd <- vapply(seq_along(theta), function(j) { e <- replace(numeric(length(theta)), j, 1e-5); (ll(theta + e) - ll(theta - e)) / 2e-5 }, 0)
  data.frame(spec = n, n_parameters = f$n_parameters, coordinates = sum(.lta_identified_coordinates(f)), length_theta = length(theta),
             max_score_vs_fd = max(abs(sc - fd)), ll_roundtrip = ll(theta) - f$log_likelihood)
}))
print(tab, row.names = FALSE)

cat("\n--- Second order, occasion transitions, empty inner occasion under FIML\n")
de <- d; de[de$time == 3, vars] <- NA
w <- list()
f <- withCallingHandlers(lta(de, vars, "id", 2, "time", order = 2, transitions = "occasion", missing = "fiml", n_starts = 2, seed = 1),
  warning = function(c) { w[[length(w) + 1L]] <<- c; invokeRestart("muffleWarning") })
for (c in w) cat(class(c)[1], ":", conditionMessage(c), "\n")
print(subset(.lta_coefficient_table(f, "transition2", .95, "observed"), term %in% c("occasion_3","occasion_4"))[, c("from","to","term","estimate","standard_error","estimated")])
cat("n_parameters", f$n_parameters, "; counted move-into-3 logits (occasion_3 terms):", sum(grepl("transition2.*occasion_3", names(coef(f)))), "\n")
print(tryCatch(head(parameter_inference(f)[, c("block","term","standard_error","estimated")], 4), error = function(e) paste("refused:", class(e)[1])))

cat("\n--- False boundary from one extreme covariate row\n")
set.seed(11)
n <- 300; T <- 3
pid <- rep(seq_len(n), each = T); occ <- rep(seq_len(T), n)
z <- rnorm(n * T); z[7] <- -5.5   # a single person-row at -5.5 SD
s <- matrix(0L, n, T); s[, 1] <- sample(2, n, TRUE)
for (t in 2:T) { zz <- z[occ == t]; p <- plogis(-1.5 + 2 * zz); move <- runif(n) < p; s[, t] <- ifelse(move, 3L - s[, t-1], s[, t-1]) }
dz <- data.frame(id = pid, time = occ, z = z, y1 = rnorm(n * T, c(0, 3)[as.vector(t(s))]), y2 = rnorm(n * T, c(0, 3)[as.vector(t(s))]))
w <- character()
fz <- withCallingHandlers(lta(dz, c("y1","y2"), "id", 2, "time", transition_covariates = "z", n_starts = 3, seed = 1),
  warning = function(c) { w <<- c(w, class(c)[1]); invokeRestart("muffleWarning") })
cat("warnings:", unique(w), " boundary:", fz$boundary, " probability_boundary:", fz$probability_boundary, "\n")
print(fz$transition_coefficients)
lp <- exp(.lta_log_transition(fz$designs$transition[[1]], array(.lta_parameters(fz)$transition[,,,1], dim(fz$transition_coefficients)[1:3])))
cat("min move prob at occ 2:", min(lp), " second smallest row:", sort(lp[, 1, 2])[1:3], "\n")
cat("row 2 share: fitted probabilities averaged:\n"); print(get_results(fz, "transitions"))
print(tryCatch(subset(parameter_inference(fz), term == "z")[, c("from","to","estimate","standard_error")], error = function(e) paste("refused:", class(e)[1], conditionMessage(e))))
dz2 <- dz; dz2$z[7] <- 0
fz2 <- suppressWarnings(lta(dz2, c("y1","y2"), "id", 2, "time", transition_covariates = "z", n_starts = 3, seed = 1))
cat("without the extreme row: boundary", fz2$boundary, "\n")
print(tryCatch(subset(parameter_inference(fz2), term == "z")[, c("from","to","estimate","standard_error")], error = function(e) paste("refused:", class(e)[1])))

cat("\n--- Homogeneous engine, a move that never happens\n")
set.seed(4); n <- 200
s <- cbind(sample(2, n, TRUE, c(.5,.5))); s <- cbind(s, ifelse(s[,1]==1, 1L, sample(1:2, n, TRUE, c(.3,.7))), NA); s[,3] <- ifelse(s[,2]==1, 1L, sample(1:2, n, TRUE, c(.3,.7)))
dn <- data.frame(id = rep(1:n, 3), time = rep(1:3, each = n), y1 = rnorm(3*n, c(0,5)[as.vector(s)]), y2 = rnorm(3*n, c(0,5)[as.vector(s)]))
w <- character(); fh <- withCallingHandlers(lta(dn, c("y1","y2"), "id", 2, "time", n_starts = 2, seed = 1), warning = function(c) { w <<- c(w, class(c)[1]); invokeRestart("muffleWarning") })
cat("homog warnings:", unique(w), " boundary", fh$boundary, "\n"); print(round(fh$transition_probabilities[,,1], 8))
print(tryCatch(head(parameter_inference(fh)), error = function(e) paste("refused:", class(e)[1])))
dn$u <- 1
w <- character(); fg <- withCallingHandlers(lta(dn, c("y1","y2"), "id", 2, "time", n_starts = 2, seed = 1, weights = "u"), warning = function(c) { w <<- c(w, class(c)[1]); invokeRestart("muffleWarning") })
cat("general warnings:", unique(w), " boundary", fg$boundary, " LL hom", fh$log_likelihood, " gen", fg$log_likelihood, " npar", fh$n_parameters, fg$n_parameters, "\n")

cat("\n--- term names with dots\n")
dd <- d; dd$my.var <- dd$z; dd$title <- factor(ifelse(dd$id %% 2 == 0, "Dr.", "Mr"))
fd <- suppressWarnings(lta(dd, vars, "id", 2, "time", transition_covariates = c("my.var","title"), n_starts = 2, seed = 1))
print(.lta_coefficient_table(fd, "transition", .95, "observed")[, c("from","to","term","estimate","standard_error","estimated")])
print(table(.lta_identified_coordinates(fd))); cat("n_par", fd$n_parameters, "\n")
