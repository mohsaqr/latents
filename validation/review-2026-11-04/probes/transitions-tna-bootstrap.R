suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
sim <- function(n, T, seed) {
  set.seed(seed)
  P <- list(matrix(c(.9,.1,.1,.9), 2, byrow = TRUE), matrix(c(.4,.6,.6,.4), 2, byrow = TRUE))
  do.call(rbind, lapply(seq_len(n), function(i) {
    h <- sample(2, 1, prob = c(.6,.4)); s <- integer(T); s[1] <- sample(2, 1, prob = if (h == 1) c(.7,.3) else c(.3,.7))
    for (t in 2:T) s[t] <- sample(2, 1, prob = P[[h]][s[t-1], ])
    data.frame(id = i, time = seq_len(T), y1 = rnorm(T, c(0, 3)[s]), y2 = rnorm(T, c(1, -2)[s], .8), z = rnorm(T), w = rep(rnorm(1), T))
  }))
}
vars <- c("y1","y2"); d <- sim(200, 5, 1)
f <- suppressWarnings(lta(d, vars, "id", 2, "time", n_group_classes = 2, n_starts = 5, seed = 1))
cat("--- get_tna\n")
m <- get_tna(f); print(class(m)); w <- m$weights; print(round(unclass(w)[1:2,1:2], 6))
counts <- apply(f$transition_counts, 1:2, sum); man <- counts / rowSums(counts)
cat("max diff vs counts normalised by hand:", max(abs(unclass(w)[1:2,1:2] - unname(man))), "\n")
print(round(m$inits, 6)); hand_init <- as.vector(f$group_probabilities %*% f$initial_probabilities)
cat("max diff inits:", max(abs(unname(m$inits) - hand_init)), "\n")
# Model-implied marginal over classes: P(s_t=l | s_{t-1}=k) averaged over occasions by occupancy (independent of posteriors)
hand_marg <- function(f, T = 5) {
  K <- 2; H <- 2; num <- matrix(0, K, K)
  for (h in 1:H) { pr <- f$initial_probabilities[h, ]; for (t in 2:T) {
     num <- num + f$group_probabilities[h] * diag(pr) %*% f$transition_probabilities[, , h]
     pr <- as.vector(pr %*% f$transition_probabilities[, , h]) } }
  num / rowSums(num)
}
cat("model-implied marginal (no data):\n"); print(round(hand_marg(f), 4))
gm <- get_group_tna(f); print(class(gm)); print(names(gm))
cat("group_tna class matrices max diff:", max(abs(unname(unclass(gm[[1]]$weights)[1:2,1:2]) - unname(f$transition_probabilities[,,1])), abs(unname(unclass(gm[[2]]$weights)[1:2,1:2]) - unname(f$transition_probabilities[,,2]))), "\n")
cat("group_tna inits diff:", max(abs(unname(gm[[1]]$inits) - unname(f$initial_probabilities[1, ])), abs(unname(gm[[2]]$inits) - unname(f$initial_probabilities[2, ]))), "\n")

cat("\n--- resampling units\n")
draws <- c(3L, 3L, 7L)
rs <- .multilpa_resample_groups(d, split(seq_len(nrow(d)), match(d$id, unique(d$id))), "id", draws)
print(table(table(rs$id))); print(nrow(rs)); print(length(unique(rs$id)))
print(rs[rs$time == 1, c("id", "time", "y1")]); print(subset(d, id %in% c(3, 7) & time == 1)[, c("id","y1")])

cat("\n--- relabelling: probabilities must be permuted, not changed\n")
K3 <- within(sim(150, 4, 3), { y1 <- y1 + rep(sample(c(0, 3, 6), 150, TRUE), each = 4) })
g <- suppressWarnings(lta(K3, vars, "id", 3, "time", transition_covariates = "z", initial_covariates = "w", n_starts = 2, seed = 1, max_iter = 30))
est <- .lta_bootstrap_estimates(g); order3 <- c(3L, 1L, 2L)
rel <- .lta_relabel_profiles(est, order3)
softmax <- function(e) { e <- e - max(e); exp(e)/sum(exp(e)) }
xrow <- c(1, .7); zrow <- c(1, -1.3)
initial_of <- function(values) { v <- values[grepl("^initial", names(values))]
  vapply(1:2, function(k) sum(v[sprintf("initial.group_class_1.profile_%d.%s", k, c("(Intercept)", "w"))] * xrow), 0) }
p_old <- softmax(c(initial_of(est), 0)); p_new <- softmax(c(initial_of(rel), 0))
cat("initial: new[j] vs old[order[j]] max diff", max(abs(p_new - p_old[order3])), "\n")
move_of <- function(values, a) { e <- numeric(3); for (b in setdiff(1:3, a)) e[b] <- sum(values[sprintf("transition.group_class_1.profile_%d->profile_%d.%s", a, b, c("(Intercept)", "z"))] * zrow); softmax(e) }
Told <- t(vapply(1:3, function(a) move_of(est, a), numeric(3))); Tnew <- t(vapply(1:3, function(a) move_of(rel, a), numeric(3)))
cat("transition: max diff", max(abs(Tnew - Told[order3, order3])), "\n")
g2 <- suppressWarnings(lta(K3, vars, "id", 3, "time", order = 2, n_starts = 1, seed = 1, max_iter = 20))
est2 <- .lta_bootstrap_estimates(g2); rel2 <- .lta_relabel_profiles(est2, order3)
second_of <- function(values) { arr <- array(0, c(3,3,3)); for (i in 1:3) for (j in 1:3) { e <- numeric(3); for (l in setdiff(1:3, j)) e[l] <- values[sprintf("transition2.group_class_1.profile_%d->profile_%d->profile_%d.(Intercept)", i, j, l)]; arr[i,j,] <- softmax(e) }; arr }
cat("second order: max diff", max(abs(second_of(rel2) - second_of(est2)[order3, order3, order3])), "\n")

cat("\n--- bootstrap vs Wald, 40 resamples, homogeneous and general\n")
d1 <- sim(120, 4, 5)
h <- suppressWarnings(lta(d1, vars, "id", 2, "time", n_starts = 3, seed = 1))
bh <- suppressWarnings(parameter_inference(h, data = d1, method = "bootstrap", iter = 40, n_starts = 2, seed = 1))
wh <- parameter_inference(h)
print(cbind(subset(wh, parameter == "transition_probability")[, c("outcome","term","estimate")], wald = subset(wh, parameter == "transition_probability")$standard_error, boot = subset(bh, parameter == "transition_probability")$standard_error))
g1 <- suppressWarnings(lta(d1, vars, "id", 2, "time", transition_covariates = "z", n_starts = 3, seed = 1))
bg <- suppressWarnings(parameter_inference(g1, data = d1, method = "bootstrap", iter = 40, n_starts = 2, seed = 1))
wg <- parameter_inference(g1)
print(transform(wg[, c("block","from","to","term","estimate")], wald = wg$standard_error, boot = bg$standard_error))
cat("valid:", attr(bg, "n_valid"), "/", attr(bg, "iter"), "\n")

cat("\n--- bootstrap with one long sequence under occasion transitions\n")
d3 <- sim(60, 3, 6); extra <- transform(subset(sim(1, 5, 7), time > 3), id = 1L, y1 = y1); d3 <- rbind(d3, extra)
o <- suppressWarnings(lta(d3, vars, "id", 2, "time", transitions = "occasion", n_starts = 2, seed = 1))
cat("n_occasions", o$n_occasions, "\n")
bo <- withCallingHandlers(parameter_inference(o, data = d3, method = "bootstrap", iter = 20, n_starts = 1, seed = 1),
   warning = function(c) { cat("warning:", class(c)[1], conditionMessage(c), "\n"); invokeRestart("muffleWarning") })
print(table(attr(bo, "messages"), useNA = "ifany"))
print(bo[, c("block","from","to","term","estimate","standard_error")])
