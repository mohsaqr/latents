suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
# Null: iid profile membership, independent of group. What do cross-level criteria select?
draw <- function(seed, J, n, d) { set.seed(seed)
  k <- sample(1:2, J * n, TRUE); pm <- rbind(c(-1.5, -1), c(1.5, 1))[, seq_len(d), drop = FALSE]
  out <- data.frame(team = rep(sprintf("t%03d", seq_len(J)), each = n))
  out[paste0("y", seq_len(d))] <- lapply(seq_len(d), function(r) rnorm(J * n, pm[k, r], .6))
  out }
scenario <- function(label, J, n, d, vm, bv, seeds = 1:3) {
  vars <- paste0("y", seq_len(d))
  rows <- lapply(seeds, function(seed) {
    D <- draw(seed, J, n, d)
    one <- function(fam, H) suppressWarnings(multilpa(D, vars, "team", n_profiles = 2, n_group_classes = H, family = fam, seed = 1, n_starts = 3, variance_model = vm, between_variance = bv))
    f1 <- one("full_cross_level", 1); f2 <- one("full_cross_level", 2); r2 <- one("restricted_cross_level", 2)
    p1 <- suppressWarnings(multilpa(D, vars, "team", n_profiles = 2, n_group_classes = 1, variance_model = vm, seed = 1, n_starts = 3))
    p2 <- suppressWarnings(multilpa(D, vars, "team", n_profiles = 2, n_group_classes = 2, variance_model = vm, seed = 1, n_starts = 3))
    data.frame(scenario = label, seed = seed,
      rest_dBIC_H2_H1 = r2$bic - f1$bic, full_dBIC_H2_H1 = f2$bic - f1$bic,
      full_minus_rest_LL_H2 = f2$log_likelihood - r2$log_likelihood, extra_par = f2$n_parameters - r2$n_parameters,
      full_minus_rest_dBIC_H2 = f2$bic - r2$bic,
      profile_family_dBIC_H2_H1 = p2$bic - p1$bic)
  })
  do.call(rbind, rows)
}
print(rbind(
  scenario("d=2 n=6 equal/equal", 100, 6, 2, "equal", "equal"),
  scenario("d=1 n=6 equal/equal", 100, 6, 1, "equal", "equal"),
  scenario("d=2 n=6 varying/varying", 100, 6, 2, "varying", "varying"),
  scenario("d=2 n=20 equal/equal", 60, 20, 2, "equal", "equal")), digits = 4)
