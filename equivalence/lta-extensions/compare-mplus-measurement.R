# Mplus 9 (demo) two-occasion LTA with occasion-specific thresholds (measurement
# not invariant): 1200 persons, 3 binary items per occasion, 2 classes.
# Reference: JStats validation/r-reference/mplus-demo-lta-classvarying-measurement.*
source(file.path("equivalence", "lta-extensions", "common.R"))
raw <- read.table(file.path(jstats, "validation", "r-reference",
                            "mplus-demo-lta-classvarying-measurement.dat"))
long <- do.call(rbind, lapply(1:2, function(t) {
  items <- raw[, (t - 1L) * 3L + 1:3]
  names(items) <- paste0("u", 1:3)
  cbind(person = seq_len(nrow(raw)), occasion = t, items)
}))
long[paste0("u", 1:3)] <- lapply(long[paste0("u", 1:3)], factor)
fit <- lta(long, paste0("u", 1:3), "person", n_profiles = 2, time = "occasion",
           categorical = paste0("u", 1:3), measurement = "occasion", n_starts = 10,
           seed = 1, tol = 1e-10, max_iter = 3000)
thresholds <- function(block) {
  vapply(block$response_probabilities, function(p) log(p[, 1L] / p[, 2L]), numeric(2))
}
mplus_t1 <- rbind(c(-2.616, -1.975, -1.965), c(2.080, 1.723, 1.951))
mplus_t2 <- rbind(c(-1.392, -2.938, -1.014), c(1.386, 2.720, 0.857))
perm1 <- align(thresholds(fit$measurement[[1L]]), mplus_t1)
perm2 <- align(thresholds(fit$measurement[[2L]]), mplus_t2)
table <- get_results(fit, "transitions")
ours <- matrix(0, 2, 2)
ours[cbind(match(table$from, c("profile_1", "profile_2")),
           match(table$to, c("profile_1", "profile_2")))] <- table$probability
result <- rbind(
  row_result("mplus_measurement", "log likelihood", fit$log_likelihood, -4015.054, 5e-3),
  row_result("mplus_measurement", "parameter count", fit$n_parameters, 15, 0),
  row_result("mplus_measurement", "max |threshold diff|",
             max(abs(thresholds(fit$measurement[[1L]])[perm1, ] - mplus_t1),
                 abs(thresholds(fit$measurement[[2L]])[perm2, ] - mplus_t2)), 0, 2e-3),
  row_result("mplus_measurement", "max |transition diff|",
             max(abs(ours[perm1, perm2] - rbind(c(0.842, 0.158), c(0.265, 0.735)))),
             0, 2e-3))
print(result, row.names = FALSE, digits = 8)
saveRDS(result, file.path("equivalence", "lta-extensions", "result-mplus-measurement.rds"))
