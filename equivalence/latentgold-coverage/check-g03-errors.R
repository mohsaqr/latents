# Independent check of the g03 standard errors (random-intercept growth
# mixture, one random-intercept variance shared by both classes, a residual
# variance per class). The likelihood is written here from its definition and
# its observed information taken by numDeriv, in Latent GOLD's natural
# parameters, at the latents estimate. At a maximum the observed-information
# errors do not depend on the parameterization, so the three can be compared.
# Run from the project root:
#   Rscript equivalence/latentgold-coverage/check-g03-errors.R

suppressMessages(pkgload::load_all(".", quiet = TRUE))
directory <- file.path("equivalence", "latentgold-coverage")
source(file.path(directory, "compare-functions.R"))
case <- readRDS(file.path(directory, "kit", "cases.rds"))$g03

g03_log_likelihood <- function(par, wide, times) {
  beta <- matrix(par[1:4], 2L, 2L)               # columns: classes; rows: intercept, slope
  residual <- par[5:6]
  random <- par[[7L]]
  share <- c(1 - stats::plogis(par[[8L]]), stats::plogis(par[[8L]]))
  if (any(residual <= 0) || random <= 0) return(-Inf)
  per_class <- vapply(1:2, function(k) {
    mu <- beta[1L, k] + beta[2L, k] * times
    sigma <- random + diag(residual[[k]], length(times))
    root <- chol(sigma)
    centred <- sweep(wide, 2L, mu)
    z <- backsolve(root, t(centred), transpose = TRUE)
    log(share[[k]]) - sum(log(diag(root))) - 0.5 * length(times) * log(2 * pi) -
      0.5 * colSums(z^2)
  }, numeric(nrow(wide)))
  top <- apply(per_class, 1L, max)
  sum(top + log(rowSums(exp(per_class - top))))
}

fit <- fit_case(case)
theta <- stats::coef(fit)
natural <- c(theta[["coefficient.class_1.(Intercept)"]], theta[["coefficient.class_1.time"]],
             theta[["coefficient.class_2.(Intercept)"]], theta[["coefficient.class_2.time"]],
             exp(2 * theta[["log_sigma.class_1"]]), exp(2 * theta[["log_sigma.class_2"]]),
             exp(2 * theta[["random.log_cholesky.shared.(Intercept):(Intercept)"]]),
             theta[["membership.class_2.(Intercept)"]])
wide <- as.matrix(stats::reshape(case$data, idvar = "id", timevar = "time",
                                 direction = "wide", sep = "")[, -1L])
times <- sort(unique(case$data$time))
stopifnot("the independent likelihood reproduces latents' maximum" =
            isTRUE(all.equal(g03_log_likelihood(natural, wide, times),
                             as.numeric(stats::logLik(fit)), tolerance = 1e-10)))
hessian <- numDeriv::hessian(g03_log_likelihood, natural, wide = wide, times = times)
independent <- sqrt(diag(solve(-hessian)))

#' The three sets of errors, term by term
#' @noRd
g03_error_table <- function(listing, terms, independent) {
  labels <- c("F1<-1|Class(1)", "F2<-1|Class(1)", "F1<-1|Class(2)", "F2<-1|Class(2)",
              "y1|Class(1)", "y1|Class(2)", "F1", "Class(2)<-1")
  lg <- listing$parameters$standard_error[match(labels, listing$parameters$term)]
  ours <- terms$standard_error[match(labels, terms$term)]
  data.frame(term = labels, latent_gold = lg, latents = round(ours, 6),
             independent = round(independent, 6),
             latents_vs_independent = signif(ours / independent - 1, 3),
             latent_gold_vs_independent = signif(lg / independent - 1, 3))
}

print(g03_error_table(read_listing_any(file.path(directory, "returned", "g03.lst")),
                      delta_terms(gmm_map(1:2, "intercept", 5L), theta, stats::vcov(fit)),
                      independent), row.names = FALSE)
