Sys.setenv(NOT_CRAN = "true", VECLIB_MAXIMUM_THREADS = "1", OMP_NUM_THREADS = "1")
pkgload::load_all(".", quiet = TRUE)
set.seed(251105)
cat("Version", as.character(utils::packageVersion("latents")), "\n")
z <- matrix(rnorm(900), 300, 3)
z[, 3] <- .8*z[, 1] + .3*z[, 2] + .5*z[, 3]
design <- cbind(1, z)
coefficients <- matrix(c(.2, .4, -.3, .2, -.5, .1, .4, -.2), 4, 2)
affine <- diag(c(1, 3, -.4, 11))
affine[1, -1] <- c(1e8, -2e7, 3e6)
other_design <- design %*% affine
other_coefficients <- backsolve(affine, coefficients)
check_design <- function(x) {
  basis <- .r3step_design_basis(x)
  whitened <- t(backsolve(basis, t(x), transpose = TRUE))
  cat("Orthonormal Gram error:", max(abs(crossprod(whitened)/nrow(x) - diag(ncol(x)))), "\n")
  stopifnot(max(abs(crossprod(whitened)/nrow(x) - diag(ncol(x)))) < 1e-7)
}
check_design(design)
check_design(other_design)
reference <- .multilpa_logit_prior(coefficients, design)
actual <- .multilpa_logit_prior(other_coefficients, other_design)
cat("Affine predicted probability error", max(abs(actual-reference)), "\n")
stopifnot(max(abs(actual-reference)) < 1e-7)
.r3step_require_interior(coefficients, design)
.r3step_require_interior(other_coefficients, other_design)
pieces <- list(n_classes=3L, group_index=rep(seq_len(60),each=5L), modal=sample.int(3L,300,TRUE))
pieces$posteriors <- matrix(.1,300,3)
pieces$posteriors[cbind(seq_len(300),pieces$modal)] <- .8
x <- list(n_group_classes=3L, profile_probabilities=matrix(c(.5,.3,.2,.2,.6,.2,.1,.3,.6),3,3,byrow=TRUE),
          group_probabilities=c(.3,.4,.3))
errors <- diag(3)*.8+.2/3
model <- .multilpa_r3step_group_model(x,z,pieces,errors)
shifted_model <- .multilpa_r3step_group_model(x,other_design[,-1],pieces,errors)
theta <- model$start + seq_along(model$start)*.007
natural <- model$transform %*% theta
natural2 <- natural
nalpha <- 6L
nbeta <- 6L
beta <- matrix(natural[nalpha+seq_len(nbeta)],3,2)
back <- solve(affine[-1,-1], beta)
natural2[nalpha+seq_len(nbeta)] <- as.vector(back)
natural2[seq_len(nalpha)] <- natural[seq_len(nalpha)] -
  rep(as.vector(affine[1,-1] %*% back),each=3)
theta2 <- as.vector(solve(shifted_model$transform,natural2,tol=0))
cat("Group likelihood affine error", model$objective(theta)-shifted_model$objective(theta2), "\n")
stopifnot(abs(model$objective(theta)-shifted_model$objective(theta2)) < 1e-6)

g <- rep(seq_len(45),each=6)
w <- rep(c(0,1,3),length.out=45)
d <- data.frame(g=g,w=w[g],a=rnorm(length(g)),b=rnorm(length(g),.7))
d$a[seq(7,nrow(d),11)] <- NA_real_
d$b[seq(8,nrow(d),13)] <- NA_real_
str(d)
fit <- multilpa(d,c("a","b"),"g",n_profiles=1L,n_group_classes=1L,
                weights="w",centering="grand",scaling="grand",missing="fiml",
                n_starts=1L,seed=1L,tol=1e-11)
nuisance <- .multilpa_grand_transform(fit)
cat("Grand weighted mean equations", colSums(nuisance$contributions), "\n")
cat("Grand transform replay error", max(abs(nuisance$transform(nuisance$eta)-fit$indicator_data),na.rm=TRUE), "\n")
stopifnot(max(abs(colSums(nuisance$contributions))) < 1e-8,
          max(abs(nuisance$transform(nuisance$eta)-fit$indicator_data),na.rm=TRUE) < 1e-10)
cat("Shared helper zeroSE result (existing path):\n")
print(.inference_wald(1,0,.95))
cat("Independent algebra and nuisance-equation checks passed.\n")
