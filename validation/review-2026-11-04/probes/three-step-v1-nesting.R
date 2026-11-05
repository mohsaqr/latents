suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
# Nested data, one group class: step three is cluster-robust by default, while V1 is the
# observed-information covariance of a likelihood that treats observations as independent.
ns <- asNamespace("latents")
orig <- get(".multilpa_step_one_correction", ns)
src <- paste(deparse(orig, width.cutoff=500L), collapse="\n")
pat <- "stats::vcov(x, scale = \"unconstrained\")"
stopifnot(grepl(pat, src, fixed=TRUE))
robust <- eval(parse(text=sub(pat, "attr(parameter_inference(x, vcov_type = \"robust\"), \"covariance_unconstrained\")", src, fixed=TRUE)), envir=ns)
environment(robust) <- environment(orig)
one <- function(seed, sep=.8, G=60, m=15, rho=1.2) {
  set.seed(seed)
  g <- rep(seq_len(G),each=m); u <- rep(rnorm(G,0,rho),each=m); n <- G*m
  tr <- 1L+(runif(n)<plogis(u))
  # indicators with a group-level shift as well: real nesting of the measurement
  e_a <- rep(rnorm(G,0,.5),each=m); e_b <- rep(rnorm(G,0,.5),each=m)
  d <- data.frame(g=g, a=rnorm(n,c(-sep,sep)[tr])+e_a, b=rnorm(n,c(-sep,sep)[tr])+e_b, y=rnorm(n,c(0,1)[tr])+rep(rnorm(G,0,.4),each=m))
  f <- multilpa(d,c("a","b"),"g",n_profiles=2,n_group_classes=1,n_starts=5,seed=1,tol=1e-10)
  V_o <- vcov(f,scale="unconstrained"); V_r <- attr(parameter_inference(f,vcov_type="robust"),"covariance_unconstrained")
  a <- three_step(f,d,"y"); b <- three_step(f,d,"y",step_one="estimated")
  unlockBinding(".multilpa_step_one_correction", ns); assign(".multilpa_step_one_correction", robust, ns)
  c2 <- three_step(f,d,"y",step_one="estimated")
  assign(".multilpa_step_one_correction", orig, ns); lockBinding(".multilpa_step_one_correction", ns)
  cat(sprintf("seed %d: SE fixed %s | estimated (observed V1) %s | estimated (robust V1) %s | median sqrt(diag ratio) V_r/V_o %.2f\n",seed,
      paste(format(a$standard_error,digits=4),collapse=" "),paste(format(b$standard_error,digits=4),collapse=" "),paste(format(c2$standard_error,digits=4),collapse=" "),
      median(sqrt(diag(V_r)/diag(V_o)))))
}
invisible(lapply(1:3,one))
invisible(lapply(4:5,one,sep=1.2))
# Numbers an attribute claim needs
d <- data.frame(a=c(rnorm(150,-2),rnorm(150,2)),b=c(rnorm(150,-2),rnorm(150,2)),y=rnorm(300))
f <- lpa(d,c("a","b"),n_profiles=2,n_starts=2,seed=1)
cat("variance attr none/pairs/overall:", attr(three_step(f,d,"y",method="ml"),"variance"),"|",attr(three_step(f,d,"y",method="ml",contrast="pairs"),"variance"),"|",attr(three_step(f,d,"y",method="ml",contrast="overall"),"variance"),"\n")
