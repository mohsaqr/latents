suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
w <- function(expr) { warns <- character(); v <- withCallingHandlers(tryCatch(expr,error=function(e) { warns <<- c(warns,paste0("ERROR:",class(e)[1])); NULL }),warning=function(w){ warns <<- c(warns,class(w)[1]); invokeRestart("muffleWarning") }); list(v=v,w=warns) }
n <- 400; tr <- 1L+(runif(n)<.5)
d <- data.frame(a=rnorm(n,c(-2.5,2.5)[tr]),b=rnorm(n,c(-2.5,2.5)[tr]))
f <- lpa(d,c("a","b"),n_profiles=2,n_starts=3,seed=1)
# covariate that (almost) determines the class: strong and exact
d$x_strong <- ifelse(f$subject_profiles==1,1,-1)*runif(n,.2,2)
d$x_noise <- rnorm(n)
print(.multilpa_error_matrix(.multilpa_level_assignments(f,"individuals")))
a <- w(r3step(f,d,"x_strong")); cat("warnings/errors:",a$w,"\n"); if (!is.null(a$v)) print(as.data.frame(a$v)[,c("term","estimate","standard_error","p_value")])
d$x_mod <- ifelse(f$subject_profiles==1,1,-1)*.6+rnorm(n,0,.15)
a <- w(r3step(f,d,"x_mod")); cat("near separated, warnings/errors:",a$w,"\n"); if (!is.null(a$v)) print(as.data.frame(a$v)[,c("term","estimate","standard_error","p_value")])
b <- w(multilpa(transform(d,g=seq_len(n)),c("a","b"),"g",n_profiles=2,n_group_classes=1,profile_covariates="x_strong",n_starts=3,seed=1)); cat("one-step warnings:",b$w,"\n")
# Chain seed 6 of the separation probe, inspected
one <- function(seed, n=80, K=3L, TT=4, sep=1.2) {
  set.seed(seed)
  s <- matrix(0L,n,TT); s[,1] <- sample.int(K,n,TRUE); z <- matrix(rnorm(n*TT),n,TT)
  invisible(lapply(2:TT,function(t){ mv <- runif(n)<plogis(-1.3+.8*z[,t]); s[,t] <<- ifelse(mv, sample.int(K,n,TRUE), s[,t-1]) }))
  mus <- rbind(c(-sep,-sep),c(sep,-sep),c(0,sep)); sv <- as.vector(s)
  d <- data.frame(id=rep(seq_len(n),TT),time=rep(seq_len(TT),each=n),z=as.vector(z),a=rnorm(n*TT,mus[sv,1]),b=rnorm(n*TT,mus[sv,2]))
  f <- lta(d,c("a","b"),id="id",time="time",n_profiles=K,n_starts=3,seed=1)
  r <- w(r3step(f,d,transition_covariates="z")); print(as.data.frame(r$v)[,c("block","from","to","term","estimate","standard_error","p_value","p_adjusted")]); cat("warnings:",r$w,"\n")
  ob <- w(r3step(f,d,transition_covariates="z",vcov_type="observed")); print(as.data.frame(ob$v)[,c("from","to","term","estimate","standard_error","p_value")]); cat("warnings:",ob$w,"\n")
  inp <- .lta_three_step_inputs(f,d,character(),"z"); print(round(inp$errors,3))
  M <- matrix(inp$modal,n,TT); cat("modal move counts\n"); print(table(from=as.vector(M[,-TT]),to=as.vector(M[,-1])))
  print(as.data.frame(get_results(f,"transitions"))[,c("from","to","probability")])
  os <- w(lta(d,c("a","b"),id="id",time="time",n_profiles=K,transition_covariates="z",n_starts=3,seed=1))
  cat("one-step warnings:",os$w,"\n")
  if (!is.null(os$v)) { v <- w(parameter_inference(os$v)); cat("one-step inference warnings:",v$w,"\n") }
}
one(6); one(4)
