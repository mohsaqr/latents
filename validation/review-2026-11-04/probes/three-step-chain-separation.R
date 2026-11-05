suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
# Small chains: when a transition cell is empty, what does the table report?
one <- function(seed, n=80, K=3L, TT=4, sep=1.2) {
  set.seed(seed)
  s <- matrix(0L,n,TT); s[,1] <- sample.int(K,n,TRUE); z <- matrix(rnorm(n*TT),n,TT)
  invisible(lapply(2:TT,function(t){ mv <- runif(n)<plogis(-1.3+.8*z[,t]); s[,t] <<- ifelse(mv, sample.int(K,n,TRUE), s[,t-1]) }))
  mus <- rbind(c(-sep,-sep),c(sep,-sep),c(0,sep)); sv <- as.vector(s)
  d <- data.frame(id=rep(seq_len(n),TT),time=rep(seq_len(TT),each=n),z=as.vector(z),a=rnorm(n*TT,mus[sv,1]),b=rnorm(n*TT,mus[sv,2]))
  f <- lta(d,c("a","b"),id="id",time="time",n_profiles=K,n_starts=3,seed=1)
  warns <- character()
  r <- withCallingHandlers(tryCatch(r3step(f,d,transition_covariates="z"),latents_no_converge=function(e) NULL),warning=function(w){ warns <<- c(warns,class(w)[1]); invokeRestart("muffleWarning") })
  if (is.null(r)) { cat(sprintf("seed %2d refused no_converge\n",seed)); return(invisible(NULL)) }
  big <- which(abs(r$estimate)>10)
  cat(sprintf("seed %2d max|est| %7.2f  SE at it %8.3f  min p among big %s | min SE %.3f | warnings: %s\n",seed,max(abs(r$estimate)),r$standard_error[which.max(abs(r$estimate))],
     if(length(big)) format(min(r$p_value[big]),digits=3) else "-",min(r$standard_error,na.rm=TRUE),paste(unique(warns),collapse=",")))
  if (length(big)) print(as.data.frame(r)[big,c("block","from","to","term","estimate","standard_error","p_value")])
}
invisible(lapply(1:14,one))
