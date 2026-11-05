suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
# r3step chain and by_group_class run BFGS from one start. Do other starts find a higher maximum?
multi <- function(obj, gr, start, sd=1, B=12) {
  ends <- vapply(seq_len(B), function(i) {
    o <- tryCatch(optim(start+rnorm(length(start),0,sd), obj, gr, method="BFGS", control=list(maxit=2000,reltol=1e-12)), error=function(e) list(value=NA))
    -o$value }, 0)
  ends
}
# Chain, K=3, weak separation
chain <- function(seed, sep, K=3L, n=150, TT=4) {
  set.seed(seed)
  s <- matrix(0L,n,TT); s[,1] <- sample.int(K,n,TRUE)
  z <- matrix(rnorm(n*TT),n,TT)
  invisible(lapply(2:TT,function(t){ mv <- runif(n)<plogis(-1+.8*z[,t]); s[,t] <<- ifelse(mv, sample.int(K,n,TRUE), s[,t-1]) }))
  mus <- rbind(c(-sep,-sep),c(sep,-sep),c(0,sep))
  sv <- as.vector(s)
  d <- data.frame(id=rep(seq_len(n),TT),time=rep(seq_len(TT),each=n),z=as.vector(z),a=rnorm(n*TT,mus[sv,1]),b=rnorm(n*TT,mus[sv,2]))
  f <- lta(d,c("a","b"),id="id",time="time",n_profiles=K,n_starts=3,seed=1)
  r <- tryCatch(r3step(f,d,transition_covariates="z"), latents_no_converge=function(e) { cat(sprintf("chain seed %d sep %.1f REFUSED: %s\n",seed,sep,conditionMessage(e))); NULL })
  if (is.null(r)) {
    tp <- get_results(f,"transitions")
    print(as.data.frame(tp)[,intersect(names(tp),c("from","to","probability","estimate"))]); return(invisible(NULL))
  }
  inp <- .lta_three_step_inputs(f,d,character(),"z")
  stacked <- do.call(rbind,inp$designs$transition)
  ins <- unlist(lapply(seq_along(inp$designs$transition),function(t) inp$layout$within[,t+1L]))
  T1 <- .lta_three_step_basis(stacked,which(ins))
  st <- inp; st$designs <- list(initial=inp$designs$initial, transition=lapply(inp$designs$transition,function(x) x%*%T1))
  obj <- function(th) -sum(.lta_three_step_evaluate(th,st)$log_likelihood)
  gr <- function(th) -colSums(.lta_three_step_evaluate(th,st)$scores)
  np <- (K-1)+2*(K-1)*K
  ends <- multi(obj,gr,rep(0,np),1,10)
  cat(sprintf("chain seed %d sep %.1f: reported LL %.4f; D-diag %s; random-start ends max %.4f (%d of %d within 1e-3 of reported, %d above)\n",seed,sep,attr(r,"log_likelihood"),paste(format(diag(attr(r,"classification_errors")),digits=2),collapse=" "),max(ends,na.rm=TRUE),sum(abs(ends-attr(r,"log_likelihood"))<1e-3,na.rm=TRUE),length(ends),sum(ends>attr(r,"log_likelihood")+1e-3,na.rm=TRUE)))
  if (any(ends>attr(r,"log_likelihood")+1e-3,na.rm=TRUE)) print(as.data.frame(r)[,c("block","from","to","term","estimate","standard_error")])
}
invisible(lapply(1:3,chain,sep=1.5)); invisible(lapply(4:6,chain,sep=1))
# by_group_class, K=3, H=2 and H=3
bg <- function(seed, K=3L, H=2L, G=60, m=10, sep=1.2) {
  set.seed(seed)
  g <- rep(seq_len(G),each=m); N <- G*m; gc <- rep(sample.int(H,G,TRUE),each=m); x <- rnorm(N)
  mix <- rbind(c(.7,.2,.1),c(.1,.2,.7),c(.3,.4,.3))[seq_len(H),seq_len(K),drop=FALSE]
  tr <- vapply(seq_len(N),function(i) sample.int(K,1,prob=mix[gc[i],]),1L)
  mus <- rbind(c(-sep,-sep),c(sep,-sep),c(0,sep))
  d <- data.frame(g=g,x=x,a=rnorm(N,mus[tr,1]),b=rnorm(N,mus[tr,2]))
  f <- multilpa(d,c("a","b"),"g",n_profiles=K,n_group_classes=H,n_starts=5,seed=1)
  r <- r3step(f,d,"x",by_group_class=TRUE)
  pc <- .multilpa_level_assignments(f,"individuals")
  mod <- .multilpa_r3step_group_model(f,as.matrix(d["x"]),pc,.multilpa_error_matrix(pc))
  ends <- multi(mod$objective,mod$gradient,mod$start,1.5,10)
  cat(sprintf("by_group seed %d K=%d H=%d: reported LL %.4f; random ends max %.4f (%d of %d within 1e-3, %d above)\n",seed,K,H,attr(r,"log_likelihood"),max(ends,na.rm=TRUE),sum(abs(ends-attr(r,"log_likelihood"))<1e-3,na.rm=TRUE),length(ends),sum(ends>attr(r,"log_likelihood")+1e-3,na.rm=TRUE)))
}
invisible(lapply(1:4,bg)); invisible(lapply(5:6,bg,H=3L))
