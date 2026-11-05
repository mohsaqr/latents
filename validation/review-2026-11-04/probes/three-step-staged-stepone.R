suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
G <- 100; m <- 10; n <- G*m
gc <- rep(1:2,length.out=G); g <- rep(seq_len(G),each=m)
tr <- 1L+(runif(n)<c(.3,.7)[gc[g]])
# poorly separated, realistic: entropy about 0.6
d <- data.frame(g=g, a=rnorm(n,c(-.8,.8)[tr]), b=rnorm(n,c(-.7,.7)[tr]), c=rnorm(n,c(.7,-.7)[tr]), y=rnorm(n,c(0,1)[tr]))
full <- multilpa(d,c("a","b","c"),"g",n_profiles=2,n_group_classes=2,n_starts=5,seed=1,tol=1e-10)
st <- fit_staged(d,c("a","b","c"),"g",n_profiles=2,n_group_classes=2,n_starts=5,seed=1)
cat("entropy full:", get_results(full,"entropy")$entropy %||% NA, "\n")
print(class(st)); print(sum(full$n_parameters), quote=FALSE); print(st$n_parameters)
show <- function(fit,label) { invisible(lapply(c("bch","ml"),function(m){
  a <- three_step(fit,d,"y",method=m); b <- three_step(fit,d,"y",method=m,step_one="estimated")
  cat(sprintf("%-7s %-3s means %s | SE fixed %s | SE estimated %s | ratio %s | step_one attr=%s\n",label,m,
    paste(format(a$estimate,digits=5),collapse=" "),paste(format(a$standard_error,digits=5),collapse=" "),paste(format(b$standard_error,digits=5),collapse=" "),
    paste(format(b$standard_error/a$standard_error,digits=4),collapse=" "),attr(b,"step_one")))
})) }
show(full,"full"); show(st,"staged")
cat("dim V1 full", dim(vcov(full,scale="unconstrained")), " staged", dim(vcov(st,scale="unconstrained")), "\n")
print(rownames(vcov(st,scale="unconstrained")))
