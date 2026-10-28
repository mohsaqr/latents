library(latents)
set.seed(912)
g <- 60; n <- g*6
c <- sample(1:2,g,TRUE)
k <- vapply(rep(c,each=6),function(h)sample(1:2,1,prob=if(h==1)c(.8,.2) else c(.2,.8)),integer(1))
d <- data.frame(id=rep(1:g,each=6),x=rnorm(n),y=rnorm(n,c(-1.2,1.2)[k]),w=rep(runif(g,.1,4),each=6))
x <- enumerate_regressions(y~x,d,n_classes=2,n_group_classes=2,id='id',n_starts=2,seed=3,weights='w')
f <- get_results(x,'model',n_classes=2,n_group_classes=2)
quality <- latents:::.latents_grid_quality_table(x)
p <- f$expectation$rho; w <- f$spec$sampling_weights
expected <- 1+sum(w*rowSums(ifelse(p>0,p*log(p),0)))/(sum(w)*log(ncol(p)))
cat('grid group entropy:',quality$group_entropy,'direct weighted group entropy:',expected,'fit group entropy:',get_results(f,'fit')$group_entropy,'\n')
