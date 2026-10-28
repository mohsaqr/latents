library(latents)
set.seed(21)
d <- do.call(rbind,lapply(1:80,function(i) {
 states <- Reduce(function(previous,t) if(runif(1)<.7) previous else 3L-previous,1:3,init=sample.int(2,1),accumulate=TRUE)
 data.frame(id=i,time=1:4,y1=rnorm(4,c(0,2)[states],.7),y2=rnorm(4,c(0,1.5)[states],.7))
}))
tail <- d[d$time==4 & d$id<=20,]; tail$time <- 5; tail[c('y1','y2')] <- NA
d <- rbind(d,tail)
f <- lta(d,c('y1','y2'),'id',2,'time',n_starts=2,seed=1,missing='fiml',transitions='occasion')
b <- parameter_inference(f,data=d,method='bootstrap',iter=4,n_starts=1,seed=45)
print(b)
print(diag(attr(b,'covariance')))
cat('unidentified covariance rows:'); print(attr(b,'covariance')[!b$estimated,,drop=FALSE])
saveRDS(list(f=f,b=b,d=d),'/tmp/astra-lta-bootstrap-probe.rds')
