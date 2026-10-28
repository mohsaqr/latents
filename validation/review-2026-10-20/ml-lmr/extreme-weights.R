pkgload::load_all(quiet=TRUE)
b <- c(1e8,1e9,1e10,1e11,1e12)
res <- data.frame(large_negative=b,
 actual=vapply(b,function(size).latents_imhof_upper(0,c(1,-size)),numeric(1)),
 exact=stats::pf(b,1,1,lower.tail=FALSE))
res$absolute_error<-abs(res$actual-res$exact)
print(res,digits=12,row.names=FALSE)
