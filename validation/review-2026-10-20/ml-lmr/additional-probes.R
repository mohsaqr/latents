pkgload::load_all(quiet=TRUE)
z <- readRDS('validation/review-2026-10-20/ml-lmr/models.rds')
x <- get_candidate(z$base,2)
baseline <- three_step(x,z$data,'distal',method='ml',vcov_type='independent')
res <- lapply(c(1e6,1e8,1e9,1e10,1e11,1e12),function(offset) {
 ans <- tryCatch(three_step(x,transform(z$data,distal=distal+offset),'distal',method='ml',vcov_type='independent'),error=function(e)e)
 if(inherits(ans,'error'))return(data.frame(offset,error=conditionMessage(ans),max_se_ratio=NA_real_,max_mean_error=NA_real_))
 data.frame(offset,error=NA_character_,max_se_ratio=max(ans$standard_error/baseline$standard_error), max_mean_error=max(abs(ans$estimate-offset-baseline$estimate)))
})
print(do.call(rbind,res),row.names=FALSE)
# Pure-coordinate repair demonstrates the scaled matrix is regular and gives
# the same reference distribution; no source changed.
a <- .latents_vuong_parts(get_candidate(z$scaled,2))
b <- .latents_vuong_parts(get_candidate(z$scaled,1))
rescale <- function(p) {
 s <- 1/sqrt(diag(p$information))
 p$scores <- sweep(p$scores,2,s,'*')
 p$information <- p$information*tcrossprod(s)
 p
}
cat('\nScaled-coordinate comparison (same public fits):\n')
print(.latents_vuong_row(rescale(a),rescale(b)))
# Public example whose second converged fit is a local maximum worse than
# its simpler neighbour: generate overfitted rare/skewed components using one
# start and retain only cases with usable robust information.
# Avoid forging fit likelihoods or scores.
find_negative <- function(seed) {
 set.seed(seed)
 n<-120
 cl<-sample(1:3,n,TRUE,prob=c(.45,.45,.1))
 dat<-data.frame(a=rnorm(n,c(-1,1,3)[cl]),b=rnorm(n,c(-1,1,-3)[cl]))
 g<-suppressWarnings(enumerate_lpa(dat,c('a','b'),n_profiles=2:3,model='VVI',n_starts=1,seed=seed,tol=1e-10))
 tab<-as.data.frame(g)
 if(!all(tab$converged)||diff(tab$log_likelihood)>=-1e-6)return(NULL)
 scored<-tryCatch(lapply(1:2,function(i)suppressWarnings(.latents_vuong_parts(get_candidate(g,i)))),error=function(e)NULL)
 if(is.null(scored))return(NULL)
 list(seed=seed,grid=g,parts=scored,table=tab)
}
candidates<-lapply(1:50,find_negative)
found<-Filter(Negate(is.null),candidates)
cat('\nNegative LR public candidate cases:',length(found),'\n')
if(length(found)) {
 one<-found[[1]]
 print(one$table)
 print(tryCatch(get_results(one$grid,'lmr'),error=function(e)conditionMessage(e)))
 saveRDS(one,'validation/review-2026-10-20/ml-lmr/negative-lr-model.rds')
}
