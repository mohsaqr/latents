library(latents, lib.loc='/private/tmp/astra-review/lib')
set.seed(1027)
d <- data.frame(id=rep(1:20,each=5),x=sample(0:1,100,TRUE),y=sample(0:1,100,TRUE))
calls <- 0L
next_seed <- function() {calls <<- calls+1L; calls}
f <- suppressMessages(multilca(d,vars=c('x','y'),id=NULL,n_profiles=1,n_starts=1,seed=next_seed()))
print(c(seed_expression_evaluations=calls, seed_used=f$seed))
# Standard error/condition call stays bounded for the originally repaired path.
wrong <- tryCatch(multilca(d,vars='absent',id=NULL,n_profiles=1),error=identity)
print(list(condition_class=class(wrong),call=conditionCall(wrong)))
