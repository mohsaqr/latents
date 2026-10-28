pkgload::load_all('.',quiet=TRUE)
set.seed(12)
g <- rep(1:30,each=6)
z <- rnorm(length(g))
d <- data.frame(id=g,z=z,y=100+rep(rnorm(30,sd=3),each=6)+ifelse(z>0,3,-3)+rnorm(length(g),sd=.5))
fit <- function(centering) withCallingHandlers(multilpa(d,'y','id',2,n_group_classes=1,
                              profile_covariates='z',centering=centering,n_starts=1,seed=2,tol=1e-10),
                              warning=function(w) invokeRestart('muffleWarning'))
none <- fit('none')
person <- fit('person')
grand <- fit('grand')
print(list(none_means=none$means,person_means=person$means,grand_means=grand$means,
           same_none_person=identical(none$means,person$means),same_none_grand=identical(none$means,grand$means),
           person_centering_field=person$centering,grand_centering_field=grand$centering))
manual <- d
manual$y <- d$y-ave(d$y,d$id)
manual_fit <- withCallingHandlers(multilpa(manual,'y','id',2,n_group_classes=1,profile_covariates='z',
                          n_starts=1,seed=2,tol=1e-10), warning=function(w) invokeRestart('muffleWarning'))
print(list(manual_person_means=manual_fit$means, advertised_person_ll=person$log_likelihood,
           manually_centered_ll=manual_fit$log_likelihood))
saveRDS(list(none=none,person=person,grand=grand,manual=manual_fit),'validation/audit-2026-10-13/model-core/dispatch-probes.rds')
