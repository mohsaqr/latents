devtools::load_all(quiet=TRUE)
Sys.setenv(NOT_CRAN='true')
library(testthat)
set.seed(1310)
d <- data.frame(g=rep(1:20,each=10), x=rnorm(200), z=rnorm(200))
d$x <- d$x + rep(c(-2,2),each=100)
noise <- lpa(d,c('x','z'),2,noise=TRUE,n_starts=1,seed=1)
test_that('audit confirms named diagnostics bypass noise refusal', {
  expect_error(plot(noise,what='entropy'),class='latents_unsupported_noise')
  entropy <- plot_diagnostics(noise,what='entropy')
  expect_s3_class(entropy,'ggplot')
  expect_gt(max(entropy$data$entropy),1)
  expect_error(plot(noise,what='avepp'),class='latents_unsupported_noise')
  expect_s3_class(plot_diagnostics(noise,what='avepp'),'ggplot')
})
w <- d;w$w <- c(rep(1,100),rep(10,100))
grid <- enumerate_lpa(w,c('x','z'),1:2,model='EEI',weights='w',n_starts=1,seed=1)
test_that('audit confirms weighted tree disagrees with fitted shares', {
  fit <- get_candidate(grid,n_profiles=2)
  tree <- plot_tree(grid)$layers[[2]]$data
  actual <- tree$share[tree$n_profiles==2]
  expected <- latents:::.gg_profile_key(fit)$share
  expect_gt(max(abs(sort(actual)-sort(expected))),0.3)
})
test_that('audit confirms former abbreviated enumeration views fail', {
  expect_error(plot(grid,what='enum'),class='latents_unknown_criterion')
  expect_error(plot(grid,what='tr'),class='latents_unknown_criterion')
})
b <- data.frame(x=rnorm(160),n=rep(c(10L,20L),80))
b$success<-rbinom(160,b$n,plogis(.5+b$x));b$failure<-b$n-b$success
bin <- mixture_regression(cbind(success,failure)~x,b,1,family='binomial',n_starts=1,seed=1)
test_that('audit confirms binomial training/newdata prediction scale changes', {
  expect_equal(predict(bin)$fitted/predict(bin,b)$fitted,as.numeric(b$n))
})
test_that('audit confirms NULL seed does not restore RNG as consolidated docs promise', {
  before<-.Random.seed
  invisible(simulate(bin,seed=NULL))
  expect_false(identical(before,.Random.seed))
  before<-.Random.seed
  invisible(simulate(bin,seed=1))
  expect_identical(before,.Random.seed)
})
test_that('audit confirms plot help aliases were removed', {
  rd<-tools::parse_Rd('man/latents-plots.Rd')
  aliases<-vapply(Filter(function(x)identical(attr(x,'Rd_tag'),'\\alias'),rd),paste0,character(1),collapse='')
  expect_false('plot.multilpa' %in% aliases)
  expect_false('plot.latents_growth_mixture' %in% aliases)
})
