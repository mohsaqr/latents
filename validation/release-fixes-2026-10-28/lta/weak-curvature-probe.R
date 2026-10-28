pkgload::load_all('.',quiet=TRUE)
# A near-boundary initial/transition Bernoulli logit's exact score and
# curvature. Inverse score norm is huge despite a local logistic length of 1.
eta <- 16
score <- function(x) stats::plogis(-x)
parameter_scale <- 1 / score(eta)
total <- function(displacement) score(eta + displacement * parameter_scale) * parameter_scale
unbounded <- latents:::.inference_score_information(total,0,1e-4,scheme='five_point',h_rule='absolute')[1,1]/parameter_scale^2
bounded <- latents:::.inference_score_information(total,0,1e-4,scheme='five_point',h_rule='scaled',scale=1/parameter_scale)[1,1]/parameter_scale^2
expected <- stats::plogis(eta)*stats::plogis(-eta)
print(data.frame(eta,parameter_scale,unbounded,bounded,expected))
cat('Relative errors: unbounded',abs(unbounded/expected-1),'bounded',abs(bounded/expected-1),'\n')
stopifnot(abs(bounded/expected-1)<1e-8,abs(unbounded/expected-1)>.9)
# Final test-file runtime, with fixture model refits included.
Sys.setenv(NOT_CRAN='true')
timing <- system.time(result <- testthat::test_file('tests/testthat/test-release-audit-lta.R',reporter='summary'))
print(timing)
stopifnot(!any(as.data.frame(result)$failed>0),!any(as.data.frame(result)$error))
