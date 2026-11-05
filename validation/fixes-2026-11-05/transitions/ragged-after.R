pkgload::load_all('.', quiet = TRUE)
set.seed(110505)
data <- data.frame(id = rep(1:50, each = 3), time = rep(1:3, 50))
state <- sample(0:1, nrow(data), TRUE)
data$a <- rnorm(nrow(data), 4 * state, .7)
data$b <- rnorm(nrow(data), -3 * state, .7)
data <- rbind(data, data.frame(id = 1L, time = 4:5, a = c(.2, 3.7), b = c(.1, -2.8)))
str(data); print(head(data))
reference <- readRDS('validation/review-2026-11-05/engines/results.rds')$transitions
warnings <- character()
bootstrap <- withCallingHandlers(parameter_inference(reference, data = data,
  method = 'bootstrap', iter = 15, n_starts = 1, max_iter = 1000, seed = 48),
  warning = function(w) {
    warnings <<- c(warnings, paste(class(w)[1], conditionMessage(w)))
    invokeRestart('muffleWarning')
  })
print(warnings)
print(attributes(bootstrap)[c('n_valid','messages')])
stopifnot(inherits(bootstrap,'latents_bootstrap'), attr(bootstrap,'n_valid') >= 2,
  attr(bootstrap,'n_valid') < 15,
  any(grepl('occasion or design term is absent', attr(bootstrap,'messages')),na.rm=TRUE))
saveRDS(bootstrap,'validation/fixes-2026-11-05/transitions/ragged-after.rds')
