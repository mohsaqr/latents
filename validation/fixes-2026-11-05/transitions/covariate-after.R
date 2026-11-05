pkgload::load_all('.', quiet = TRUE)
generate <- function(n, slope, seed) {
  set.seed(seed)
  occasions <- 3L
  id <- rep(seq_len(n), each = occasions)
  time <- rep(seq_len(occasions), n)
  z <- rnorm(n * occasions)
  states <- matrix(0L, n, occasions)
  states[, 1L] <- sample(2L, n, TRUE)
  invisible(lapply(2:3, function(t) {
    chance <- plogis(-1.5 + slope * z[time == t])
    states[, t] <<- ifelse(runif(n) < chance, 3L - states[, t - 1L], states[, t - 1L])
  }))
  data.frame(id = id, time = time, z = z,
    y1 = rnorm(n * occasions, c(0, 3)[as.vector(t(states))]),
    y2 = rnorm(n * occasions, c(0, 3)[as.vector(t(states))]))
}
fit <- function(data, covariate = 'z', starts = 3L) {
  warnings <- character()
  result <- withCallingHandlers(lta(data, c('y1','y2'), 'id', 2, 'time',
    transition_covariates = covariate, n_starts = starts, seed = 1),
    warning = function(w) {warnings <<- c(warnings, class(w)[1]); invokeRestart('muffleWarning')})
  cat('boundary',result$boundary,'warnings',paste(warnings,collapse=','),'\n')
  result
}
data <- generate(300L, 2, 11L)
str(data);print(head(data))
base <- fit(data)
data$z[which(data$time == 2L)[3L]] <- -6.5
outlier <- fit(data)
table <- subset(parameter_inference(outlier), term == 'z')
print(table[,c('from','to','estimate','standard_error')])
stopifnot(!outlier$boundary, all(is.finite(table$standard_error)))
dotted <- generate(150L, 1, 3L)
dotted$title <- factor(ifelse(dotted$id %% 2L == 0L, 'Aa', 'Zed.'))
dotfit <- fit(dotted, 'title', 2L)
dottable <- get_results(dotfit, 'transition_coefficients')
print(dottable[,c('from','to','term','standard_error','estimated')])
stopifnot(any(dottable$term == 'titleZed.'),
  all(dottable$estimated[dottable$term == 'titleZed.']),
  all(is.finite(dottable$standard_error[dottable$term == 'titleZed.'])),
  sum(.lta_identified_coordinates(dotfit)) == dotfit$n_parameters)
saveRDS(list(base = base, outlier = outlier, dotted = dotfit),
 'validation/fixes-2026-11-05/transitions/covariate-after.rds')
