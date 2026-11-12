# What does this lcmm return for subject-specific predictions?
suppressPackageStartupMessages(library(lcmm))
cat("lcmm", as.character(utils::packageVersion("lcmm")), "\n")
data <- data.frame(id = rep(1:40, each = 4), time = rep(0:3, 40))
set.seed(1)
data$y <- 1 + 0.5 * data$time + rep(stats::rnorm(40), each = 4) + stats::rnorm(160)
one <- lcmm::hlme(y ~ time, random = ~ 1, subject = "id", data = data, verbose = FALSE)
two <- lcmm::hlme(y ~ time, mixture = ~ time, random = ~ 1, subject = "id", ng = 2, data = data,
                  B = one, verbose = FALSE)
new <- data[data$id <= 3, ]
marginal <- lcmm::predictY(two, newdata = new, var.time = "time")
subject <- lcmm::predictY(two, newdata = new, var.time = "time", marg = FALSE, subject = "id")
str(marginal$pred); str(subject$pred); print(utils::head(subject$pred, 3)); print(names(subject))
