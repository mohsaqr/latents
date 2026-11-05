suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(5)
n <- 300; kk <- rbinom(n,1,.5)+1
dp <- data.frame(a = ifelse(kk==1, "x", ifelse(runif(n)<.5,"x","y")), b = ifelse(kk==1,"u", ifelse(runif(n)<.3,"u","v")),
                 c = ifelse(kk==1, ifelse(runif(n)<.8,"r","s"), ifelse(runif(n)<.2,"r","s")), d = ifelse(kk==1, ifelse(runif(n)<.9,"p","q"), ifelse(runif(n)<.3,"p","q")),
                 e = ifelse(runif(n)<c(.7,.2)[kk],"m","n"))
l0 <- lca(dp, c("a","b","c","d","e"), n_classes = 2, n_starts = 6, seed = 1)
l1 <- lca(dp, c("a","b","c","d","e"), n_classes = 2, n_starts = 6, seed = 1, local_dependence = "c:d")
floor0 <- .multilpa_probability_floor(l0); floor1 <- .multilpa_probability_floor(l1)
cat("null: any probability at floor", any(unlist(l0$response_probabilities) <= floor0), " boundary flag", l0$boundary, "\n")
cat("alt : any probability at floor", any(unlist(l1$response_probabilities) <= floor1), " boundary flag", l1$boundary, "\n")
print(l0$response_probabilities$a)
warned <- character()
tb <- withCallingHandlers(as.data.frame(compare_models(free = l0, joint = l1, test = "lrt"))[,c("model","n_parameters","log_likelihood","lr_statistic","lr_df","lr_p_value")],
                          warning = function(w) { warned <<- c(warned, class(w)[1]); invokeRestart("muffleWarning") })
print(tb); cat("warnings raised by the LRT:", if (length(warned)) warned else "none", "\n")
q <- diagnostics(l0, plots = FALSE); ck <- get_results(q, "checks"); print(ck[ck$check == "estimates at a bound", c("check","display","status","detail")])
cat(capture.output(print(q))[1:2], sep = "\n")
pi0 <- tryCatch(parameter_inference(l0), error = function(e) e); cat(class(pi0)[1], ":", conditionMessage(pi0), "\n")
# count mean and NB dispersion on their bound
cat("Number of grid rows called boundary:\n")
g <- enumerate_lca(dp, c("a","b","c","d","e"), n_classes = 1:3, n_starts = 6, seed = 1)
cat(sprintf("%s: boundary %s; probability at floor %s\n", names(g$fits), g$table$boundary,
            vapply(g$fits, function(f) any(unlist(f$response_probabilities) <= .multilpa_probability_floor(f)), TRUE)), sep = "")
