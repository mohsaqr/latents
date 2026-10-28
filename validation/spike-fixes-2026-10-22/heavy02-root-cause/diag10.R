# How many distinct partitions do 200 random-centre k-means runs give on heavy_02?
suppressMessages(pkgload::load_all(".", quiet = TRUE))
invisible(lapply(c("helper-montecarlo.R", "helper-mplus-classification.R"), \(file) {
  sys.source(file.path("equivalence", file), envir = globalenv())
}))
root <- file.path("equivalence", "features", "direct-effects")
sys.source(file.path(root, "helpers.R"), envir = globalenv())
case <- de_cases(root)$heavy_02
y <- scale(as.matrix(subset(case$data, select = case$spec$vars)))
set.seed(1)
signature <- vapply(seq_len(200), function(run) {
  cluster <- .multilpa_start_kmeans(y, 2L)$cluster
  # Labels are arbitrary: name the cluster of row 1 "1".
  paste(as.integer(cluster == cluster[[1]]), collapse = "")
}, character(1))
counts <- sort(table(signature), decreasing = TRUE)
cat("distinct partitions in 200 runs:", length(counts), "; runs per partition:", unname(counts), "\n")
