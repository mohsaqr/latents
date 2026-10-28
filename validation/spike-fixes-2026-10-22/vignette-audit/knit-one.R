# Knit one vignette in a fresh process; save .md, figures and the knit env.
# Usage: Rscript knit-one.R <vignette.Rmd> <outdir>
args <- commandArgs(trailingOnly = TRUE)
input <- normalizePath(args[[1]])
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
stem <- tools::file_path_sans_ext(basename(input))
work <- file.path(out_dir, stem)
dir.create(work, showWarnings = FALSE)
file.copy(input, work, overwrite = TRUE)
setwd(work)
env <- new.env(parent = globalenv())
t0 <- proc.time()[["elapsed"]]
status <- tryCatch({
  knitr::opts_knit$set(root.dir = work)
  knitr::knit(basename(input), output = paste0(stem, ".md"), envir = env,
              quiet = TRUE)
  "ok"
}, error = \(e) paste("ERROR:", conditionMessage(e)))
elapsed <- proc.time()[["elapsed"]] - t0
objs <- mget(ls(env, all.names = TRUE), envir = env)
objs <- Filter(\(o) !is.function(o) && !inherits(o, "environment"), objs)
saveRDS(objs, file.path(work, "env.rds"))
writeLines(c(sprintf("status\t%s", status), sprintf("elapsed\t%.1f", elapsed),
             sprintf("latents\t%s", find.package("latents"))),
           file.path(work, "status.tsv"))
