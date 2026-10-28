library(testthat)
destination <- "validation/fixes-2026-10-20"
suite_names <- c("unit", "equivalence")
summary <- do.call(rbind, lapply(suite_names, function(kind) {
  table <- as.data.frame(readRDS(file.path(destination, paste0(kind, "-tests.rds"))))
  if (kind == "unit") {
    corrected <- as.data.frame(readRDS(file.path(destination, "catalogue-unit-tests.rds")))
    affected <- table$file %in% corrected$file
    stopifnot(sum(affected) == nrow(corrected),
              identical(table$test[affected], corrected$test),
              identical(table$nb[affected], corrected$nb),
              sum(table$failed[affected]) == 1L,
              !any(table$failed[!affected]), !any(table$error),
              !any(corrected$failed), !any(corrected$error))
    table[affected, ] <- corrected
    saveRDS(table, file.path(destination, "unit-final-cases.rds"))
    write.csv(table[, setdiff(names(table), "result")],
      file.path(destination, "unit-final-cases.csv"), row.names = FALSE)
    cat("Unit aggregation replaces the complete condition-catalogue test file",
        "after its documentation-only correction. Original run is retained.\n")
  }
  data.frame(suite = kind, cases = nrow(table),
    as.list(colSums(table[, c("nb", "failed", "error", "warning", "skipped", "passed")])))
}))
write.csv(summary, file.path(destination, "validation-summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
original <- read.csv(file.path(destination, "pre-catalogue-source-manifest.csv"))
current_md5 <- unname(tools::md5sum(original$file))
stopifnot(identical(original$file[current_md5 != original$md5], "R/conditions.R"))
before <- readRDS(file.path(destination, "code-before-catalogue.rds"))
after <- setNames(lapply(names(before), parse, keep.source = FALSE), names(before))
stopifnot(identical(before, after))
cat("All production executable expressions match the complete suite run;",
    "only condition-catalogue comments changed.\n")
write.csv(data.frame(file = original$file, md5 = current_md5),
  file.path(destination, "final-source-manifest.csv"), row.names = FALSE)
paths <- c(list.files("R", "[.]R$", full.names = TRUE),
  list.files("tests/testthat", "[.]R$", full.names = TRUE),
  list.files("tests/equivalence", "[.]R$", full.names = TRUE),
  list.files("man", "[.]Rd$", full.names = TRUE),
  list.files("vignettes", "[.]Rmd$", full.names = TRUE),
  "DESCRIPTION", "NAMESPACE", "NEWS.md")
write.csv(data.frame(file = paths, md5 = unname(tools::md5sum(paths))),
  file.path(destination, "final-delivery-manifest.csv"), row.names = FALSE)
comparison <- read.csv(file.path(destination, "latentgold", "comparison.csv"))
cat("Retained Latent GOLD comparison statuses:\n")
print(table(comparison$status))
stopifnot(!any(comparison$status == "disagree"),
          all(summary$failed == 0), all(summary$error == 0))
