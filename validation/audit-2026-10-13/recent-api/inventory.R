files <- c("R/enumeration.R", "R/selection.R", "R/candidates.R", "R/compare-models.R", "R/get-data.R", "R/accessors.R", "R/additive-compare.R", "R/lta-compare.R")
tests <- c("enumeration", "enumeration", "enumeration, enumerate-verbs, compat-0921", "enumeration", "get-data, tidy-classes", "get-data, tidy-classes", "additive-compare", "lta-compare")
findings <- list(
  enumerate_classes = "RA-09 (preexisting grid side); RA-10 (preexisting)",
  compare_models = "RA-02; RA-03",
  .multilpa_compare_fits = "RA-02; RA-07 (contract gap); RA-09",
  .latents_comparison_classification = "RA-01; RA-05",
  .gg_selection_smallest = "RA-04",
  .latents_kneedle = "RA-08",
  .multilpa_catalogue = "RA-05",
  get_results.latents_family_enumeration = "RA-06",
  get_results.latents_transition_enumeration = "RA-06"
)
inventory <- do.call(rbind, lapply(seq_along(files), function(i) {
  expressions <- parse(files[i], keep.source = TRUE)
  refs <- attr(expressions, "srcref")
  do.call(rbind, lapply(seq_along(expressions), function(j) {
    e <- expressions[[j]]
    if (!is.call(e) || !is.symbol(e[[1]]) || !as.character(e[[1]]) %in% c("<-", "=") || length(e) != 3L || !is.call(e[[3]]) || !identical(as.character(e[[3]][[1]]), "function")) return(NULL)
    name <- as.character(e[[2]])
    data.frame(file = files[i], function_name = name, start_line = as.integer(refs[[j]][1]),
      end_line = as.integer(refs[[j]][3]), review_status = "static body and call-contract reviewed",
      related_passing_tests = tests[i], direct_runtime_probe = if(name %in% names(findings)) "see findings and probes.R" else "not asserted individually",
      finding = findings[[name]] %||% "no actionable finding found in this audit", stringsAsFactors = FALSE)
  }))
}))
write.csv(inventory, "validation/audit-2026-10-13/recent-api/function-inventory.csv", row.names = FALSE)
print(table(inventory$file)); cat("Total top-level functions: ", nrow(inventory), "\n")
