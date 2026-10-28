help_aliases <- function(doc) {
  blocks <- Filter(function(x) identical(attr(x, "Rd_tag"), "\\alias"), doc)
  vapply(blocks, function(x) paste(x, collapse = ""), character(1))
}
old_help <- function(p) {
  text <- system2("git", c("show", paste0("a10a919:man/", p)), stdout = TRUE)
  connection <- textConnection(text)
  on.exit(close(connection))
  tools::parse_Rd(connection)
}
old_files <- system2("git", c("ls-tree", "--name-only", "a10a919:man"), stdout = TRUE)
old_files <- old_files[grepl("[.]Rd$", old_files)]
old_aliases <- unique(unlist(lapply(old_files, function(p) help_aliases(old_help(p)))))
current_files <- list.files("man", pattern = "[.]Rd$", full.names = TRUE)
current <- lapply(current_files, tools::parse_Rd)
new_aliases <- unique(unlist(lapply(current, help_aliases)))
lost <- setdiff(old_aliases, new_aliases)
cat("Previous help aliases:", length(old_aliases), "missing:", length(lost), "\n")
print(lost)
utils::write.csv(data.frame(alias = old_aliases, preserved = old_aliases %in% new_aliases),
                 "validation/fixes-2026-10-13/help-alias-check.csv", row.names = FALSE)
stopifnot(length(lost) == 0L)
plot_doc <- paste(readLines("man/latents-plots.Rd"), collapse = "\n")
stopifnot(all(vapply(c("multilpa", "multilpa_enumeration", "multilpa_covariates",
                       "multilpa_transitions", "multilpa_lta", "multilpa_additive",
                       "multilpa_cross_level", "latents_family_enumeration",
                       "latents_transition_enumeration", "latents_comparison",
                       "latents_growth_mixture", "latents_mixture_regression",
                       "latents_pooled", "multilpa_bootstrap_lrt"),
                     function(class) grepl(paste0("\\method{plot}{", class, "}"),
                                           plot_doc, fixed = TRUE), logical(1))))
pkgload::load_all(quiet = TRUE)
pkgdown::check_pkgdown()
cat("All released aliases, plot method usages and pkgdown references pass\n")
