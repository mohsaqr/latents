# Every exported plot verb and wild-card table verb has a help page of its
# own, named after it, so `?plot_avepp` opens plot_avepp and nothing else
# (the combined pages latents-plots, wildcard-plots, wildcard-tables and
# wildcard-grid-tables were split in 1.1.15).
.latents_wild_verbs <- c("get_checks", "get_entropy", "get_classification",
                         "get_avepp", "get_residuals", "get_replication",
                         "get_candidates", "get_best", "get_deltas",
                         "get_knees", "get_tree")

test_that("each exported plot and table verb has its own help page", {
  man <- test_path("..", "..", "man")
  skip_if_not(dir.exists(man), "package source `man/` not available")
  exports <- getNamespaceExports("latents")
  verbs <- c(grep("^plot_", exports, value = TRUE),
             intersect(.latents_wild_verbs, exports))
  expect_gt(length(verbs), 30L)
  pages <- file.path(man, paste0(verbs, ".Rd"))
  expect_true(all(file.exists(pages)), info = paste(verbs[!file.exists(pages)], collapse = ", "))
  aliases <- lapply(pages, \(page) {
    rd <- tools::parse_Rd(page)
    tags <- vapply(rd, \(x) attr(x, "Rd_tag"), character(1))
    unlist(lapply(rd[tags == "\\alias"], as.character))
  })
  # One function per page: no page carries a second plot verb as an alias.
  expect_identical(vapply(aliases, length, integer(1)), rep(1L, length(verbs)))
  expect_identical(unlist(aliases), verbs)
  bundles <- c("latents-plots", "wildcard-plots", "wildcard-tables",
               "wildcard-grid-tables")
  expect_false(any(file.exists(file.path(man, paste0(bundles, ".Rd")))))
})

test_that("the plot and table pages stay short", {
  man <- test_path("..", "..", "man")
  skip_if_not(dir.exists(man), "package source `man/` not available")
  pages <- c(list.files(man, pattern = "^plot_.*\\.Rd$", full.names = TRUE),
             file.path(man, paste0(.latents_wild_verbs, ".Rd")))
  lines <- vapply(pages, \(page) length(readLines(page)), integer(1))
  expect_lte(max(lines), 100L)
})

# The long single-function pages were rewritten for readers in 1.1.15
# (multilpa 778 -> 448 Rd lines, get_results 768 -> 483); no page may grow
# back past 600 lines.
test_that("no help page exceeds 600 Rd lines", {
  man <- test_path("..", "..", "man")
  skip_if_not(dir.exists(man), "package source `man/` not available")
  pages <- list.files(man, pattern = "\\.Rd$", full.names = TRUE)
  lines <- vapply(pages, \(page) length(readLines(page)), integer(1))
  names(lines) <- basename(pages)
  expect_lte(max(lines), 600L, label = paste(names(which.max(lines)), max(lines)))
})
