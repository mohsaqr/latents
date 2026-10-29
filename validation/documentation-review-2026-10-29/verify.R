# Verify that the editorial revision preserves executable package/example code.
base_text <- function(path) {
  system2("git", c("show", shQuote(paste0("HEAD:", path))), stdout = TRUE)
}
paths <- list.files("R", pattern = "\\.R$", full.names = TRUE)
unchanged <- vapply(paths, function(path) {
  identical(parse(text = base_text(path), keep.source = FALSE),
            parse(path, keep.source = FALSE))
}, logical(1))
stopifnot(all(unchanged))
stopifnot(identical(base_text("NAMESPACE"), readLines("NAMESPACE")))
cat(length(paths), "R files have unchanged executable code; NAMESPACE unchanged.\n")

rd_paths <- list.files("man", pattern = "\\.Rd$", full.names = TRUE)
issues <- lapply(rd_paths, function(path) {
  tools::checkRd(tools::parse_Rd(path))
})
stopifnot(all(lengths(issues) == 0L))
cat(length(rd_paths), "help topics parse without Rd problems.\n")

vignette_paths <- list.files("vignettes", pattern = "\\.Rmd$", full.names = TRUE)
code <- function(text) {
  out <- tempfile(fileext = ".R")
  on.exit(unlink(out))
  knitr::purl(text = text, output = out, quiet = TRUE, documentation = 0)
  parse(out, keep.source = FALSE)
}
stopifnot(all(vapply(vignette_paths, function(path) {
  baseline <- base_text(path)
  # Approved removal of the advanced helper from the introductory catalogue.
  if (identical(path, "vignettes/lpa.Rmd")) {
    baseline <- baseline[!grepl("^starting_values\\(fit\\)", baseline)]
  }
  identical(code(baseline), code(readLines(path)))
}, logical(1))))
cat(length(vignette_paths), "bundled vignettes match, allowing removal of the introductory starting_values call.\n")

# Documentation examples are part of the source-comment comparison, so check
# their extracted R expressions independently of the prose and Rd formatting.
examples <- function(path, baseline = FALSE) {
  text <- if (baseline) base_text(path) else readLines(path)
  # The condition-handling example has one deliberately revised message.
  if (baseline && identical(path, "man/latents-conditions.Rd")) {
    text <- sub("this fit carries no ordering",
                "No occasion ordering was specified for this fit.", text,
                fixed = TRUE)
  }
  rd <- tools::parse_Rd(textConnection(text))
  ex <- Filter(function(x) identical(attr(x, "Rd_tag"), "\\examples"), rd)
  if (!length(ex)) return(NULL)
  parse(text = paste(unlist(ex), collapse = ""), keep.source = FALSE)
}
old_rd <- system2("git", c("ls-tree", "-r", "--name-only", "HEAD", "man"),
                  stdout = TRUE)
old_rd <- grep("\\.Rd$", old_rd, value = TRUE)
stopifnot(all(vapply(old_rd, function(path) {
  identical(examples(path, TRUE), examples(path))
}, logical(1))))
cat(length(old_rd), "existing help examples match, allowing the revised condition message.\n")
