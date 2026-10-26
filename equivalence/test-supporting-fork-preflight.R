# Pure safety/protocol checks. Historical kits, fits and studies are never run.

support_fork_guard <- function() {
  scope <- new.env(parent = baseenv())
  sys.source(file.path("..", "validation", "fork-runtime.R"), envir = scope)
  scope$.validation_fork_cores
}

support_fork_source <- function(x) {
  is.call(x) && identical(x[[1L]], as.name("source")) &&
    grepl("fork-runtime.R", paste(deparse(x), collapse = ""), fixed = TRUE)
}

support_fork_early <- function(x) {
  is.call(x) && identical(x[[1L]], as.name("invisible")) &&
    is.call(x[[2L]]) && identical(x[[2L]][[1L]], as.name(".validation_fork_cores"))
}

support_fork_unwrap <- function(x) {
  if (!is.call(x)) return(x)
  if (identical(x[[1L]], as.name(".validation_fork_cores"))) return(x[[2L]])
  as.call(lapply(as.list(x), support_fork_unwrap))
}

test_that("supporting fork guards preserve the complete scientific AST", {
  baseline <- c(
    "features/ml-three-step/make-kit.R" = "4f0c32e6b041944bd23dabf86c27d666",
    "search/blrt-emem.R" = "bf7f8082f1e095eb2ec1f1c0cafb5c22",
    "search/pilot.R" = "b02a264a913f2f76bf360ce5b93d7b6f",
    "montecarlo/audit-maxima.R" = "2d16900b84091f29b78eb8bbae017544")
  expected_calls <- c(2L, 1L, 1L, 1L)
  invisible(lapply(seq_along(baseline), function(i) {
    expressions <- parse(names(baseline)[i], keep.source = FALSE)
    sources <- vapply(expressions, support_fork_source, logical(1))
    preflight <- vapply(expressions, support_fork_early, logical(1))
    expect_identical(sum(sources), 1L)
    cleaned <- lapply(expressions[!sources & !preflight], support_fork_unwrap)
    text <- vapply(cleaned, function(x) paste(deparse(x, width.cutoff = 500),
                                              collapse = "\n"), "")
    path <- tempfile()
    on.exit(unlink(path), add = TRUE)
    writeLines(text, path)
    expect_identical(unname(tools::md5sum(path)), unname(baseline[i]))
    calls <- list()
    walk <- function(x) {
      if (!is.call(x)) return(invisible(NULL))
      head <- x[[1L]]
      bare <- identical(head, as.name("mclapply"))
      qualified <- is.call(head) && identical(head[[1L]], as.name("::")) &&
        identical(head[[2L]], as.name("parallel")) &&
        identical(head[[3L]], as.name("mclapply"))
      if (bare || qualified) calls[[length(calls) + 1L]] <<- x
      invisible(lapply(as.list(x)[-1L], walk))
    }
    invisible(lapply(expressions, walk))
    expect_length(calls, expected_calls[i])
    invisible(lapply(calls, function(call) {
      expect_identical(call$mc.cores[[1L]], as.name(".validation_fork_cores"))
      expect_length(call$mc.cores, 2L)
      expect_identical(call$mc.preschedule, FALSE)
    }))
  }))
})

test_that("supporting preflight refuses unsafe forks and preserves safe requests", {
  guard <- support_fork_guard()
  withr::local_seed(77)
  state <- .Random.seed
  expect_identical(guard(2L, "unix", "Linux", "reference BLAS"), 2L)
  expect_identical(guard(2, "unix", "Darwin", "libRblas.dylib"), 2)
  expect_identical(guard(1L, "windows", "Windows", "BLAS"), 1L)
  expect_error(guard(2L, "windows", "Windows", "BLAS"),
                 class = "latents_validation_unsafe_fork")
  expect_error(guard(2L, "unix", "Darwin", "Accelerate/vecLib"),
                 class = "latents_validation_unsafe_fork")
  expect_identical(.Random.seed, state)
})

test_that("the historical kit refuses unsafe forks before deleting its data", {
  expressions <- parse("features/ml-three-step/make-kit.R", keep.source = FALSE)
  preflight <- which(vapply(expressions, support_fork_early, logical(1)))
  deletion <- which(vapply(expressions, function(x) {
    is.call(x) && identical(x[[1L]], as.name("unlink"))
  }, logical(1)))
  expect_length(preflight, 1L)
  expect_lt(preflight, deletion)
  guard <- support_fork_guard()
  scope <- new.env(parent = baseenv())
  scope$.validation_fork_cores <- function(cores) {
    guard(cores, "unix", "Darwin", "Accelerate")
  }
  expect_error(eval(expressions[[preflight]], scope),
                 class = "latents_validation_unsafe_fork")
})

test_that("search scripts locate the common guard from their script location", {
  root <- normalizePath("..")
  target <- normalizePath(file.path(root, "validation", "fork-runtime.R"))
  invisible(lapply(c("search/blrt-emem.R", "search/pilot.R"), function(path) {
    expressions <- parse(path, keep.source = FALSE)
    scope <- new.env(parent = baseenv())
    scope$commandArgs <- function(trailing_only) {
      c("Rscript", paste0("--file=", file.path(root, "equivalence", path)))
    }
    scope$source <- function(file) normalizePath(file)
    expect_identical(eval(expressions[[1L]], scope), target)
    preflight <- which(vapply(expressions, support_fork_early, logical(1)))
    creation <- which(vapply(expressions, function(x) {
      is.call(x) && identical(x[[1L]], as.name("dir.create"))
    }, logical(1)))
    expect_lt(preflight, creation)
  }))
})
