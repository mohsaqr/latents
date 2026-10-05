# Equivalence tests share the helpers of tests/testthat.
helper_files <- list.files(test_path("..", "testthat"), "^helper-.*[.]R$", full.names = TRUE)
invisible(lapply(helper_files, sys.source, envir = environment()))
