Sys.setenv(NOT_CRAN='true')
pkgload::load_all('.',quiet=TRUE)
testthat::test_dir('tests/testthat',filter='^(release-audit-lta|lta-general|review-lta|review-1026-lta-open|review-1026-lta-r2|release-transitions)$',reporter='summary',stop_on_failure=TRUE)
