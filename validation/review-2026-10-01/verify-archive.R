# Verify the built package contains the reviewed R, Rd, tests and NEWS sources.
# R CMD build legitimately adds DESCRIPTION metadata, so check it separately.
evidence <- file.path('validation', 'review-2026-10-01')
archive <- file.path('tmp', 'review-package', 'latents_0.9.7.tar.gz')
manifest <- read.csv(file.path(evidence, 'SOURCE.csv'))
selected <- grepl('^(R|man|tests)/', manifest$path) |
  manifest$path %in% c('NAMESPACE', 'NEWS.md')
selected <- selected & basename(manifest$path) != 'Rplots.pdf'
manifest <- manifest[selected, ]
extracted <- file.path('tmp', 'review-package', 'source-snapshot')
dir.create(extracted, recursive = TRUE, showWarnings = FALSE)
utils::untar(archive, exdir = extracted)
built_paths <- file.path(extracted, 'latents', manifest$path)
actual <- unname(tools::md5sum(built_paths))
stopifnot(identical(actual, manifest$md5))
write.csv(data.frame(path = manifest$path, md5 = actual),
          file.path(evidence, 'ARCHIVE-SOURCE.csv'), row.names = FALSE)
writeLines(c(paste('Archive MD5:', unname(tools::md5sum(archive))),
             paste('Matching source files:', nrow(manifest)),
             'DESCRIPTION build metadata is excluded from byte comparison.'),
           file.path(evidence, 'archive-identity.log'))
cat('Built package matches', nrow(manifest), 'reviewed source files.\n')
