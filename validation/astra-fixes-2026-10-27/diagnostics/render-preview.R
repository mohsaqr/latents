pkgload::load_all('.', quiet=TRUE)
rmarkdown::render('tmp/astra-diagnostics.Rmd', output_file='astra-diagnostics.html',
                  knit_root_dir=getwd(), quiet=TRUE)
stopifnot(file.exists('tmp/astra-diagnostics.html'),
          file.info('tmp/astra-diagnostics.html')$size > 10000)
cat('Rendered tmp/astra-diagnostics.html\n')
