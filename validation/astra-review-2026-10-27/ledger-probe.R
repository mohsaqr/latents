pkgload::load_all('.',quiet=TRUE)
source('equivalence/harness/registry.R')
inv <- data.frame(output='measure$estimate',table='measure',column='estimate',n_values=3L)
cmp <- compare_values('one of three',1,1,tolerance=1e-8)
cmp$output <- 'measure$estimate'; cmp$n_values <- 1L
ledger <- output_ledger('partial-probe',inv,cmp)
print(ledger)
cat('run.R uncovered gate detects:', nrow(subset(ledger,status=='not covered')), '\n')
