#!/bin/bash
# Usage: run-all.sh <base|cur>
which=$1
S=/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/9bfae39a-cc34-4432-8b4f-db4045409538/scratchpad
export VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export R_LIBS=$S/lib/$which
for f in $S/src/cur/vignettes/*.Rmd $S/src/cur/vignettes/articles/*.Rmd; do
  echo "$(date +%T) start $(basename $f)"
  nice -n 10 Rscript $S/knit-one.R "$f" $S/out/$which > $S/out/$which-$(basename $f .Rmd).log 2>&1
  echo "$(date +%T) done  $(basename $f) rc=$?"
done
echo ALLDONE
