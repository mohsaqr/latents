#!/bin/sh
# Run every starts-study input through Mplus, JOBS at a time; an input that
# does not terminate normally is reported. From the project root:
#   sh equivalence/starts-study/mplus/run.sh
here="$(cd "$(dirname "$0")" && pwd)"
mplus="${MPLUS:-/Applications/MplusDemo/mpdemo}"
cd "$here"
ls *.inp | xargs -P "${JOBS:-2}" -I{} sh -c '
  stem="${1%.inp}"
  rm -f "$stem.out" "$stem-results.dat"
  "$2" "$1" "$stem.out" > /dev/null 2>&1
  grep -q "THE MODEL ESTIMATION TERMINATED NORMALLY" "$stem.out" && echo "$stem done" || echo "$stem did not terminate normally" >&2
' _ {} "$mplus"
