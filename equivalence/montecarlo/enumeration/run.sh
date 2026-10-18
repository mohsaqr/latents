#!/bin/sh
# Run every enumeration input through Mplus, in this folder; an input that does
# not terminate normally stops the run. From the project root:
#   sh equivalence/montecarlo/enumeration/run.sh
set -e
here="$(cd "$(dirname "$0")" && pwd)"
mplus="${MPLUS:-/Applications/MplusDemo/mpdemo}"
cd "$here"
for input in *.inp; do
  stem="${input%.inp}"
  rm -f "$stem.out" "$stem-results.dat"
  "$mplus" "$input" "$stem.out" > /dev/null 2>&1
  grep -q "THE MODEL ESTIMATION TERMINATED NORMALLY" "$stem.out" || { echo "$stem did not terminate normally" >&2; exit 1; }
  echo "$stem done"
done
