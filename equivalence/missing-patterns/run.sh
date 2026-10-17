#!/bin/sh
# Run every missing-pattern case through Mplus and Latent GOLD.
# From the project root: sh equivalence/missing-patterns/run.sh [mplus|latentgold]
# Outputs are written to equivalence/missing-patterns/returned/.
set -e
here="$(cd "$(dirname "$0")" && pwd)"
kit="$here/kit"
returned="$here/returned"
mkdir -p "$returned"
which_program="${1:-all}"

if [ "$which_program" = "all" ] || [ "$which_program" = "mplus" ]; then
  mplus="${MPLUS:-/Applications/MplusDemo/mpdemo}"
  for input in "$kit"/m*.inp; do
    case_name="$(basename "$input" .inp)"
    rm -f "$returned/$case_name.out" "$returned/$case_name-posteriors.dat" \
      "$returned/$case_name-results.dat"
    cp "$input" "$kit/$case_name.dat" "$returned/"
    (cd "$returned" && "$mplus" "$case_name.inp" "$case_name.out" > "$case_name-console.txt" 2>&1)
    if ! grep -q "THE MODEL ESTIMATION TERMINATED NORMALLY" "$returned/$case_name.out"; then
      echo "Mplus $case_name did not terminate normally" >&2
      exit 1
    fi
    echo "Mplus $case_name done"
  done
fi

if [ "$which_program" = "all" ] || [ "$which_program" = "latentgold" ]; then
  wine="${WINE:-$HOME/Applications/Wine Devel.app/Contents/Resources/wine/bin/wine}"
  prefix="${WINEPREFIX:-$HOME/.wine-latentgold}"
  work="$prefix/drive_c/lgmissing"
  mkdir -p "$work"
  for input in "$kit"/m*.lgs; do
    case_name="$(basename "$input" .lgs)"
    # Outputs of an earlier run are removed first, so a failed run cannot leave
    # them behind to be copied as if new.
    rm -f "$work/$case_name.lst" "$work/${case_name}_posteriors.txt" \
      "$returned/$case_name.lst" "$returned/${case_name}_posteriors.txt"
    cp "$input" "$kit/$case_name-lg.dat" "$work/"
    # Latent GOLD's exit status is not reliable under Wine; the fresh files are.
    (cd "$work" && WINEPREFIX="$prefix" "$wine" "C:\\Program Files\\LatentGOLD6.1\\lg61.exe" \
      "$case_name.lgs" /b /o "$case_name.lst" > /dev/null 2>&1) || true
    for produced in "$case_name.lst" "${case_name}_posteriors.txt"; do
      if [ ! -s "$work/$produced" ]; then
        echo "Latent GOLD $case_name: no $produced was written" >&2
        exit 1
      fi
      cp "$work/$produced" "$returned/"
    done
    echo "Latent GOLD $case_name done"
  done
fi
