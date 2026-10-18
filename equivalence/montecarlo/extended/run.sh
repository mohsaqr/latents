#!/bin/sh
# Run every Monte Carlo case through Mplus and Latent GOLD.
# From the project root: sh equivalence/montecarlo/run.sh [mplus|latentgold] [pattern]
# Outputs land in equivalence/montecarlo/returned/. Earlier outputs of a case
# are removed before it runs, and a case that does not terminate normally, or
# writes no output, stops the run.
set -e
here="$(cd "$(dirname "$0")" && pwd)"
kit="$here/kit"
returned="$here/returned"
mkdir -p "$returned"
which_program="${1:-all}"
pattern="${2:-}"

run_mplus_case() {
  input="$1"; returned="$2"; kit="$3"; mplus="$4"
  case_name="$(basename "$input" .inp)"
  rm -f "$returned/$case_name.out" "$returned/$case_name-posteriors.dat" "$returned/$case_name-results.dat"
  cp "$input" "$kit/$case_name.dat" "$returned/"
  (cd "$returned" && "$mplus" "$case_name.inp" "$case_name.out" > "$case_name-console.txt" 2>&1)
  if ! grep -q "THE MODEL ESTIMATION TERMINATED NORMALLY" "$returned/$case_name.out"; then
    echo "Mplus $case_name did not terminate normally" >&2
    exit 1
  fi
  echo "Mplus $case_name done"
}

if [ "$which_program" = "all" ] || [ "$which_program" = "mplus" ]; then
  mplus="${MPLUS:-/Applications/MplusDemo/mpdemo}"
  export -f run_mplus_case 2>/dev/null || true
  ls "$kit"/*.inp | grep "${pattern}" | while read -r input; do
    printf '%s\n' "$input"
  done | xargs -P "${JOBS:-4}" -I{} sh -c "$(declare -f run_mplus_case 2>/dev/null || type run_mplus_case | tail -n +2); run_mplus_case \"\$1\" \"$returned\" \"$kit\" \"$mplus\"" _ {}
fi

if [ "$which_program" = "all" ] || [ "$which_program" = "latentgold" ]; then
  wine="${WINE:-$HOME/Applications/Wine Devel.app/Contents/Resources/wine/bin/wine}"
  prefix="${WINEPREFIX:-$HOME/.wine-latentgold}"
  work="$prefix/drive_c/lgextended"
  mkdir -p "$work"
  ls "$kit"/*.lgs | grep "${pattern}" | while read -r input; do
    case_name="$(basename "$input" .lgs)"
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
