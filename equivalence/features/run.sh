#!/usr/bin/env bash
# Run one feature kit through Latent GOLD 6.1 (under Wine) and the Mplus 9 Demo.
#
# From the project root:
#   sh equivalence/features/run.sh <feature> [all|latentgold|mplus] [pattern]
# e.g. sh equivalence/features/run.sh direct-effects latentgold small
#
# A kit is equivalence/features/<feature>/kit/. Each case is a set of files
# named after it: <case>.lgs and/or <case>.inp, and the data files they read,
# all named <case>.* or <case>-*. Outputs land in <feature>/returned/. Earlier
# outputs of a case are removed before it runs, and a case that writes no
# listing stops the run. JOBS (default 4) cases run at once.
[ -n "${BASH_VERSION:-}" ] || exec bash "$0" "$@"
set -e
here="$(cd "$(dirname "$0")" && pwd)"
feature="${1:?usage: run.sh <feature> [all|latentgold|mplus] [pattern]}"
which_program="${2:-all}"
pattern="${3:-}"
kit="$here/$feature/kit"
returned="$here/$feature/returned"
[ -d "$kit" ] || { echo "no kit at $kit" >&2; exit 1; }
mkdir -p "$returned"
jobs="${JOBS:-4}"

run_latentgold_case() {
  input="$1"; kit="$2"; returned="$3"; feature="$4"
  wine="${WINE:-$HOME/Applications/Wine Devel.app/Contents/Resources/wine/bin/wine}"
  prefix="${WINEPREFIX:-$HOME/.wine-latentgold}"
  case_name="$(basename "$input" .lgs)"
  work="$prefix/drive_c/lgfeatures/$feature/$case_name"
  rm -rf "$work"
  mkdir -p "$work"
  rm -f "$returned/$case_name.lst" "$returned/${case_name}_"*.txt
  for file in "$kit/$case_name".* "$kit/$case_name"-*; do
    [ -e "$file" ] && cp "$file" "$work/"
  done
  # Latent GOLD's exit status is not reliable under Wine; the fresh files are.
  (cd "$work" && WINEPREFIX="$prefix" "$wine" "C:\\Program Files\\LatentGOLD6.1\\lg61.exe" \
    "$case_name.lgs" /b /o "$case_name.lst" > /dev/null 2>&1) || true
  if [ ! -s "$work/$case_name.lst" ] ||
     ! grep -q "Log-likelihood (LL)" "$work/$case_name.lst"; then
    echo "Latent GOLD $case_name: no fitted likelihood listing was written" >&2
    exit 1
  fi
  cp "$work/$case_name.lst" "$returned/"
  for produced in "$work/${case_name}_"*.txt; do
    [ -e "$produced" ] && cp "$produced" "$returned/"
  done
  echo "Latent GOLD $case_name listing collected; validate convergence and estimates downstream"
}

run_mplus_case() {
  input="$1"; kit="$2"; returned="$3"
  mplus="${MPLUS:-/Applications/MplusDemo/mpdemo}"
  case_name="$(basename "$input" .inp)"
  rm -f "$returned/$case_name.out" "$returned/$case_name"-*.dat "$returned/$case_name-console.txt"
  for file in "$kit/$case_name".inp "$kit/$case_name".dat; do
    [ -e "$file" ] && cp "$file" "$returned/"
  done
  (cd "$returned" && "$mplus" "$case_name.inp" "$case_name.out" > "$case_name-console.txt" 2>&1) || true
  if [ ! -s "$returned/$case_name.out" ] ||
     ! grep -q "THE MODEL ESTIMATION TERMINATED NORMALLY" "$returned/$case_name.out"; then
    echo "Mplus $case_name: no normally terminated output was written" >&2
    exit 1
  fi
  echo "Mplus $case_name done"
}

export -f run_latentgold_case run_mplus_case

if [ "$which_program" = "all" ] || [ "$which_program" = "latentgold" ]; then
  for input in "$kit"/*.lgs; do
    [ -e "$input" ] || continue
    case "$(basename "$input")" in *"$pattern"*) printf '%s\n' "$input" ;; esac
  done | xargs -P "$jobs" -I{} bash -c 'run_latentgold_case "$1" "$2" "$3" "$4"' _ {} "$kit" "$returned" "$feature"
fi

if [ "$which_program" = "all" ] || [ "$which_program" = "mplus" ]; then
  for input in "$kit"/*.inp; do
    [ -e "$input" ] || continue
    case "$(basename "$input")" in *"$pattern"*) printf '%s\n' "$input" ;; esac
  done | xargs -P "$jobs" -I{} bash -c 'run_mplus_case "$1" "$2" "$3"' _ {} "$kit" "$returned"
fi
