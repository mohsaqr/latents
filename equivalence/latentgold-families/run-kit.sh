#!/bin/sh
# Runs every kit case through Latent GOLD 6.1 under Wine and collects the
# output. Run from the project root after make-kit.R.
set -e
WINE="${WINE:-$HOME/Applications/Wine Devel.app/Contents/Resources/wine/bin/wine}"
WORK="${LATENTGOLD_WORK:-$HOME/.wine-latentgold/drive_c/lgfamilies}"
RETURNED="${LATENTGOLD_RETURNED:-$(pwd)/equivalence/latentgold-families/returned}"
mkdir -p "$RETURNED"
cd "$WORK"
for f in g*.lgs; do
  case_name="${f%.lgs}"
  echo "running $case_name"
  rm -f "$case_name.lst" "${case_name}_posteriors.txt"
  WINEPREFIX="$HOME/.wine-latentgold" "$WINE" "C:\\Program Files\\LatentGOLD6.1\\lg61.exe" "$f" /b /o "$case_name.lst" > "$case_name.wine.log" 2>&1
  if [ ! -s "$case_name.lst" ] || [ ! -s "${case_name}_posteriors.txt" ] ||
     ! grep -q "Log-likelihood (LL)" "$case_name.lst"; then
    echo "Latent GOLD $case_name: missing fresh fitted listing or posteriors" >&2
    exit 1
  fi
  cp "$case_name.lst" "${case_name}_posteriors.txt" "$RETURNED/"
done
echo done
