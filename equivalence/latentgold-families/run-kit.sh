#!/bin/sh
# Runs every kit case through Latent GOLD 6.1 under Wine and collects the
# output. Run from the project root after make-kit.R.
set -e
WINE="$HOME/Applications/Wine Devel.app/Contents/Resources/wine/bin/wine"
WORK="$HOME/.wine-latentgold/drive_c/lgfamilies"
RETURNED="$(pwd)/equivalence/latentgold-families/returned"
cd "$WORK"
for f in g*.lgs; do
  case_name="${f%.lgs}"
  echo "running $case_name"
  WINEPREFIX="$HOME/.wine-latentgold" "$WINE" "C:\\Program Files\\LatentGOLD6.1\\lg61.exe" "$f" /b /o "$case_name.lst" > "$case_name.wine.log" 2>&1
  cp "$case_name.lst" "${case_name}_posteriors.txt" "$RETURNED/"
done
echo done
