#!/bin/sh
# Runs every kit case through Latent GOLD 6.1 under Wine and collects the
# listings. Run from the project root after make-kit.R. Pass case names to run
# only those.
set -e
WINE="${WINE:-$HOME/Applications/Wine Devel.app/Contents/Resources/wine/bin/wine}"
WORK="${LATENTGOLD_WORK:-$HOME/.wine-latentgold/drive_c/lgdistal}"
KIT="$(pwd)/equivalence/latentgold-categorical-distal/kit"
RETURNED="${LATENTGOLD_RETURNED:-$(pwd)/equivalence/latentgold-categorical-distal/returned}"
mkdir -p "$RETURNED" "$WORK"
cp "$KIT"/*.dat "$KIT"/*.lgs "$WORK"/
cd "$WORK"
if [ "$#" -gt 0 ]; then cases="$*"; else cases=$(ls *.lgs | sed 's/[.]lgs$//'); fi
for case_name in $cases; do
  echo "running $case_name"
  rm -f "$case_name.lst"
  WINEPREFIX="$HOME/.wine-latentgold" "$WINE" "C:\\Program Files\\LatentGOLD6.1\\lg61.exe" \
    "$case_name.lgs" /b /o "$case_name.lst" > "$case_name.wine.log" 2>&1 || true
  if [ ! -s "$case_name.lst" ] || ! grep -q "Log-likelihood (LL)" "$case_name.lst"; then
    echo "Latent GOLD $case_name: no fitted listing" >&2
    continue
  fi
  cp "$case_name.lst" "$RETURNED/"
done
echo done
