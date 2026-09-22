# Retained Latent GOLD output

Produced 2026-09-22 by Latent GOLD 6.1 (Basic+Advanced/Syntax+Choice), free academic single-user licence,
running under Wine Devel 11.17 (Gcenx build) on macOS 26.3.1, arm64; prefix ~/.wine-latentgold.

Each model was estimated in batch as `lg61.exe <case>.lgs /b /o <case>.lst`, from a directory inside the
Wine drive. The kit that produced it was built by `make-kit.R` under multilpa
0.5.0.

## Why it is kept

Latent GOLD is commercial and Windows-only, and its licence is per user. The
listings and posterior files here are the evidence, so `compare.R` runs
offline and nobody needs Latent GOLD again unless the kit's data or syntax
change. `INPUTS.md5` fingerprints exactly those inputs; `compare.R` warns with
`multilpa_stale_latentgold_output` when the kit no longer matches.

## Files

- 26 listings and posterior files (`OUTPUTS.md5`)
- from 24 input files (`INPUTS.md5`)

These are Latent GOLD's own output, kept for verification only. They are not
redistributed as part of the package: `validation/` is in `.Rbuildignore`.
