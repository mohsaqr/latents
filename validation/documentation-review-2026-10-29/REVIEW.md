# Documentation revision — 2026-10-29

The revision covers the CRAN reference manual and function help: all 37
existing topics, a new package overview, DESCRIPTION, and the four bundled
vignettes. Roxygen comments in 24 R files were revised and the corresponding
Rd topics regenerated. The package overview is maintained directly in
`man/latents-package.Rd`. Website-only articles and historical review records
were outside this revision's scope.

Changes replace conversational wording and development-history narrative
with descriptions of the models, arguments, outputs, assumptions and
limitations. Scientific corrections include the current start-selection tie
rule, conditional uncertainty in staged and three-step analyses, appropriate
interpretation of entropy and BIC weights, local rather than guaranteed
global convergence, and accurate descriptions of simulated datasets. The
vignettes distinguish model selection from classification quality and avoid
interpreting a nonsignificant coefficient as evidence of no association.

Verification:

- All 103 R files have identical parsed executable code relative to HEAD
  (`96c4b58`); NAMESPACE is byte-identical.
- All 38 Rd topics pass parsing and `tools::checkRd`.
- All four vignettes retain identical executable chunks. Existing help
  examples retain their executable expressions except for one deliberately
  revised string returned by the condition-handling example.
- Synthetic-data smoke checks covered single-level profiles, mixture
  regression, staged fitting, diagnostics, parameter inference, BCH results,
  prediction, starting values, and descriptive summaries. An initial staged
  example with independent Gaussian indicators produced singular observed
  information; this is a model-identification limitation, not a regression.
- Roxygen regeneration completed without warnings. `R CMD build` successfully
  rendered all four vignettes. The 111-page manual was generated with
  `R CMD Rd2pdf`; its first page and a representative function-help page were
  visually inspected. The PDF is `docs/latents-manual.pdf`.
- Claim/attribution checks: 247 passed, no failures, warnings or skips.
  Thirteen external-software agreement statements were removed from the
  public prose and claim registry; three retained conventions were reanchored.
  Backing tests and result files were preserved.
- Verification-matrix safeguards: 148 assertions passed; four Python tests
  passed. Nineteen existing result files were re-attested after verifying
  source-code identity, without rerunning statistical benchmarks. Previous
  sidecars are archived in
  `validation/documentation-prose-attestation-2026-10-29/`.
- The matrix's regenerated text, HTML and CSV outputs are unchanged. Its
  optional Chrome PDF rerender failed in the restricted environment; the
  existing, unchanged matrix PDF was retained. This is separate from the
  newly generated R reference manual, which built successfully.

The first temporary archive omitted `tests/fixtures`, causing 14 missing-file
test errors; all documentation checks passed. The archive was rebuilt with
all fixtures (verified byte-for-byte), and the complete package check was
repeated. The final test run reports 4,295 passes, 336 skips, no failures and
no warnings. Final `R CMD check --no-build-vignettes` status: **OK** (0 errors, 0 warnings,
0 notes). All vignette scripts and the PDF-manual check passed. Full output
is recorded in `check.log`. The
checked archive is `tmp/documentation-review/repaired/latents_1.1.1.tar.gz`, with
SHA-256 `f08576e357f3a45abd760446713361c4169215f226179ffef3beaf52724e1251`.
The version remains 1.1.1; this is an editorial revision. No release,
installation into the user's library, commit or push was performed.
