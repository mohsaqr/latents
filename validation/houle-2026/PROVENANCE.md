# Supplemental Table S1 numerical fixture

`table-s1.csv` transcribes only the likelihood, free-parameter count, and four
information criteria in the 24 rows of Table S1 of the user-supplied DOCX
supplement to Houle et al. (2026), DOI `10.1177/10944281261469432`.
The supplied supplement's filename is
`sj-docx-1-orm-10.1177_10944281261469432.docx`.
Its SHA-256 is
`ff92b677092062f9dcc02bd03627bed59e17729a599283b2bbcfe2ca71a3fbb8`.

The fixture was extracted from the first and only `w:tbl` in
`word/document.xml`. Row order and printed numerical precision are preserved.
`sabic` is the package's name for the source table's ABIC. The enumerated class
count denotes L2 classes in `additive_dispersion_l2`, L1 profiles in
`single_level_l1`, and L2 classes with five L1 profiles in
`cross_level_fixed_l1`. These block identifiers are concise descriptions,
not assertions that the package implements the fitted source models.

Run from the repository root, using only base R:

```sh
Rscript validation/houle-2026/check-fit-indices.R
```

To additionally check every fixture value against the supplied original DOCX:

```sh
Rscript validation/houle-2026/check-fit-indices.R /path/to/sj-docx-1-orm-10.1177_10944281261469432.docx
```

The script infers the sample-size convention from rounded BIC minus AIC,
checks joint rounding intervals for the likelihood and all four criteria,
and compares the current repository's `information_criteria()` output with
independent formulas. It writes `fit-index-results.csv` and
`fit-index-results.txt` in this directory.

The package call uses an explicitly artificial metadata adapter containing
published LL and parameter counts. Its placeholder posterior matrices permit
the diagnostics API to run; they are not estimates. Only the four likelihood
penalties under the individual convention are checked. No likelihood fitting,
posterior, entropy, standard-error, model-comparison-test, or whole-model
replication is claimed. A group count of 100 in the adapter is an arbitrary
sentinel that distinguishes conventions; it is not inferred from Table S1.

Every reported number has three decimal places. The check allows a half-unit
of rounding per printed criterion and twice that amount for the deviance
computed from printed LL. It requires a common unrounded deviance to satisfy
all five intervals, rather than treating a discrepancy of approximately
0.001 as an error. The source document itself is not redistributed here.
