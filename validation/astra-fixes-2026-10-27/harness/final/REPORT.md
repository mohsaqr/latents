# Final combined full-output validation

This is a combined validation matrix, not a second uninterrupted full run.
The original full run executed all 15 suites. Four affected suites were rerun:
`depmixs4` and `hard-depmixs4` after completing mixed residual references,
`invariance` after completing excluded-observation references, and `refusals`
after restoring the classed empty-data error. The other 11 suites retain their
original measured comparisons. Each CSV row identifies its validation_run.

The strict acceptance gate passes; all 15 registered suites are present.

- suites:      15
- quantities:    9068
- values: 3143431
- disagreements:       0
- incomplete_outputs:       0
- warnings:       4

Original and rerun REPORT.md/SESSION.txt files remain in the sibling folders.
No reference value, tolerance, status or exemption was changed during combination.
