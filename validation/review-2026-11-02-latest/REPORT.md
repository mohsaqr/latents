# Review of the latest uncommitted changes (2026-11-02)

Reviewer: Codex (static reading; its sandbox could not start R, so it ran nothing and
could not write this file). Every finding below was then reproduced, fixed and
re-run in the main session; the "observed" lines are from those runs.

Scope: `print.summary_multilpa_enumeration()` and `.multilpa_print_enumeration_warnings()`
(R/enumeration.R); `get_candidates()`, `get_best()`, `get_deltas()`, `get_knees()`,
`get_tree()` (R/wildcards.R); compatibility with CRAN 1.1.9 (commit 8518dffe).

## 1. Failed candidates moved the grid's edges (P2) - fixed

`R/enumeration.R`, edge logic. The edges came from the converged candidates. With
counts 1:3 requested, 3 failed and every minimum at 2, the print warned that the
minimum was at the edge; with only 3 converged and the minimum at 3, it said nothing.

- Before: edge warning TRUE for the first case, FALSE for the second.
- After: FALSE, TRUE. The edges are those of the requested grid.

## 2. `get_candidates()` needed a buildable tree (P2) - fixed

`R/wildcards.R`, `.latents_wild_grid_table()`. Availability was read from
`get_results(x, "all")`, which builds every table. A regression grid with one usable
candidate cannot build its tree (`latents_no_tree`), so `get_candidates()` failed.

- Before: `[latents_no_tree] A tree needs at least two numbers of classes for one model.`
- After: returns the table. The table list is read from each grid's catalogue without
  building anything.

## 3. `get_tree()` refused family and transition grids (P2) - fixed

Same cause: those grids' `"all"` holds only `candidates` and `best`, though `"tree"`
is a table of theirs.

- Before: `[latents_bad_argument] get_tree() has no table for an enumeration across model families.`
- After: `identical(get_tree(g), get_results(g, "tree"))` is TRUE.

## 4. The minima table ignored `rows` (P3) - fixed

`head(shown, max(rows, nrow(shown)))` printed every minimum whatever `rows` was.

- Before: 14 rows for `rows = 0`, `1` and the default.
- After: default 14 (every criterion), `rows = 1` one row, `rows = 0` headings only.

## Checked and found correct (Codex, static)

- Constant dimensions, single-level grids, a lower bound above one; NA winners filtered.
- No printer error identified for all-failed, LCA or model-structure grids.
- No `what =` advice in the rewritten print; `digits` reaches the minima; argument
  validation intact.
- `get_candidates()` keeps the retained values and drops only all-NA or all-empty columns.
- The six fixes of validation/review-2026-11-02-wildcards/REPORT.md are present, each
  with a regression test.
- Against `git archive 8518dffe`: no removed exports or S3 registrations, no moved
  arguments, no `get_results()` column removed.

## Tests after the fixes (main session)

test-wildcards.R 207/0 (three new blocks, one per finding group), test-tidy-classes.R
61/0, test-enumeration.R 315/0, test-api-compat-0921.R 44/0.
