# Audit 2026-11-01 — slice `plots`

Files: R/plot.R, R/plot-gg.R, R/plot-named.R, R/plot-docs.R.
94 functions (plus the palette/shape constants, the `utils::globalVariables()`
block and the `latents-plots` documentation object in R/plot-docs.R) checked;
4 defects fixed, 8 items logged, 1 item for another slice (an existing test
that asserts the old unclassed message). No estimate,
default, signature, table layout or printed table changed. Visual output
changed for the pairs view (all fits) and for the classification-certainty
and -uncertainty views of weighted fits only; rendered before/after at
`tmp/audit-1101-plots.html`.

Scratch evidence scripts (not part of the package):
`/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/f561d1fc-3262-4e98-94dd-fd3bab536edf/scratchpad/plots/`
(`fits.R` builds nine fits: two-level LPA, LPA, row-weighted LPA with zero
weights, group-weighted multilpa with zero-weight groups, LCA, weighted LCA,
mixed LPA/LCA, LTA, noise fit; `c3.R`/`c3b.R` strip means and entropy vs
tables; `c4.R` profile/bars/heatmap vs `get_results`; `c5.R` responses and
whiskers vs `parameter_inference()`; `c6.R` key shares, prevalence, group
views; `c7.R`/`c8.R` transitions and sequences; `c9.R`/`c10.R` pairs ellipses
and raincloud quartiles; `s1.R` categorical mean/median/mode by hand;
`e1.R`–`e4.R` enumeration criteria, minimum ring, quality views; `b1.R` error
classes; `p1.R` `what = "all"` on every fit; `audit-1101-plots.Rmd` the
render).

## Independent checks (layer data against tables; no defect unless noted)

- **Profile views** (two-level, weighted, mixed, LTA fits): line/point means
  == `get_results(x, "profiles")$mean`; whisker half-width / qnorm(.975) ==
  `mean_standard_error`; bar heights (layer 2 `y`) == means; heatmap tiles ==
  `get_results(x, "profiles", scale = "standardized")$mean`. All exact.
- **Categorical summaries**: mean/median/mode on numeric labels {1,2,4,5} and
  {2,3,10} recomputed by hand from `response_probabilities`; equal.
- **Responses**: `category = "last"`, `"first"`, `2` on LCA, weighted LCA,
  mixed fit == `get_results(x, "responses")$probability` for that category;
  whiskers == `conf_low`/`conf_high` of `parameter_inference(boundary = "fix")`.
- **Profile key**: shares == `get_results(x, "counts")$effective_proportion`
  (also with a noise component in the denominator); sizes view bar lengths and
  labels == `effective_count`/`effective_proportion` at both levels; total in
  the subtitle == table total (zero-weight units excluded: 1,259 of 1,422 rows).
- **AvePP**: tiles == `average_posteriors` table, weighted; "n = " ==
  `n_assigned` (zero-weight rows and groups excluded) at both levels.
- **Prevalence**: points == `profile_probabilities` mapped through the display
  order; class labels == `group_class_probability`.
- **Transitions**: tiles == `get_results(x, "transitions")$probability`.
- **Sequences**: every drawn cell (row -> group ordered by class then id,
  column -> sorted occasion) == `get_results(x, "sequences")$profile`;
  1,422/1,422 cells, 0 mismatches, on an LTA and a `multilpa(time = )` fit.
- **Pairs**: ellipse centres == `fit$means` to 3e-16 (diagonal, mixed and full
  covariance); ellipse second moment == covariance × chi2(.95, 2)/2.
- **Raincloud**: quartiles == quantiles of the indicator by assigned profile
  from `get_results(x, "assignments")`.
- **Enumeration**: every plotted criterion == `as.data.frame(grid)` column
  over converged candidates (single- and two-level grids, defaults and
  `"all"`); minimum ring at the column minimum; level-duplicate criteria
  collapsed only when equal; quality view == `.latents_grid_quality_table()`.
- **Palette**: only Okabe-Ito colours, the house diverging pair and greys.
  `globalVariables()` covers every aesthetic name used.

## Function table

| function | file:line | checked how | verdict |
|---|---|---|---|
| plot.multilpa | R/plot.R:3 | ran synthetic (9 fits, every view, `what = "all"`); b1.R | ok |
| .multilpa_check_plot_statistic | R/plot.R:32 | ran (`plot(what = "bars", statistic = "median")` -> latents_bad_argument) | ok |
| .multilpa_check_plot_arguments | R/plot.R:43 | ran b1.R | fixed (F2) |
| .multilpa_draw_view | R/plot.R:62 | ran every view | ok |
| .gg_every_view | R/plot.R:108 | ran `plot(what = "all")` on 9 fits | ok |
| .multilpa_require_case_posteriors | R/plot.R:148 | ran (1-profile, no group classes) | ok |
| .multilpa_plot_observations | R/plot.R:184 | ran raincloud/pairs/parallel; quartiles vs assignments | ok |
| plot.multilpa_enumeration | R/plot.R:211 | read; covered by test-plot.R, test-api-compat-0921.R | ok |
| .latents_tree_arguments | R/plot.R:270 | read | ok |
| .multilpa_enumeration_plot_advice | R/plot.R:277 | read | ok |
| .latents_grid_plot_advice | R/plot.R:285 | read | ok |
| .latents_grid_criterion | R/plot.R:303 | read; covered by test-grid-contract.R | ok |
| plot_enumeration | R/plot.R:320 | ran e2–e4.R, b1.R | fixed (F2) |
| .multilpa_distinct_criteria | R/plot.R:453 | ran e4.R on one- and two-level grids | ok |
| plot.multilpa_covariates | R/plot.R:481 | read; covered by test-plot.R | ok |
| .multilpa_default_view | R/plot.R:535 | ran (LCA default draws responses) | ok |
| .multilpa_scored_profiles | R/plot.R:547 | ran s1.R | ok |
| plot_views | R/plot.R:555 | ran | ok |
| .multilpa_mean_error_matrix | R/plot.R:595 | ran c4.R (SE == table) | ok |
| .multilpa_response_error_table | R/plot.R:627 | ran c5.R (bounds == inference) | ok |
| .multilpa_plot_named_view | R/plot-named.R:19 | ran (named verbs on 9 fits) | ok |
| .latents_plot_kind | R/plot-named.R:46 | read; covered by test-plot-named.R | ok |
| .latents_named_plot_registry | R/plot-named.R:75 | read; covered by test-plot-named.R | ok |
| .latents_named_plot_lookup | R/plot-named.R:113 | read | ok |
| .latents_plot_offered_views | R/plot-named.R:125 | read | ok |
| .latents_refuse_named_plot | R/plot-named.R:151 | ran `plot_transitions(lpa fit)` | ok |
| .latents_draw_named | R/plot-named.R:190 | ran | ok |
| .multilpa_plot_method | R/plot-named.R:212 | read | ok |
| .multilpa_supplied_arguments | R/plot-named.R:226 | ran | ok |
| .gg_require_fit | R/plot-named.R:241 | ran (`plot_profiles(data.frame())` -> latents_not_a_fit) | ok |
| .multilpa_plot_named_views | R/plot-named.R:275 | ran `plot_profiles(what = "all")` | ok |
| .multilpa_plot_grouped | R/plot-named.R:314 | ran | ok |
| plot_profiles | R/plot-named.R:324 | ran b1.R; test-audit-1101-plots.R | fixed (F2) |
| plot_responses | R/plot-named.R:350 | ran c5.R | ok |
| plot_sequences | R/plot-named.R:360 | ran c8.R | ok |
| plot_prevalence | R/plot-named.R:369 | ran c6.R | ok |
| plot_diagnostics | R/plot-named.R:377 | ran c3b.R, c6.R, b1.R | fixed (F1, F2) |
| plot_transitions | R/plot-named.R:422 | ran c7.R | ok |
| plot_tree | R/plot-named.R:430 | ran (`level = "foo"` classed); drawing in R/tree.R | ok |
| plot_selection | R/plot-named.R:451 | ran b1.R | fixed (F2) |
| .gg_require | R/plot-gg.R:34 | read | ok |
| .gg_draw_on_call | R/plot-gg.R:47 | read | ok |
| .gg_available | R/plot-gg.R:57 | read | ok |
| .gg_theme | R/plot-gg.R:67 | read | ok |
| .gg_top_legend | R/plot-gg.R:94 | read | ok |
| .gg_percent | R/plot-gg.R:103 | read | ok |
| .gg_profile_key | R/plot-gg.R:115 | ran c6.R (shares == counts table, noise) | ok |
| .gg_scales | R/plot-gg.R:141 | ran | ok |
| .gg_adjusted_note | R/plot-gg.R:164 | read | ok |
| .gg_titles | R/plot-gg.R:173 | read | ok |
| .gg_spread | R/plot-gg.R:186 | read (order kept, min gap) | ok |
| .gg_ink | R/plot-gg.R:199 | read | ok |
| .gg_probability_fill | R/plot-gg.R:208 | read | ok |
| .gg_probability_scale | R/plot-gg.R:217 | read | ok |
| .gg_tile_theme | R/plot-gg.R:227 | read | ok |
| .gg_group_class_text | R/plot-gg.R:237 | read | ok |
| .gg_mean_frame | R/plot-gg.R:252 | ran c4.R | ok |
| .gg_series_plot | R/plot-gg.R:293 | ran (layer data) | ok |
| .gg_category_scores | R/plot-gg.R:356 | ran s1.R | ok |
| .gg_label_scores | R/plot-gg.R:374 | ran s1.R | ok |
| .gg_scorable_categories | R/plot-gg.R:391 | ran | ok |
| .gg_categorical_summary_frame | R/plot-gg.R:397 | ran s1.R (by hand) | ok |
| .gg_view_profiles | R/plot-gg.R:423 | ran c4.R, s1.R | ok |
| .gg_view_bars | R/plot-gg.R:464 | ran c4.R | ok |
| .gg_view_heatmap | R/plot-gg.R:503 | ran c4.R | ok |
| .gg_view_response_heatmap | R/plot-gg.R:538 | ran c5.R | ok |
| .gg_view_polytomous_heatmap | R/plot-gg.R:600 | ran c5.R | ok |
| .gg_view_raincloud | R/plot-gg.R:634 | ran c9.R | logged (L3) |
| .gg_view_responses | R/plot-gg.R:739 | ran c5.R, b1.R | fixed (F3) |
| .gg_match_bounds | R/plot-gg.R:815 | ran c5.R | ok |
| .gg_undrawn_caption | R/plot-gg.R:839 | ran (mixed fit) | ok |
| .gg_view_probabilities | R/plot-gg.R:859 | ran c6.R | logged (L1) |
| .gg_view_sequences | R/plot-gg.R:901 | ran c8.R | logged (L4) |
| .gg_view_transitions | R/plot-gg.R:1005 | ran c7.R | ok |
| .gg_level | R/plot-gg.R:1064 | ran c3b.R | fixed (F1) |
| .gg_has_group_classes | R/plot-gg.R:1113 | read | ok |
| .gg_view_sizes | R/plot-gg.R:1119 | ran c2.R, c6.R | ok |
| .gg_capitalize | R/plot-gg.R:1151 | read | ok |
| .gg_case_values | R/plot-gg.R:1162 | ran c3b.R | fixed (F1) |
| .gg_case_strip | R/plot-gg.R:1180 | ran c3b.R | fixed (F1) |
| .gg_view_posteriors | R/plot-gg.R:1241 | ran c3b.R | fixed (F1) |
| .gg_view_entropy | R/plot-gg.R:1261 | ran c3b.R | fixed (F1) |
| .gg_view_avepp | R/plot-gg.R:1282 | ran c2.R, c6.R | ok (now reads F1's weights) |
| .gg_observed_cases | R/plot-gg.R:1341 | ran | ok |
| .gg_view_parallel | R/plot-gg.R:1355 | ran | ok |
| .gg_ellipse | R/plot-gg.R:1413 | ran c9.R/c10.R | ok |
| .gg_pair_covariance | R/plot-gg.R:1429 | ran (diagonal and full) | ok |
| .gg_view_pairs | R/plot-gg.R:1438 | ran c10.R, p1.R | fixed (F4); logged (L5) |
| .gg_enumeration_frame | R/plot-gg.R:1570 | ran e2.R | ok |
| .gg_view_enumeration | R/plot-gg.R:1596 | ran e2–e4.R | logged (L6) |
| .gg_flow_band | R/plot-gg.R:1704 | read (used by R/tree.R) | ok |
| .gg_plots | R/plot-gg.R:1718 | read | ok |
| print.latents_plots | R/plot-gg.R:1724 | read (restores devAskNewPage on exit) | ok |
| .gg_view_enumeration_quality | R/plot-gg.R:1750 | ran e2.R | logged (L7) |
| constants `.gg_okabe_ito`, `.gg_reserved_grey`, `.gg_class_palette`, `.gg_shapes`; `utils::globalVariables()` | R/plot-gg.R:13-30 | read | ok |
| `latents-plots` (documentation object) | R/plot-docs.R:1-408 | read against behaviour; @return covers every export; @examples tidy | ok |

## Fixed

**F1. Weighted fits: the certainty and uncertainty views did not show the
table's weighted numbers.** `plot_diagnostics(fit, "posteriors")` labels each
class "mean x.xx" and its subtitle says that mean "is its average posterior
probability"; `"entropy"` prints "relative entropy = 1 - mean = x.xxx". Both
were plain means over the classified units, while the `classification` and
`entropy` tables weight by the sampling weights (review 2026-10-26 S-3/S-4).
Reproduction (group-weighted `multilpa`, a third of students weight 0, rest
1 or 3): group strip "mean 0.99 / mean 0.97" vs table average_posterior
0.983 / 0.976; group entropy subtitle 0.905 vs table 0.918; individual
entropy 0.813 vs 0.819; row-weighted `lpa`: entropy 0.698 vs 0.699, profile-2
strip "mean 0.78" vs 0.785. Zero-weight units were already excluded (the
views use the informative rows). Change: `.gg_level()` now returns each
posterior row's weight from `.multilpa_level_weights()` (the tables' own
helper; R/plot-gg.R:1064); `.gg_case_values()` carries it (R/plot-gg.R:1162);
`.gg_case_strip()` takes a weighted mean under weights and the plain mean
otherwise, so unweighted output is bit-identical (R/plot-gg.R:1180);
`.gg_view_entropy()` prints `.multilpa_relative_entropy(posteriors, weights)`,
the table's value (R/plot-gg.R:1261); `.gg_view_avepp()` reuses the same
weights. After: strip means equal `average_posterior` to 1e-12 at both levels
on all four fits, subtitles equal `relative_entropy` to the printed digits.
Guard: test-audit-1101-plots.R, "weighted strip means and entropy use the
sampling weights (individuals)" (hand-computed weighted means and entropy from
the posteriors and weights) and "(groups)"; "unweighted strip means are the
plain means, as before". Sabotage (unweighted means/entropy restored): 9
failures.

**F2. Bad argument values raised unclassed errors.** `plot_profiles(fit,
what = "foo")`, `statistic = "foo"`, `scale = "foo"`, `plot_diagnostics(what
= "foo")`, `plot_selection(what = "foo")` used `match.arg()`;
`plot(fit, labels = NA)`, `intervals = NA`, `category = c(1, 2)`,
`plot_enumeration(combine = NA)`/`labels`/`mark_minimum` and `cutoff = 2`
used `stopifnot()`. All returned `simpleError` (b1.R output), where the
documented contract is `latents_bad_argument`. Change: `.latents_match_arg()`
with the formals' choices (keeps the full-vector default and unique
abbreviations, as `match.arg()` did) in R/plot-named.R:331-334, 383, 456;
`.latents_require()` (the package's classed `stopifnot()`, same messages) in
R/plot.R:49, 339, 368, 395. Guard: "the named plot verbs and plot() refuse bad
arguments by class" (also checks `what = "bar"`, `"avep"` still match).
Sabotage (`stopifnot`/`match.arg` restored): the test errors.

**F3. A fractional `category` index was truncated.** `plot(lca_fit, what =
"responses", category = 2.7)` drew category 2 with no message
(`as.integer(2.7)`). Change: a numeric `category` must be whole (within
`all.equal()`), else it names no category and raises
`latents_unknown_category` (R/plot-gg.R:748-757). Guard: "a fractional
category index is refused, a whole one is used" (also checks `category = 2`
draws column 2 of every block). Sabotage: 1 failure.

**F4. Pairs view encoded the profiles by colour alone.** Cases were points
coloured by profile (shape 19), ellipses coloured paths; no shape, linetype or
label paired. Change (R/plot-gg.R:1438): cases take the profile's marker
shape and fill (the shapes every other profile view uses, from
`.gg_profile_key()`); each ellipse is centred on a marker of that shape at the
profile's mean; one merged colour/fill/shape legend; subtitle says so. Layer
data: centre markers == `fit$means`; one shape per profile == key shape.
Guard: "pairs marks each profile by shape and each ellipse by its mean".
Sabotage (shape mapping removed): 1 failure. Rendered before/after:
`tmp/audit-1101-plots.html`.

## Logged, not fixed

- L1. `.gg_view_probabilities()` subtitle "over %d observed groups" uses
  `x$n_groups`, which counts zero-weight and all-empty groups (106 vs 71
  carrying groups in the probe). It equals `get_results(x, "model")$n_groups`;
  changing it would change a printed number on valid input. Owner decision.
- L2. Categorical "median" uses `which(cumsum(p) >= 0.5)[1]`; a cumulative
  probability within 1e-16 of 0.5 could pick the neighbouring category.
  Floating-point edge only.
- L3. Raincloud (and the pairs diagonal densities) summarise the drawn cases
  unweighted, zero-weight units included (their rows are kept, per the owner
  rule). They describe the drawn observations, not a table; no table to agree
  with.
- L4. Sequence plot: on a dense grid (> 60 groups or > 40 occasions) the
  profile digits are not printed, so profiles are told apart by fill colour
  alone. Inherent to an index plot at that density.
- L5. Pairs diagonal densities remain colour-coded only (their profile is
  readable from the points and centre markers in the same column).
- L6. `.gg_view_enumeration()`: a converged candidate with a non-finite
  criterion is drawn as an "x" on the floor, the mark the subtitle calls "did
  not converge"; the subtitle count itself is correct. None seen on realistic
  grids.
- L7. `.gg_view_enumeration_quality()` puts a profile grid with a single
  number of profiles and several group-class counts in one x column, where the
  criteria view switches to a group-class axis. Readable, numbers correct.
- L8. `plot_enumeration(what = "aic", cutoff = "a")` does not refuse `cutoff`,
  which only the quality and selection views use.

## Needs another slice

- tests/testthat/test-plot-named.R:142 asserts `match.arg()`'s message:
  `expect_error(plot_profiles(fits$lpa, statistic = "average"), "should be one of")`.
  After F2 the refusal is `latents_bad_argument` ("must be one of ..."), so
  this one expectation errors. Not edited (outside the assigned files).
  Replacement: `expect_error(plot_profiles(fits$lpa, statistic = "average"),
  class = "latents_bad_argument")`.
- equivalence/engine-golden/compare.R: besides the 13 known
  `replicate_rows` -> `replicate_identity` cases, 22 cases differ, all at
  `probes/tables/value/group_classes@display` ("length 4 vs 5"): the 21
  `reg_*` mixture-regression cases and `growth_lcga`. That is a table probe of
  mixture regressions/growth mixtures, not a plot; R/plot*.R build no tables
  and these files were also modified by other slices
  (R/mixture-regression-*.R, R/get-data.R). Owner of that table to check.

## Test runs

Baseline before any change (test-plot.R ... test-review-1015-plots.R): TOTAL
pass 753 fail 0 error 0 skip 0.

After the fixes, `NOT_CRAN=true`, one R process, `testthat::test_file()` per
file:

```
                                  file pass fail error skip warn
               test-audit-1101-plots.R   33    0     0    0    0
                           test-plot.R  161    0     0    0    0
                        test-plot-gg.R   56    0     0    0    0
                     test-plot-named.R  135    0     1    0    0
                    test-plot-polish.R  109    0     0    0    0
     test-plot-categorical-summaries.R  134    0     0    0    0
        test-audit-plots-predictions.R   86    0     0    0    0
      test-audit-drawing-integration.R    6    0     0    0    0
                test-math-plot-audit.R   14    0     0    0    0
              test-review-1015-plots.R   49    0     0    0    0
                test-api-compat-0921.R   44    0     0    0    0
                    test-compat-0921.R   10    0     0    0    0
        test-audit-extension-repairs.R   43    0     0    0    0
        test-audit-reporting-repairs.R   25    0     0    0    0
                    test-categorical.R   75    0     0    0    0
 test-empty-classification-consumers.R   39    0     0    0    0
                    test-enumeration.R  315    0     0    0    0
                      test-fit-verbs.R  255    0     0    0    0
                       test-get-data.R  230    0     0    0    0
                  test-grid-contract.R  237    0     0    0    0
               test-informative-rows.R   57    0     0    0    0
             test-mixture-regression.R   73    0     0    0    0
         test-post-submission-review.R   60    0     0    0    0
             test-posterior-rounding.R    7    0     0    0    0
                test-release-surface.R   88    0     0    0    0
                 test-report-surface.R   40    0     0    0    0
               test-review-1015-lead.R   19    0     0    0    0
        test-review-1015-transitions.R   38    0     0    0    0
       test-review-1026-selection-r2.R   41    0     0    0    0
     test-review-1027-error-contract.R  216    0     0    0    0
                       test-tna-fits.R   11    0     0    0    0
                            test-tna.R   62    0     0    0    0
                           test-tree.R  131    0     0    0    0
        test-review-1027-diagnostics.R  158    0     0    0    0
TOTAL pass 3057 fail 0 error 1 skip 0
```

The one error is test-plot-named.R:143 (the message assertion above):
`<latents_bad_argument> `statistic` must be one of "mean", "median", "mode", not "average".`

Sabotage (each fix broken once, test-audit-1101-plots.R): weights -> pass 24
fail 9; classes -> pass 19 error 1; category -> pass 32 fail 1; pairs -> pass
32 fail 1. All restored.

`devtools::document()`: no change to NAMESPACE or man/ from this slice
(roxygen edits were in `@noRd` helpers only; the 8 man files already modified
by other slices showed the same diffstat before and after).
`codetools::checkUsage()` over the 94 functions: every non-local variable is
declared in `utils::globalVariables()`.
