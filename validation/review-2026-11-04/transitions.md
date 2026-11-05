# Review 2026-11-04: transitions

Probes are in `validation/review-2026-11-04/probes/transitions-*.R`. Nothing in R/, man/, tests/ or the session files was changed. The mutation checks override functions in memory only.

## 1. What was reviewed

| File | Depth |
|---|---|
| R/lta-engine.R | Read in full. Probed: likelihood, posteriors, parameter counts, flags. |
| R/transitions.R | Read in full. Probed: likelihood, expected counts, prevalence, tables. |
| R/transition-inference.R | Read in full. Probed: scores, SEs, natural-scale Jacobian. |
| R/lta-inference.R | Read in full. Probed: scores, SEs, identification, boundary rule. |
| R/lta-bootstrap.R | Read in full. Probed: units, relabelling, resamples that change shape. |
| R/lta-compare.R | Read lines 1-470. `.lta_simulate()` and `.lta_nested()` by reading only; no new simulation recovery. |
| R/lta-methods.R | Read in full. Tables probed. Plot method not checked. |
| R/tna.R | Read. `get_tna()` and `get_group_tna()` on homogeneous fits probed. Modal-sequence path read only. |
| R/sequences.R | Read only, no probe of my own. |
| R/lta-starts.R | NOT read. It was exercised by every fit, but no line-by-line check. |

Also not reviewed:
- `.lta_enumerate()` beyond its argument checks, and the enumeration tables.
- The `.inference_sandwich()` formula and the robust and OPG paths, beyond confirming they run.
- `.lta_bootstrap_lrt()` beyond reading it. No LRT size work, which is settled.
- Ordinal and count measurement beyond score and count checks.
- R/lta-three-step.R, which is out of scope.
- External software agreement, which the verification matrix marks strong for `lta()`. I did not repeat it.

Limit on the rank results: they cover designs without covariates. With covariates the same models may be identified, and I did not test that.

## 2. Findings

Counts: CRITICAL 0, MAJOR 3, MINOR 6, EDGE 0.

### F1. MAJOR, CONFIRMED. Effectively empty origins are reported as estimates and counted as parameters

**Where:** R/lta-engine.R:1370 (`!(totals > 0)`), R/transitions.R:972 (`<= 0`), docs at R/transitions.R:656-660.

**What is wrong:** the TRANS-1 repair and the documentation treat an origin as empty only when its expected moves are exactly zero in floating point. An origin with 1e-9 or 1e-245 expected moves has no information, yet it is reported as estimated. It is also counted in `n_parameters`, AIC and BIC. The documented scenario, "a profile that no unit occupies before the unit's last occasion", therefore behaves as documented only when the numbers underflow.

**Evidence, realistic design:** `transitions-near-empty-realistic.R`. K = 3, T = 4, n = 300, 5 indicators, and a state reached only at the final wave. Three seeds gave identical results.

Homogeneous engine, row out of the top profile:

| to | probability | expected_count | estimated |
|---|---|---|---|
| 1 | 1e-10 | 9.3e-21 | TRUE |
| 2 | 0.4928 | 1.2e-09 | TRUE |
| 3 | 0.5072 | 1.3e-09 | TRUE |

- `n_parameters` is 38. Without that row it would be 36, so BIC carries 2 log(300) = 11.4 too much.
- There were no warnings at all.

General engine with `transitions = "occasion"`:
- The row is 0.1255 / 0.1375 / 0.7370 at occasions 2, 3 and 4 alike, with `estimated = TRUE`.
- Expected moves are 2.2e-21, 2.9e-21 and 5.9e-21.
- Identical values at three occasions mean the logits never moved. They also fit the start formula `0.7 I + 0.3 pi`, with pi = 0.418, 0.458, 0.123 summing to 0.999. This is the TRANS-1 symptom, start values reported as estimates, surviving the exact-zero repair.
- I inferred the start-formula match by arithmetic. I did not print the start itself.
- `n_parameters` is 50 against 44, so BIC carries 6 log(300) = 34.2 too much.
- A `latents_boundary` warning does appear, but it concerns the zero initial probability and says nothing about this row.

**Evidence, separation sweep:** `transitions-near-empty.R`. A state appears only at T = 2, two indicators, homogeneous engine.

| separation | row out of that state | expected count | n_parameters | estimated |
|---|---|---|---|---|
| 6 | 1e-10 / 1 | 1.9e-10 | 11 | TRUE |
| 12 | 1e-10 / 1 | 2e-34 | 11 | TRUE |
| 30 | 1e-10 / 1 | 1e-245 | 11 | TRUE |
| 60 | 0.5 / 0.5 | 0 | 10 | FALSE, warned |

A reader of the first three rows sees an absorbing state. By my reading of the M-step the row is set by whichever person has the largest tiny weight, which I did not separately trace.

**Mitigation:** Wald inference is refused on these fits by the probability floor, confirmed in the sweep. The table, the criteria and `get_tna()` are not protected.

**Who is affected:** any fit with a state that appears only late. PLAUSIBLE but not probed: a group class that never visits a state, and second-order origin pairs, where near-empty cells should be more common.

**Smallest fix:** compare the origin total with a small tolerance instead of zero, in both engines, and use the existing empty-row path.

**Changes reported numbers:** yes, on affected fits only: those rows, `n_parameters`, AIC, BIC, and possibly enumeration choices. The owner must choose the tolerance. Probed rows sit at 1e-9 or below, so a tolerance around 1e-6 expected moves would catch them, but that value is a judgement and has not been tested on other data.

### F2. MAJOR, CONFIRMED. Unidentified specifications are fitted and reported without a warning

**Where:** no check exists in `.lta_setup()` or `.lta_setup_general()`. R/transitions.R:872 tests only for singleton groups. Counts are at R/transitions.R:977 and R/lta-engine.R:1472.

**What is wrong:** several specifications have more free parameters than the latent path distribution can identify. Their shares and matrices are arbitrary, yet they are reported as estimates and counted in full. `n_best_replicated` reads as full replication, because the likelihood is flat across the solutions.

**Evidence, rank:** `transitions-identification.R` is plain mathematics with no package code. It gives the rank of d p(path) / d logits at random points, against the count `lta()` uses. Cases with excess:

| K | T | specification | counted | rank | excess |
|---|---|---|---|---|---|
| 2 | 2 | 2 group classes | 7 | 3 | 4 |
| 3 | 2 | 2 group classes | 17 | 8 | 9 |
| 2 | 2 | stayer | 5 | 3 | 2 |
| 3 | 2 | stayer | 11 | 8 | 3 |
| 2 | 3 | 2 group classes, occasion transitions | 11 | 7 | 4 |
| 3 | 3 | 2 group classes, occasion transitions | 29 | 23 | 6 |
| 2 | 3 | 2 group classes, order 2 | 15 | 7 | 8 |
| 2 | 3 | stayer, order 2 | 9 | 7 | 2 |

No excess at: K=2 T=3 H=2 homogeneous, K=3 T=3 H=2, T=4 H=2 or H=3, T=5 H=2, T=3 stayer, T=3 or 4 stayer with occasion transitions, T=4 H=2 with occasion transitions, T=3 order 2, T=4 H=2 order 2, T=4 stayer order 2, T=4 order 2 with occasion transitions.

**Evidence, fits:** `transitions-identification-fits.R`, n = 300.

| specification | logLik | n_parameters | replicated | warnings | Wald |
|---|---|---|---|---|---|
| T=2, H=1 | -2843.09933 | 15 | 6/6 | none | given |
| T=2, H=2 | -2843.09933 | 19 | 6/6 | none | refused |
| T=2, stayer | -2843.09933 | 17 | 6/6 | none | **given, SE 297 to 760** |
| T=3, H=1, order 2 | -4153.33832 | 19 | 6/6 | none | given |
| T=3, H=2, order 2 | -4153.33832 | 27 | 6/6 | none | refused |
| T=3, stayer, order 2 | -4153.33832 | 21 | 6/6 | none | refused |
| T=3, H=2, occasion transitions | -4153.33832 | 23 | 6/6 | none | refused |

**Evidence that the stayer share is arbitrary:** `transitions-stayer-t2.R`. The fit reports a stayer share of 0.381. I built other parameter points from its implied latent joint distribution and evaluated them with the package's own `.lta_expectation()`. Stayer shares of 0.05, 0.20 and 0.40, each with stayer initial 0.3 and 0.6, all give a log likelihood of -2843.09932940, a difference of 0.00e+00. On pre/post data a reader would report "38% stayers", and the data say nothing about that number.

**Singularity check:** at the T = 2 stayer fit the information eigenvalues run from 1.0e-6 to 497. That passes `min_le_eps_max`, so Wald SEs are returned instead of `latents_singular_information`. The H = 2 case is refused. Which one you get depends on numerical noise.

**Effect on selection:** the over-count alone did not flip a choice here. The likelihoods are equal, so BIC favours the smaller model either way (5771.76 against 5794.57 in `enumerate_classes()`). The harm is the arbitrary estimates, not the selection.

**Smallest fix:** a classed `latents_unidentified` warning at setup, using the fixed list above or a rank computed at a random point. A simple `count > K^T - 1` test would miss K=3 T=3 H=2 occasion, where 29 is compared with 26 but the rank is 23. It should apply only without covariates, because I did not test identification with them.

**Changes reported numbers:** a warning changes none. Correcting `n_parameters` would change AIC and BIC, and that is the owner's call.

### F3. MAJOR, CONFIRMED. Bootstrap inference crashes when a resample changes the coefficient set

No wrong number is reported, so this is outside the literal MAJOR wording. I rate it MAJOR because a supported path is unusable on ordinary unbalanced data and fails with an internal, unclassed error, possibly after a long run.

**Where:** R/lta-bootstrap.R:184 and :197 (alignment), reached from R/simulation-service.R:151, where `extract(fit)` sits outside the `tryCatch`. The intended refusal at R/lta-bootstrap.R:166, "could not be matched", runs after alignment, so it is not reached in these cases.

**Evidence:** `transitions-boot-crash.R`, 60 persons.

| fit | result |
|---|---|
| Occasion transitions, one person longer than the rest, iter 20 and 60 | `simpleError`: "the two signatures must have the same shape" |
| Occasion measurement, one longer person, iter 30 | `simpleError`: "`blocks` and `scales` must be lists of the same length" |
| Transition covariate, factor with a one-person level, iter 30 | `simpleError`: "the two signatures must have the same shape" |
| Initial covariate, same factor, iter 30 | `simpleError`: "the two signatures must have the same shape" |
| Homogeneous engine, one longer person, iter 30 | completes, 30 of 30 valid |

A single longest person is absent from about 36% of resamples, so with the default 199 iterations a crash is close to certain whenever few persons hold the longest sequence or a factor level. Short runs can pass by chance: in `transitions-misc2.R` iter 15 with seed 1 completed on the same kind of data, and the first 12 draws I traced all contained the long person.

**Smallest fix:** compare the names before aligning, so the replicate takes the existing failure path.

**Changes reported numbers:** none today, because this path errors. But dropping those replicates conditions the bootstrap on the rare person or level being present, so `occasion_5`-type intervals would rest on a selected subset. The owner should decide between that and refusing such coefficients.

### F4. MINOR, CONFIRMED. `second_order_transitions` averages over occasion-specific matrices without saying so

**Where:** R/lta-methods.R:148 onward; man/get_results.Rd:283-284 says "averaged over the groups".

**Evidence:** `transitions-misc.R`, `order = 2, transitions = "occasion"`. The table has 8 rows and no `occasion` column. Row 1->1->2 reads 0.1948, while the coefficients are -1.322 at occasion 3 and -1.523 at occasion 4. Those give 0.2105 and 0.1790, whose mean is the reported value. First-order tables do carry `occasion`.

**Fix:** add an `occasion` column, or state the averaging in the documentation.

**Changes reported numbers:** an added column changes the table layout. A documentation fix changes nothing.

### F5. MINOR, CONFIRMED. Under `order = 2`, the `"transitions"` table covers only the first move but does not say so for homogeneous fits

**Where:** R/lta-methods.R:103 and the `transition_occasions` line after it.

**Evidence:** `transitions-misc.R`. A homogeneous second-order fit reports `occasion = NA`, the same label a homogeneous first-order fit uses to mean every occasion. By the code these logits govern only the move into occasion 2. With `transitions = "occasion"` the same table correctly says occasion 2. `plot.multilpa_lta()` draws this table under the title "Transition probabilities", by reading only.

**Fix:** report `occasion = 2`, or document it.

**Changes reported numbers:** the label only.

### F6. MINOR, CONFIRMED. The warning and the help page misstate second-order behaviour at an empty inner occasion

**Where:** warning text at R/lta-engine.R:987-997; R/lta-engine.R:1316 claims the warning states the upper bound; help at R/transitions.R:661-664.

**Evidence:** `transitions-counts-and-flags.R`, `order = 2, transitions = "occasion"`, FIML, occasion 3 empty.
- The warning says the move-into logits "are not counted in n_parameters, AIC or BIC and have no standard error".
- In fact `n_parameters` is 23, those logits are counted, and they are flagged `estimated = TRUE`. The second-order move into occasion 3 still has its coordinates.
- The information is singular, so inference is refused for all coefficients.
- The upper-bound statement promised in the code comment is not in the warning.

**Fix:** a separate second-order wording.

**Changes reported numbers:** none.

### F7. MINOR, CONFIRMED. A probability boundary is flagged differently by the two engines on one model

**Where:** R/transitions.R:999 (variance only) against R/lta-inference.R:770-796.

**Evidence:** `transitions-counts-and-flags.R`, a move that never occurs. Both engines give log likelihood -1921.699 and 11 parameters. The homogeneous fit has `boundary = FALSE` and no warning. The general fit, with weights of 1, has `boundary = TRUE` and warns `latents_boundary`. Wald inference is refused on the homogeneous fit by the floor; I did not call it on the general fit, whose code refuses a boundary fit.

**Fix:** set the flag and warn in the homogeneous engine.

**Changes reported numbers:** the `boundary` field and the `boundary` column of the tables and grids. No estimates change.

### F8. MINOR, CONFIRMED. One row at an extreme covariate value marks an interior fit as a boundary fit and withholds every SE

**Where:** R/lta-inference.R:796, which takes the minimum over person-level rows.

**Evidence:** `transitions-boundary-false.R`, slope about 2 per SD, n = 300.
- Without the extreme row: `boundary = FALSE`, slope SEs 0.255 and 0.274, smallest row probability 2.8e-4.
- With one destination row at -6.5 SD: smallest row probability 3.0e-7, `boundary = TRUE`, all Wald SEs refused, and the warning says a probability "is estimated at zero (a move or start that never occurs)".
- The coefficients barely move (-1.287, 2.106 to -1.276, 2.115) and the fit is interior.
- Strong slopes without an outlier did not trigger it: slope 3 with n = 1000 returned SEs, and slope 3 with n = 300 had a smallest probability of 3.1e-6. This needs an outlier or a very strong effect.

**Fix:** judge divergence on the logit coefficients in the fitting chart rather than on a person-row probability.

**Changes reported numbers:** SEs would appear on such fits.

### F9. MINOR, CONFIRMED. A design term ending in a dot loses its identification and its SE

**Where:** R/lta-inference.R:618 and R/lta-methods.R:38. `strsplit(..., fixed = TRUE)` drops a trailing empty piece.

**Evidence:** `transitions-boundary-false.R`, factor level `"Zed."`.
- The term appears as `titleZed`, with `standard_error = NA` and `estimated = FALSE`.
- There are 11 identified coordinates against `n_parameters` of 13, with no warning.
- Internal dots are fine: `my.var` was correct.
- Levels such as `"Dr."` are ordinary data.

**Fix:** split with a parser that keeps the trailing piece.

**Changes reported numbers:** SEs would appear on such fits. The likelihood and `n_parameters` would not change.

### Noted, not counted

- A stayer row has `estimated = TRUE` although its identity transitions are fixed by the model (R/lta-methods.R, `.lta_logits_identified()`, deliberate per its comment). Seen in `transitions-identification-fits.R`. The help does not say what `estimated` means for stayers. A one-line documentation point.
- Under both layouts, occasion-specific parameters count from each group's own first row, so a late entrant's moves are read against earlier occasion labels. This is documented in `lta()` (`occasions`), and I do not report it as a defect. Likewise grid steps assume equal spacing and only positions seen in the data exist, which the help states.

## 3. Verified sound

**Likelihood and posteriors.** `transitions-brute.R` enumerates every path in plain base R, with covariate designs built by hand from the data. Gaps and late entry are included.

| specification | LL difference | max posterior difference | max count difference |
|---|---|---|---|
| homogeneous, grid gaps, H=1 | 0 | 1.3e-15 | 2.1e-14 |
| homogeneous, observed, H=1 | 5.7e-14 | 1.0e-15 | 1.4e-14 |
| homogeneous, grid gaps, H=2 | 0 | 1.2e-15 | 1.4e-14 |
| homogeneous, observed, H=2 | 1.1e-13 | 1.2e-15 | 5.3e-15 |
| general, occasion transitions, grid gaps | 0 | 1.3e-15 | not compared |
| order 2, grid gaps, H=2 | 5.7e-14 | 1.1e-15 | not compared |
| stayer, occasion measurement, initial covariate, grid gaps | 0 | 1.6e-15 | not compared |
| K=3, transition covariates read at the destination, factor levels c/a/b | 2.8e-14 | 2.7e-15 | not compared |

A skipped grid position consumes a transition with no emission, as documented. Class posteriors differ by at most 1.3e-15, and prevalence for H = 1 by at most 2.2e-16.

**Reference coding.** Factor levels c/a/b give terms `fa` and `fb`, with the first level as reference. Staying is the reference for transitions and the last profile for initial logits; the brute-force check depends on both.

**Parameter counts.** `n_parameters` equals the identified coordinates in all 29 specifications of `transitions-counts-and-flags.R`, including all 14 covariance structures. Hand counts agree: 38 (K=3 covariates), 18, 19, 23, 51, 17, 57 and 23.

**General scores.** Against central differences, the maximum difference is at most 1.7e-7 across those 29 specifications.

**Homogeneous scores and SEs.** `transitions-engines.R`, seven specifications: continuous H=1 and H=2, equal variance, full covariance, mixed categorical, FIML, grid gaps.
- Score against differences of the likelihood: at most 1.4e-7.
- Package SE against an SE from a likelihood-only Hessian: at most 1.2e-5 relative.
- Natural-scale Jacobian against differences: at most 1.0e-10.

**General SEs.** `transitions-misc2.R`, against `optimHess` with analytic gradient: 9.5e-8 relative on first order with occasion transitions and covariates, and 9.5e-8 on order 2 with occasion transitions and covariates. This Hessian uses the package's own analytic gradient, whose accuracy the score check above establishes. Those fits have one group class; the probe heading says H = 2, but the calls do not set it.

**Intervals and p-values.** The 90% interval half-width over SE is 1.644854, equal to `qnorm(.95)`. Two-sided p-values match exactly. Natural-scale intervals stayed inside [0, 1] on the fit printed.

**Engine agreement.**
- T = 4, general engine via weights of 1: LL -2715.8937024 in both, 15 parameters, BIC 5506.9469, individual BIC 5527.7413, `nobs` 150 and `BIC()` 5506.9469 in both. Transition difference at most 1.6e-8.
- T = 2, `transitions = "occasion"`: LL -1350.3194154 in both, BIC equal, logit SE 0.34934 printed in both for the same row, and 0.34116 for the other row on the general fit. The delta-method SEs of the general logits, 0.049841 and 0.060142, match the homogeneous natural SEs 0.04984 and 0.06014.
- Covariance structures EII, VEI, EVE, VEV: LL differences at most 1.4e-6 (EVE; others about 1e-12), parameter counts equal.

**Optimum.** `transitions-misc.R`: the largest score at the reported optimum is at most 5.8e-5 across 10 fits (EVV with weights and occasion measurement is the largest; diagonal occasion measurement 2.9e-5). The set covers EM-only structures with occasion measurement, FIML and weights.

**Sample size.** BIC uses persons, `bic_individual` uses rows. Weights of 2 against duplicated persons: LL -2715.893702 and -5431.787405, a factor of 2 to the printed digits, with BIC using log(150) and log(300) as the convention says.

**Bootstrap units.** `.multilpa_resample_groups()` keeps each person's rows whole and gives a person drawn twice new ids (draws 3, 3, 7 gave 3 ids of 5 rows each, with matching values). Wald and bootstrap SEs agree at 40 resamples (0.0317 and 0.0322; 0.0384 and 0.0372; general logits 0.168 and 0.180).

**Relabelling.** Permuting labels leaves probabilities permuted, not changed: initial with covariate 0, transition with covariate 1.1e-16, second order 1.1e-16.

**Tna export.** `get_tna()` weights equal counts normalised by hand (difference 0) and initial probabilities match (1.1e-16). The result is close to the model-implied marginal computed without data (0.7018 against 0.7021). `get_group_tna()` carries the class matrices (1.1e-16) and initials (0). The help describes this rule correctly.

## 4. Test gaps

Method: one function overridden in memory, then related test files run. The full suite was not run, as the brief asks, so "survived" means survived the files named.

| Mutant | Files | Result |
|---|---|---|
| Natural Jacobian: last transition rows x 1.1 | test-transition-inference, test-transitions, test-transitions-numerics | **Survived** |
| Homogeneous `estimated` column always TRUE | test-transitions, test-review-1015-transitions, test-release-transitions, test-tna, test-review-1026-lta-r2 | **Survived** |
| Prevalence columns reversed | test-transitions, test-release-transitions survived; caught after adding test-transitions-numerics, test-review-1015-transitions, test-lta-bootstrap | Caught, 2 failures, narrow |
| Chart Jacobian set to identity | test-transition-inference, test-review-1026-lta-r2, test-lta-general survived; caught by test-review-1026-lta-open, test-audit-inference-repairs, test-review-1101-bootstrap-tables | Caught, 2 failures, narrow |
| Expected counts x 1.01 | test-transitions, test-review-1015-transitions, test-release-transitions survived; caught once test-tna and test-transitions-numerics were added | Caught, 1 failure, narrow |
| Stayer initial logits not counted | test-lta-general, test-review-lta survived; caught by test-review-1026-lta-r2, test-review-1015-transitions, test-release-units, test-compare-lrt set | Caught, narrow |
| get_tna using class means | test-tna, test-tna-fits | Caught, 2 failures |
| Relabel without rebasing | test-lta-bootstrap, test-review-1026-lta | Caught, 11 failures |
| Second-order previous and from swapped | test-lta-general, test-review-lta | Caught, 1 failure |
| Second-order count K^2 replaced by K | test-lta-general, test-review-lta | Caught, 1 failure |

Central formulas whose tests would still pass if broken:
1. `.multilpa_transition_jacobian()`, the delta method behind every homogeneous transition SE. The test in test-transition-inference.R checks that two entries of a row share one SE, which scaling a whole row leaves true. It needs a comparison with numerical differences.
2. `.multilpa_transition_frame()`'s `estimated` column on a fit with an empty row. No value-level assertion in the five files named.

Missing tests for the findings:
- Near-empty origins: the tests I read exercise exact zero; I found none for tiny positive counts, though I did not read every test file.
- Unidentified T = 2 and T = 3 specifications.
- Bootstrap on unbalanced data with occasion-specific parameters, or with rare factor levels. The `boot_lta_data()` fixture in test-lta-bootstrap is balanced, which would explain why F3 went unseen.
- Skipped grid positions under second order or a stayer class. tests/equivalence/test-review-lta.R:95 has a path sum with skipped positions for occasion transitions with two group classes, and its other specifications run without skipped positions. This review's probe `transitions-brute.R` covers the two missing cases, and the suite could adopt it.
