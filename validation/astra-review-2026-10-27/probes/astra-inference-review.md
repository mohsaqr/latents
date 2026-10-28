# Astra subreview: inference, classification, criteria and bootstrap

Reviewed release baseline `24773e6` through `5596a15`, plus the current uncommitted diagnostics/reporting additions. No production source, package tests or repository documents were changed by this subreview. Live synthetic reproductions ran both through `pkgload::load_all()` and the coordinator's fresh `/private/tmp/astra-review/lib` install, under R 4.5.2 with one BLAS thread. Full-package tests belong to the coordinator's separate run. Findings below remain unfixed.

## Confirmed findings

### INF-A1 — P2: unidentified LTA bootstrap coefficients retain finite covariance rows

**Committed, introduced as an incomplete identification repair in `d5d4ae8`; still present in the working tree.** Location: `R/lta-bootstrap.R:65–71`.

The bootstrap correctly marks a transition into an entirely unobserved trailing occasion `estimated = FALSE` and replaces its standard error and interval with `NA`. It then calculates `covariance <- stats::cov(kept)` without applying the same mask. The returned covariance attribute therefore treats arbitrary, unidentified transition logits as estimated parameters, including cross-covariances with the identified parameters. This is internally inconsistent and can contaminate downstream contrast or delta-method calculations. Their variability reflects the bootstrap fits' arbitrary values in a direction the data do not identify.

Reproduction: 80 persons with four observed occasions and two Gaussian indicators; append an all-missing fifth occasion for 20 persons; fit `lta(..., transitions = "occasion", missing = "fiml", n_starts = 2, seed = 1)`, then `parameter_inference(..., method = "bootstrap", iter = 4, n_starts = 1, seed = 45)`. All four draws were usable. Both fifth-occasion logits report `estimated = FALSE`, `standard_error = NA`, but their covariance diagonals are **0.04989773 and 0.06113789**, with finite cross-covariances. The number of draws is intentionally small: this tests masking, not bootstrap calibration.

Evidence: `/tmp/astra-lta-bootstrap-probe.R`, its `.log`, `-installed.R` and `-installed.log`; saved fixture `/tmp/astra-lta-bootstrap-probe.rds`. The installed probe exactly reproduced the issue.

Recommendation: mask corresponding rows **and columns** of the returned covariance with `NA`; leave raw replicate estimates available as raw evidence if desired, with identification metadata. Add a regression asserting consistency between `estimated`, standard errors and the covariance attribute for trailing and interior missing occasions. Do not assert zero variance for an unidentified coefficient.

### INF-A2 — P2: new cross-level diagnostics put zero-weight groups back into class counts

**Uncommitted addition.** Location: `R/diagnostics-checks.R:180–197`, especially 185, 188 and 196.

The posterior-based diagnostic helper drops only rows with missing posterior entries. A zero sampling weight does not make a posterior missing. The helper then rescales weights to `nrow(posteriors)`, calculates modal counts over every row, and multiplies estimated shares by that same full row count. This reverses the committed rule that zero-weight units are absent from sample sizes and classification.

Real fit: 60 groups of six observations, restricted cross-level model, 2 profiles and 2 group classes; ten groups have weight zero. `nobs(fit) = 50`. The fit's `group_classes` counts are **30 and 20**. New diagnostics produce expected counts **36 and 24**, modal counts **30 and 30**, and print the smallest class as **40.0% (n = 24)**. The expected count is inflated by 60/50 = 1.2. The source model's existing `n_assigned` table also includes the zero-weight groups, but the newly created diagnostics expected-n inflation is introduced by this helper.

Evidence: `/tmp/astra-inference-probes.R`, `.log`, `-installed.R`, `-installed.log`. Installed and loaded-source runs agree.

Recommendation: subset posterior rows and aligned weights to positive-weight contributing units before any modal count, average posterior or effective-count calculation; retain the package's normalized positive-sample count convention. Add a zero-weight cross-level fixture; the current weighted cross-level test uses positive weights only and repeats the helper's all-row normalization, so it cannot expose this defect.

### INF-A3 — P2: new general-transition diagnostics classify all-missing rows despite `drop_empty = TRUE`

**Uncommitted addition.** Locations: `R/diagnostics-checks.R:180–197`, `R/diagnostics-checks.R:315–325`; dispatch at `R/report.R:123–146`.

`multilpa_lta` (general transitions: occasion-varying, covariate or second order) is accepted by the new diagnostics API but is not among `.multilpa_any_fit()`'s three classes. It therefore uses the posterior-only fallback. That fallback tests whether posterior values are missing, rather than whether the original row is classified under `drop_empty`. A row with no observed indicators has a finite smoothed posterior and is therefore incorrectly retained in entropy, class counts and average posteriors.

On the real fit used for INF-A1, the model stores **340 rows, 320 informative**, with `drop_empty = TRUE`. Correct expected counts over informative rows are **151.8988 and 168.1012**; diagnostics totals 340 instead of 320 and reports the smallest count as approximately 162. The direct informative-row relative entropy is **0.8888741**; diagnostics reports **0.8543787**. Appending unclassified occasions therefore changes the new classification assessment.

Evidence: `/tmp/astra-lta-diagnostics-probe.R` and `.log`, using the fresh installed library and saved actual fit. The direct entropy calculation uses only `-sum(p log p)/(n log K)`, not the diagnostic implementation.

Recommendation: give general LTA an explicit eligibility/weight adapter using its stored observed-row mask, `drop_empty` policy, sampling weights and `group_index`; use it consistently for counts, modal assignments, AvePP and entropy. The existing "every family" diagnostics fixture includes only homogeneous `lta()` and misses this separate class. Test both homogeneous and general engines, with empty rows, zero weights and both `drop_empty` policies. This and INF-A2 share an eligibility-design problem, but have separate user-visible triggers and may be consolidated into one finding in the final report.

### INF-A4 — P3: newly accepted families do not fulfil the documented diagnostics table/coercion contract

**Uncommitted addition.** Location: `R/report.R:123–146`, plus the accessor's table selection.

The new help promises transition fits' entropy/classification/average-posterior tables and says `as.data.frame(diagnostics)` returns entropy. For the actual general-LTA fit above, only `replication` and `checks` are stored; `get_results(dg, "entropy")` raises `latents_bad_argument`, and `as.data.frame(dg)` silently returns the **replication table** (`n_starts`, `n_screened`, ...). The same omission follows for newly accepted additive, cross-level, regression and growth objects. Those families did not previously support diagnostics, so this is a newly advertised contract mismatch rather than a regression of their former successful coercion.

Evidence: last two outputs of `/tmp/astra-lta-diagnostics-probe.log`.

Recommendation: either provide the documented primary entropy table for all supported families, with explicit n/a rows, or precisely document and implement a consistent alternative primary table for these families. Avoid data-frame coercion changing meaning merely because only a subset of result tables was populated.

## Review of changes since release

| Area | Review assessment and release impact |
|---|---|
| Weighted sample convention | `R/weights.R` now normalizes over raw-positive units, preserves tiny positive weights whose ratios underflow, and tests within-unit constancy relatively. This intentionally changes likelihoods/criteria when old analyses contained zero-weight units. The same eligibility convention must reach new reporting paths; INF-A2/A3 show it is not yet complete. |
| Complete-data entropy and criteria | `R/classification-entropy.R` carries group entropy plus conditional profile entropy from the E-step; `R/diagnostics.R` uses the same complete-data entropy under both sample-size conventions for nontransition fits. This fixes earlier group-criterion behaviour and intentionally changes ICL/AWE/CLC/model rankings. Single-group models no longer have ICL_groups equal BIC merely because group entropy is zero. Transition-level entropy is an explicitly separate convention. Old multigroup fits lacking stored conditional entropy return NA rather than reconstructing an unsupported result; users need to refit. |
| Classification summaries | Weighted entropy, AvePP, classification error matrices and comparison/enumeration summaries are now consistent for ordinary profile/class models. Modal counts remain actual row counts by design; population shares and posterior averages are weighted. Avoid describing modal-versus-estimated gaps as solely classification error when sampling weights differ. |
| Empty rows and three-step | `.multilpa_level_assignments()` subsets posterior/assignment matrices and keeps the original row/group positions for outcome/design alignment. The error matrix, BCH and step-three data now use classified units. Weighted three-step estimation remains explicitly refused even though weighted error-matrix reporting exists. |
| Centring/scaling inference | `R/inference-service.R` stacks estimated grand transform equations into covariance, while person transforms default to robust covariance and refuse observed/OPG. Covariate inference and pooling forward the same rule. These are intentional SE and API changes, not output-neutral repairs. Existing calibration evidence and regression files are substantial, but this subreview did not repeat Monte Carlo calibration. |
| Numerical conditioning | Inference differences use coordinate scales and caps; additive covariance is inverted in correlation coordinates; LTA covariate inference uses the fitting chart and maps scores/information/covariance back. These changes target translated/scaled predictors and small-unit outcomes. Reviewed source plumbing; did not rerun all family-specific reference-package comparisons. |
| Category bounds/labels | Boundary handling holds only the appropriate response coordinates; complements at the upper bound have unavailable SEs. Shared escaping helpers disambiguate labels containing `:`, `%`, commas and brackets across coefficients, covariance names and table lookups. Consumers parsing old raw paired labels must adapt to escaped terms. |
| ML three-step search | Equal-variance distal ML adds BFGS from quantile placements, accepts only finite converged stationary candidates with positive-definite conditioned information, and preserves the previous winner under ties. Finite searches are explicitly not global-maximization guarantees. The Box-Muller cache limitation is documented rather than claiming complete RNG restoration. |
| General LTA identification | Flat trailing logits and empty measurement blocks are excluded from Wald inversion and counted dimensions. Interior empty occasions retain conditional inference for the exit move, with a warning that only the product is identified. This is an intentional conditional convention, not an unconditional SE guarantee. Bootstrap metadata repair is incomplete (INF-A1). |
| BLRT | P-values now use valid replicates and report invalid counts; fewer than half valid gives NA. The documentation acknowledges conditioning on successful refits. A 50% validity threshold is an operational rule, not evidence of nominal-size calibration when failure is outcome-dependent; the historical bootstrap-calibration claims should not be generalized to a large invalid fraction. |
| Comparison/selection | Data signatures are compared exactly except allowed normalization roundoff and provable forward centring; nominal versus ordinal encodings are compared on observed labels. Historical model labels equal to newly introduced filter argument names are restored. Replication-unavailable candidates no longer display as replicated. |
| Diagnostics redesign | Separates reported quantities from factual failures, keeps screened starts out of unfinished counts, and explicitly reports unsupported tests. This improves transparency. New family dispatch needs eligibility and table-contract fixes above before claiming full-family support. |

## Evidence quality and boundaries

- Existing report `validation/review-2026-10-26/round-2/inference.md` documents numeric-score checks, robust/OPG comparisons, bound invariance and Monte Carlo calibration. These are **historical project evidence**, not independently rerun in this subreview.
- Existing tests cover category-order invariance, partial paired items, weighted selection, zero-weight likelihoods, centring calibration, LTA rank/flat directions and error contracts. The newly discovered gaps are combinations of previously repaired conventions and new adapters.
- `equivalence/harness/known-open.csv` is empty. This means the registered harness currently has no listed exception; it cannot establish absence of the defects above, which need their own coverage.
- Focused reproductions used synthetic data only. No external commercial statistical program was run. Bootstrap probe assesses output consistency, not coverage.
- Release recommendation for this area: fix INF-A1 plus the new diagnostics eligibility defects, resolve/document INF-A4, add targeted regression cases, then rerun relevant tests/harness and the final source-archive check. Keep the explicit behavioural-change notes for centring SEs, complete-data criteria, zero-weight sample counts, partial joint-item FIML, bootstrap validity policy and saved-fit compatibility.

## Follow-up of coordinator's regression-grid lead

### INF-A5 — P2: regression enumeration quality view still discards sampling weights

**Committed mismatch after round-3 weighting changes, still present in the working tree.** Location: `R/grid-views.R:264–267`, affecting lines 285–288 and the public enumeration quality plot via `R/plot-gg.R:1682`.

`level_weights()` explicitly returns `NULL` for every regression grid. The fit's own tables were changed in round 3 to weighted classification, but the quality-view adapter remains unweighted and reads full posterior matrices (also retaining zero-weight units). Therefore the displayed group entropy and AvePP can disagree with the corresponding candidate's own tables. The primary individual entropy is copied from the grid and remains weighted, so different panels silently use different conventions.

Two real-fit reproductions on the fresh installed library:

1. Weighted single-level `enumerate_regressions(y ~ x, ..., n_classes = 2, n_starts = 2, seed = 3)`: quality-view minimum AvePP **0.6876206**, candidate classification-table minimum **0.6902239**.
2. Weighted two-level enumeration (60 groups of six, 2 individual classes and 2 group classes): quality-view group entropy **0.6751874**, direct weighted entropy **0.6669745**, candidate fit table **0.6669745**.

Evidence: `/tmp/astra-regression-grid-probe.R` and `.log`; `/tmp/astra-regression-group-grid-probe.R` and `.log`. These probe the numeric table consumed by the plot, without modifying or rendering any plot. Fix by obtaining the appropriate regression/growth row or group weights and eligibility mask, or consuming the candidate's own classification summaries directly. Add an end-to-end assertion that every plotted quality measure equals the corresponding weighted candidate table, covering positive unequal weights and added zero-weight units. The coordinator identified the source lead; this subreview independently reproduced its two numeric consequences.
