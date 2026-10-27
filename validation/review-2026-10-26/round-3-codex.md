# Adversarial review and fixes — round 3, Codex

2026-10-26. Starting tree identified by the owner as d5d4ae8, branch `features`.
All processes ran from the repository root. No git command was run. The released
0.9.21 sources at a10a919 were read directly from the stored object files using
Python's standard library, without changing repository metadata.

## Result and scope

**Eleven confirmed defects fixed**, including three families affected by one
missed sampling-weight convention. The requested regression file is
`tests/testthat/test-review-1026-codex.R`: **15 test blocks, 73 passing assertions,
zero failures, errors, warnings or skips**. Every block fails with the relevant
pre-fix definitions restored together: 35 failed assertions, seven errors, one
pre-existing boundary warning, eleven passing guards. Thus all fifteen blocks
were demonstrated to fail before the fixes; the passing guards are not counted
as defect discoveries.

I read both review rounds and the harness instructions before investigating.
This review concentrated on interactions and propagation, rather than repeating
those rounds' mathematical derivations. No full test suite or full harness run
was started. Neither `equivalence/harness/outputs.R` nor `registry.R` was edited;
no harness suite was edited. No plot implementation or palette changed.

Evidence directories:

- `/private/tmp/latents-codex-round3/`: root reproductions, release extraction,
  original source backups, combined pre-fix run and final regression results.
- `/private/tmp/codex-kernel-r3/`: counts, extreme weights, fixed/FIML/SQUAREM,
  centring and bounded-prior checks.
- `/private/tmp/codex-criteria/`: weighted summaries, joint-entropy enumeration,
  paired covariate inference and its independent likelihood/Hessian.
- `/private/tmp/latents-codex-growth/`: Gaussian floors, singleton identification,
  Brownian units, predictor maps and original-source checks.

The permanent minimal reproductions are the named tests below. Scratch scripts
preserve additional measurements and independent checks.

## Confirmed defects, fixes and regressions

### C3-1 — Ordinary category labels interpreted as partial joint observations

**Locations:** `R/categorical.R:33`, `:137`, `:249`; `R/predict.R:117`.
Supporting provenance preservation: `R/accessors.R:215`,
`R/local-dependence.R:184`, `R/simulation-service.R:105`.

A plain factor with categories `a:x`, `b:x`, `a::` was mistaken for a joined
item with a missing second response. With category counts 40, 30, 30 and a
balanced independent binary item, the one-class log likelihood was
**−130.401148261484**, whereas the independent multinomial maximum is
`40 log(.4) + 60 log(.3) + 100 log(.5)` = **−178.204715590517**: an excess of
**47.803567329033**. The reported item probabilities were .7/.3, with the third
category lost. Prediction also silently accepted an unseen `a::` response on an
ordinary item whose fitted categories were `a:x` and `b:x`.

The half-observed interpretation now requires explicit pair provenance.
Typed category values, reconstructed data and bootstrap row subsets retain
that provenance; older saved fits recover it from their recorded
`local_dependence`. Real partial pairs still integrate their missing item out.
A bootstrap resample retains the full declared joint grid even when a cell is
unobserved. The existing rare-cell bootstrap test's expectation of a dropped
replicate was updated accordingly; the affected file passes all 61 assertions.

**Reproduction:** `labels-before.R`, `labels-before-corrected.log` in the root
evidence directory. The initial exploratory hand formula indexed sorted factor
levels incorrectly; that was corrected before accepting the defect. The numbers
above are from the corrected run.

**Tests:** “ordinary category labels cannot become partial joint observations”
and “prediction accepts partial markers only for declared joint items”
(lines 3 and 18). Guards cover actual pair marginal densities, rebuilding,
resampling and inference on saved fits without the new value attribute.

### C3-2 — New comparison controls captured released model labels

**Location:** `R/compare-models.R:105`.

Release 0.9.21 accepted arbitrary model labels in `compare_models(...)`.
`compare_models(n_profiles = regression_fit, other = regression_fit)` returned a
comparison in the released implementation, but raised `latents_bad_argument` in
the reviewed tree: the added `n_profiles` filter consumed the fit. The same
collision applies to `n_group_classes` and `cutoff`.

Fit-valued controls are now restored as labelled models, in their original call
order. Numeric grid filters and numeric cutoff retain their current behavior.
The four-model regression example returns the four original labels and equal
BIC weights .25; numeric grid selection still selects the requested counts.

**Reproduction:** `api-collision.R`, `api-label-before.log`; independent reference:
the actual extracted a10a919 `compare_models` implementation.

**Test:** “compare_models retains released model labels that now name controls”
(line 53); invalid numeric filtering without a grid is checked by condition class.

### C3-3 — Observed-indicator counts still included absent rows

**Location:** `R/multilpa-starts.R:224`.

Appending zero-weight groups to a mixed continuous/categorical FIML fit changed
`n_observed_by_indicator` from **46/45/46** to **57/57/58**, despite unchanged
likelihood and estimates. Direct nonmissing counts on positive-weight rows give
46/45/46 under either `drop_empty` convention.

The stored count now masks out zero-weight rows.

**Reproduction:** `probes.R` in the kernel evidence directory.
**Test:** “zero-weight rows do not increase observed indicator counts under FIML”
(line 71), including both `drop_empty` values.

### C3-4 — A positive weight could underflow into an absent unit

**Location:** `R/weights.R:57`.

Raw weights `(1e-300, 1e300, 1e300)` normalized to `(0,1,1)`. Although all three
units were sampled, `nobs` became two and the other two weights were normalized
to the wrong total. At a supplied standard-normal start the reported LL/BIC were
**−2.8378770664093453 / 7.0620484939385815**. The independent weighted normal
formula, with three raw-positive units and dominant weights 1.5 each, gives
**−4.2568155996140185 / 10.7108557765642569**.

Normalization now uses the raw-positive count and preserves an underflowed
positive ratio at the smallest representable positive double. Its negligible
mass no longer changes the count or scaling of all other units. Finite weights
whose sum overflows, `(1e308,1e308,1e308)`, still normalize to ones.

**Reproduction:** `weight-extremes.R` in the kernel evidence directory.
**Test:** “positive sampling weights stay present when their ratio underflows”
(line 99), using an independently evaluated fixed normal likelihood.

### C3-5 — Weighted classification convention missed three families

**Locations:** `R/mixture-regression-methods.R:101`, `:209`, `:248`, `:349`;
`R/growth-mixture-methods.R:315`, `:375`, `:391`;
`R/additive-methods.R:63`, `:130`, `:184`, `:198`.

The profile family's revised convention had not reached regression, growth or
additive reporting. Their likelihoods/shares were weighted, while entropy,
AvePP and classification tables remained unweighted. Independent formulas were
`sum_i w_i [-sum_k p_ik log p_ik]` and, within a modal class,
`sum_i w_i p_ik / sum_i w_i`.

Measured old → hand/fixed values:

| Family | Quantity | Old | Hand/fixed |
|---|---|---:|---:|
| Regression | ICL | 819.5279784 | 819.5910811 |
| Regression | Relative entropy | .3352921 | .3350392 |
| Regression | Class-two AvePP | .7659512 | .7710580 |
| Growth | ICL | 811.2367281 | 808.8380133 |
| Growth | Relative entropy | .8517356 | .8783558 |
| Growth | Class-one AvePP | .9677288 | .9811282 |
| Additive | Relative entropy | .5700079 | .5528055 |
| Additive | Class-one AvePP | .8625565 | .8724203 |

All three now apply sampling weights consistently, using shared entropy helpers
and weighted class means. Modal assignment counts remain counts.

**Reproduction:** `probe.R`, `pre-fix.out`, `post-fix.out` in the criteria directory.
**Tests:** “regression classification estimates follow the sampling weights”,
“growth classification estimates follow the sampling weights”, “additive
classification estimates follow the sampling weights” (lines 151, 168, 182).

### C3-6 — Two-level regression ICL omitted joint classification entropy

**Location:** `R/mixture-regression-methods.R:217`.

ICL charged only marginal row-class entropy, whereas the revised shared
convention charges the group class and all row classes jointly. Enumerating
all `2^3` row-class paths and both group classes independently, using fitted
normal densities and class priors, gives joint entropy **47.36013**. Old ICL was
**491.08191**; the hand result and corrected ICL are **544.07627**.

New `.mixture_complete_entropy()` adds group-class entropy to expected row
entropy conditional on the group class, with the same sampling weights as the
likelihood. This is an identity at any parameter point. The selected start in
this fixture reports nonconvergence; no claim that it reached the maximum is
needed for this exact entropy check.

**Reproduction:** `probe.R`, `twolevel.rds` in the criteria directory.
**Test:** “two-level regression ICL charges the brute-force joint entropy”
(line 199). Its independent reference uses neither stored posteriors nor a
package entropy helper.

### C3-7 — Covariate inference lost half-observed paired responses

**Locations:** `R/covariates.R:648`, `R/covariate-inference.R:425`.

On a 300-row membership-covariate fit, all **34** half-observed pair rows were
rebuilt as fully missing. `vcov(fit)` raised `latents_bad_inference_data`;
`vcov(fit, data = original_data)` also failed because the original items were
never joined on this inference path.

Stored-data reconstruction now restores the observed-item markers and pair
metadata. Explicit original data goes through the shared pair resolver. All
34 observations are recovered; supplied/default covariance matrices agree
exactly. A hand likelihood using normal densities, logistic membership and
sums of compatible joint probabilities agrees with LL within **1.136868e-13**.
The inverse numDeriv Hessian of that hand likelihood gives SEs within
**8.690182e-7** relative. Changed observed responses remain refused by class.
Older saved covariate fits without the new typed-value attribute also work.

**Reproduction:** root `partial-cov.R`; criteria `partial-cov-oracle.R` and
`partial-cov-oracle.out`.
**Test:** “paired FIML covariate inference retains each half-observed item”
(line 245).

### C3-8 — Absent outcomes set Gaussian floors and growth coordinates

**Locations:** `R/mixture-regression.R:378`, `:494`; `R/growth-mixture.R:354`.

Gaussian floors were computed before weights were attached; growth transforms
also used all outcome/design rows. Forty persons with four observations fit at
LL **−275.7766** (regression) / **−172.5993** (random-intercept growth). Appending
four zero-weight persons with outcomes ±1e5 made both raise
`latents_no_valid_start`. Independent OLS gives residual ML variance
**1.8391998**; deleting zero-weight rows is an exact likelihood identity.

Floors and growth coordinates now use carrying rows. Both padded fits recover
the original likelihood and floor **2.038943e-6**; residual variances are
**1.8391998 / .2139618**. A constant outcome among carrying rows is refused as
`latents_bad_data` even when absent rows vary.

**Reproduction:** growth `zero-floor.R`, `zero-floor-after.R`, `zero-design.R`.
**Tests:** “zero-weight outcomes cannot change Gaussian residual floors” and
“a constant Gaussian outcome on carrying rows is refused by class”
(lines 278, 305).

### C3-9 — Relative Gaussian floor became absolute in tiny outcome units

**Locations:** `R/mixture-regression.R:514`, `:735`.

`min_variance * max(var(y), double.eps)` compared a unitless constant with an
outcome variance. Multiplying the preceding data's outcomes by 1e-12 made
regression raise `latents_no_valid_start` and growth raise `latents_em_decrease`
(with a reported decrease of 31.4). The independent density Jacobian requires
LL **4145.1868 / 4248.3641**, obtained by subtracting `160 log(1e-12)` from the
original LL.

The positive numerical guard is now `double.xmin`. Weighted/unweighted
regression and growth satisfy the hand likelihood and variance transformations.

**Reproduction:** growth `tiny-outcome.R`.
**Test:** “relative regression variance floors remain relative below machine
epsilon” (line 314).

### C3-10 — Singleton growth persons counted an unidentified variance

**Location:** `R/mixture-regression.R:383`.

With 30 one-observation persons, one class and `random = ~1`, the model returned
LL **−41.76077**, `converged = TRUE`, no random boundary, random-intercept
variance **.1143557**, residual variance **.8332413**, and four parameters.
OLS gives the identical likelihood with three identified parameters. The hand
marginal covariance is `G00 + sigma²`; adding delta to one variance and
subtracting it from the other changes no likelihood.

A random-intercept growth fit with every carrying person a singleton now raises
`latents_not_identified`. Repeated zero-weight persons cannot bypass it.
Correlated residual processes do not identify the intercept/nugget split either.

**Reproduction:** growth `singleton.R`.
**Test:** “singleton persons cannot identify a random intercept and residual”
(line 336), checking conditions by class and the independent variance identity.

### C3-11 — Brownian diffusion bounds and starts ignored time units

**Location:** `R/growth-mixture-correlated.R:214`, `:252`.

Fifty persons with five observations, one class, a random intercept and Brownian
residuals: multiplying time by 1e8 changed LL from **−296.3683934339222** to
**−756.4263322851039**, losing **460.057938851182**. Diffusion SD mapped to the
original units changed from **.7769658603014684** to **13.82497668613537**;
residual variance changed from **.05374561871848357** to
**1.911299803721864e-6**. The dense normal likelihood independently reproduces
both values. Since `W_ij = diffusion_sd² min(t_i,t_j)`, diffusion SD must map as
`sd / sqrt(c)` when time maps to `c*time`.

The Brownian floor/start now refer to process variance at the largest carrying
time; process time scales and AR gap calculations exclude absent persons.
Both units now give LL **−296.36839343**, mapped SD **.77696586**, residual
variance **.05374562**, and no residual boundary. Absent persons with times up
to 4e12 leave LL and diffusion SD unchanged.

**Reproduction:** growth `brownian.R`, `measure-old.R`, `measure-old.log`;
`zero-process.R` checks absent-time interactions.
**Test:** “Brownian growth residuals preserve their model when time units change”
(line 356), with independent dense likelihood and conditional prediction checks.

## Attempts that held

| Probe and independent check | Produced result |
|---|---|
| Released surface: all exported verb formals and S3 registrations against a10a919 | 28 exported verbs and 160 S3 registrations examined; none removed, no old formal name removed. The dots/control collision above was the confirmed break. |
| Released result fields/table names/column names on Gaussian, categorical, covariate, additive, full cross-level, staged, transition, regression and growth fits | No removed field, table or column on these nine representative fits. Default/control additions were recorded in `api.log`; this is not a claim about every possible optional result. |
| Grand/person centring and scaling × FIML × zero-weight padding × `drop_empty = FALSE` | LL, BIC, means and variances exactly unchanged. |
| Fixed means × weighted two-level FIML × SQUAREM | Hand normal-marginal/group-logsumexp LL error ≤1.136868e-13; tight plain EM versus SQUAREM LL difference 8.299139e-12; fixed-mean error zero; minimum history gain zero; both converged. |
| Bounded-prior VEI/VEV M-steps, floors .01, 1, 10, against independent shape optimization | Objective discrepancy ≤9.24e-11. |
| All-missing continuous FIML column | Refused with `latents_bad_data`. |
| One-class FIML items with 25 and 30 levels | LL versus multinomial formula −3.092282e-11; parameter count 53, as required. |
| Sparse 22×22 joined item | 484 categories retained, probabilities sum to one; direct LL discrepancy zero; parameter count 483. |
| Ordinary two-level model, every group a singleton | Refused with `latents_unidentified`. |
| 100,000 rows, supplied Gaussian parameters, tiny class expected membership one | Independent logsumexp LL discrepancy zero; largest posterior discrepancy 2.514629e-22; thin-profile warning raised. |
| Identical Gaussian class starts with SQUAREM | Posterior discrepancy from .5 exactly zero; LL versus the one-normal MLE 8.731149e-11; converged. This checks density arithmetic, not identification of duplicated components. |
| Corrected regression criteria through `compare_models()` and weighted `enumerate_regressions()` | ICL/relative-entropy discrepancies zero. |
| Brownian time multiplied by 1e-8 after the fix | LL difference −5.684342e-14; mapped diffusion-SD difference −1.593170e-13. |
| Noise/prior with unsupported weights/missing indicators; staged incompatible measurement settings | Existing refusals held; no additional defect confirmed. |

Scripts: root `api.R`, `api-methods.R`, `api-fields.R`, `api-families.R`,
`extremes.R`, `large-fixed.R`; kernel `interaction-probes.R`, `fixed-squarem.R`;
criteria `probe-grid.R`; growth `held-probes.R`.

## Verification

All test-file runs used `NOT_CRAN=true OMP_NUM_THREADS=1
VECLIB_MAXIMUM_THREADS=1`, `pkgload::load_all(".", quiet=TRUE)`,
`library(testthat)` and `test_file(...)`. The combined pre-fix proof loads the
current package, restores backed-up definitions in both its namespace and
attached exports, and runs the same final regression file. Its script is
`/private/tmp/latents-codex-round3/all-before.R`; its table asserts that every
block has a failed assertion or error. `final-tests.rds` retains the fixed result.

There were **23 distinct targeted test files, 1,226 passing assertions** across
the successful runs in this session. This is a deduplicated per-file inventory,
not a full-suite run or a single simultaneous run; the final consolidated
regression file and directly affected compatibility checks were rerun after
consolidation. The one stale rare-cell bootstrap assertion was corrected as
explained under C3-1. Initial commands mentioning nonexistent test paths were
corrected; they are command-path errors, not package test failures.
| File under tests/testthat/ | Passing assertions |
|---|---:|
| test-review-1026-codex.R | 73 |
| test-weights.R | 76 |
| test-missing-covariance.R | 50 |
| test-review-1026-zero-weight-likelihood.R | 19 |
| test-categorical.R | 75 |
| test-review-1026-categorical-r2.R | 28 |
| test-local-dependence.R | 61 |
| test-features-joint-items.R | 61 |
| test-predict.R | 19 |
| test-review-1026-selection-r2.R | 41 |
| test-bootstrap-inference.R | 66 |
| test-mixture-regression.R | 78 |
| test-growth-mixture.R | 122 |
| test-multilevel-growth.R | 33 |
| test-additive-methods.R | 82 |
| test-review-1026-mixreg-r2.R | 84 |
| test-review-1026-growth-r2.R | 71 |
| test-review-1026-additive-r2.R | 31 |
| test-covariate-inference.R | 40 |
| test-covariate-missing.R | 24 |
| test-covariate-measurement.R | 46 |
| test-growth-correlated.R | 25 |
| test-growth-predict.R | 21 |

Logs are in the four evidence directories listed above. No new harness suite was
needed or changed, so no harness-suite invocation was required. The owner's
background full run uses an installed snapshot and does not certify these
new source changes.

## Unresolved or deliberately unclaimed

- AR residuals with multilevel growth K=2/H=2, random intercept/slope, and
  `time -> 1000*time + 2000` agree on returned LL within 8.526513e-13,
  coefficients within 1.161295e-7 after mapping, and posteriors within
  2.767307e-8. Both fits report nonconvergence. This verifies the coordinate
  map at the returned points; an independent maximal-fit comparison remains
  unresolved. No defect is asserted from this probe.
- `.growth_start()` still reads absent persons for some starting covariance
  calculations. The tested one-class fits recover the same optimum after the
  confirmed fixes. A separate multi-start counterexample was not independently
  confirmed, so that code was not changed.
- Existing documented low-priority observations (SQUAREM steplength units,
  conservative BVR p-values, near-duplicate additive diagnostics) were not
  reclassified or re-reported as new defects.
- No fresh full suite, complete external harness, live parallel-worker audit,
  exhaustive optional-result compatibility matrix, or new Monte Carlo coverage
  study was run. No claim of exhaustive correctness follows from these checks.
