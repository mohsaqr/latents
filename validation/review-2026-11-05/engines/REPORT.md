# Estimation engines and data transformations review

Eight outstanding issues were reproduced on the current 1.1.16 working tree: five can change reported estimates, parameter counts or selection, two accept unidentified models, and one crashes bootstrap inference on unbalanced panels. The Gaussian likelihood and small-unit checks passed. These findings were already present in 1.1.9; they should not be described as regressions introduced by the changes after that release.

## Scope and evidence

Reviewed the relevant history from `a10a919b` (0.9.21) and the changes from `8518dffe` (1.1.9) through HEAD `a5551081` plus the existing working-tree changes. The detailed source review focused on centering/scaling, regression design construction and assignment tables, growth identifiability, transition parameter identification and bootstrap alignment.

All experiments use newly generated synthetic data with new seeds, R 4.5.2, `VECLIB_MAXIMUM_THREADS=1`, `OMP_NUM_THREADS=1`, `NOT_CRAN=true`, and one fitting core. No package code, package tests or generated documentation was edited. Root runs the package-wide suite separately.

| Finding | Impact | Current reproduction |
|---|---|---|
| E1 Partial singleton centering | A false zero profile can win BIC | Confirmed |
| E2 Weighted grand transformations | Integer weights disagree with duplicated units, including fitted raw means | Confirmed |
| E3 Unused factor levels | Extra parameters, inflated BIC, unavailable standard errors | Confirmed |
| E4 ID named `class` | Assignment/group tables lose the original IDs | Confirmed |
| E5 Two-wave random slope growth | Arbitrary covariance decomposition reported as a converged fit | Confirmed |
| E6 Ragged transition bootstrap | Internal alignment error aborts the bootstrap | Confirmed |
| E7 Two-wave stayer model | Arbitrary stayer share reported without identification warning | Confirmed |
| E8 Near-empty transition origins | Effectively unsupported transition rows reported as estimated and counted | Confirmed |

The main six experiments and their assertions are in [probes.R](probes.R), with output in [probes.log](probes.log). Transition-specific checks are [transition-ridge.R](transition-ridge.R) and [near-empty.R](near-empty.R). Retained RDS files contain the fits and tables.

## Findings

### E1 Person centering creates an artificial profile from singleton units

**Location:** `R/centering.R:203–210` subtracts every unit's own mean; `R/centering.R:315–322` refuses the result only when the indicator is flat in every carrying unit.

With 75 persons, including 10 with one observation, person centering produces 10 exact all-zero vectors. A three-profile fit puts a profile on those vectors: its means are zero to numerical precision and its three variances are at the floor, `9.27e-7`. The only warning is the generic variance-boundary warning.

| Sample | BIC for 2 profiles | BIC for 3 profiles | Preferred count |
|---|---:|---:|---:|
| Including singleton persons | 3228.997 | 2919.770 | 3 |
| Excluding singleton persons | 3155.763 | 3157.454 | 2 |

The BIC values from different samples should not be compared with one another; each row demonstrates the within-sample selection. The singleton deviations carry no within-person variation, yet the likelihood treats them as observations that support an almost point-mass profile.

**Repair direction:** flag partially singleton indicators explicitly. Treat their centered values as missing, or refuse the affected call, according to the intended missing-data contract. This changes affected fits and needs a regression test with a mixture of singleton and repeated persons.

### E2 Grand centering and scaling ignore positive sampling weights

**Location:** `R/centering.R:116–140` excludes zero-weight rows but computes unweighted standard deviations; `R/centering.R:192` uses unweighted `colMeans`.

The independent experiment uses 320 observations with integer weights 1 and 3. For an EII model without transformations, weighted and duplicated-row fits agree in raw profile means to `2.03e-8`. With `centering="grand", scaling="grand"`, their raw means differ by `0.291447`; their reported transformed means differ by `0.642443`.

The first indicator's centering offset is `1.477410` for the weighted fit and `0.773712` for duplicated rows. Its scale is `1.750772` versus `1.599845`. For a spherical covariance model, changing the relative indicator scales changes the model itself, so this is more than a reporting-origin difference.

**Repair direction:** define weighted offsets and scales consistently with the package's weight-normalization convention, then test weighted fits against duplicated units under both scale-invariant and non-invariant covariance structures. If the intended convention is an unweighted transformation followed by a weighted fit, document the exception to the claim that integer weights act as repeated units.

### E3 Unused membership factor levels inflate parameter counts

**Location:** `R/mixture-regression.R:891` constructs the membership model frame without dropping unused levels; `R/mixture-regression.R:690` constructs the random design directly. The response regression frame already uses `drop.unused.levels=TRUE` at `R/mixture-regression.R:606–607`.

For 70 persons, a membership factor with observed levels `a`, `b` and unused level `unused` creates a zero `grpunused` design column. The same data after `droplevels()` give:

| Membership factor | Parameters | Log likelihood | BIC |
|---|---:|---:|---:|
| Unused level retained | 9 | -330.325 | 698.8865 |
| Unused level dropped | 8 | -330.325 | 694.6380 |

The BIC excess is `log(70)=4.2485`. The retained-level fit raises `latents_singular_information`, blaming unidentified parameters, empty classes, boundaries or separation; the unused design column is the actual cause in this example.

**Repair direction:** drop unused levels when constructing membership and random model frames, retaining the resulting coding for prediction. The membership case was executed here; the identical random-design issue is a source-level confirmation, supported by the previous review's independent random-design experiment.

### E4 An ID column named class is overwritten in assignment tables

**Location:** `R/mixture-regression-methods.R:283–285` writes the ID and then overwrites `table$class`; the grouped table repeats this collision at `R/mixture-regression-methods.R:310–311`.

A valid grouped regression with `id="class"` and IDs 1001 through 1070 fits without warnings. Its group table has one `class` column containing only the estimated labels 1 and 2; every original ID is lost. The stored fit retains the correct group structure, but a caller cannot join the returned table to the original IDs through its advertised ID column.

**Repair direction:** reserve names at input validation or retain IDs under an explicit non-colliding name. Test both per-row and per-group tables; the growth-method table code has related collisions and should be covered by the same contract.

### E5 Two-wave random intercept and slope growth has an exact covariance ridge

**Location:** `R/mixture-regression.R:315–324` checks only the all-singleton random-intercept case; it does not check the rank of the covariance design.

For 120 persons at times 0 and 1, the model returns `converged=TRUE`, `random_boundary=FALSE`, six parameters and no warning. Changing the random covariance and residual variance according to

`g00 -> g00-d, g01 -> g01+d, g11 -> g11-2d, residual_variance -> residual_variance+d`

leaves every modeled two-wave covariance unchanged. The independently evaluated likelihood remains `-367.7593` for `d=-0.05`, `0`, and `0.025`, while all random covariance matrices remain positive definite. The residual variance ranges from `0.411961` to `0.486961`.

These are distinct admissible parameter vectors with the same data distribution. Convergence and absence of a boundary do not establish identification.

**Repair direction:** check the covariance-design rank before fitting, and warn or refuse a deficient specification. Tightening a numerical Hessian threshold alone cannot reliably establish structural identification.

### E6 Occasion-specific transition bootstrap crashes when the longest person is absent

**Location:** `R/lta-bootstrap.R:184` and `R/lta-bootstrap.R:197` align signatures before checking the coefficient names at line 166. `R/simulation-service.R:151` extracts estimates outside the refit error handler.

The synthetic panel has 50 persons, 49 observed at three occasions and one observed at five. Fitting occasion-specific transitions succeeds. A 15-replicate `parameter_inference(method="bootstrap")` aborts with a plain `simpleError`: `the two signatures must have the same shape`.

When the unique longest person is omitted from a resample, the refit has fewer occasions and therefore a different coefficient/signature shape. The intended failed-replicate path is reached after the alignment that already throws.

**Repair direction:** validate matching shapes/parameter sets before alignment and catch extraction errors. Decide explicitly how to report intervals for coefficients supported only by resamples that contain the rare person or level; silently discarding other replicates conditions their distribution.

### E7 A two-wave stayer share is not identified

**Location:** transition setup has no guard for this design; the general count in `R/lta-engine.R:1472` counts the stayer coordinates.

For 220 two-wave persons, a two-state stayer model reports a stayer share of `0.429915`, `converged=TRUE`, three replicated best starts and no warning. Reconstructing its implied two-wave latent joint distribution and decomposing it using stayer shares `0.01`, `0.10` and `0.20` gives exactly the same likelihood, `-1853.50921685`, in every case.

The two-wave distribution determines the total mass in each latent state pair, but cannot distinguish staying within the mover component from membership in a separate stayer component.

**Repair direction:** warn or refuse unsupported short-panel specifications before fitting. Parameter counts and Wald checks do not substitute for the identification check. See [transition-ridge.log](transition-ridge.log).

### E8 Near-empty transition origins remain counted and labeled estimated

**Location:** `R/transitions.R:972` treats an origin as empty only when its expected outgoing total is zero; `R/lta-engine.R:1370` uses the same exact-zero rule.

The new simulation has three states and four waves, with state 3 appearing only at the final wave. The homogeneous fit has expected outgoing count `2.922e-7` for state 3, but reports all its outgoing transitions as estimated, counts 38 parameters, sets `boundary=FALSE` and raises no warning. Its reported state-3 self-transition probability is `0.99999999813`.

The occasion-specific engine has corresponding outgoing totals `2.789e-17`, `4.508e-24` and `2.064e-19`. It counts 50 parameters and marks all these rows estimated. Their reported probabilities vary widely despite practically no expected outgoing moves. This fit does issue a general probability-boundary warning, which does not mark these rows as unsupported.

**Repair direction:** specify a practical expected-count threshold or a separate weak-information diagnostic. The threshold is a modeling decision; these experiments establish that the current exact-zero test leaves numerically negligible origins in the reported estimate set. See [near-empty.log](near-empty.log).

## Chronology and checks that passed

The six experiments in `probes.R` were rerun against isolated source from `8518dffe` (1.1.9). [baseline-119-probes.log](baseline-119-probes.log) is byte-for-byte identical to the current-tree output. The stayer and near-empty checks also have baseline logs. This separates outstanding older defects from the latest changes.

Source inspection of `a10a919b` shows the singleton-centering, unweighted-grand-centering, membership-frame, ID-overwrite, stayer and bootstrap-alignment mechanisms already in 0.9.21. Grand scaling was added afterward, in the 0.9.38 period; its weight interaction predates 1.1.9. The two-wave growth path also exists in 0.9.21. These statements concern source lineage; the complete scripts were executed against 1.1.9 and current code only.

[sanity.R](sanity.R) independently verifies the Gaussian calculations on a newly simulated sample. The VVI observed-data log likelihood agrees with a direct product-mixture calculation exactly; multiplying both indicators by `1e-6` gives the expected likelihood Jacobian shift within `5.46e-12`; a one-component Gaussian regression agrees with `lm`'s ML log likelihood within `1.14e-13`. Output is retained in [sanity.log](sanity.log).

The full set of covariance structures, regression bootstrap likelihood-ratio tests, growth bootstrap, live parallel worker recovery, and all existing minor findings were not independently re-exercised in this slice. The retained November 4 reports cover broader checks; their results should be attributed to that review, rather than counted as fresh verification here.
