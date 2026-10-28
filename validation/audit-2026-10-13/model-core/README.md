# Model core audit — 2026-10-13

This review covers 31 estimation, measurement and inference files, containing 402 top-level function definitions and approximately 15,916 source lines. Production files and the owner's pending changes were not edited. The audit artifacts are separate from the shipped tests.

## Findings

### MC-1 · P2 · Covariate dispatch silently ignores requested centering

`multilpa()` matches `centering` at R/multilpa.R:1225, but the covariate dispatch at lines 1280–1292 omits it. Neither `.multilpa_covariate_model()` nor `.multilpa_fit_covariates()` accepts it. Consequently, a call combining `centering = "person"` or `"grand"` with profile/group covariates fits the original indicator values without a refusal or notice.

Verified public reproduction: 30 persons with six occasions each, a continuous indicator near 100 and a profile membership covariate. The three requested modes produce exactly identical profile means (96.42578, 103.06085) and likelihood -476.2394. The returned covariate model has no centering field. Fitting the explicitly person-centered data instead produces profile means -2.391015 and 2.555912, likelihood -288.3874. This changes the model being estimated, beyond a cosmetic labeling issue.

Evidence: `dispatch-probes.R`, `dispatch-probes.log`, `dispatch-probes.rds`. Relevant functions: `multilpa`, `.multilpa_covariate_model`, `.multilpa_fit_covariates`, `.multilpa_cov_designs`. Corrective direction: carry centering and its recorded offsets through the covariate estimator and its downstream data/inference contract, or refuse unsupported combinations before fitting. Provenance: the lost-argument dispatch predates 0.9.12; its existing forwarding path originates in 504ae1f4 with later additions through the 0.9.8 review (0d14b540).

### MC-2 · P2 · Variance floors violate constrained covariance models

R/covariance-structure.R:175 applies componentwise `pmax()` after the constrained diagonal update. R/covariance-structure.R:403 similarly clips eigenvalues separately after the constrained full-covariance update. A floor applied independently to one profile can break constraints tying profiles together.

Verified internal examples, effective profile weights `(10, 10)`:

- EVI, scatter diagonals `(0.001, 10)` and `(1, 1)`, floor 0.1: returned covariance determinants are 0.10 and 0.01. EVI requires equal determinants.
- VEI, scatter diagonals `(0.01, 0.1)` and `(1, 10)`, floor 0.1: the first-to-second-profile variance ratios are 0.4805062 and 0.1519494. VEI requires the same ratio across dimensions.
- EVV with diagonal scatter matrices `(0.001, 10)` and `(1, 1)`, floor 0.1: returned determinants are 0.55 and 0.01. EVV requires equal determinants.

A public synthetic two-profile example with `min_variance = 0.05` fails every start for EVI, VEI, EVV and VEV with `latents_all_starts_failed` caused by an EM likelihood decrease. The input and floor are valid under the public API. The simple unconstrained VVI/EEE/VVV flooring conventions are not implicated by this finding.

Evidence: `covariance-probes.R`, `covariance-probes.log`, `covariance-probes.rds`; the VEI and EVV direct examples were additionally executed during review. Corrective direction: solve the floor constraint inside the structured optimization, preserving volume/shape/orientation constraints and ascent. Repeated clipping cannot supply the promised constrained M-step.

Provenance: the diagonal clipping line originates in commit 504ae1f4, 2026-09-21; it predates the 0.9.12–0.9.23 work. This is a newly discovered existing defect, not evidence of a recent regression.

### MC-3 · P2, extreme-scale numerical limitation · Growth quadratic loses residual information

R/growth-mixture-engine.R:148 evaluates the Woodbury quadratic by subtracting two large terms: `(r'r - u'A^-1u / sigma2) / sigma2`. The 0.9.15 direct residual fallback fixes cancellation in forming `r'r`, but does not fix this subsequent cancellation when random-effect variance is huge compared with residual variance.

An independent scalar random-intercept reference evaluates the same exact Gaussian density using

`sum((y - mean(y))^2) / sigma2 + n * mean(y)^2 / (sigma2 + n * G)`.

For `y = c(1e8 - 1, 1e8, 1e8 + 1)`, regression mean zero, `G = 1e16`, and `sigma2 = 1`, `.growth_class_terms()` returns log density -21.72680 instead of -23.22680: its implied quadratic is zero instead of three. At amplitude 1e7 the likelihood error is 0.03125; ordinary-scale controls at 1e3/1e4 agree closely.

Evidence: `probes.R`, `probes.log`, `probes.rds`. This is an internal likelihood evaluation reproduction, not a demonstrated bad public fitted model at ordinary scales. The public variance floor is relative to outcome variance and reduces exposure at default settings; unusually small user floors or strongly ill-scaled random-effect designs deserve particular attention. A stable fallback based on posterior residuals/orthogonal projection could preserve the low-dimensional speed while avoiding subtraction of nearly equal quadratic terms. The covariance/residual score calculations should be checked alongside it.

Provenance: this expression originates in the 0.9.12 growth release, commit 6ef582f7. The subsequent 0.9.15 stabilization leaves it present.

### MC-4 · P3 documentation · Prior arguments have intentional, undocumented model exceptions

The source introduction in R/prior.R:8–24 explicitly preserves mclust EM conventions that differ from textbook conjugate updates: EII omits prior scale, EEI does not shrink means, and VEV ignores degrees of freedom. This is intentional compatibility behavior and is covered by equivalence tests; it should not be characterized as a newly incorrect estimator implementation.

The public help nonetheless presents `shrinkage`, `mean`, `dof` and `scale` as generic effective conjugate-prior controls and requests MAP estimation without exposing those structure-specific exceptions. A user supplying strong prior information can reasonably expect a control to affect the advertised parameter.

Verified knob probes on one-profile synthetic data: with prior mean 20 and shrinkage 100, EEI means remain 0.004908 and 0.089744, whereas VVI means become 10.00245 and 10.04487. Changing EII scale from 1 to 1e6 leaves its variance exactly unchanged. Changing VEV dof from 4 to 400 leaves its variances exactly unchanged.

Evidence: `prior-probes.R`, `prior-probes.log`, `prior-probes.rds`. Corrective direction: explain the structure-specific conventions in `prior_control()` help and qualify the MAP wording where the implementation deliberately follows the reference program's departures. Source provenance: commit b876c259, 2026-09-27, before 0.9.12.

## Coverage and evidence

`function-inventory.csv` contains every top-level function in the assigned files, exact source line bounds, formals, review status, confirmed finding references, and static usage diagnostics. All 402 definitions were parsed/indexed and checked with `codetools::checkUsage()` in the actual loaded package namespace. The source review considered estimator dispatch, likelihood definitions, missing-data integration, constrained updates, starts, pack/unpack charts, scores, boundary behavior, sampling weights and table-label transforms across the assigned families.

There are direct name references in test source for 140 of the 402 functions. The remaining 262 cannot be called "untested" on that basis: many internal functions are exercised through exported estimators, S3 dispatch, and helper calls. Conversely, a string reference does not prove runtime execution or adequate edge-case coverage. The full root-run test and equivalence suite reported 6,349 successful expectations across 958 cases; this is integration evidence, not a proof for every mathematical branch.

The inventory deliberately distinguishes static review from targeted synthetic runtime probes. The four desired invariants in `test-findings.R` currently fail as expected (two covariance constraints, the growth quadratic and covariate centering); the final saved regression log has no harness errors. Only the functions involved in the saved synthetic cases receive targeted probe status. No assertion of individual runtime execution for all 402 functions is made. Formal identification proofs, a complete parameter-space numerical sweep, real-world datasets and a fresh independent optimizer reference for every covariance/prior structure are outside this pass.

### Areas reviewed

| Area | Review focus and observations |
|---|---|
| Shared EM and latent structures | Relative stopping rules, decrease guards, SQUAREM fallback, start selection, quasi-Newton retry, three regression nestings, weighted score aggregation and forward-backward handling. No additional confirmed defect beyond the findings above. |
| Gaussian and categorical measurement | Conditional sufficient statistics by missingness pattern, observed marginal densities, held means, covariance bounding, exact bounded multinomial shares, deterministic category encoding and dimensions. Conditional covariance arithmetic and categorical KKT updates are coherent in reviewed paths. |
| Covariance structures | Fourteen structure codes, parameter counts, starting-value projection, diagonal/full M-steps, common-orientation partial solves and floors. MC-2 is the main unsupported constraint/ascent edge. |
| Ordinal/count indicators | Adjacent-category ordinal charts, scores, location identification, Poisson/NB density functions, large-size dispersion asymptotics, boundary updates, typed draws, names and tables. Recent stable NB density/derivative changes are present. |
| Profile/class and covariate estimators | Single-level wrappers, argument refusals, category/start alignment, fixed blocks, group indexing, weighted posteriors and counts, inference data identity. MC-1 identifies a lost preprocessing argument. |
| Transition engines and inference | First/second-order augmentation, mover/stayer constraints, unobserved occasions, occasion coordinate removal, full-covariance chart packing, inference refusals and bootstrap relabeling. Covariance-chart support for general models is deliberately separated from their restricted Wald score path. |
| Growth mixtures | Per-person cross-products, pattern sharing, trajectory starts, multilevel structure, covariance scales/permutations, EM updates, analytic scores and quasi-Newton. MC-3 concerns extreme random-effect/residual ratios; the recent trajectory-start and location-residual remedies are present. |
| Regression mixtures | Outcome support, finite evaluated designs, Gaussian/binomial/Poisson/NB/ordinal blocks, Newton steps, screening, coefficient coordinates and boundary conditional inference. Recent single-class binary and stable count-density corrections are present. |
| Additive and cross-level families | Sufficient statistics, exact random-intercept density, weighted M-steps, zero-between KKT polishing, Newton finish, IC parameter counts and group-score inference. Cross-level inference is correctly refused because the working likelihood uses ratings and their derived group means. |
| Weights, noise, prior and staged measurement | Overflow-resistant normalization, independent-unit weighting, weighted-only robust inference policy, noise restrictions/splitting, held-block consistency, stage compatibility and explicit reference-prior conventions. MC-4 concerns public explanation of deliberate prior exceptions. |

## Recent-work assessment

The 0.9.12/0.9.13 additions materially expanded the engine surface: random effects, multilevel growth, NB and ordinal regression, shared EM/structure services, and new derivative/coordinate paths. The 0.9.15 fixes address several real numerical and inference issues, and the corresponding source logic remains present. MC-3 is a remaining numerical edge in the recently introduced growth engine. MC-1 and MC-2 are earlier defects exposed by combination/constraint probes rather than failing broad suite tests. The latest uncommitted enumeration/named-plot/documentation work was preserved and reviewed by other audit tracks; this core review is not a separate release/API compatibility certification.

Known open issue from HANDOFF/LEARNINGS: regression specification formulas can capture fitting frames and inflate serialized fits. This was not presented as a newly discovered issue. Additional unverified concerns, such as approximation from 12-significant-digit growth design-pattern keys and extended-LTA iteration/history bookkeeping, were not promoted to confirmed findings.

## Reproduction

From the package root, run the individual `*-probes.R` scripts with `Rscript`; they load the current local package through pkgload. Their saved logs and RDS files retain actual synthetic outputs. `inventory.R` rebuilds the function inventory and compact comment-free source extracts used to inspect the larger files. `test-findings.R` records currently failing desired invariants for the three estimator issues and is intentionally outside `tests/testthat/`.
