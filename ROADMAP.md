# latents development roadmap

**Planning baseline:** repository source inspected on 23 September 2026;
`DESCRIPTION` reports version 0.7.0.  
**Status:** proposed development plan. Features listed as planned are
not commitments or implemented capabilities.  
**Objective:** provide a coherent R workflow for multilevel latent
profiles, latent classes, and latent transitions, with flexible
specification, dependable inference, and independently validated
results.

## 1. Direction and scope

The immediate objective is to complete the analytical workflows already
supported by the package. The next objective is to support richer
longitudinal and measurement models. Continuous random effects are a
subsequent extension that requires substantial statistical development.

Success means that a researcher can specify the model required by the
study, evaluate its assumptions, quantify uncertainty, and reproduce the
complete analysis through a consistent interface. Each release should
deliver a usable scientific workflow, including documentation and
evidence for its inferential claims.

Full Mplus equivalence is outside the scope of this roadmap. Growth
mixtures, general factor mixtures, survival models, and general
structural equation modeling remain possible later directions. Their
inclusion should follow an explicit decision to broaden the package’s
remit and account for existing R implementations.

## 2. Existing foundation

The following capabilities already exist and should be preserved:

| Area | Current implementation | Important boundary |
|----|----|----|
| Measurement | Gaussian, categorical, and mixed indicators | Categorical measurement uses unrestricted response probabilities |
| Multilevel structure | Observation profiles nested within discrete latent group classes | Continuous group random effects are deferred |
| Gaussian covariance | Fourteen covariance structures; enumeration across covariance models | Wald inference covers EEI, VVI, EEE, and VVV; bootstrap inference covers all fourteen |
| Missing indicators | Observed-data likelihood in base models and LTA | The one-step membership-covariate path rejects FIML |
| Covariates | Profile- and group-membership regressions | Complete numeric predictors; restricted combinations with other features |
| Auxiliary analysis | BCH means/contrasts and R3STEP at observation or group level | First-stage classification-error uncertainty is treated as fixed |
| Staged estimation | Measurement estimated first, then held fixed | First-stage uncertainty is not propagated |
| Transitions | Gaussian/categorical/mixed LTA, optional latent group classes, observed/grid occasion handling | Homogeneous first-order transitions, invariant measurement, no transition predictors or parameter inference |
| Evaluation | Multiple starts, diagnostics, information criteria, sensitivity analysis, parametric bootstrap class comparisons for supported base models | LTA enumeration and likelihood-ratio testing are not implemented |
| Independent checks | Retained external comparisons and simulation infrastructure | Coverage and difficult-data validation need expansion for new workflows |

Baseline source: [main
estimator](https://pak.dynasite.org/latents/R/multilpa.R),
[transitions](https://pak.dynasite.org/latents/R/transitions.R),
[inference](https://pak.dynasite.org/latents/R/inference.R), [three-step
analysis](https://pak.dynasite.org/latents/R/three-step.R), [deferred
features](https://pak.dynasite.org/latents/future/README.md), and
[equivalence
infrastructure](https://pak.dynasite.org/latents/equivalence/README.md).

The current enumeration interface uses `model =` for covariance models.
New examples and tests must follow the current API rather than earlier
`structure =` examples.

## 3. Phases, effort, and costs

Estimates assume one experienced statistical R developer, with
statistical review available. They include implementation, tests,
documentation, and independent validation. A person-week means 40 hours.
Costs use an **illustrative rate of €100/hour**; this is a budgeting
assumption, not a market quotation.

| Phase | Deliverable | Dependencies | Effort | Illustrative cost |
|----|----|----|---:|---:|
| 0 | Capability contracts and validation design | None | 1–2 person-weeks | €4k–8k |
| 1 | LTA bootstrap uncertainty; FIML with membership covariates; predictor formulas | Phase 0 | 6–12 person-weeks | €24k–48k |
| 2 | Initial/transition predictors; time-varying transitions; LTA candidate comparison | Shared design infrastructure and LTA validation from Phase 1 | 8–14 person-weeks | €32k–56k |
| 3 | Known-group/partial-invariance models; uncertainty-aware auxiliary workflows | Parameter design from Phases 0–2 | 12–22 person-weeks | €48k–88k |
| 4A, optional | Poisson and negative-binomial measurement | Stable measurement and inference contracts | 4–8 person-weeks | €16k–32k |
| 4B, optional | Continuous random-intercept mixtures with inference | Integration design and independent reference models | 8–16 person-weeks | €32k–64k |
| 4C, optional | Random-intercept LTA | Phase 2 plus a validated random-effects design | 12–24 person-weeks | €48k–96k |

**Core programme, Phases 0–3:** approximately **27–50 person-weeks
(€108k–200k)** before contingency. The first practical milestone, Phases
0–1, is approximately **7–14 person-weeks (€28k–56k)**.

Phase estimates assume shared infrastructure and are not sums of
independent feature quotations. Optional phases should be re-estimated
after design work; overlap with earlier phases is uncertain. Allow a
separate 25% planning reserve for numerical or statistical redesign.
Compute, external software licences, procurement, and review waiting
time are excluded. Person-weeks describe effort, not guaranteed calendar
completion; parallel staffing will not reduce time proportionally.

## 4. Phase 0 — Establish contracts and validation targets

**Scientific outcome:** a precise statement of which model is estimated
by each feature combination and what evidence is required for a release.

### Work

Record a capability matrix crossing model family, indicator type,
covariance model, missingness, covariates, fixed parameters, and
inference method.

Assign each combination a status: supported, planned, or intentionally
unsupported. Verify supported combinations from executable examples.

Define parameter storage and naming for initial logits, transition
logits, predictor terms, and equality constraints.

Specify how model objects retain factor levels, contrasts,
transformations, reference categories, and fitted row order.

Design bootstrap replicate records, label alignment, failure handling,
and seed management. Reuse the existing bootstrap infrastructure where
appropriate.

Choose benchmark models and record their exact parameterization,
identification constraints, missingness handling, and sample-size
convention.

Define simulation scenarios and acceptance criteria before examining
results.

### Acceptance gate

Every Phase 1 task has a written statistical specification, an
identified source implementation, and at least one independently
checkable reference case. Existing supported behaviour is captured by
regression checks. Proposed public arguments remain provisional until
the design review is complete.

## 5. Phase 1 — Complete current analysis workflows

### 1A. Bootstrap inference for latent transitions

**Question enabled:** how precisely are initial and transition
probabilities estimated?

Resample complete observed groups or people, retaining their sequences,
indicator missingness, and occasion grid semantics.

Give duplicated sampled groups distinct bootstrap identifiers.

Refit the same model specification using reproducible multiple starts.

Align both profile and latent group-class labels. Apply profile
permutations consistently to both axes of each transition matrix.

Return parameter standard errors, percentile intervals, and
requested/completed/failed replicate counts.

Distinguish nonconvergence, collapsed classes, boundary estimates, and
ambiguous label matches.

Define an explicit failure policy. Never silently discard failures and
report unconditional validity for the surviving replicates.

Flag transition rows without estimable outgoing information; do not give
fallback probabilities ordinary inferential status.

**Acceptance:** recovery and interval coverage are assessed for one and
multiple group classes, unequal sequence lengths, and each supported
measurement family. Label permutations leave results unchanged.
Published limitations specify where reliable intervals could not be
established.

### 1B. Missing indicators in membership-covariate models

**Question enabled:** how do predictors relate to membership when some
indicators are missing?

Extend the covariate expectation step to the observed-data measurement
likelihood.

Use the appropriate conditional moments for partially observed Gaussian
vectors and omit unobserved categorical responses from their likelihood
contributions.

Extend parameter inference consistently with the new likelihood; guard
combinations whose inference remains unsupported.

Preserve the complete-data result as a limiting case.

Specify behaviour for fully missing indicator rows and groups and
confirm identifiability checks.

**Boundary:** missing predictors require a separate modeling or
imputation strategy. Indicator FIML does not automatically handle them.

**Acceptance:** estimates and inferential quantities agree with
independent reference fits within declared numerical tolerances.
Simulation scenarios include complete data, MCAR, and explicitly
specified MAR mechanisms; unsupported mechanisms are not described as
solved.

### 1C. Predictor formulas and design matrices

**Question enabled:** how do categorical predictors, interactions, and
specified transformations relate to membership?

Add a formula-based interface while retaining existing character-vector
calls.

Construct and store design matrices using familiar R conventions.

Preserve contrasts and reference levels across refits and bootstrap
samples.

Detect rank deficiency, missing predictors, and invalid group-level
variation before estimation.

Test loss of a factor level in a resample without silently changing the
parameterization.

**Acceptance:** formula and manually constructed numeric-design fits
agree. Output identifies predictor contrasts and membership reference
categories. A worked vignette combines factor predictors, missing
indicators, and interpretable membership effects.

### Phase 1 release gate

The three workflows have documented capability boundaries, successful
external comparisons, simulation evidence, and executable tutorials. No
claim of support is based solely on an argument being accepted by a
function.

## 6. Phase 2 — Model predictors and changes in latent dynamics

**Scientific outcome:** estimate how initial membership and transitions
depend on covariates and occasion.

### Work

Add separate specifications for initial-profile membership, transitions,
and, where supported, latent group-class membership.

Define how a time-varying predictor is aligned: origin occasion,
destination occasion, or an explicitly specified lag. Record that choice
in output.

Estimate multinomial transition regressions conditional on the origin
profile.

Begin with shared predictor slopes and clearly identified intercepts;
add class-specific slopes only with a justified identification and
validation plan.

Permit transitions to differ across predefined occasions or periods,
alongside the homogeneous model.

Reuse Phase 1 bootstrap inference for the new parameters after
validating label alignment and sparse transition behaviour.

Extend candidate enumeration to LTA with common diagnostics, parameter
counting, information criteria, and restart records.

Report predicted transition probabilities at specified covariate values
as well as logit coefficients.

### Acceptance gate

The no-predictor model reproduces existing LTA, and equal time-specific
parameters reproduce homogeneous transitions. Independent fits agree for
at least one transition-covariate and one time-varying specification.
Simulations examine rare transitions, sparse covariate patterns, and
short sequences. Information-criterion conventions are explicit.

LTA likelihood-ratio tests require a separate calibration task.
Enumeration does not imply that an ordinary chi-square class-count test
is valid. Analytic and sandwich standard errors can follow once the full
transition score has been derived and checked; they are not assumed to
be included in this phase’s bootstrap deliverable.

## 7. Phase 3 — Evaluate measurement and external-variable models

### 3A. Known groups and partial measurement invariance

**Question enabled:** do profiles have comparable meanings across
populations or occasions?

Introduce observed known groups separately from the nesting identifier
and from latent group classes.

Implement a parameter map for fixed values, free values, and selected
equality constraints.

Support fully invariant measurement and selected departures across
observed groups or occasions.

Require sufficient anchors and document identification and
label-comparability assumptions.

Add direct covariate effects on selected indicators as a separately
validated specification.

Provide comparisons for supported nested invariance models with fixed
class counts and appropriate regularity checks.

**Acceptance:** the fully constrained model reproduces the current
invariant model. Constraint accounting, free-parameter counts, and
inference agree with independent references. Simulated noninvariance is
distinguishable from changes in membership within the tested scenarios.

A model allowing arbitrary differences across every group, class, and
occasion is not an initial deliverable. Each added dimension needs
identification and data-support checks.

### 3B. First-stage uncertainty in staged and auxiliary analysis

**Question enabled:** how uncertain are external-variable associations
when the measurement solution was itself estimated?

Bootstrap the full sequence: measurement fit, label alignment,
classification-error calculation, and auxiliary model.

Resample at the independent group level and preserve the required
within-group structure.

Apply the same principle to staged estimation with estimated first-stage
measurement parameters.

Report conditional and full-pipeline uncertainty as distinct inferential
options.

Evaluate sensitivity to weak separation and unstable
classification-error matrices.

**Acceptance:** simulation reports bias, coverage, interval width, and
failed replications. It compares the current conditional approach with
the proposed full-pipeline procedure; bootstrap support is released only
for combinations with defensible evidence.

### 3C. Broader external outcome models

Add covariate-adjusted continuous outcome models.

Add binary/categorical outcome models with an explicitly justified
classification-error correction.

Specify handling of missing external outcomes separately from missing
measurement indicators.

Validate each outcome family independently; do not assume that generic
weighted regression is valid for every BCH weight configuration.

Count distal outcomes and arbitrary secondary SEM are subsequent
extensions, requiring revised effort estimates.

### Phase 3 release gate

A complete vignette demonstrates an observed-group or longitudinal
invariance analysis and an external-variable analysis with propagated
measurement uncertainty. Documentation distinguishes known groups,
latent group classes, and sampling clusters throughout.

## 8. Phase 4 — Optional expansion of the model family

### 4A. Count measurement

Start with Poisson and negative-binomial indicator likelihoods,
exposure/offset handling where appropriate, and mixtures of count and
existing indicator families. Validate overdispersion, zero frequencies,
boundaries, missing indicators, and parameter inference.
Zero-inflated/hurdle models are a subsequent decision and are not
included in the initial estimate.

### 4B. Continuous random-intercept mixtures

Use the deferred code as a prototype to audit. Specify where the random
effect enters: indicator means, membership logits, or both. These are
different models and should have separate contracts.

Before release, establish multiple-profile recovery against an
independent implementation, integration accuracy across quadrature
settings, identifiable parameterizations, and an explicit
variance-boundary inference strategy. The 8–16-week estimate concerns a
bounded random-intercept specification; it does not include general
random slopes or arbitrary correlated random effects.

### 4C. Random-intercept LTA

Develop a stated model separating persistent person heterogeneity from
occasion-specific dynamics. The deferred static random-intercept code
does not by itself provide RI-LTA. Validate the
measurement/random-effect relationship, transition parameterization,
numerical integration, predictor effects, and uncertainty against
matched RI-LTA references.

### Later decisions

Random slopes and general continuous random effects were provisionally
estimated at 16–32+ person-weeks (€64k–128k+). Growth/factor mixtures
have a similarly broad preliminary range. Survey weights, strata,
additional observed nesting levels, and crossed classifications require
their own design studies. None is included in the core programme budget.

A VLMR test is lower priority than these workflows because a parametric
bootstrap comparison already exists for supported base models. Restore
the deferred test only after independently validating its reference
distribution and p-values; a conventional chi-square substitution is not
a solution.

## 9. Validation and release policy

Every release should supply four types of evidence:

1.  **Likelihood and implementation checks:** limiting cases, row/order
    invariance where appropriate, label permutations, parameter counts,
    and derivatives checked away from an optimum.
2.  **Independent numerical comparisons:** matched identification and
    constraints; retained inputs and outputs; justified tolerances for
    likelihoods, parameters, predictions, and standard errors where
    available.
3.  **Monte Carlo evidence:** bias, empirical standard deviation,
    estimated standard errors, interval coverage and width, convergence,
    boundaries, and label-alignment failures.
4.  **User workflow checks:** executable vignettes, stable extraction
    and plotting interfaces, informative unsupported-combination errors,
    and package checks on supported platforms.

The simulation design should vary the dimensions relevant to each
feature: number and size of groups, unequal group sizes, separation,
rare classes, covariance structure, indicator type, missingness,
sequence length, and transition sparsity. Use a planned scenario design
rather than an unnecessarily large full Cartesian grid.

For each coverage estimate, report Monte Carlo uncertainty. With 1,000
successful replications and coverage near 95%, the approximate Monte
Carlo standard error is 0.7 percentage points. Pilot runs may use fewer
replications; final replication counts should follow the desired
precision and computational burden. Report failed runs separately and
explain whether coverage is conditional on successful estimation.

Establish acceptable performance and numerical tolerances before
inspecting final results. Failures must lead to correction, an explicit
support restriction, or a documented research limitation. Unit-test
success alone does not establish calibrated inference.

## 10. First development cycle

The first two person-weeks should produce reviewable groundwork:

| Task | Output | Completion condition |
|----|----|----|
| Audit current support | Capability matrix | Every claimed combination linked to an executable check |
| Design LTA resampling | Bootstrap specification | Group resampling, grid semantics, label alignment, and failures defined |
| Audit covariate likelihood | FIML implementation design | Required expectation, maximization, and inference changes identified |
| Specify predictor interface | Formula/design-matrix proposal | Backward compatibility and resampling behaviour documented |
| Define evidence | Reference and simulation plan | Models, scenarios, measures, tolerances, and release decisions stated |

Implementation can then begin with LTA bootstrap inference, followed by
FIML with covariates and formula support. Independent reference
preparation can proceed alongside implementation. Statistical
definitions and parameter contracts should be settled before separate
contributors build dependent components.

## 11. Positioning and reference material

The package’s contribution is the combination of multilevel measurement,
dynamics, external variables, and evaluation in a consistent workflow.
Existing tools already cover parts of this space, so development
priorities should be reviewed against their capabilities before
commissioning a new model family.

Reference material used in the preceding capability assessment:

- [Mplus multilevel mixture
  modeling](https://www.statmodel.com/HTML_UG/chapter10V8.htm).
- [Mplus latent transition analysis and
  covariates](https://www.statmodel.com/Webtalk2.shtml).
- [Mplus random-intercept LTA models and
  examples](https://www.statmodel.com/RI-LTA.shtml).
- [Mplus BCH auxiliary-variable technical
  note](https://www.statmodel.com/examples/webnotes/webnote21.pdf).
- [glca
  documentation](https://search.r-project.org/CRAN/refmans/glca/html/glca.html).
- [multilevLCA
  documentation](https://search.r-project.org/CRAN/refmans/multilevLCA/html/multiLCA.html).
- [LMest mixed latent Markov
  models](https://search.r-project.org/CRAN/refmans/LMest/html/lmestMixed.html).
- [lcmm documentation](https://cecileproust-lima.github.io/lcmm/).
- [OpenMx model examples](https://openmx.ssri.psu.edu/node/4480) and
  [tidySEM](https://cjvanlissa.github.io/tidySEM/).

Reassess estimates and priorities after each phase using actual
development effort, benchmark results, and researcher needs. Assign
release versions only after the scope and acceptance evidence are
established.
