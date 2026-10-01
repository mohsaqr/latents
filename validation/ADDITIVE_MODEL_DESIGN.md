# Additive multilevel model specification and validation design

Draft implementation contract, 30 September 2026, against latents 0.8.8.
The proposed first implementation estimates group classes from continuous raw
ratings, with a Gaussian distribution of group means inside each class and a
shared within-group residual variance. It includes exact marginal likelihood,
multiple starts, group classification, and inference for regular interior fits.
It is a new model family; this document and its numerical checks do not add a
shipped estimator or establish parameter recovery or interval coverage.

## Source and scientific scope

Houle et al. distinguish additive models with shared within-group variances
and either shared or class-specific between-group variances. Their additive
climate syntax places raw ratings at both levels; their contextual syntax
uses manifest group means as separate between-level indicators. These are
different reference specifications. See the [author prepublication](https://www.statmodel.com/download/ML-LPA%20FINAL%20-%20Prepub.pdf),
printed p. 8 and supplemental PDF pp. 35–36. The
[published article](https://doi.org/10.1177/10944281261469432) identifies the six
families. The retained [package comparison](../equivalence/houle-2026/REPORT.md)
records the limitations of the original study materials.

**Proposed choice:** implement the raw-rating Gaussian model below first,
with both between-variance restrictions. Its equations and tests are derived
independently here. This is a candidate structural match to the climate syntax;
a fresh matched external fit must establish numerical equivalence. Interpreting
ratings as contextual or climate measurements remains a study-design decision.
No automatic centering, manifest aggregation, or measurement-error correction
is implied by the model name.

## Model definition

Let j index J independent observed groups, i index n_j observations in a group,
r index d continuous indicators, and h index H latent group classes.
Conditional on the supplied group sizes, the proposed generative model is

```text
C_j ~ Categorical(omega_1, ..., omega_H)
B_j | C_j = h ~ Normal_d(mu_h, diag(T_h))
Y_ij | B_j, C_j = h ~ Normal_d(B_j, diag(W))
```

W is a positive d-vector of within-group variances shared by all classes;
T_h is a nonnegative d-vector of between-group variances. The between restriction
is either varying (a vector per class) or equal (one vector for all classes).
Each indicator has its own group intercept; there is no scalar intercept with
unit loadings shared by all indicators. Conditional residuals and the
between-level covariance matrices are diagonal. Observations are independent
given their group intercept, and groups are independent. Class weights sum to
one and do not depend on group size; informative group sizes require a separate
model or sensitivity study.

There are no individual latent profiles. Every member shares the group's class,
but has individual residual variation. Group classes may differ in both location
and between-group spread when T varies. The defining additive restriction is
the shared within-group variance, not exclusively differences in means.

## Exact likelihood and estimation

For complete data define a_jr = mean_i(Y_ijr) and
S_jr = sum_i((Y_ijr - a_jr)^2). The original-observation component log density is

```text
ell_jh = -0.5 * sum_r [
  n_j * log(2*pi) + (n_j - 1) * log(W_r)
  + log(W_r + n_j*T_hr)
  + S_jr/W_r
  + n_j*(a_jr - mu_hr)^2/(W_r + n_j*T_hr)
]

log L = sum_j logsumexp_h(log(omega_h) + ell_jh)
```

This follows from the covariance W_r*I + T_hr*11' of a group's raw observations.
It preserves the normalizing constants. A density for group means alone drops
the residual evidence and cannot be substituted when estimating W or comparing
raw-observation likelihoods. Sufficient statistics make evaluation O(JHd)
after O(Nd) preparation, without quadrature.

The E-step returns class weights R_jh and Gaussian intercept moments:

```text
R_jh = softmax_h(log(omega_h) + ell_jh)
M_jhr = mu_hr + n_j*T_hr/(W_r + n_j*T_hr) * (a_jr - mu_hr)
V_jhr = T_hr*W_r/(W_r + n_j*T_hr)
```

With Q_h = sum_j R_jh and N = sum_j n_j, the proposed interior EM updates are

```text
omega_h_new = Q_h/J
mu_hr_new = sum_j R_jh*M_jhr / Q_h
T_hr_new = sum_j R_jh*[V_jhr + (M_jhr - mu_hr_new)^2] / Q_h
W_r_new = sum_jh R_jh*[S_jr + n_j*((a_jr - M_jhr)^2 + V_jhr)] / N
```

All conditional moments come from the old parameters. For equal T, pool its
numerator across classes and divide by J. Between-level updates count groups;
within-level updates count observations. A declared positive floor on W is a
constraint, and bound-active fits must report that fact.

Initialize T strictly positive. At T_hr = 0 the latent-intercept EM update
freezes both that mean and variance. Support for exact zero therefore requires
optimization of the marginal likelihood with nonnegative T, including explicit
boundary candidates and convergence/KKT checks. Do not call an EM start at zero
an independently optimized boundary solution. Preserve multiple-start records,
caller RNG state, and all failed starts; final selection must reproduce the
likelihood at the returned parameters. Use plain EM initially. Acceleration
requires separate monotonicity and same-solution evidence.

## Identification and parameter accounting

For varying T the nominal free count is H*d + H*d + d + H - 1.
For equal T it is H*d + d + d + H - 1. Thus H=2, d=2 gives 11 or 9 parameters.
Any later fixed blocks reduce the count explicitly. A zero estimated variance
remains part of the nominal model; it must not silently reduce an information
criterion penalty as though it had been fixed before fitting.

Repeated observations identify within-group variance through their contrasts.
If all groups are singletons, W_r and T_hr are inseparable: only W_r + T_hr
is observed. Refuse that specification. Mixed singleton and repeated groups
are permitted, subject to identification diagnostics. An indicator with zero
pooled within-group scatter can drive W to its floor and must be flagged.

H must not exceed J. Duplicate components, vanishing group-class weights,
too few informative groups, and a singular information matrix can still make
a permitted input unidentified or weakly identified. These conditions need
classed diagnostics and inference refusals. A minimum group count by itself
does not prove identification. Labels are arbitrary; align class means and
between variances jointly. Include H=3 with a non-self-inverse permutation in
alignment tests.

## Proposed public interface and results

A separate fitting verb avoids giving the existing required `n_profiles`
argument an artificial meaning. The provisional interface is

```r
additive_lpa(data, vars, id, n_group_classes,
             between_variance = c("varying", "equal"),
             n_starts = 20, max_iter = 1000, tol = 1e-8,
             min_variance = 1e-6, seed = NULL)
```

This name and defaults are proposals, not available functions. Resolve naming
before exporting. Store an explicit family tag and a distinct S3 class so
cross-sectional accessors cannot accidentally manufacture individual profiles.
The same engine can later implement dispersion and additive-dispersion by
changing which mean and variance blocks are shared.

The first public increment uses complete numeric ratings and known group IDs,
diagonal covariance at both levels, and H fixed by the caller. It returns
group-class means, within and between variances, class weights, group posteriors,
posterior intercept means/variances, convergence and boundary diagnostics,
restart records, parameter counts, log likelihood, and information criteria.
Return class-conditional intercept moments with an explicit class index;
also offer marginal moments including uncertainty about class membership:

```text
intercept_mean_jr = sum_h R_jh*M_jhr
intercept_variance_jr = sum_h R_jh*[V_jhr + (M_jhr - intercept_mean_jr)^2]
```

Distinguish uncertainty about a group's intercept conditional on fitted
parameters from uncertainty in estimated population parameters.

Tables use the existing `group_class`, `indicator`, `estimate` vocabulary;
variance rows carry `level = "within"` or `"between"`. Group posterior and
assignment tables have one row per observed group. A join to member rows must
identify those values as inherited group assignments. Entropy and class shares
use groups as their denominator. Define `nobs()` as J and expose N separately.
Use J for the default BIC, retain N as an
explicit alternative convention, and report both sample sizes. Never mix a
raw-observation likelihood with an incompatible comparison likelihood.

Existing `get_results()`, `coef()`, `vcov()`, `confint()`, `logLik()`, `nobs()`,
`summary()`, and ggplot results should dispatch on the new class. Define
`names(coef()) == rownames(vcov())` on the free estimation scale; tidy parameter
tables use natural units through a checked delta-method map. Reuse extraction
conventions, not an old estimator class whose fields imply individual profiles.

FIML, correlated residuals or group intercepts, covariates, sampling weights,
staged fitting, auxiliary outcomes, transition dynamics, and bootstrap class
tests are later increments. Requests for these combinations must raise classed
refusals. Grand-mean translation can be reversible; group-mean centering removes
the location information required by this specification and is not an optional
preprocessing switch in the first increment.

## Inference contract

For regular interior solutions encode means, log W, log T, and H-1 baseline
class logits. Differentiate the exact per-group likelihood to obtain scores;
derive observed information and the group-score sandwich, with an independent
numerical derivative check away from the optimum. Retain the package's CR0
convention and distinguish observed, robust, and OPG covariance estimates.
Log-scale intervals keep variances positive; baseline-logit and softmax maps
provide natural-scale class-weight uncertainty.

Refuse ordinary Wald inference for zero or bound-active variances, collapsed
classes, unconverged estimates, nonstationary solutions, and singular
information. Group-score covariance needs sufficient rank, not just finite
entries. Bootstrap boundary inference is a separate validation task. Comparing
H and H+1 must not use an ordinary chi-square reference distribution.

## Independent reference cases

1. Compare the compact density with Cholesky evaluation of the full group
   covariance at identical supplied parameters, including unequal group sizes.
2. For one indicator, integrate the conditional Gaussian density over its
   random intercept using an independent numerical integrator.
3. For H=1 and balanced n>1, check the independent interior MLE:
   mu = grand mean, W = sum_j S_j/[J*(n-1)],
   T = sum_j(a_j-mu)^2/J - W/n. If this T is negative, the boundary MLE has
   T=0 and W=[sum_j S_j + n*sum_j(a_j-mu)^2]/(J*n).
4. With T fixed to zero, compare with a shared-residual Gaussian mixture whose
   class is assigned to the whole group, including posterior probabilities.
5. Prepare fresh Mplus H=1, H=2 varying T, and H=2 equal T runs on identical
   synthetic raw data. Use the climate syntax as a starting point, explicit
   zero covariances, retained TECH1 constraints and parameter counts, and ML
   first. Compare MLR separately with the matching sandwich convention.

Retain data, input/output, seeds, software versions, checksums and label maps
under `equivalence/`. Compare likelihoods at identical parameters before
comparing optimizer maxima. Proposed full-precision tolerances: per-group
likelihood discrepancy <=1e-8, posterior discrepancy <=1e-7, aligned mean
discrepancy <=1e-4 in observed indicator-SD units, aligned variance discrepancy
<=1e-4 in squared observed indicator-SD units, class-weight discrepancy <=1e-6
in absolute probability units, and regular-fit SE discrepancy <=1%
relative. Round-limited external outputs use explicit printing intervals;
looser precision must not hide an optimizer disagreement. Any changed
tolerance needs justification before examining the final benchmark.

## Simulation design and release decisions

Prepare a selected scenario design rather than a full Cartesian grid. First
pilot with 100 attempted datasets per cell, then use at least 1,000 attempted
datasets for each regular inference cell. Do not replace failed datasets until
1,000 successes are accumulated.

| Scenario | Proposed generating conditions | Primary purpose |
|---|---|---|
| Regular reference | H=2, d=2, J=200, n=10; equal class weights; means (-1,-1) and (1,1); W=(1,1); varying T=(0.25,0.5)/(0.5,0.25) or shared T=(0.4,0.4) | Recovery, observed and robust intervals |
| Small group sample | Same model, J=50 | Finite-group limitations |
| Unequal sizes | J=200, sizes cycled through 2,5,10,20; sizes independent of class | Correct likelihood and variance weighting |
| Singleton admixture | J=200, 25% of groups size 1 and others size 10, independent of class | Identification with incomplete replication |
| Three classes | J=300, equal weights, means (-1.5,-1.5),(0,0),(1.5,1.5); W=(1,1), shared T=(0.3,0.3) | Alignment and parameter accounting |
| Rare class and weak separation | J=200, weights 0.9/0.1; means (-0.4,-0.4)/(0.4,0.4); W=(1,1), shared T=(0.4,0.4) | Collapse and misleading inference |
| Variance boundary | J=200, means (-1,-1)/(1,1), W=(1,1), shared T=(0,0) or (0.01,0.01) | Boundary optimization and refusal |
| Misspecification | Regular reference with correlated Gaussian residuals (correlation 0.4), or variance-scaled t(5) residuals | Limits of the diagonal model and sandwich claims |

The first six cells vary precision and identification within the stated model;
the boundary cell requires valid optimization and refusals, not ordinary Wald
coverage. Misspecification cells describe sensitivity and cannot certify the
Gaussian model for general nonnormal data. Record the exact parameter choice
and seed for each subcell in a machine-readable registry before running it.

For each parameter report bias, empirical SD, mean estimated SE, 95% coverage,
interval width, Monte Carlo SE, and effective successful count. Report
unconditional fit/interval availability and coverage conditional on an interval
being available. Failed fits and withheld intervals remain separate categories.

**Proposed release criteria for the regular reference and unequal-size cells:**
at least 98% fit-and-interval
availability; absolute bias <=0.1 empirical SD; estimated-to-empirical SE ratio
0.9–1.1; 95% coverage between 92.5% and 97.5% with Monte Carlo uncertainty
reported. These are engineering acceptance choices, not statistical guarantees.
Cells failing them require correction, a stated support limitation, or more
evidence; selecting only favorable parameters is not a remedy. Difficult cells
must reveal their failures reliably, even where calibrated inference cannot
be established. No claim of coverage is based on likelihood identity checks.

## Work sequence and current evidence

1. Settle the public name and the raw-rating contract; preserve existing model
   behavior. Implement the exact density, sufficient statistics and group
   posterior calculations with independent tests.
2. Implement both T restrictions, starts, EM and marginal boundary polishing.
   Establish the balanced H=1 oracle and fresh external comparison before
   treating the fitter as validated.
3. Implement checked derivatives, inference, tables and ggplot dispatch;
   exercise the complete user workflow and unsupported-combination refusals.
4. Run the predeclared simulations, document supported scenarios and failures,
   add a runnable tutorial, then run package tests and R CMD check. Release
   scope depends on the evidence rather than the prototype's runtime.

The [design check script](additive-design-check.R), run on synthetic data,
currently verifies density identity, Gaussian intercept moments, interior EM
monotonicity, row and label invariance, the T=0 density limit, the singleton
variance ridge, and balanced one-class interior/boundary MLEs against a separate
dense-likelihood optimizer. On 30 September 2026 the compact/dense discrepancy was
7.11e-15, the integration discrepancy 3.55e-15, and the one-step EM likelihood
increase 5.2448744158. The retained [check output](additive-design-check.txt)
records the numerical results. These checks validate the proposed formulas only.
External fits, a shipped estimator, parameter recovery, and coverage remain
outstanding. Step 1's internal engine (R/additive.R: sufficient statistics,
exact density, group posteriors and intercept moments) is implemented and
tested in tests/testthat/test-additive-engine.R; it is not exported and the
public name is unresolved. Step 2 ships as the experimental
`multilpa(family = "additive", between_variance = ...)` (owner decision,
30 September 2026, replacing the proposed separate verb): multistart EM,
exact zero-between updates and a KKT boundary pass, tested against the
balanced H=1 closed forms, dense-likelihood stationarity and a group-level
mixture (tests/testthat/test-additive-fit.R). Mplus references, accessors,
inference and simulations remain outstanding. Step 3 (tables, S3 methods,
analytic-score inference, plots, vignette) and step 4's predeclared
simulation are done; see ADDITIVE_SIMULATION.md (31/34 reference-cell
parameters meet the gates; between-variance ML bias is a stated limitation;
rare weakly separated classes are an open issue). Mplus references remain
outstanding (Mplus unavailable here). Full delivery exceeds the backlog's one-to-two-day prototype
estimate; re-estimate after the engine and external reference pass.
