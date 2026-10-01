# Group-class families: simulation evidence

Predeclared study of `multilpa(family = "additive")`, run 30 September –
1 October 2026 on R 4.5.2, latents 0.8.8 plus the uncommitted additive family.
Design and release gates: [ADDITIVE_MODEL_DESIGN.md](ADDITIVE_MODEL_DESIGN.md),
"Simulation design and release decisions". Script:
[additive-simulation.R](additive-simulation.R). The cell registry
([additive-simulation-registry.csv](additive-simulation-registry.csv)) was
written before the first replicate; the script refuses to run if it changes.
Every number below is read from
`additive-simulation-full-summary.csv` (per parameter) and
`additive-simulation-full.txt` (per cell); per-replicate rows are in
`additive-simulation-full-replicates.csv`.

## Procedure

- 11 cells x 1000 attempted datasets (pilot: 100). Each replicate has its own
  seed (`base_seed + replicate`), so results do not depend on scheduling.
- Each dataset: `multilpa(..., family = "additive", between_variance = <cell>,
  n_starts = 10)`; then the `"parameters"` table with `vcov_type = "observed"`
  and `"robust"`. Classes are aligned to the generating classes by the
  permutation minimizing the distance between class means.
- Every warning and message is recorded per replicate, not suppressed.
- Failed fits and withheld intervals are separate categories; availability is
  counted over attempted datasets.

## The pilot found a defect, fixed before the full run

The first pilot ([additive-simulation-pilot-before-newton.txt](additive-simulation-pilot-before-newton.txt))
recorded `latents_unconverged` from inference (scaled score above 1e-2) in
13/100 datasets of `regular_varying`, 29 of `unequal_sizes`, 87 of
`three_classes` and 97 of `rare_weak`, although every EM run had met its
likelihood-change rule. EM had stopped short of the maximum. Interior fits are
now finished by safeguarded Newton steps (`.additive_newton()`); the rerun
pilot and the full run recorded no such warning in any cell. Regression test:
"Newton polishing reaches the maximum EM stops short of".

## Cell status (full run)

| Cell | Fitted | Boundary (SEs withheld) | Intervals available |
|---|---|---|---|
| regular_varying | 1000 | 0 | 1000 |
| regular_shared | 1000 | 0 | 1000 |
| small_groups (J = 50) | 1000 | 6 | 994 |
| unequal_sizes | 1000 | 0 | 1000 |
| singleton_admixture | 1000 | 0 | 1000 |
| three_classes | 1000 | 0 | 1000 |
| rare_weak | 1000 | 0 | 995 |
| boundary_zero | 1000 | 773 | 217 |
| boundary_small | 1000 | 355 | 634 |
| misspecified_correlated | 1000 | 0 | 1000 |
| misspecified_t5 | 1000 | 0 | 1000 |

No fit failed and no replicate raised a script error.

## Per-cell summary

Pass = parameters meeting all four predeclared gates (availability >= 98%,
|bias| <= 0.1 empirical SD, observed SE / empirical SD in 0.9-1.1, observed
95% coverage in 92.5-97.5%). Coverage Monte Carlo SE is about 0.007 at 1000
datasets.

| Cell | Pass | max abs(bias)/SD | SE ratio (observed) | Coverage observed | Coverage robust |
|---|---|---|---|---|---|
| regular_varying | 12/12 | 0.076 | 0.92-1.01 | 0.926-0.955 | 0.912-0.951 |
| regular_shared | 8/10 | 0.146 | 0.95-1.02 | 0.930-0.956 | 0.924-0.956 |
| unequal_sizes | 11/12 | 0.119 | 0.93-1.03 | 0.925-0.953 | 0.911-0.954 |
| singleton_admixture | 12/12 | 0.091 | 0.97-1.01 | 0.930-0.961 | 0.929-0.951 |
| three_classes | 12/13 | 0.108 | 0.96-1.03 | 0.938-0.959 | 0.932-0.961 |
| small_groups | 6/12 | 0.149 | 0.87-1.05 | 0.911-0.972 | 0.896-0.958 |
| misspecified_correlated | 9/12 | 0.136 | 0.96-1.02 | 0.932-0.955 | 0.937-0.951 |
| misspecified_t5 | 10/12 | 0.072 | 0.56-1.03 | 0.730-0.966 | 0.923-0.967 |
| rare_weak | 2/10 | 0.743 | 0.39-1.01 | 0.544-0.950 | 0.544-0.949 |
| boundary_zero | 0/10 | 0.659 | 0.95-1.68 | 0.000-0.954 | 0.000-0.959 |
| boundary_small | 0/10 | 0.090 | 0.98-1.38 | 0.907-0.962 | 0.905-0.961 |

## Release gates: regular reference and unequal sizes

The gates apply to `regular_varying`, `regular_shared` and `unequal_sizes`.
31 of 34 parameters pass. The three that do not are all between-group
variances, and all fail on bias alone:

| Cell | Parameter | bias / SD | SE ratio | Coverage (observed) |
|---|---|---|---|---|
| regular_shared | between.shared.y1 | -0.102 | 0.978 | 0.947 |
| regular_shared | between.shared.y2 | -0.146 | 0.951 | 0.930 |
| unequal_sizes | between.1.y2 | -0.119 | 0.932 | 0.925 |

Every between-group variance in every regular cell is estimated low (20 of
20). At 200-300 groups the bias is 0.8-2.8% of the true value; at 50 groups
(`small_groups`) it is 3.3-6.4%, roughly four times larger. That is the size and the
1/J rate of the ordinary small-sample downward bias of a maximum-likelihood
variance component. **Stated support limitation:** between-group variances
carry ML small-sample bias of about 0.1-0.15 empirical SD at 200 groups; their
interval coverage stays within the gate (0.925-0.955). A REML-type correction
for the mixture is not implemented. Per the design, the gates were not changed
after seeing these results, and the family stays experimental.

## Difficult cells

- **Small groups (J = 50).** Means cover 0.911-0.929 and between-variance
  robust SEs are 14-17% too small (robust coverage 0.896-0.937). Observed SEs
  are closer (0.87-0.92 for between variances). Prefer `vcov_type =
  "observed"` with few groups.
- **Rare, weakly separated class (0.9/0.1, means +-0.4).** Inference is
  misleading and *not* refused: intervals are available in 99.5% of datasets,
  but the large class's means are biased by 0.40-0.43 SD (the rare class's by
  0.14-0.15 SD), mean SEs are 0.41-0.54 of the empirical SD (coverage
  0.640-0.785), and the weight SEs are 0.39 of the empirical SD with
  intervals covering 0.544. Addressed by the weak-class diagnostic below.
- **True between variance zero.** 77.3% of datasets estimate it at zero and
  standard errors are withheld (correct refusals). In the other 21.7% the
  estimate is positive and its log-scale interval cannot contain zero, so
  coverage of the between variance is 0.000 by construction; means, within
  variances and weights still cover 0.931-0.954. A variance whose true value may
  be zero needs a boundary test, not a Wald interval.
- **Between variance 0.01.** 35.5% on the boundary; observed SEs for the
  between variance are 1.35-1.37 times the empirical SD (conservative).
- **Non-normal residuals, t(5).** The observed SEs of the within variances are
  0.56-0.57 of the empirical SD (coverage 0.730-0.733); the robust SEs are 0.94-0.97
  (coverage 0.923-0.945). Use `vcov_type = "robust"` when residual tails may
  be heavy. Correlated residuals (0.4) left inference calibrated (coverage
  0.932-0.955); three parameters missed the bias gate (max 0.136 SD).

## Weak-class diagnostic

Protocol written before any value was computed:
[additive-diagnostic-protocol.md](additive-diagnostic-protocol.md); script
[additive-diagnostic.R](additive-diagnostic.R); output
[additive-diagnostic.txt](additive-diagnostic.txt), per-dataset values in
`additive-diagnostic-datasets.csv` (9840 datasets with intervals; the refits
reproduce the stored estimates, e.g. 5.6e-16 on `rare_weak` replicate 7).

Statistic: effective groups per class, `w^2 (1 - w) / Var(w)`, the minimum
over classes. Rule on replicates 1-500: largest grid value with at most 1%
flags in the five reference cells, giving **50**. Held-out replicates 501-1000:

| | Effective groups < 50 (shipped) | Entropy < 0.80 | Smallest AvePP < 0.85 |
|---|---|---|---|
| Flag rate, reference cells | 0.0080 | 0.0048 | 0.0004 |
| Flag rate, rare_weak | 1.000 | 0.722 | 0.881 |
| Coverage of class means and weights, unflagged | 0.9477 (MCSE 0.0014) | 0.9406 | 0.9398 |
| Coverage, flagged | 0.7998 (MCSE 0.0051) | 0.6360 | 0.6890 |
| All three declared criteria | met | detection fails | met |

The effective-groups warning also flags every `small_groups` dataset (50
groups; no class can exceed 25 effective groups) and 1.2% of
`singleton_admixture`, 2.2% of `three_classes` and 0.6% of `unequal_sizes`.
It flags `small_groups` for sample size, not separation; there
class-parameter coverage was 0.927 (MCSE 0.005), a real but milder shortfall.
Median effective groups in `rare_weak` were 2.7, against 20.3 in
`small_groups` and 69-95 elsewhere. The threshold was not revisited after the
held-out results. Shipped as the `latents_weak_class` warning with
`effective_groups` in the `group_classes` table and `min_effective_groups` in
the `fit` table; the message reports both the expected class size and the
share of information kept, so the reader can tell small from poorly separated.
It does not refuse inference.

## Dispersion and additive-dispersion families

Script [families-simulation.R](families-simulation.R), registry
[families-simulation-registry.csv](families-simulation-registry.csv) written
before the first replicate; the same procedure and gates as above, 1000
datasets per cell, 200 groups. Output: `families-simulation.txt`,
`families-simulation-summary.csv`, `families-simulation-replicates.csv`.
No fit failed, every fit converged, none was on a boundary and none raised
`latents_weak_class`.

| Cell | Pass | max abs(bias)/SD | SE ratio (observed) | Coverage observed | Coverage robust |
|---|---|---|---|---|---|
| dispersion_reference (within 0.5 vs 2.0, shared means and between) | 10/10 | 0.064 | 0.97-1.04 | 0.937-0.964 | 0.936-0.958 |
| dispersion_unequal_sizes (sizes 2/5/10/20) | 10/10 | 0.100 | 0.95-1.02 | 0.938-0.960 | 0.928-0.960 |
| additive_dispersion_varying | 14/14 | 0.079 | 0.97-1.03 | 0.933-0.962 | 0.934-0.960 |
| additive_dispersion_equal | 11/12 | 0.103 | 0.96-1.01 | 0.931-0.959 | 0.926-0.959 |

The one miss is the shared between-group variance of y2 in
`additive_dispersion_equal` (bias -0.103 SD, SE ratio 0.977, coverage 0.941):
the same ML small-sample downward bias of a between-group variance as in the
additive cells, under the same stated limitation.

## Not established

- External equivalence with Mplus: Mplus is not installed on this machine.
  Latent GOLD 6.1 serves as the external reference instead
  (`equivalence/latentgold-families/`).
- Simulation evidence for the cross-level families, which report no standard
  errors; their likelihood is checked against the published formula only.
- Coverage for designs outside these cells; any claim about between-variance
  intervals when the true value is near zero.
