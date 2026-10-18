# Monte Carlo equivalence

Simulated datasets with known parameters, each fitted by genuine **Mplus
VERSION 9 DEMO (Mac)**, **Latent GOLD 6.1** (under Wine) and this package, all
on 2026-10-17 and 2026-10-18. Part of Equivalence 2.0 (roadmap, Priority E).
Every comparison reads the package's reported values, matches classes once
from the posteriors and reuses that matching for every class-level quantity,
checks ids for completeness, and takes its tolerance from what the reference
prints plus the measured posterior difference `d` of the two optimizers.

## Designs

**Main kit** (`kit/`, `make-kit.R`, `run.sh`): 8 configurations x 15
replicates, within the Mplus Demo limits (6 dependent, 2 independent, 2
between variables). Separation is moderate (entropy 0.55 to 0.87).

| Configuration | Model | Data |
|---|---|---|
| lpa_eei | 3 profiles, EEI | n = 500, 4 indicators (simulab) |
| lpa_vvi | 3 profiles, VVI | n = 500, 4 indicators (simulab) |
| lpa_vvv | 2 profiles, full covariance | n = 500, 3 correlated indicators (simulab) |
| lca_binary | 3 classes | n = 800, 6 binary items (simulab) |
| ml_lpa | 2 profiles x 2 group classes, VVI | 60 groups of 15 (simulab) |
| ml_lca | 2 classes x 2 group classes | 60 groups of 15, 5 binary items |
| lpa_missing | 3 profiles, VVI, FIML | MAR on 3 indicators (simulab) and 5% empty rows |
| lpa_covariates | 2 profiles, VVI, `c ON x1 x2` | n = 500 |

**Extended kit** (`extended/`): 4 configurations x 10 replicates: two-level
EEE and VVV against both programs; latent transition models (2 states, 4
occasions) with missing values and with shared variances against Latent GOLD.

**Enumeration kit** (`enumeration/`): `lpa_vvi_01` at k = 1 to 5 and
`lca_binary_01` at k = 1 to 4, through Mplus.

## Tests and results

| Test | What is compared | Result |
|---|---|---|
| `test-montecarlo.R` | per case, against Mplus: log likelihood (8 significant digits), parameters, N, BIC, row-by-row posteriors, entropy (7 decimals single level), AvePP and the whole AvePP matrix, OCC, class counts, profile means implied by Mplus's posteriors, group-level posteriors, entropy and AvePP (Mplus's row averaging), every criterion (AIC, KIC, BIC, SABIC, CAIC, ICL, AWE, CLC) from Mplus's own quantities, `predict()`; against Latent GOLD: log likelihood, parameters, posteriors of every row; `drop_empty = FALSE` counts against Latent GOLD | 3,461 expectations pass |
| `test-montecarlo-extended.R` | the extended kit, against both programs | 403 pass |
| `test-enumeration-mplus.R` | every candidate of an enumeration grid against Mplus at its k; the BIC minimum; `compare_models()` BIC weights; `cores = 2` identical to `cores = 1` | 49 pass |
| `audit-maxima.R` (retained `audit-maxima.csv`, checked by `test-montecarlo.R`) | every case from 5 seeds with 10 and 20 starts against the best of Mplus and Latent GOLD | 1,200 of 1,200 fits reach it; worst gap -4.8e-5 |

Largest differences on the main kit: log likelihood 4.0e-5, BIC 3.5e-5,
posteriors 1.8e-4, single-level entropy 3.6e-6.

Recovery of the generating parameters is reported in
`validation/MONTECARLO_RECOVERY.md`.

## Findings

- **AvePP depends on the modal class, which jumps.** In `lpa_covariates_05`
  two rows sit within 2 d of a tie, so the two programs assign them to
  different classes and AvePP differs by 0.0017 although the posteriors
  differ by 1.6e-4. The tolerance adds the number of such rows over the
  smallest class (`modal_slack()`).
- **Over-extracted models.** For `lpa_vvi_01` at 4 profiles Mplus's best
  (-3105.390, not replicated) is degenerate (a class of 2 observations with
  a variance of 0.000); its second (-3106.613) was reached twice. This
  package reaches -3108.709 from 30 to 300 starts and -3106.051 from 1000:
  its random starts are converged k-means partitions, which explore such a
  surface less widely than Mplus's screened random starts. At the true
  number of classes every fit reaches the reference maximum (the audit
  above). Nothing is asserted where Mplus warns that its solution is not
  replicated.
- **`predict()` refuses membership-covariate models**
  (`latents_unsupported_prediction`), which the test checks.

## Rerunning

From the project root: `Rscript equivalence/montecarlo/make-kit.R`, then
`sh equivalence/montecarlo/run.sh` (Mplus four at a time, Latent GOLD in
turn; `JOBS=` sets the parallelism). Each run removes a case's earlier
outputs first and stops on a run that does not terminate normally or writes
nothing. The extended and enumeration kits have their own `make-kit.R` and
`run.sh`. The tests run offline from the retained `returned/` outputs.
