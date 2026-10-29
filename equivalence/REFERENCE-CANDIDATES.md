# Reference candidates for every major function

Compiled 2026-10-28 from three searches: the JStats collection (`../JStats`)
and the published literature. **None of these is wired into the harness yet.**
The numbers below are the published values, not latents results: nothing was
fitted to collect them.

The goal is at least three independent references per major function. A
published table with obtainable data is preferred; genuine output of another
program comes next.

**Status.** ✓ means the row's log-likelihood (or its first listed number) was
found in the saved text of the paper, manual or output file. ○ means it was
read online during the search but not re-checked. Other values in a row come
from the search notes; check each against its source before encoding it as a
reference.

## Current coverage

From `validation/release-fixes-2026-10-28/equivalence/report.csv` and the
testthat files listed in `README.md`.

| Function | References already in use |
|---|---|
| `lpa()` | mclust 6.1.2; Scrucca et al. (2023) mclust book, Example 3.2; mclust vignettes; Mplus UG 7.9/7.10; Latent GOLD 6.1; Houle et al. (2026) |
| `multilpa()` | Latent GOLD 6.1 (`twolevel-lpa`); Mplus two-level runs |
| `lca()` | poLCA 1.6.0.1; Linzer & Lewis (2011), Agresti & Lang (1993), Zhou & Lange (2010), Goodman (2002), Dayton (1998), McCutcheon (1987) |
| `multilca()` | Mplus 8.8 UG 10.7; glca 1.4.2; multilevLCA 2.1.6; Lyrvall et al. (2025) |
| `lta()` | depmixS4 1.5.1; LMest 3.2.8; one Mplus occasion-specific LTA demo; Latent GOLD step-3 Markov; Muthén & Asparouhov (2022), regular-LTA row |
| `mixture_regression()` | flexmix 2.3.20; EM with MASS `glm.nb` and VGAM `vglm(acat)` M-steps; one Mplus demo |
| growth (`random =`) | lcmm 2.2.1 |
| `three_step()`, `r3step()` | tidySEM 0.2.12; depmixS4; sandwich `vcovCL`; Mplus R3STEP/DU3STEP; Latent GOLD step 3 |
| `bootstrap_lrt()`, VLMR | Mplus TECH14, TECH11; Latent GOLD VLMR |
| `pool_imputations()` | mice 3.19.0 `pool.scalar` |

The weakest are `mixture_regression()`, growth and `lta()`: no published
paper, and growth has a single program.

## Candidates by function

### `lpa()`

| Reference | Data and access | Published numbers | To match | |
|---|---|---|---|---|
| Masyn (2013), *The Oxford Handbook of Quantitative Methods* vol. 2, ch. 25, doi:10.1093/oxfordhb/9780199934898.013.0025 | Reaven–Miller diabetes divided by 10, n = 145: `https://www.stata-press.com/data/r19/gsem_lca2.dta` (needs haven) | Table 25.12: LL, parameters, BIC, CAIC, AWE and adjusted LMR p for K = 1–5 under EEI, VVI, EEE, VVV (EEI-2 −1702.55; VVV-3 −1536.64, BIC 3217.61). Table 25.14: entropy .88, AvePP, OCC. Table 25.15: means, SDs, correlations | Mplus; BIC is −2LL + p ln n; values rounded to 2 dp | ✓ |
| Stata 19 [SEM] Example 52g | same data | EEI-2 LL −1702.5542 and VVV-3 LL −1536.6409, every parameter with its SE; `lcstats` AIC, BIC, LMR for K = 2–5 | Stata's 4-class LL (−1626.83) is a local maximum; Masyn reaches −1606.30 | ✓ |
| Morin et al. (2025), StepMix, *JSS* 113(8), Table 11, doi:10.18637/jss.v113.i08 | `iris`; `mclust::banknote` | iris VVV-3: −180.185 (mclust −180.186, Mplus −180.185). banknote VVI-2 with Status on membership: −771.669 (MoEClust the same, Mplus −771.668) | skip their diabetes row: it has a distal outcome and no dataset version | ✓ |
| Scrucca, Fop, Murphy & Raftery (2016), *R Journal* 8(1), doi:10.32614/RJ-2016-021 | old `mclust::diabetes`; `gclus::wine` | 14 structures × G = 1–9 BIC grid; wine EVE-3 BIC −6873.257; bootstrap LRT statistics 361.186445 / 114.703559 / 7.437806 | VVV-3 (BIC −4760.091) is a local maximum: Mplus reaches LL −2303.492 | ✓ |
| mclust 6.1.3 vignette | current `mclust::diabetes` | VVV-3 LL −2295.118, BIC −4734.561, all parameters; bootstrap LRT statistics 369.98537 / 121.10652 / 19.90071 | the dataset differs from the 2016 paper's | ✓ |
| Lebret et al. (2015), Rmixmod, *JSS* 67(6), doi:10.18637/jss.v067.i06 | `geyser` | K = 3 LL −1131.0738, BIC 2312.5998 | equal proportions forced; useful mainly as an engine with `weight =` | ✓ |

Gaps: no published FIML table on public data (StepMix's `gaussian_*_nan` is the
software fallback). No published weighted LPA (Rmixmod `weight =`, and the
identity that integer weights equal replicated rows). The noise component has
only mclust.

### `lca()`

| Reference | Data and access | Published numbers | To match | |
|---|---|---|---|---|
| Stata 19 [SEM] Examples 50g and 51g | Stouffer–Toby role conflict, n = 216: `gsem_lca1.dta` | 2 classes: LL −504.46767, every intercept with its SE, class probability .7207539 (SE .0580926). K = 1/2/3: G² 81.084 / 2.720 / 0.387, entropy .7193 / .6110, LMR | adds SEs, entropy and LMR to the Goodman (2002) values already used | ✓ |
| Beath (2017), randomLCA, *JSS* 81(13), doi:10.18637/jss.v081.i13 | `randomLCA::myocardial`, n = 94 | BIC 524.7441 (K = 1), 402.2951 (K = 2); LL −180.7027; class probabilities .5422/.4578; bootstrap LRT p = .001 (999 replicates) | the default `penalty = 0.01` shifts the optimum (penalised LL −180.7829) | ✓ |
| Mplus 8.8 UG ex7.6, ex7.7, ex7.8 | output and data in `../JStats/validation/r-reference/` (`mplus-ug-ex7.x-full.out.txt`, `mplus-ug-ex7.x.dat`) | ex7.6, 3 categories, n = 5000: LL −19214.480. ex7.7/7.8, nominal: LL −21408.482. Estimates with SEs | ex7.6 uses STARTS = 0, so its LL may be a local maximum | ✓ |
| StepMix Table 11, covariate simulation | `data_bakk_covariate(n_samples = 2000, sep_level = 0.9, random_state = 42)`, exported from Python | LL −5192.049 in StepMix, poLCA and Mplus | 6 binary items, covariate on membership | ✓ |
| Kleinman & Horton (2011), SAS and R blog, Example 8.21 | `mosaicData::HELPrct`, n = 431 complete (homeless, cesd > 20, satreat, linkstatus) | K = 3: SAS PROC LCA LL −1032.48; poLCA −1032.484, 14 parameters, AIC 2092.967, BIC 2149.893, shares .2162/.0785/.7053 (in `lca-packages/`) | df = 1; some probabilities on the boundary; current randomLCA refuses it | ✓ |

Also in JStats, unused: Mplus UG ex7.4 (binary, LL −1325.213) and ex7.12
(covariate with a class-invariant direct effect `u4 ON x`, LL −1255.396; only if
`lca()` supports that constraint).

### `multilca()`

Already four independent references. No new public source was found:
Vermunt (2003) and Henry & Muthén (2010) use private data, Mplus ex10.6 has a
parametric random effect, and Di Mari et al. (2023) shares data and software
with Lyrvall et al. Add the singleton-group identity below.

### `multilpa()`

No reproducible published reference: Mäkikangas et al. (2018) and Vermunt
(2008) use private data. Options are an Mplus `TWOLEVEL MIXTURE` run with
continuous indicators (adapt UG ex10.7) and the singleton-group identity,
which carries the `lpa()` references over.

### `lta()`

| Reference | Data and access | Published numbers | To match | |
|---|---|---|---|---|
| Bartolucci, Pandolfi & Pennoni (2017), LMest, *JSS* 81(4), doi:10.18637/jss.v081.i04 | bundled with LMest | `RLMSdat` (one ordinal item, T = 7, n = 1718), k = 3, homogeneous transitions: LL −13557.21, 20 parameters, BIC 27263.39, transition matrix and SEs. `data_SRHS_long`, covariates on initial and transition probabilities, k = 5: LL −62426.58. `data_criminal_sim` females, 2 clusters × 2 states: LL −18347.08, BIC 36925.18 | categories coded from 0; some transition logits near −30, so compare LL; LMest 3.x renamed `est_lm_*` to `lmest()` | ✓ |
| IMMERSE lca-bookdown, `github.com/immerse-ucsb/lca-bookdown` | LSAY, 2 waves × 5 binary items, N = 3092 (`data/lsay_lta.csv`, `lta/lta_model/`) | Mplus 8.11, 4 classes. Measurement invariant: LL −14542.768, BIC 29366.816. Occasion-specific: LL −14531.997, BIC 29506.005 | — | ✓ |
| Mplus 8.8 UG ch. 8 and 10, statmodel.com | posted with each example | ex8.12 hidden Markov: LL −13427.811. ex8.15 mover-stayer: LL −19504.367, mover/stayer .487/.513. ex10.13 two-level LTA with a between-level class: LL −20716.373 | ex10.13 needs no numerical integration, so it matches group classes | ○ |
| Pennoni, Pandolfi & Bartolucci (2024), LMest, *R Journal* RJ-2024-036, section 3.3 | `data_market_sim` (200 customers, T = 5; `brand` and `price` as two categorical items) | `lmestSearch()`, homogeneous transitions, k = 2/3/4: LL −2520.131 / −2445.538 / −2434.262, np 17 / 29 / 43, AIC 5074.263 / 4949.075 / 4954.524, BIC 5130.334 / 5044.727 / 5096.352 (in `lta-packages/`). Its covariate model (LL −2444.437, np 12) holds Psi fixed and uses LMest's `difflogit` slopes, which `lta()` does not parametrize the same way | k = 1 is not a transition model | ✓ |

No published weighted LTA: integer weights should equal LMest's
frequency-aggregated (`yv`) fit.

### `mixture_regression()`

| Reference | Data and access | Published numbers | To match | |
|---|---|---|---|---|
| Stata 19 [FMM] Examples 1a, 1b, 2 | `https://www.stata-press.com/data/r19/mus03sub.dta`; `gsem_mixture.dta` (needs haven) | 1a, `fmm 3: regress lmedexp income c.age##c.age totchr i.sex`: LL −4727.6738, N = 2955, class probabilities .2216/.7085/.0699, full coefficient and SE tables (1 and 2 classes: −4807.39, −4758.18). 1b, `lcprob(totchr)`: −4712.3871. 2, Poisson: −11502.686 | concomitant model for 1b | ✓ |
| SAS/STAT 14.2 PROC FMM, Example 40.3 | Ames salmonella (Margolin et al. 1981), 18 rows printed in the manual | 2-class Poisson, `dose` and `logd` slopes shared across classes: −2LL 121.8; intercepts 1.9097 (.2654), 2.4770 (.2731); mixing .8173 | shared slopes through `common =`; −2LL printed to one decimal | ✓ |
| Grün & Leisch (2008), FlexMix 2, *JSS* 28(4), doi:10.18637/jss.v028.i04 | bioChemists; betablocker (Aitkin 1999) | Poisson with kid5/mar/ment shared and `~fem` on membership: BIC 3182.328. Binomial with classes by Center, k = 3, common Treatment: LL −159.3605 | flexmix output, but fixed published anchors; betablocker exercises group-level classes | ✓ |
| Benaglia et al. (2009), mixtools, *JSS* 32(6), doi:10.18637/jss.v032.i06 | `mixtools::CO2data`, n = 28 | 2-class Gaussian: λ .7549/.2451, σ 2.0493/0.8094, LL −66.93977 | started from user values, possibly a local maximum | ✓ |
| Deb & Trivedi (1997), *J Appl Econometrics* 12(3) | `AER::NMES1988` | 2-class negative binomial | NB1 vs NB2; tables not opened | ○ |

No published ordinal mixture regression with open data was found. Also in
JStats: an Mplus 2-class Poisson with offset (`mplus-demo-offset-poisson-classspecific`,
LL −2045.335) and MoEClust ML fixtures (`tests/fixtures/latent/mixture-experts-moeclust-ref.json`).

### Growth: `mixture_regression(random =)`

| Reference | Data and access | Published numbers | To match | |
|---|---|---|---|---|
| Proust-Lima, Philipps & Liquet (2017), lcmm, *JSS* 78(2), doi:10.18637/jss.v078.i02 | `lcmm::paquid`, `normMMSE` | fixed `poly(age65, 2, raw) + CEP`, random and mixture on `poly(age65, 2)`, equal random covariance. ng = 2: LL −8899.228, 15 parameters, BIC 17891.67. ng = 3: −8891.351 from the default start, −8890.148 from the grid search | the ng = 3 pair is a ready local-maximum test | ✓ |
| Mplus 9 Demo LCGA and GMM, in `../JStats/tests/fixtures/mplus/` (`lcga-ml/`, `runner-lcga/`, `gmm-ml/`, `runner-gmm-growth/`) | `.dat` beside each output | LCGA K = 2: −2091.416 (300 × 5). GMM K = 2: −2727.520, K = 3: −2715.247. Estimates with SEs, BIC, entropy | residual variance class-specific and equal across waves, as latents fits; GMM matches `random_covariance = "varying"` | ✓ |
| Mplus 8.8 UG ex8.9 | posted | binary LCGA: LL −1272.307, classes .462/.538 | binomial, subject-level classes | ○ |
| IMMERSE LSAY science GMM | `gmm/gmm_enum` in the repo above, N = 5866, 6 waves | 3 classes: LL −63254.481, 31 parameters | its residual variances vary by occasion | ○ |

Mplus UG ex8.1 (LL −3183.703) leaves residual variances free by occasion and is
out of scope unless constrained. AR and Brownian residuals have no published
GMM numbers: use one-class identities against `nlme::lme` with `corAR1` and
`corCAR1`.

### `three_step()`, `r3step()`

| Reference | Data and access | Published numbers | To match | |
|---|---|---|---|---|
| StepMix (2025), *JSS* 113(8), Table 9 | GSS 1976–77 through `gssr` | BCH income means 26.71 (.74) / 36.81 (1.64) / 44.68 (4.26); ML 21.05 / 44.73 / 61.26; ML model LL −18754.72, 38 parameters | FIML on income | ✓ |
| Bakk & Kuha (2021), *BJMSP* 74(2), open access PMC8247311 | same GSS data (ICPSR) | Table 4: BCH 26.65 (.60) / 37.06 (1.35) / 44.70 (3.09); ML unequal variances 20.63 / 44.94 / 65.67; Table 6, covariates | Latent GOLD syntax in Appendix A | ○ |
| Mplus Web Talk 8, in `../JStats/validation/r-reference/mplus-webtalk8-bch.out.txt` with `mplus-webtalk8-distal.dat` | N = 1000, missing coded 999 | 3-class LCA LL −3634.865, BIC 7407.885; y1 BCH means −0.830 (.137), 0.337 (.140), 1.173 (.101); overall Wald 132.564; BCHC for y2, y3 | — | ✓ |
| IMMERSE replication of Nylund-Gibson, Grimm & Masyn (2019), doi:10.1080/10705511.2019.1590146 | `distals/data/lsay_df.csv`, N = 2675 | BCH step 3: LL −12744.043. ML step 3: LL −12742.252 | step-1 files saved in the repo | ✓ |

Also in JStats: a manual ML three-step distal regression
(`mplus-demo-lca-three-step-distal-regression`, LL −4272.523).

### `bootstrap_lrt()`, VLMR and diagnostics

- mclust bootstrap LRT statistics (2016 paper and the current vignette, above).
  The statistics are deterministic; the p-values are not. Both papers' VVV fits
  are partly at local maxima. ✓
- randomLCA `myocardial`: bootstrap LRT p = .001. ✓
- Adjusted LMR p-values: Masyn Table 25.12 (Mplus TECH11) and Stata
  `lcstats`. Masyn's row K tests K against K + 1; Stata's tests K − 1 against K. ✓
- Mplus 9 Demo runs in JStats: `validation/r-reference/mplus-lca-tech11-s1/s2`,
  `mplus-lca-tech14-s1/s2`, `mplus-lca-cprob` (data in the `lca_2class` and
  `lca_3class` rows of `tests/fixtures/latent/mplus-ref.json`), and
  `tests/fixtures/mplus/lca-boundary-tech11/` (12 seeds). ✓
- Classification diagnostics: Masyn Table 25.14 (VVV-3: entropy .88, AvePP
  .958/.918/.973, OCC 21.74/41.86/94.06). ✓

### `pool_imputations()`

No published pooled LCA or LPA table exists. Add `mitml::testEstimates` and
`mitools::MIcombine` as implementations independent of mice. mice's
`pool.scalar(n = Inf)` fmi is (r + 2/(df + 3))/(r + 1): check which fmi latents
reports. The Mplus imputation demos in JStats are factor models and check only
the pooling arithmetic.

## Identity checks that need no external reference

- Every group of size one: `multilca()` must equal `lca()`, and `multilpa()`
  must equal `lpa()`, for any number of group classes. This carries every
  single-level reference to the two-level functions.
- Integer weights equal replicated rows, in every family.
- One class with AR or Brownian residuals equals `nlme::lme` with `corAR1` or
  `corCAR1`.

## Before using these

- Three different "diabetes" datasets circulate: (a) old `mclust::diabetes`,
  identical to UCLA `EXAM1.dat`; (b) mclust ≥ 6, relabelled fpg/glucose/insulin
  with `glucose[104] <- 455`; (c) Reaven–Miller divided by 10, used by Masyn and
  Stata.
- Where a published fit is a local maximum (mclust 2016 VVV, Stata 52g with
  four classes, lcmm ng = 3 from its default start), assert that latents
  reaches at least that LL rather than equality.
- Data outside base R: Stata `.dta` files need haven; GSS needs `gssr`; the
  StepMix simulation must be exported from Python.
- The `jstats-mplus/PROVENANCE.md` scan covered only
  `../JStats/validation/r-reference/` and missed `tests/fixtures/mplus/`. The UG
  data files there are named `mplus-ug-ex7.x.dat`, while the inputs name
  `ex7.x.dat`.
- flexmix's Gaussian M-step divides by n − rank, so its fixed point is not the
  ML optimum. The flexmix suite already replaces that step
  (`harness/suites/flexmix.R`, lines 55–63).

## Sources

- Masyn (2013): https://www.statmodel.com/download/Masyn_2013.pdf
- Stata [SEM] examples: https://www.stata.com/manuals/semexample50g.pdf, `51g`, `52g`
- Stata [FMM] examples: https://www.stata.com/manuals/fmmexample1a.pdf, `1b`, `2`
- SAS PROC FMM 14.2: https://support.sas.com/documentation/onlinedoc/stat/142/fmm.pdf
- StepMix: https://www.jstatsoft.org/article/view/v113i08
- mclust 5 (2016): https://journal.r-project.org/archive/2016/RJ-2016-021/RJ-2016-021.pdf
- mclust vignette: https://cran.r-project.org/web/packages/mclust/vignettes/mclust.html
- Rmixmod: https://www.jstatsoft.org/article/view/v067i06
- randomLCA: https://www.jstatsoft.org/article/view/v081i13
- SAS and R Example 8.21: https://sas-and-r.blogspot.com/2011/01/example-821-latent-class-analysis.html
- LMest (2017): https://www.jstatsoft.org/article/view/v081i04
- LMest (2024): https://journal.r-project.org/articles/RJ-2024-036/
- Mplus User's Guide: https://www.statmodel.com/usersguide/chapter7.shtml, `chapter8.shtml`, `chapter10.shtml`
- IMMERSE lca-bookdown: https://github.com/immerse-ucsb/lca-bookdown
- FlexMix 2: https://www.jstatsoft.org/article/view/v028i04
- mixtools: https://www.jstatsoft.org/article/view/v032i06
- lcmm: https://www.jstatsoft.org/article/view/v078i02
- Bakk & Kuha (2021): https://pmc.ncbi.nlm.nih.gov/articles/PMC8247311/
- `AER::NMES1988`: https://rdrr.io/pkg/AER/man/NMES1988.html
- Vermunt (2003): https://jeroenvermunt.nl/sm2003.pdf
- mice `pool.scalar`: https://amices.org/mice/reference/pool.scalar.html
