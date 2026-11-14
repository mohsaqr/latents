# Mplus User's Guide examples: what latents fits, and what it would need

Every example of the Mplus 8.8 User's Guide chapters 7, 8 and 10
(statmodel.com/usersguide, `chapter7.zip`, `chapter8.zip`, `chapter10.zip`;
each posts input, data and Mplus's output). Mplus's log likelihood and
parameter count are from its output. "In this kit" examples are compared in
`compare-functions.R`; the rest name the function or the feature latents
would need. **Before writing a new function in the mixture family, look here
(and in chapters 11-13 when they apply) for an Mplus example to compare it
with.**

Status: **kit** (this kit) · **elsewhere** (compared in another file) ·
**needs** (a latents feature that would make it comparable) · **out of
scope** (not a mixture model latents means to fit).

## Chapter 7: mixtures with cross-sectional data

| Example | Model | Mplus LL (p) | Status |
|---|---|---|---|
| 7.1 | Mixture regression, common and class-specific slopes, x on membership | -836.899 (9) | kit: `mixture_regression(common =, membership =)` |
| 7.2 | Zero-inflated Poisson mixture regression | -902.299 (10) | needs: zero-inflated count family in `mixture_regression()` |
| 7.3 | LCA, binary indicators, automatic starts | -965.244 (9) | kit; also `test-jstats-mplus.R` |
| 7.4 / 7.5 | LCA, binary indicators, user-specified starts | -1325.213 (9) | kit (7.4; 7.5 is the same model and data) |
| 7.6 | LCA, three-category indicators | -19214.480 (17) | kit; also the JStats fixtures |
| 7.7 / 7.8 | LCA, nominal indicators | -21408.482 (17) | kit (7.7; 7.8 the same model) |
| 7.9 | LPA, class-invariant variances (EEI) | -3177.162 (13) | kit: `lpa(model = "EEI")` |
| 7.10 | LPA, class-varying variances (VVI) | -3174.564 (17) | kit: `lpa(model = "VVI")` |
| 7.11 | LCA with binary, censored, nominal and count indicators | -4348.434 (15) | needs: censored indicators (binary, nominal and count exist in `multilpa()`) |
| 7.12 | LCA with a covariate and a direct effect on a binary indicator | -1255.396 (11) | needs: `direct_effects` on categorical indicators (today continuous only) |
| 7.13 | LCA with parameter constraints | -2208.165 (5) | needs: equality/fixed-parameter constraints in `lca()` |
| 7.14 | Confirmatory LCA, two categorical latent variables | -9141.129 (29) | needs: several latent class variables in one model |
| 7.15 | Loglinear model, three-way table | -1380.635 (5) | out of scope |
| 7.16 | LCA with a local dependence in one class | -2598.111 (10) | needs: class-specific `local_dependence` (today in every class) |
| 7.17 | Mixture CFA | -2643.223 (17) | needs: factor mixture models |
| 7.18 | LCA with a second-order factor (twins) | -20142.956 (15) | out of scope |
| 7.19 | SEM, categorical latent variable on a continuous one | -5332.536 (18) | out of scope |
| 7.20 | SEM mixture | -7076.192 (23) | needs: factor/SEM mixtures |
| 7.21 | Mixture with known classes (multiple groups) | -6273.830 (19) | needs: a `known_class` argument |
| 7.22 | LPA, indicators correlated within class (EEE, 3 classes) | -3054.963 (24) | kit: `lpa(model = "EEE")` |
| 7.23 | (not in `chapter7.zip`) | — | — |
| 7.24 | CACE with missing data on the class indicator | -920.329 (8) | out of scope |
| 7.25 | Zero-inflated Poisson as a two-class model | -724.663 (6) | needs: a class with a fixed (zero) count |
| 7.26 | CFA, non-parametric non-normal factor | -2183.081 (18) | needs: factor mixture models |
| 7.27 | Factor (IRT) mixture | -8193.989 (33) | needs: factor mixture models |
| 7.28 / 7.29 | Twin models | -9927.558 (4) / -24258.023 (11) | out of scope |
| 7.30 | Continuous-time survival (Cox) | -98.989 (3) | out of scope |

## Chapter 8: mixtures with longitudinal data

| Example | Model | Mplus LL (p) | Status |
|---|---|---|---|
| 8.1 / 8.2 | GMM, continuous outcome, covariate on growth and class | -3183.703 (15) / -3337.654 (15) | needs: occasion-specific residual variances in `mixture_regression(random =)` (Mplus's default) |
| 8.3 | GMM, censored outcome | -3211.438 (15) | needs: censored family |
| 8.4 | GMM, categorical outcome (random effects by integration) | -3463.892 (14) | needs: random effects for binomial growth mixtures |
| 8.5 | GMM, zero-inflated Poisson / negative binomial outcome | -3881.206 (14) / -4422.166 (19) | needs: zero-inflated family; random effects for count growth mixtures |
| 8.6 | GMM with a categorical distal outcome | -3693.340 (17) | needs: occasion-specific residual variances (as 8.1) |
| 8.7 | Sequential process GMM, two latent class variables | -6606.707 (33) | needs: several latent class variables |
| 8.8 | GMM with known classes | -7051.578 (21) | needs: `known_class` |
| 8.9 | LCGA, binary outcome | -1272.307 (5) | kit |
| 8.10 | LCGA, three-category outcome | -1706.048 (6) | needs: thresholds shared by the classes in the ordinal family (latents: class-specific, 7 parameters) |
| 8.11 | LCGA, zero-inflated Poisson | -3107.216 (7) | needs: zero-inflated family |
| 8.12 | Hidden Markov, four occasions | -13427.811 (5) | kit |
| 8.13 | LTA, binary covariate as known class | -14596.240 (32) | needs: `known_class` (Mplus models the group's proportion); comparable by adding its marginal log likelihood |
| 8.14 | LTA, continuous covariate on states and transitions | -12902.140 (31) | kit |
| 8.15 | Mover-stayer LTA, three occasions | -19504.367 (32) | kit |
| 8.16 | Discrete-time survival mixture | (input only, no output posted) | out of scope |
| 8.17 | Continuous-time survival mixture (Cox) | -3025.792 (15) | out of scope |

## Chapter 10: multilevel mixtures

| Example | Model | Mplus LL (p) | Status |
|---|---|---|---|
| 10.1 | Two-level mixture regression, random intercept | -1752.471 (13) | needs: a between-level random intercept in `mixture_regression()` |
| 10.2 / 10.3 | Two-level mixture regression, between-level class | -1574.798 (11) / -1944.220 (21) | needs: a between-level random intercept with group classes |
| 10.4 | Two-level CFA mixture | -5334.228 (23) | needs: factor mixture models |
| 10.5 | Two-level IRT mixture | -9004.551 (43) | needs: factor mixture models |
| 10.6 | Two-level LCA, random class logits, covariates | -3888.076 (25) | needs: continuous random effects on class logits |
| 10.7 | Two-level LCA, between-level class | -98736.127 (99) | elsewhere: `harness/suites/mplus-twolevel-lca.R` (`multilca()`) |
| 10.8 / 10.9 / 10.10 | Two-level growth / GMM (three-level) | -13129.013 (18) / -5919.509 (20) / -11854.500 (23) | needs: three-level growth mixtures |
| 10.11 | Two-level LCGA, three-category outcome | -3317.328 (7) | needs: shared ordinal thresholds (as 8.10) and two levels |
| 10.12 | Two-level LTA, random effects | -10250.454 (18) | needs: continuous random effects on transition logits |
| 10.13 | Two-level LTA, between-level class | -20716.373 (18) | needs: group-class effects on part of the initial and transition logits (`lta(n_group_classes =)` frees them all) |

## Not yet surveyed

Chapters 11 (missing data), 12 (Monte Carlo) and 13 (special features)
hold a few more mixture examples; survey them when a matching feature is
added.
