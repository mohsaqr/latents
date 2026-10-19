# Latent GOLD runs for the newer features

Latent GOLD 6.1 run under Wine on this machine (see `../latentgold/README.md`
for how). Each folder keeps the syntax, data and listing of the runs the
package's tests are pinned to, so every number in those tests can be traced to
Latent GOLD's own output.

## `vlmr/` — Vuong-Lo-Mendell-Rubin tests (`get_results(x, "lmr")`)

| Run | Model | Latent GOLD | latents | Test |
|---|---|---|---|---|
| `c09_vlmr` | `c09.dat`, 2 continuous items, class-specific variances, 1:3 classes | 1 vs 2: LR 143.4874, aLR 139.1372, reference mean 5.3059, SD 3.2909 | 143.48737, 139.13724, 5.3059102, 3.2909135 | checked by hand (p = 0.0000 printed) |
| `null2_equal` | `null2.dat`, 4 continuous items, one variance per item, 2:3 classes | LR 8.3027, aLR 8.0439, mean 4.1216, SD 2.8763, VLMR p 0.0848, aLMR p 0.0941 | 8.30273, 8.04386, 4.12162, 2.87633, 0.0847534, 0.0940538 | `tests/testthat/test-lmr.R` |

The 2-vs-3 comparisons with class-specific variances (`c09_vlmr`, and
`null2.lgs` before the variances were made equal) end at a 3-class solution
with active variance constraints, which is why the pinned case uses equal
variances. The adjustment `LR / (1 + 1 / (d log n))` is confirmed by Latent
GOLD's own ratio: 26.3304 / 25.5321 = 1 + 1 / (5 log 600).

## `local-dependence/` — joint items (`local_dependence = "a:b"`)

| Run | Model | Latent GOLD | latents |
|---|---|---|---|
| `ld_plain` | 4 binary items, 2 classes | LL -3645.2184, 9 parameters | -3645.2184, 9 |
| `ld_joint` | as `ld_plain` with `a <-> b \| Cluster` | LL -3614.1388, 11 parameters, BIC 7308.7229, sizes 0.5093 / 0.4907, P(a = yes) 0.6914 / 0.2184, P(b = yes) 0.6827 / 0.2750 | the same to every printed digit (`tests/testthat/test-local-dependence.R`) |

## `local-dependence-continuous/` — a covariance for a continuous pair (`local_dependence = "y1:y2"`)

| Run | Model | Latent GOLD | latents |
|---|---|---|---|
| `cov_pair` | 4 continuous items, 2 classes, class-specific variances, `y1 <-> y2 \| Cluster` | LL -7272.5017, 19 parameters, BIC 14679.7150, covariances -0.4362 (s.e. 0.0442) and 0.5968 (s.e. 0.0580) | the same to every printed digit, classes in the other order (`tests/testthat/test-local-dependence.R`) |

The ML three-step cases are part of the main kit (`../latentgold/`, case
`c09`).
