# Ordinal and count indicators against Latent GOLD 6.1

Latent GOLD's `ordinal` dependent variable is an adjacent-category logit with
fixed equally spaced scores ("Ord-Fixed"), category intercepts and one class
effect, which is exactly latents' `ordinal =` model; `poisson` is latents'
`count =`. The kit fits five simulated datasets in both programs.

| case | model | latents LL | Latent GOLD LL | Npar | max parameter gap |
|---|---|---|---|---|---|
| o01 | 2 classes: continuous + ordinal (4) + count | -3772.866800 | -3772.8668 | 11 | 4.5e-05 |
| o02 | 3 classes: ordinal (3), ordinal (5), count | -4007.057264 | -4007.0573 | 15 | 4.9e-05 |
| o03 | 2 classes: two counts | -2212.099624 | -2212.0996 | 5 | 4.4e-05 |
| o04 | 2 classes: three ordinal (4, 3, 5) | -2649.637659 | -2649.6377 | 13 | 5.0e-05 |
| o05 | two-level, 2 group classes x 2 profiles: continuous + ordinal + count | -3446.194290 | -3446.1943 | 13 | 5.0e-05 |
| o06 | 2 classes: continuous + negative binomial, dispersion by class | -4837.397234 | -4837.3972 | 9 | 5.5e-05 |
| o07 | 2 classes: ordinal + negative binomial, shared dispersion | -3624.950104 | -3624.9501 | 8 | 9.0e-05 |

Latent GOLD's `poisson overdispersed` is NB2 (variance mu + sigma^2 mu^2),
with the dispersion by class (`k | Cluster;`) or shared (`k;`), which is
latents' `count_model = "negative_binomial"` with `count_dispersion =
"varying"` or `"equal"`. o07 is weakly separated: EM converges at a linear
rate, so latents is run at tol 1e-10 (1e-12 is not reached in 5000
iterations); its likelihood agrees to 4e-6.

Parameters compared: class sizes (single-level), every ordinal category
probability per class, count and continuous means per class, after matching
class labels. Latent GOLD prints four decimals, so 5e-05 is agreement.

An earlier o04 (two ordinal indicators, three classes) was saturated, df = 0,
with no unique maximum; Latent GOLD flagged estimation warnings and stopped
on the ridge 0.22 below latents. It was replaced by an identified case.

Run (from the package root): `make-kit.R`, the Latent GOLD batch runs (see the
memory note / `equivalence/latentgold-families/README.md` for the Wine
command), copy the `.lst` files to `returned/`, then `compare.R`;
`make-fixture.R` stores the results as `tests/testthat/fixtures/latentgold-ordinal.rds`.
