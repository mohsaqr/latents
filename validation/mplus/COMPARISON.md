# Mplus comparison — 2026-09-17

The native R `multilpa` estimator agrees with Mplus in six retained comparisons
across three datasets. Two comparisons reproduce **published** single-level
results; four compare against **new, actual Mplus runs** of our exact two-level
model. The fitting engine was unchanged.

## Published results reproduced

Official Mplus 8.8 User's Guide outputs from 2022-04-19:

| Model | Observations / indicators | Mplus log likelihood (published) | R log likelihood | Absolute difference |
|---|---|---:|---:|---:|
| Example 7.9, equal variances | 500 / 4 | -3177.162 | -3177.161923046 | 0.000076954 |
| Example 7.10, varying variances | 500 / 4 | -3174.564 | -3174.563712388 | 0.000287612 |

These examples use **identical data**, confirmed by MD5
`dd304e84be81745eabd9669fc859a432`; they are two models, not two independent
datasets. The supplied generating-class column is excluded from fitting.

Means, variances, AIC, BIC, and entropy match within half of the published
0.001 unit. Class proportions match within 0.000005. Modal class counts match
exactly. Effective membership counts require an explicit 0.001-person
optimization allowance; this is not claimed to be just printed rounding.
Both R fits used 20 independent starts and replicated their best likelihood.
These comparisons validate only the single-level limit.

Original sources:

- Example 7.9: [output](https://www.statmodel.com/usersguide/chap7/ex7.9.html),
  [data](https://www.statmodel.com/usersguide/chap7/ex7.9.dat),
  [input](https://www.statmodel.com/usersguide/chap7/ex7.9.inp).
- Example 7.10: [output](https://www.statmodel.com/usersguide/chap7/ex7.10.html),
  [data](https://www.statmodel.com/usersguide/chap7/ex7.10.dat),
  [input](https://www.statmodel.com/usersguide/chap7/ex7.10.inp).

The downloaded originals are retained in `public/`. Parsed references and
source hashes are in `tests/fixtures/mplus/public-*.rds` relative to the project
root.

## Direct two-level comparisons

Mplus VERSION 9 DEMO (Mac) was actually executed locally. Every model has two
individual profiles and two discrete group classes. Each Mplus analysis used
100 initial starts and 20 final optimizations, with independently initialized
R estimation. Each Mplus run terminated normally and replicated its best
likelihood. No R-fitted values were supplied as Mplus starting estimates.

| Data and variance constraint | Individuals / groups / indicators | Parameters | Mplus saved log likelihood | R log likelihood | Absolute difference |
|---|---|---:|---:|---:|---:|
| Project simulation, varying | 1200 / 60 / 2 | 11 | -3040.5252 | -3040.525244426 | 0.000044426 |
| Project simulation, equal | 1200 / 60 / 2 | 9 | -3041.9169 | -3041.916897387 | 0.000002613 |
| Public Example 10.4 data, varying | 1000 / 110 / 5 | 23 | -6968.1494 | -6968.149408163 | 0.000008163 |
| Public Example 10.4 data, equal | 1000 / 110 / 5 | 18 | -6968.8336 | -6968.833585853 | 0.000014147 |

The Mplus saved-results file uses eight significant digits. Likelihood
differences above are within that output precision; printed `.out` likelihoods
have only three decimals. Posterior probabilities were separately saved with
`FORMAT = F20.12`.

Maximum absolute differences after aligning class labels and subject IDs:

| Quantity | Simulation, varying | Simulation, equal | Public data, varying | Public data, equal |
|---|---:|---:|---:|---:|
| Profile means | 5.05e-8 | 3.01e-8 | 1.85e-6 | 1.78e-6 |
| Profile variances | 4.81e-8 | 1.32e-8 | 1.44e-6 | 1.91e-7 |
| Profile proportions within group classes | 5.12e-9 | 3.31e-9 | 3.84e-6 | 1.12e-6 |
| Group-class proportions | 2.56e-10 | 9.26e-12 | 1.02e-5 | 3.44e-6 |
| Individual posterior probabilities | 1.16e-10 | 1.69e-10 | 1.15e-5 | 7.37e-6 |
| Group posterior probabilities | 5.09e-12 | 7.88e-13 | 2.95e-5 | 1.06e-5 |

The public dataset has less decisive group memberships than our deliberately
separated simulation. Its differences reflect iterative stopping criteria as
well as output rounding. The retained test tolerance is **1e-4 absolute** for
its parameters/posteriors, not an assertion of machine-precision equivalence.
Native R starts/seeds differ between the comparison scripts and regression
tests, so last digits can vary while all checks continue to pass. Both public
data models give identical marginal modal classifications in R and Mplus.

### Public data model distinction

The [original Example 10.4](https://www.statmodel.com/usersguide/chap10/ex10.4.html)
fits a **two-level CFA mixture**, with continuous factors and a different
likelihood. Our package cannot reproduce that specification. Its original
reported log likelihood, -5334.228, is deliberately **not** our comparison target.

We downloaded the [original data](https://www.statmodel.com/usersguide/chap10/ex10.4.dat)
and [input](https://www.statmodel.com/usersguide/chap10/ex10.4.inp), preserved all
five indicator values and the 110 observed group IDs, excluded the generating
class label, and fitted a **new identical LPA specification in both programs**.
The new Mplus input/output/results are retained in `public-refit/`. This is
public simulation/example data, not a newly collected empirical dataset.

### Exact model matching

The input files and TECH1 parameter maps confirm:

- `WITHIN = y1-yD` assigns the Gaussian indicators to the individual level.
- `CLASSES = cb(2) cw(2)` and `BETWEEN = cb` define both discrete levels.
- `cw#1 ON cb` lets the individual profile proportions differ across group classes.
- `MODEL cw` estimates profile means and diagonal variances shared across group
  classes. Reused variance labels impose equality across individual profiles
  in the equal-variance variants.
- There are no continuous random effects, covariates, missing values, or
  residual indicator correlations. Parameter counts match exactly.

Mplus prints warnings that the indicators are uncorrelated within class;
these describe our intentional diagonal-covariance specification, not an
ignored fitting failure. The warnings are retained in the original outputs.
The R variance bound is inactive in every compared fit; these checks do not
establish equivalence when that constraint is active.

Mplus saves **joint** class-pattern posteriors. We sum across group classes for
individual marginals and across individual profiles for group marginals. We
align saved records on explicit subject IDs, since Mplus can reorder them.
We do not equate joint-modal labels with marginal-modal labels.

**BIC convention:** Mplus's BIC uses the individual count here. Comparisons use
`fit$bic_individual`; the package's default `fit$bic` uses the group count and
is intentionally different. Standard errors, MLR corrections, and confidence
intervals are not compared because the R package does not implement them.

## Reproduce and inspect

Run all commands from the project root unless stated otherwise:

```sh
Rscript validation/mplus/compare-public.R
Rscript validation/mplus/generated/compare-generated.R
Rscript validation/mplus/compare-public-refit.R
Rscript -e 'pkgload::load_all("."); testthat::test_dir("tests/testthat")'
```

The comparison scripts use retained data/results and do not require Mplus.
To regenerate the genuine Mplus outputs, run the installed `mpdemo` executable
with `varying.inp varying.out` and `equal.inp equal.out` from `generated/`, and
with `public-varying.inp public-varying.out` and
`public-equal.inp public-equal.out` from `public-refit/`. On this machine the
executable is `/Applications/MplusDemo/mpdemo`. Models satisfy its demo limits.
Regenerating Mplus files changes output timestamps and therefore provenance
hashes; rerun comparison scripts to regenerate the reference fixtures.

Actual estimates and comparisons are in each directory's `comparison.rds`.
Readable logs are `generated/comparison-log.txt` and
`public-refit/comparison.log`. Original `.inp`, `.out`, `*-results.dat`, and
`*-posteriors.dat` files remain available alongside them. Offline fixtures in
`tests/fixtures/mplus/` contain only genuine Mplus-derived targets, source
data, version information, and MD5 provenance.

## Scope of the conclusion

This establishes close numerical agreement with Mplus for the tested
two-profile/two-group-class Gaussian models, with either variance constraint.
It does not establish full Mplus feature parity, robust standard errors,
missing-data equivalence, high-class-count behavior, or reliable global
optimization under all weak-separation/rare-class conditions.

A search for an exact existing published two-level oracle also inspected
[Houle et al. (2026)](https://smslabstats.weebly.com/uploads/1/0/0/6/100647486/ml-lpa_final_-_prepub.pdf)
and [Collie et al. (2020)](https://www.frontiersin.org/journals/psychology/articles/10.3389/fpsyg.2020.00626/full).
The former's provided workflow fixes measurement parameters from a previous
fit; the latter uses factor scores and sampling weights. Those reported
analyses were not represented as exact matches to our unweighted, jointly
estimated model.
