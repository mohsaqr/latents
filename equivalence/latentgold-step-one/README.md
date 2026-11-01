# Step-one uncertainty in `three_step()` against Latent GOLD 6.1

`three_step(step_one = "estimated")` adds the sampling variance of the
step-one estimates to the step-three covariance: V3 + G V1 G', with V1 the
step-one covariance and G the derivative of the step-three estimates with
respect to the step-one parameters, the modal assignments held fixed (the
first-order correction of Bakk, Oberski and Vermunt 2014, *Political
Analysis* 22, 520-540). This kit checks it against Latent GOLD's
implementation of the same correction (Syntax Guide 4.3.12, "Reading fixed
values as well as their variance-covariance matrix from a text file").

Run from the project root:

```sh
Rscript equivalence/latentgold-step-one/make-kit.R
sh equivalence/latentgold-step-one/run-kit.sh      # Latent GOLD 6.1 under Wine
Rscript equivalence/latentgold-step-one/compare.R
```

`equivalence/test-latentgold-step-one.R` reruns the comparison offline from
`returned/`. Last run 2026-11-01: 156 rows agree, 2 not comparable, none differ.

## Cases

| case | step one | n | classes | indicators |
|---|---|---|---|---|
| s01 | `lpa(model = "VVI")` | 600 | 2 | 2 Gaussian |
| s02 | `lpa(model = "VVI")` | 900 | 3 | 3 Gaussian |
| s03 | `lca()` | 1000 | 3 | 6 binary |

Each has a numeric outcome `z` and a three-category outcome `cat`. The
step-one fits are converged to `tol = 1e-12`; at the default tolerance a
3-class fit stops 2-4e-5 below the maximum and its corrected standard errors
differ from the converged ones by about 1e-4 (0.1%).

## How Latent GOLD is set up

Step three is written by hand rather than with the `step3` option, so that
the classification errors can be read as fixed parameters:

```
dependent W nominal, z continuous;          # W = latents' modal assignment
latent Cluster nominal K coding=first;
Cluster <- (s) 1;                           # size logits log(pi_k / pi_1)
W <- (e~err) 1 | Cluster;                   # log D[t, s] / D[t, t], s != t
z <- 1 + Cluster;  z | Cluster;
e s = 's01_fixed.txt';
```

The fixed values are those latents holds fixed (the error matrix and the
average posteriors of its own fit). The file's `Step1Jacobian` (the
derivative of e and s with respect to the step-one parameters) and
`Step1Hessian` (the step-one covariance) come from a likelihood written in
`common.R`, independently of the package, so the comparison tests latents'
whole step-one path. Its maximum equals latents' to 8 decimals in every case.

Three things learned about Latent GOLD's file, from runs on s01:

- **`Step1Variances` is required for any correction.** With parameters,
  Jacobian and Hessian only (the manual calls Variances optional), the
  standard errors equal those of plain fixed values. The kit writes the same
  covariance under both keywords.
- **`Step1Gradients` (per-case gradients) switches to a second-order
  correction** that also accounts for the covariance between the steps. Here
  it moved the standard errors by up to ±3% relative to first order, in both
  directions (variant `second` in the kit, reported but not compared:
  latents implements first order).
- The `~err` logits in Latent GOLD's own `writeclassificationerrors` file are
  ordered by true class, then by assigned class without the true one.

## What is compared

- **Latent GOLD's own step one** (`writeclassificationerrors`): the error
  logits and their covariance J H J' against the independent likelihood's,
  after matching class labels. Invariant to the step-one parameterization.
  s01 and s03 agree (covariance to 3e-8 for s01). For **s02 Latent GOLD's
  search ends at another maximum** (LL -4277.02 against -4278.50: a class of
  0.5% of cases with a variance of 0.0045, a degenerate maximum of the
  heteroscedastic Gaussian mixture), so its file is reported as not
  comparable.
- **Step three**, `method = "ml"` with `vcov_type = "independent"` against
  Latent GOLD's robust standard errors: the class means or category
  probabilities, their standard errors and the overall Wald test, with
  `step_one = "fixed"` against the plain values and `step_one = "estimated"`
  against the file with its covariance. Tolerance: half a unit in Latent
  GOLD's printed place plus 1e-4 of the value. The correction raises every
  standard error, by a median factor of about 1.15-1.3 per case.

## Not covered here

BCH and proportional weighting have no Latent GOLD counterpart for this
correction; `tests/testthat/test-three-step-step-one.R` checks them against
an independent implementation. Group-level analyses (`level = "groups"`) and
the cluster-robust variance are covered by the unit tests only.
