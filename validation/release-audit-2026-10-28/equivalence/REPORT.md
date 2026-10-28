# Independent equivalence release audit — 2026-10-28

Reviewed the 0.98.0 working tree at e907bff. No production source was changed. All fresh R checks used `pkgload::load_all('.')`, so results concern current source rather than an older installed package. BLAS/OpenMP thread counts were one.

## Release blocker: shared Gaussian and negative-binomial dispersions disappear

**P1 — new regression introduced by e907bff.** `R/mixture-regression-starts.R:395` now names `params$sigma2` with class names. In `R/mixture-regression-inference.R:31` and `:34`, `c(log_sigma.all = named_scalar)` / `c(log_dispersion.all = named_scalar)` consequently produce names such as `log_sigma.all.class_1`, rather than the public parameter name expected downstream. `R/mixture-regression-methods.R:108-125` still indexes `log_sigma.all` / `log_dispersion.all` and returns missing values.

Affected: `mixture_regression(..., variance = 'equal')` for Gaussian and negative-binomial families, including **one class**, with or without shared regression coefficients. `coef()`/`vcov()` names change unexpectedly; the reported class SD/dispersion and their standard errors become `NA` even for a converged, interior model. The likelihood and fitted coefficient values in the probes are unchanged.

Evidence:

- `probe-flexmix.R` / `.log`: the existing Gaussian common-coefficient, equal-variance flexmix fixture converges with sigma² = 0.7130096; class SDs should be about 0.844399 but both are missing. The harness errors `every package estimate must be found in its tables`, reproducing the current CI failure.
- `probe-gaussian-one.R` / `.log`: even one-class Gaussian regression loses its residual SD (expected 0.9290562).
- `probe-equal-dispersion.R` / `.log`: one-class negative-binomial fit loses dispersion 0.2821447 and SE 0.06338804. Both reappear when an in-memory packer wrapper removes incidental names, with identical coefficient values.
- `probe-flexmix-patched.R` / `flexmix-patched.log`: the **whole flexmix suite passes 433/433 quantities and 29,933 values** after that in-memory-only correction, with no uncovered outputs or warnings. This is causal diagnostic evidence, **not a pass of unchanged release source**. Comparison/suite CSVs are retained beside the log.

Recommended repair: preserve the intended shared parameter names by using an unnamed scalar (`[[1L]]` or `unname`) in both shared-dispersion pack branches, then add regression tests for finite public class SD/dispersion and stable coef/vcov names for one and multiple classes. Rerun flexmix unpatched and the package suite on repaired source.

## Equivalence evidence and harness review

- The fresh wider Latent GOLD test passes: **172/172 parameter, standard-error, likelihood and parameter-count comparisons**, represented by four testthat assertions (`latentgold-wider.log`). It covers varying Gaussian variances, so this pass does not contradict the shared-variance blocker.
- `probe-lg-provenance.R` / `.log`: all 16 retained input fingerprints and 16 output fingerprints match, and all eight `cases.rds` data frames agree with the corresponding hashed `.dat` inputs. A deliberately injected 0.01 parameter error is rejected by `compare_case`; the 5.2e-5 parameter/SE tolerance is explicit.
- Harness self-tests completed successfully (`self-tests.log`). Reviewed numeric/text comparisons, NA/Inf handling, row joining, coverage/exemption accounting, known-open behavior, source fingerprinting and required-suite behavior. The known-open register is empty. No additional confirmed comparator false-positive defect was found.
- Latest CI log supplied by the parent records 12,035 agreeing quantities (3,302,066 values) from the completed suites and a failed flexmix suite. Those numbers do **not** mean a passing complete run. CI also skips its later test steps after the failed harness step.
- The tracked `equivalence/harness/REPORT.md` belongs to **0.97.2 and 15 suites**, whereas the current registry discovers 24 suites. It correctly labels its historical provenance, but must not be presented as 0.98.0 release evidence.

## Coverage gaps, not demonstrated numerical defects

1. `.github/workflows/equivalence.yaml:80-90` runs `tests/equivalence` and a narrow filter in `equivalence/run.R`; neither invokes the new `equivalence/test-latentgold-wider.R`, and no registered harness suite calls it. Add this test to CI if its 172 comparisons are to be a continuing release guarantee. Its fresh local pass is recorded above.
2. `equivalence/test-latentgold-wider.R:18-28` verifies `.dat`/`.lgs` input fingerprints but reads fitting inputs from un-fingerprinted `cases.rds`; it does not enforce RDS/DAT parity or `OUTPUTS.md5`. The retained artifacts currently agree, as independently checked above. Add these checks to prevent unnoticed future fixture drift.
3. A printed-reference comparison at 5.2e-5 establishes agreement within the declared rounding/convergence allowance, not full-double equivalence. Positive sampling-weight fixtures do not establish zero-weight semantics; those need invariant or separate reference tests. Existing coverage documentation already names remaining family/SE gaps.

## Fresh rerun commands

```sh
Rscript equivalence/harness/test-harness.R
Rscript equivalence/run.R '^latentgold-wider$'
LATENTS_EQUIVALENCE_OUTPUT=validation/release-audit-2026-10-28/equivalence/harness \
LATENTS_EQUIVALENCE_REQUIRE_ALL=true Rscript equivalence/harness/run.R \
  twolevel-lpa covariates lta-general lta-mixed ordinal-count three-step \
  predict-newdata robust-se pool-imputations
Rscript validation/release-audit-2026-10-28/equivalence/probe-flexmix.R
Rscript validation/release-audit-2026-10-28/equivalence/probe-equal-dispersion.R
Rscript validation/release-audit-2026-10-28/equivalence/probe-gaussian-one.R
Rscript validation/release-audit-2026-10-28/equivalence/probe-flexmix-patched.R
Rscript validation/release-audit-2026-10-28/equivalence/probe-lg-provenance.R
```

All commands were run with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1`. The nine-suite rerun uses unchanged source; its final status is recorded below when complete. No new Latent GOLD executable run was needed: comparisons use retained genuine listings. Other 15 suites were not rerun locally in this task; current CI evidence supplements the local checks. This audit establishes checked configurations, not exhaustive correctness across every option combination or platform.

## Final fresh nine-suite result

**Passed on unchanged 0.98.0 source: 3,715/3,715 quantities, 298,945 individual values, all nine requested suites ran, zero skips, zero warnings and zero uncovered outputs.** Exit status 0. Source fingerprints were unchanged throughout. The detailed generated evidence is in `harness/REPORT.md`, `harness/report.csv`, `harness/suites.csv`, `harness/coverage.csv` and `harness/SESSION.txt`.

The suites were twolevel-lpa (1,231), covariates (577), lta-general (102), lta-mixed (358), ordinal-count (486), three-step (271), predict-newdata (57), robust-se (541), and pool-imputations (92). This was a partial 9-of-24 run by design. The fresh harness self-test log records 42 successful blocks / 201 expectations. The release remains blocked by the independently reproduced shared-dispersion regression.
