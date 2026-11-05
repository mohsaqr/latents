# Transition parallel test diagnosis

The isolated `test-lta-parallel.R` failures arise from assumptions about the reported processor count in the tests. The sandbox returns `NA` from `parallel::detectCores()`, so the intended 1.1.13 fallback uses one core and signals `latents_unknown_cores`. The fitted-object equality checks pass. A controlled rerun that changes only the reported processor count in memory passes every executed expectation.

| Run | Passed expectations | Failed expectations | Errors | Warnings | Skipped tests |
|---|---:|---:|---:|---:|---:|
| Actual sandbox processor report | 79 | 16 | 0 | 2 | 2 |
| In-memory reported count of 8 | 99 | 0 | 0 | 0 | 2 |

Both runs use `/private/tmp/latents-review-1105-pr7czbhm/source`, its matching installed package in the adjacent `library` directory, `NOT_CRAN=true`, `_R_CHECK_LIMIT_CORES_=true`, and single-threaded BLAS environment settings. The second run replaces only `.latents_available_cores()` with `function() 8L` in that R process. It leaves the Accelerate fork exclusion and CRAN core limit active. The test paths that remain enabled use serial mapper mocks; the reported count is not permission to create eight physical workers.

The two skipped tests exercise real forked starts and larger transition grids. The first skips because the active BLAS is Apple Accelerate, for which `.latents_fork_safe()` returns `FALSE`; the second skips because of the explicit core-limit environment setting. Therefore this control validates the serialized task and scheduling contracts, not live parallel execution.

## Actual failure mechanisms

Eleven failures compare condition lists. Calls requesting two or three cores receive the additional correct unknown-core warning, while serial reference calls do not. These occur at `tests/testthat/test-lta-parallel.R:171`, `:262`, `:313`, `:346` and `:404` across the parameterized settings.

Five failures inspect task-phase counts and their costs in the test beginning at line 195. Clamping the request to one core removes the parallel candidate-setup phase. The test expects task counts `4,8` or `4,16,8`, but receives `8` or `16,8`, and subsequent cost expectations refer to the wrong phase. All corrected task-phase expectations pass when the reported core count is made deterministic in the controlled run.

The two natural-run warnings come from the lost-worker propagation test, which does not capture the additional unknown-core warning. Neither run has an R error. The assertions comparing numerical fits, likelihoods and returned objects pass in the natural environment as well as the control.

## Recommended test change

Mock the reported processor count in tests whose purpose is serialized closure contents, phase scheduling or condition replay under a specified worker budget. Continue to test unknown-core fallback separately. Keep live worker tests conditional on the actual environment. The existing fallback should not be removed to satisfy these scheduling-test assumptions.

## Related tests inspected by reading

`test-parallel-profile-pool.R` has the same environment-sensitive assumptions: its phase-count assertion at line 126 expects the parallel setup phase, its serialized-profile test compares condition lists and expects seven mapper calls, and its membership-covariate test expects a nonserial candidate mapper call. All request three cores without fixing the reported processor count. `test-parallel.R:82` also compares serial and two-core condition lists without allowing the unknown-core warning. These are plausible explanations for corresponding broad-suite failures, not independently rerun diagnoses of those files.

The intended unknown-core behavior already has direct tests in `test-review-1026-refusals.R`, under the test named “an unknown processor count runs on one core with a classed warning.” Those pass explicit `available=NA_integer_` and `available=0L` values and do not depend on the machine report.

No source files or existing tests were changed. Evidence is retained in [lta-parallel-check.log](lta-parallel-check.log), [lta-parallel-summary.log](lta-parallel-summary.log), [lta-parallel-mocked-check.log](lta-parallel-mocked-check.log), and [lta-parallel-mocked-summary.log](lta-parallel-mocked-summary.log), with structured test results in the adjacent RDS files.
