# Fit, initialization and spike-cluster review — 2026-10-22

Reviewed HEAD `12e0145` plus the previous session's uncommitted corrections.
Production code and package tests were unchanged during this review. The nine
previously fixed source/test files match the preceding validation manifest.

## Confirmed production findings

### 1. P1 — General full-covariance LTA misses an active eigenvalue bound

Location: `R/lta-engine.R:848-850`; downstream gate:
`R/lta-inference.R:473-477`.

The result assembler checks covariance **diagonals**, although a full covariance
is constrained by its smallest **eigenvalue**. A nearly collinear two-indicator
VVV model with `transitions="occasion"`, 60 subjects and three occasions,
`min_variance=.01`, converges with both covariance minimum eigenvalues exactly
`.01`. It reports `boundary=FALSE`, emits no boundary warning, and
`parameter_inference()` returns five finite Wald SEs and p-values. Changing only
the saved boundary flag to TRUE causes the intended `latents_boundary_fit`
refusal. Three independently generated noise levels reproduce this behavior.

This bypasses an existing inference restriction. It does not establish that all
returned numerical SEs are wrong; they are returned for a fit the documented
guard is supposed to refuse. LTA does not offer `boundary="exclude"`.

Correction direction: use the existing covariance eigenvalue boundary helper
on each measurement block when assembling the general LTA result.

Evidence: `lta-boundary-wald-probe.R/.log/.rds` and
`spike-diagnostic-probes.R/.log/.rds`.

### 2. P2 — Zero-weight rows control every default initialization

Location: `R/multilpa.R:1702-1706`; unweighted internal centering at line 1645;
unweighted initializer statistics and geometry at lines 181-253.

For 100 active observations from two well-separated Gaussian components, append
one zero-weight row. Changing only that excluded observation from 0 to 10000,
with the same 101-row weighted objective, changes the default reported LL from
`-152.6178709085` to `-256.4855779545`. All 10 starts, and all 100 starts in a
separate run, converge to the lower solution. Both fitted means are near zero
instead of approximately -3.026 and 3.006. `boundary=FALSE`; small/thin-profile
warnings do appear. A supplied good start and emEM recover the higher maximum
exactly on the same data.

Ward and k-means spend a component on the excluded outlier because initialization
never receives sampling weights. At a zero-weight value of `1e20`, the unweighted
centering also destroys precision in the active rows. This is an excluded-data
failure, independently of heteroscedastic spike selection.

Correction direction: compute internal origins and initial clustering from
positive-weight support, retaining full row/group mapping for subsequent EM.
Positive non-unit weights and mixed/missing indicators require explicit tests;
this finding does not prescribe a new weighted k-means algorithm.

Evidence: `zero-weight-starts.R/.log/.rds`,
`initialization-probes.R/.log`, `initialization-findings.txt`.

### 3. P2 — Non-Gaussian rows hide thin Gaussian profiles outside LTA

Locations: `R/multilpa.R:1937-1941`, `R/covariates.R:495-499`.

These paths use total classification memberships for a diagnostic about Gaussian
measurement support. A five-indicator profile with three Gaussian members
correctly warns. Appending 300 entirely empty rows with `drop_empty=FALSE`
preserves LL `-576.0075` and Gaussian memberships 3/77, but inflates classification
counts to 14.25/365.75 and removes the warning. Categorical-only rows also hide
the warning under default `drop_empty=TRUE`: ordinary counts become 303/77;
membership-covariate counts become 302.99962/77.00038, with actual Gaussian
memberships still 3/77.

Correction direction: apply the same Gaussian-observed row mask and weighted
E-step support already used in the corrected LTA paths. Preserve classification
counts, estimates and missing-data policy.

Evidence: `spike-diagnostic-probes.*`, `spike-cov-crosslevel-probes.*`.

### 4. P2 — Full cross-level thinness ignores sampling weights

Location: `R/cross-level.R:385-387`; weighted posteriors are replaced by reported
probabilities at lines 204-214 before this check.

A converged two-indicator fit with three low-weight rare groups has raw
individual memberships 6/74 but weighted memberships
0.06481231/79.93518769. Four free Gaussian parameters are estimated per profile.
No thin-profile warning appears because the check uses raw probabilities.

Correction direction: retain or reconstruct the membership weights used by
estimation for this diagnostic, consistently with other fitting paths.

Evidence: `spike-cov-crosslevel-probes.*`.

### 5. P2 — Full cross-level thinness omits the group-level Gaussian mixture

Location: `R/cross-level.R:385-387` checks only individual profiles and never
uses `between_variance` or group memberships for this warning.

With one individual profile, two group classes, two Gaussian indicators and
`between_variance="varying"`, a converged full model has group memberships 3/37
and no warning. The restricted model produces the same group means and
variances and correctly warns: three members support four Gaussian parameters.

Correction direction: assess between-group Gaussian densities separately with
group memberships and `between_variance`; identify the level in the message.

Evidence: `spike-gaps-probes.*`.

## Additional search and evidence issues

### P3 — One-final-start emEM always finalizes Ward's start

`R/search.R:61` protects start 1 unconditionally by default; the finalist builder
at `R/kernel-em.R:406-409` spends a final slot on it before ranking others. With
`n_starts=1` and a valid Ward start, screening cannot change the returned fit.
This protection is intentional and tested, but the public description that the
"best n_starts" are finalized omits the exception.

Actual public fit: `lpa(srl, names(srl), 4, model="EEE", n_starts=1,
search="emEM", n_screen=10, seed=1)` reports LL `-584.9543` from Ward. Screened
start 9 already reaches `-582.1427`. Disabling protection only in a temporary
process selects the identical pool's start 9 and converges at `-582.1419`.
Screening scores are identical between the two runs. No production mutation
was made. Resolve the one-start behavior or document the protected-slot rule.

Evidence: `root-emem-protection.R/.log/.rds`.

### P3 — Benchmark analysis classifies unmeasured finalists as non-spikes

`equivalence/search/pilot.R:85` fills missing finalist variance metrics with NA.
All 1,400 Gaussian finalists from the final no-rule run lack that metric.
`equivalence/search/analyse-pilot.R:12,22` treats NA as not a spike and admits
these unclassified candidates into its non-spike reference.

For VVI replicate 1, k4, the reported winner has LL `-1909.770` and measured
standardized variance ratio `0.005206615`, below the benchmark-only .01 cutoff.
Its duplicate finalist has ratio NA and is incorrectly admitted as non-spike.
Do not treat a missing Gaussian diagnostic as an observed interior solution;
distinguish genuinely inapplicable categorical variance metrics.

The selected-fit `compare-norule.R` analysis avoids this error. A corrected
recomputation confirms 136/180 non-spike reference attainments for emEM versus
130/180 for kmeans; 25 versus 12 selected spikes; exclusive attainments 10
versus 4. All 1,440 retained reported fits and 18,000 retained finalists
converged. These are finite-reference results under a .01 benchmark-only
criterion, not a global-optimum or universal spike-detection guarantee.

`NEWS.md:27-30` should name that criterion and qualify "best known maximum".
Its runtime claim of 3-5x also appears stale: final no-rule records have median
1.0265 versus .515 seconds (~2.0x), total 207.45 versus 143.536 seconds (~1.45x),
and family median ratios about 1.3-2.55. Earlier with-rule runs were slower;
these timings describe the retained runs, not a portable speed guarantee.

Evidence: archived `benchmark-inputs.rds`, corrected `recompute.R/.log`,
`recomputed.csv`. Run the script from the repository root.

## Known limitations reconfirmed, not new implementation defects

- A moderate-size flat component can escape every current warning. The known
  srl VVV k4 default fit (seed 1) converges at LL `-502.0884`, minimum membership
  23.67246 and standardized minimum eigenvalue .001206186, with `boundary=FALSE`
  and no warning. A floor check and a membership/parameter-count check cannot
  detect every spike. A new cutoff would be a policy change the user has
  already rejected. The benchmark .01 criterion is not a production filter.
- Many agreeing converged k-means starts can reproduce a lower local maximum.
  `equivalence/features/direct-effects/FINDINGS.md` documents heavy_02 remaining
  8.50 LL low across 600 starts; `equivalence/features/vlmr/README.md:118`
  documents nine EEI/EEE candidates 0.15-4.85 LL low. Those studies were read,
  not rerun here; they predate opt-in emEM and do not establish its current
  behavior on every case. Best-start replication is not proof of optimality.
- A larger screened pool can discard an eventual winner; emEM is not guaranteed
  to dominate default search on each dataset. Previous evidence records a .93
  LL loss on a bounded EEE example. No new selection filter was introduced.
- The old SEARCH-PLAN.md and some study comments describe abandoned plans.
  The final policy is kmeans default, emEM opt-in, no handmade spike filter.
  Final no-rule emEM BLRT calibration has not been measured.

## Validation and reproduction

Six light probe scripts were rerun serially with one BLAS thread against the
current working source via `pkgload::load_all()`. All completed successfully;
the LTA and one-start emEM probes include explicit assertions. The zero-weight
probe includes all 100 starts and good-start/emEM comparisons. Benchmark
summaries were recomputed from archived input tables without rerunning fits.

`test-review-regressions.R` expresses corrected behavior for the five production
findings. It is outside the package test suite and intentionally fails on this
snapshot: five cases, six failed expectations, seven supporting passes,
zero errors, zero warnings. Reproduce with:

```sh
Rscript validation/review-fit-search-2026-10-22/run-regressions.R
```

The runner checks that expected review failure pattern; exit 0 means the
recorded defects reproduced, not that production regressions are fixed.
Logs, CSV and result objects are retained here. The first CSV export attempt
could not serialize testthat's list column; the saved results were valid, and
the corrected exporter drops that column and was rerun successfully.

No complete unit suite, complete external feature suite, build or CRAN check
was repeated for this read-only review. The preceding successful 1320-unit-test
and 64-equivalence-test run is documented in `../fixes-2026-10-22/`; its
`tests/equivalence/` coverage does not include the larger top-level
`equivalence/` feature studies with documented failures.
