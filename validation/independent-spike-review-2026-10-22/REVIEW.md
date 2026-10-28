# Independent spike-cluster review — 2026-10-22

The proposed direction is mostly sound: fix the LTA boundary gate and excluded-data initialization, make the existing diagnostics consistent, retain kmeans as default and emEM as an option, and evaluate recovery before changing search policy. However, the present evidence does not establish that an above-cutoff fit is statistically regular, that higher likelihood improves recovery, or that mclust intrinsically prevents spikes. The action plan needs a statistical documentation correction and prominent continuous diagnostics; fixing the five reported inconsistencies alone leaves known flat components silent.

This review was read-only for the repository. I read the global instructions, current handoff/learnings, current source, archived analyses and regression scripts. I independently ran nine small mclust fits and two algebraic checks with one BLAS thread, serially. I did not rerun the five defect suites, large benchmarks, or full tests. New material is confined to this scratch directory. The preceding production fixes and git state were not changed.

## Prioritized findings

### 1. Confirmed: the full-covariance LTA boundary flag is a correctness defect; fix first

`R/lta-engine.R:848-850` inspects covariance diagonals, while a full-covariance constraint can be active in an oblique eigenvector with every diagonal well above its floor. `R/lta-inference.R:473-477` trusts that flag to refuse Wald inference. The saved probe reports minimum eigenvalues exactly .01, a false boundary flag, and finite inference. Source logic independently supports that explanation. This is an inference-contract failure, rather than evidence that every returned SE is numerically wrong.

Use the existing eigenvalue boundary helper on every estimated measurement block. Cover both invariant and occasion-specific measurements and interior controls. Review inference guards as well as warnings: changing a warning alone does not repair the refusal. The existing probe and intended-correctness regression are appropriate evidence; I did not rerun them.

### 2. Confirmed: zero-weight observations influence initialization and numerical origin; fix before search tuning

`R/multilpa.R:1645` centers using all rows; `1702-1706` passes no sampling weights to initialization; `181-253` uses all rows for imputation means, standardization, clustering and initial compositions. An excluded outlier can consume a component and destroy active-row precision. This is distinct from the legitimate possibility that stochastic initializations find different local optima.

The reported LL comparison correctly holds row count fixed: both versions contain 101 rows, including one zero-weight row whose value alone changes. Thus normalization to the number of units does not explain its approximately 104-LL gap. Comparisons against the 100-row original require the reported 101/100 rescaling; the existing script does account for that distinction.

Initialize geometry and origins from observed positive-weight support, preserve row/group mapping, and test zero-weight rows in mixed/missing-data and grouped paths. Excluding zero-weight rows does not, by itself, specify the best initialization for unequal positive weights. Do not turn this repair into an unbenchmarked wholesale change of weighted search geometry.

### 3. Confirmed statistical error: the thin-profile warning is a heuristic, not an identification test

`R/multilpa.R:2156-2159` says a profile with fewer effective members than its total density parameters cannot determine its variances and has a singular or nearly singular covariance. This is not generally true. A diagonal covariance can have positive variances in arbitrarily many dimensions from two distinct observations. For full covariance, affine support/rank matters, not the total number of means plus covariance parameters. With positive soft responsibilities, a component can retain full-rank covariance even when its summed responsibility is arbitrarily small.

The independent algebraic check gives a two-dimensional positive-definite covariance with eigenvalues 1/3 and 1/9 and responsibility mass .03, despite five Gaussian parameters. This is a direct counterexample to the documentation's implication, not proof that such a component is statistically trustworthy. Conversely, large membership does not prevent an almost planar component.

Keep a warning as an explicitly described conservative heuristic if desired, but correct its rationale and NEWS language (`NEWS.md:31-35`). Do not call its absence identification, regularity, or absence of a spike.

Normalized sampling-weight mass is also not an independent sample count. `R/weights.R:18-57` scales weights to the total number of units, including zero-weight units. Adding 90 zero-weight units to ten unit-weight observations changes normalized total mass from 10 to 100 while leaving the weighted objective's maximizers and Kish effective sample size (10) unchanged. A fixed mass-versus-parameter warning can therefore disappear solely through this normalization. The three proposed diagnostic fixes restore consistency with the package's existing convention, but do not resolve this statistical limitation. Label normalized mass explicitly; distinguish it from positive contributing units and weight concentration. Kish effective size is useful descriptive context, not a universal replacement identification test for soft mixture components.

### 4. The other three reported defects are credible implementation inconsistencies; their proposed fix needs narrower claims

Source confirms that ordinary/membership-covariate paths use overall classification counts (`R/multilpa.R:1936-1939`, `R/covariates.R:495-498`), full cross-level reporting replaces weighted counts with probabilities before its diagnostic (`R/cross-level.R:206-214,385-387`), and the full cross-level diagnostic omits the group Gaussian mixture. The saved examples plausibly reproduce each gap. Matching these paths to the existing Gaussian-observed weighted convention and assessing the between-group density separately is reasonable.

Call these fixes to coverage/consistency of an existing warning, not three proofs of unidentified Gaussian models. An any-Gaussian-observation mask is an improvement but still cannot establish covariance support under arbitrary missingness: indicator-specific observation counts and joint information patterns matter. A group-level warning should identify group classes and the density level, rather than labeling every case an individual profile.

## What the package comparison establishes

I inspected the original generating script, now copied as `original-mclust-spikes.R`. It uses the same two condition specifications, simulator calls, seed construction and indicator selection as `equivalence/search/pilot.R`; mclust uses its defaults and the shared `spike_measure()` helper on its covariance arrays. The original script also checks current latents-versus-pilot LL at replicate 1, generating K, for both conditions. The archived `mclust-benchmark-inputs.rds` contains summary rows, not datasets, fitted parameters, controls, package versions or warnings. `compare-packages.R` only verifies and joins those summaries. Archive the generating script and runtime metadata with the table; otherwise its original provenance is unnecessarily dependent on a volatile scratch directory.

My three selected cells independently reproduce archived default LL and ratios exactly, and a separate standardized-eigenvalue implementation agrees with the shared helper. This supports provenance and metric equivalence for those cells and supports the intended matching procedure. It is not a fresh verification of all 160 prior/no-prior rows. The nine runs did not rerun latents benchmarks.

The descriptive counts are supported: among 80 VVI/VVV fitted cells, default mclust has 1 below-.01 selected fit, latents kmeans 12, and emEM 25. There are only 40 independently generated datasets, each fitted at K and K+1. These are two generating conditions with 20 replicates each; the VVV fits are generated using diagonal within-component spreads. They do not establish rates across realistic correlated mixtures or other problem families.

**Do not combine the 180-cell attainment improvement with these 80-cell spike rates as if they measured the same population.** In these same 80 cells, both latents methods attain the above-.01 finite reference in **40/80** cells. emEM has 25 cutoff-small fits and 15 misses; kmeans has 12 and 28. At VVV K+1 specifically, emEM has 17/20 cutoff-small fits and kmeans 9/20; neither attains the above-cutoff reference in any of those 20 cells. The reported 136/180 versus 130/180 gain concerns the larger, mixed-family panel.

The .01 measure is a useful reproducible descriptor of a narrow covariance direction relative to marginal data variance. It is not a validated classifier of spurious components. A genuinely narrow or strongly separated generating class can satisfy it; a problematic component can lie above it. Marginal variance includes between-component separation. Standardization handles changes of individual indicator units, but does not turn the measure into a general affine-invariant or scientifically universal criterion. Use “below the benchmark's .01 cutoff” and “best observed above-cutoff reference,” not “proven spike” and “regular maximum.” Keep the cutoff strictly in benchmark analyses, consistent with user policy.

### Stopping is a demonstrated confound, not a hypothetical caveat

[mclust's official control documentation](https://mclust-org.github.io/mclust/reference/emControl.html) specifies outer LL tolerance 1e-5 and a numerical-singularity tolerance near machine precision. Current latents uses tolerance 1e-8 (`R/multilpa.R:1300`), SQUAREM by default, and an absolute variance/eigenvalue floor of 1e-6. These are different stopping/feasibility mechanisms. Matching model names does not match search, numerical path or feasible set. [mclust documents its model-based hierarchical initialization](https://mclust-org.github.io/mclust/reference/mclustBIC.html); its default tree uses VVV/SVD here, which is not latents' Ward start.

My small default-versus-tighter-control probe, holding data/model/hierarchical initialization fixed, found:

| Cell | LL at tol 1e-5 | LL at 1e-8 | LL at 1e-10 |
|---|---:|---:|---:|
| VVI replicate 1, K+1 | -1922.87916 | -1922.43980 | -1921.61250 |
| VVV replicate 1, K+1 | -2011.45052 | -2011.37591 | -2011.37571 |
| VVV replicate 12, K+1 | -2031.79829 | -2030.46713 | -2030.46666 |

The first cell stays well above the .01 cutoff. The last starts below it and becomes narrower (.001501 to about .0009006). Independently rescoring the three default parameter sets produces LLs .00848, .01481 and .02251 above the reported default LLs, consistent with the previously noted reporting lag. Compare direct likelihoods at returned parameters before interpreting small package gaps.

These nine fits demonstrate material sensitivity to stopping in selected examples. They do not estimate its contribution to the median gaps or predict the tighter-control 80-cell rate. They do rule out causal attribution of the default package differences solely to initialization breadth or inherent spike protection. A larger LL still establishes only a better in-sample objective value under compatible constraints, not a recovered or statistically regular solution.

The optional prior results answer a different estimation question. [The mclust book describes prior-based regularization](https://mclust-org.github.io/mclust-book/chapters/07_miscellanea.html). Its 0/80 count is not evidence that unregularized mclust prevents degeneracy. No prior or threshold-selection policy change is recommended here.

## Recommended revised plan

1. Fix the LTA boundary/inference contract, then zero-weight initialization/origins. Preserve unrelated working fixes and public selection policy.
2. Correct the statistical claims about thinness while making the three inconsistent diagnostics agree. State the normalized-weight convention and avoid presenting mass/parameter comparisons as rank or identification guarantees.
3. Surface per-component standardized minimum eigenvalue, marginal variance ratios, normalized Gaussian-observed mass, positive contributing support, and active-bound status prominently in model summaries/enumeration and start diagnostics. Include the component/measurement block that attains each extreme. Optional buried tables alone will leave the known moderate-size srl flat component practically silent. A neutral, always-visible continuous diagnostic needs no rejection threshold. Do not silently discard or reselect fits.
4. Resolve the emEM one-final-start contract explicitly. `R/search.R:61` and `R/kernel-em.R:406-409` make a valid protected Ward start occupy the only final slot. “Best n_starts” is therefore inaccurate. Documentation is the minimal repair; changing protection is a separate behavior decision with reproducibility consequences. The same protected-slot exception also matters for n_starts greater than one.
5. Retain kmeans default and emEM optional pending a matched recovery study. Separate an out-of-box workflow comparison from a controlled optimizer comparison. For the latter, use the same raw datasets, common independently recomputed likelihood/diagnostics, explicit compatible constraints, tight convergence plus sensitivity checks, and crossed/shared initial states where feasible. Retain failures and singularity warnings. Record generating labels/parameters; report classification and parameter recovery, selected K, held-out density where relevant, and continuous covariance/mass diagnostics alongside raw LL and finite-reference attainment. Use independent datasets as replication units and stratify K versus K+1 and covariance family. Avoid making the .01 criterion the definition of successful recovery.
6. Include shared-variance sensitivity as an explicit alternative statistical model; it reduces component-specific variance collapse but may miss genuine heteroscedasticity. Evaluate it against the same recovery targets. Automatic thresholds, extra admissibility constraints, alternative robust component families, or priors would require a separate user policy decision; none is needed for the correctness and transparency work above. Final no-rule emEM BLRT calibration remains an open empirical question.

Evidence: `check.R`, `check.log`, `mclust-small-check.csv`, and the copied original generating script in this directory; repository source locations cited above; `validation/review-fit-search-2026-10-22/{REVIEW.md,recomputed.csv,mclust-comparison.csv,test-review-regressions.R}` and retained probe objects; prior `equivalence/starts-study/codex/REVIEW.md`. The independent checks completed successfully. No production edits, full-suite claim, or universal package ranking is made.
