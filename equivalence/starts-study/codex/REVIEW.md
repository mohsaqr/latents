# Independent review of the random-starts study

The arithmetic checks out. The original interpretation does not. Ten or twenty starts usually match this study's finite reference at the generating K, but can fail almost completely on a difficult dataset. Over-extraction is much harder, yet its highest likelihood can be a singleton variance-collapse artifact. The observed tendency to select fewer classes is real in this sample; the claim that optimization error can *only* favor fewer classes is false. The proposed causal explanation about k-means versus Mplus is untested.

**Scope and verification.** I independently rebuilt the block and selection tables from the RDS start records, checked the saved summaries, parsed the Mplus outputs, performed the requested single 500-start refit, and ran three small mclust checks. There are 29 completed datasets: 16 LPA and 13 LCA, comprising 28 first replicates and one second replicate. They contain 116 fitted dataset/k cells and 58,000 start records. The planned design contains 216 datasets. This incomplete, uneven sample cannot establish reliability across the planned design. Blocks are repeated computational trials conditional on these datasets, not additional independently simulated datasets.

All reported empirical numbers below were recomputed. All writes are under this directory. No git commands, forks, or concurrent R processes were used. BLAS/OpenMP thread environment variables were set to one. The requested 500-start refit took 20.858 seconds; the three mclust fits together took 0.068 seconds. Audit computations were comfortably below the requested time budget.

**A. Code and analysis audit**

- **Block emulation:** `analyse.R` correctly combines Ward with disjoint groups of nine or nineteen k-means starts. There are 55 ten-start blocks and 26 twenty-start blocks per cell, leaving four and five k-means starts unused, respectively. I reproduced all 9,396 block rows using explicit index slicing rather than the original split/pmax implementation. Given a dataset, Ward is deterministic; the other starts use the same initialization mechanism without an adaptation based on earlier successes. This supports marginal emulation of the default start scheme. It does not prove that every block corresponds to a separately seeded real call. The comment in `simulate.R` citing `check-emulation.R` is unsupported: that file is absent.
- **Tolerance and selector:** success means a likelihood deficit strictly less than 0.001. The code uses raw `max()`, whereas `.multilpa_select_start()` permits a relative likelihood tie of `1e-10 * (1 + abs(best))`, prefers convergence within that tie, and then chooses the earliest start. Across the 116 saved cells, the selected-versus-raw-maximum difference ranged from −3.5832727e-7 to zero, immaterial here. I found no failed starts and 23 nonconverged starts; every start within 0.001 of its cell's reference was converged. Thus unconverged reference solutions do not explain these results. Boundary solutions do.
- **BIC:** the calculation `-2 * LL + n_parameters * log(n)` is correct for these complete, unweighted, single-level datasets. The saved parameter counts equal `9k−1` for four-indicator VVI LPA and `7k−1` for six-binary-item LCA. These also match Mplus in every matched cell. I independently reproduced all 2,349 selection rows and verified `summaries.rds`. Correct arithmetic does not validate ordinary BIC interpretation at singular or near-singular solutions.
- **Pairing blocks across k:** using the same block index defines one particular coupling of the per-k fits. Since the same seed is reset at each k but random-number consumption changes with k, this is not established as the distribution of real repeated enumeration calls. I independently shuffled block labels across k 1,000 times for each start budget. Every shuffle retained exactly 105 downward changes and zero upward changes for ten starts, and 47 downward changes and zero upward changes for twenty starts. This concern does not alter the observed selection result in this sample.
- **Mplus parsing:** of 71 output files, 67 terminated normally and four did not. Every normally terminated output has a usable descending final-stage ranked list whose first likelihood exactly matches the separately printed H0 likelihood. All 51 comparison CSV rows match their output and their latents raw reference. The inner join appropriately restricts comparisons to available cells, but this is a subset of both programs' completed work. The parser assumes the current heading/blank-line layout and would need defensive checks for other output formats.
- **Replication:** `mplus_replicated` is the study's strict 0.001 criterion, not necessarily Mplus's own judgment. For `lca_k4_high_n300_r1-k5`, the first two ranked values are −883.321 and −883.335; the study flag is false, although Mplus prints that its best likelihood was replicated. Keep these diagnostics distinct. Mplus's likelihoods are printed to three decimal places, limiting precision of inter-program differences.
- **Raw-record limitation:** `fits/*.rds` contain start-summary data frames, not full fitted objects. `smallest_class` and warning names are copied from the final selected fit onto every start row. They cannot identify each individual start's class size. The `boundary` column, in contrast, is recorded per start.

**B. Verdicts on the four claims**

**1. “At the true number of classes, 10 or 20 starts reliably reach the best maximum”: needs substantial qualification.** The measured event is reaching the best likelihood *observed in these same 500 starts*, not a demonstrated global optimum.

| At generating K | Ten starts | Twenty starts |
|---|---:|---:|
| LCA | 670/715 = 93.71% | 325/338 = 96.15% |
| LPA | 826/880 = 93.86% | 391/416 = 93.99% |
| Combined | 1,496/1,595 = 93.79% | 716/754 = 94.96% |

Ward alone reaches the reference on 24/29 datasets at K. All block failures at K occur in low-separation, n=300 datasets:

| Dataset | Individual starts within 0.001 | Ten-start successes | Twenty-start successes |
|---|---:|---:|---:|
| `lpa_k4_low_n300_r1` | 1/500 = 0.2% | 1/55 | 1/26 |
| `lca_k3_low_n300_r1` | 19/500 = 3.8% | 16/55 | 13/26 |
| `lca_k4_low_n300_r1` | 113/500 = 22.6% | 49/55 | 26/26 |

For the difficult LPA dataset, the worst ten- and twenty-start deficits from the latents reference are 5.557614 and 5.314898. Even the 500-start reference is 4.755592 below Mplus at K. Calling ten or twenty starts reliably sufficient for difficult cases is unjustified. Because each block helps define the finite reference, these are internal benchmark rates, not externally validated probabilities of finding an optimum.

**2. “At an over-extracted k, 10–20 starts often stop at a local maximum”: supported as a likelihood-attainment observation, with qualification about degeneracy.**

| Family and fitted k | Ten-start reference attainment | Twenty-start reference attainment |
|---|---:|---:|
| LCA, K+1 | 465/715 = 65.03% | 279/338 = 82.54% |
| LPA, K+1 | 279/880 = 31.70% | 174/416 = 41.83% |
| LCA, K+2 | 400/715 = 55.94% | 243/338 = 71.89% |
| LPA, K+2 | 189/880 = 21.48% | 119/416 = 28.61% |

The important LPA caveat is that failing to find the highest attained likelihood can mean failing to find a pathological spike. It is not automatically inferior recovery of the generating classes.

The independent Mplus comparisons are:

| Family and k | Matched cells | Range of latents minus Mplus LL |
|---|---:|---:|
| LCA, K | 13 | −0.000455 to +0.000466 |
| LCA, K+1 | 13 | −0.000797 to +0.000311 |
| LPA, K | 13 | −4.755592 to +0.000384 |
| LPA, K+1 | 12 | −12.560193 to +12.206589 |

Excluding `lpa_k4_low_n300_r1`, all matched LPA K likelihood differences have absolute value below 0.0005. The four large positive LPA K+1 differences are +8.244538, +12.206589, +8.920593, and +10.744891; all four latents selected fits have a modal singleton and boundary-active best starts.

**3. “Local-maximum error can only bias BIC toward fewer classes”: not supported as a general claim; mathematically false.** Let `e_k = LL_reference(k) − LL_short(k) ≥ 0`. Then

`BIC_short(k) = BIC_reference(k) + 2e_k`.

For a larger h and smaller l, the change in their BIC difference is `2(e_h − e_l)`. Its sign can be positive or negative. Optimization error at a smaller k can favor a larger k. The present sample showing only one direction does not establish a one-direction theorem.

Empirically, no short run chose more classes than the reference. All disagreements occur in LPA with generating K=4:

| Comparison | Ten starts | Twenty starts |
|---|---:|---:|
| Fewer classes, all LPA runs | 105/880 = 11.93% | 47/416 = 11.30% |
| Fewer classes, LPA generating K=4 | 105/275 = 38.18% | 47/130 = 36.15% |
| `lpa_k4_high_n300_r1`: five → four | 51/55 | 22/26 |
| `lpa_k4_moderate_n300_r1`: five → four | 54/55 | 25/26 |

Both selected five-class reference fits have modal singleton classes and boundary-active best starts. Every disagreement therefore selects the generating four classes rather than the spike-driven five-class reference. In the high-separation dataset, the reference BICs are 4353.073206 for four classes and 4348.322647 for five; in the moderate-separation dataset they are 4205.536519 and 4197.186100. More starts worsen class recovery in these particular cases.

As a sensitivity analysis, I removed boundary starts using the actual per-start `boundary` flag, rebuilt both the references and the short fits, and reran selection. All 2,349 short/reference comparisons then agree. In LPA K=4, the reference selects four classes on three of five datasets instead of one of five; it selects three on the other two. This isolates the observed disagreements to variance-boundary solutions. It is not a complete degeneracy screen: tiny classes and near-boundary interior solutions can remain, and this is a changed reference, not a proof of correct enumeration.

**4. “The cause is that latents' k-means starts explore less widely than Mplus's screened random starts”: not supported causally.** The schemes differ, but the comparison also changes optimization, screening, convergence rules and boundary behavior. STARTS=400 80 is not equivalent to 500 full fits. The [Mplus manual](https://www.statmodel.com/HTML_UG/chapter16V8.htm) describes its two-stage random-start procedure. This study does not isolate initialization or measure exploration breadth. The near agreement throughout LCA and mixed signs in over-extracted LPA do not establish the proposed explanation. A controlled start-strategy comparison within the same likelihood, optimizer, tolerances and admissibility policy is needed.

**C. Singleton verification and package implications**

I loaded the package with `pkgload::load_all(".")`, sourced `simulate.R`, regenerated `lpa_k4_high_n300_r1`, and fitted `lpa(..., n_profiles=5, model="VVI", n_starts=500, seed=20261035)` once. All 500 start likelihoods reproduced the saved records exactly. The regenerated indicators and Mplus's data file differed by at most 5.329071e-15.

- Selected LL: **−2048.678108984695**; 44 parameters; BIC **4348.322646854273**.
- Modal class counts: **89, 36, 114, 1, 60**.
- The singleton is profile 4, containing observation **277**; effective membership is **0.9999999999997944**.
- Its variances for **y1, y2, y3, y4 are exactly 0.000001, 0.000001, 0.000001, 0.000001**.
- Those variances are approximately **2.43154e-7, 2.63765e-7, 2.45624e-7, 3.48833e-7** times the respective sample variances.
- Its means coincide with the singleton's indicator values. Both `boundary` and `small_classes` are true; the fit raises `latents_boundary` and `latents_small_classes` warnings.

This is a spurious Gaussian likelihood spike. For a component centered on one observation, its density there grows as the reciprocal square root of the product of its variances. Shrinking the variances drives the unconstrained mixture likelihood without bound while the remaining components account for the other observations. With `min_variance=1e-6`, the actual fitted objective is constrained and finite: it is inaccurate to call that constrained likelihood itself unbounded.

I also directly rescored the fitted parameters without any additional optimization. Dividing only the singleton's four variances by ten raises LL from −2048.678108984695 to **−2044.072938798708**, a gain of **4.605170185988**. This demonstrates the spike mechanism on the actual data; the narrower parameters lie outside the original variance constraint.

Mplus is not a guaranteed regular alternative. Among the 12 matched LPA K+1 outputs, **nine** print a variance as `0.000`, **11** have a non-positive-definite derivative-product-matrix warning, and only **five** replicate the best value under the study's criterion. Printed zero is rounded, not evidence of an exact numerical floor. For the refitted dataset's Mplus five-class solution, the smallest modal class has **five** observations, effective membership **4.93324**, and displayed variances **0.138, 0.003, 0.261, 0.000**. Its LL is **−2059.423**, is unreplicated, and carries a nonidentification warning. Mplus did not report the same singleton here; it did report another suspect near-singular solution. The blanket statement that it does not report spikes is too strong.

The package should preserve its current default estimates and seeded results, while making this failure mode harder to overlook:

1. Report per-start minimum modal/effective class size, minimum variance relative to each indicator's variance, and count of active variance bounds. Keep these alongside convergence and likelihood. Do not collapse the diagnostic into one flag.
2. Add a prominent diagnostic for a tiny class combined with covariance collapse, and surface it in enumeration/model comparison. Existing boundary warnings are useful; `effective_count < 1` alone can miss an exact singleton or depend on rounding around one. Here it fired because the effective count is fractionally below one.
3. Offer explicit opt-in admissibility rules: minimum effective membership and/or a variance floor relative to the data scale, with the chosen thresholds recorded. These change the statistical problem; no universal threshold is established by this study. A modal-size filter alone is discontinuous and inadequate.
4. Offer opt-in exclusion of degenerate starts, retaining the raw best LL and exclusion reasons. Report failure if no admissible fit remains. A diagnostic filter is not the same as optimizing under a class-size constraint. Do not silently replace the default likelihood winner or claim global optimality after filtering.
5. Compare start strategies using the same admissibility policy and independent seeds. Report both raw-best and admissible-best attainment, plus recovery relative to the generating K. Otherwise more successful optimization of a spike can be mistaken for better mixture modeling.

**D. Small independent-optimizer check**

I ran mclust **6.1.2**, model `VVI`, at K+1 on three datasets, with its default hierarchical initialization, no prior, and one fit per dataset. This is a small check, not a multistart competition; [mclust documents its hierarchical initialization](https://mclust-org.github.io/mclust/reference/mclustBIC.html).

| Dataset | k | mclust LL | latents 500-start LL | Mplus LL | mclust smallest modal class |
|---|---:|---:|---:|---:|---:|
| `lpa_k2_low_n300_r1` | 3 | −1786.507704 | −1765.487411 | −1777.694 | 65 |
| `lpa_k3_high_n300_r1` | 4 | −2059.400433 | −2055.788193 | −2043.228 | 12 |
| `lpa_k4_high_n300_r1` | 5 | −2067.468105 | −2048.678109 | −2059.423 | 8 |

All three mclust likelihoods are lower than both comparison values. These particular fits avoid modal singletons, but a single initialization cannot establish a better admissible optimum or identify the cause of inter-program differences.

**Reproducibility.** `recompute.R` contains the independent reconstruction, assertions, requested refit and optimizer checks. `optional-check.R` contains the separately completed optimizer and boundary-exclusion checks; `verify-refit.R` checks all saved start likelihoods, data identity and direct density calculations. `check_pairing.py` checks output parsing and pairing sensitivity. Evidence is in `recompute-run2.log`, `optional-check.log`, `verify-refit.log`, `pairing.log`, the diagnostic CSVs, and `singleton-refit.rds`. The first audit attempt stopped at a script syntax error; the second completed the reconstruction and refit but stopped at mclust's unqualified `mclustBIC` lookup when only its namespace was loaded. Both script issues were corrected; the optimizer checks were rerun successfully with mclust attached. The expensive refit was performed only once. Final scripts parse successfully; numerical assertions and the independent follow-up checks pass.
