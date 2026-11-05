# Review of the uncommitted wild-card verbs

Date: 2026-11-02. Branch: `features`.

Reviewed `R/wildcards.R`, `tests/testthat/test-wildcards.R`, the three new Rd files, and the NAMESPACE additions, with the existing methods inspected as reference implementations. Six reproduced findings follow: five medium severity (P2), one low severity (P3). Line numbers refer to the reviewed working tree.

Run the commands below from `/Users/mohammedsaqr/Documents/Github/latents`. Every reproduction was executed with `devtools::load_all(".")`, `NOT_CRAN=true`, and both thread limits set to one. The null graphics device prevents methods that print on call from writing a PDF.

## 1. Failed regression candidates abort every grid table getter

**Severity:** P2 / medium. **Location:** `R/wildcards.R:656-658`.

Regression enumeration retains failed fits as condition objects, rather than `NULL`. Filtering only `NULL` passes these conditions to `diagnostics()`. One failed candidate consequently prevents access to the successful candidates' tables. This contradicts the promise in `man/wildcard-tables.Rd` that failed candidates have no row.

**Failing input and reproduction:** a grid with one supported and one unsupported class count.

```sh
VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 NOT_CRAN=true Rscript -e '
devtools::load_all("."); grDevices::pdf(file=NULL)
d <- data.frame(y=rep(c(-1,1),15))
g <- suppressWarnings(enumerate_regressions(y~1,d,
  n_classes=c(1,40),n_starts=1,seed=1))
print(g$table[,c("n_classes","note")])
print(lapply(g$fits,class))
print(tryCatch(get_checks(g),error=function(e)
  list(class=class(e),message=conditionMessage(e))))
'
```

**Observed:** the one-class fit succeeds; the 40-class cell records `30 rows cannot support 40 classes.` Its stored class is `c("latents_bad_data", "error", "condition")`. `get_checks(g)` aborts with `latents_not_a_fit`: `` `x` must be a fitted model of this package, but it is an object of class "latents_bad_data". `` The same unconditional `diagnostics()` call serves all six table getters.

**Expected:** return the successful candidate's checks under `classes1_groups1`, omitting the failed candidate. If no candidate supplies the requested table, raise the grid-specific `latents_bad_argument` already provided at lines 667-670. Candidate naming itself is correct here; failure filtering is the defect.

## 2. Regression posterior plots include zero-weight rows and omit the group classification

**Severity:** P2 / medium. **Location:** `R/wildcards.R:97-104`; also `R/wildcards.R:577-584` for `plot_all()`.

The dispatcher treats an offered view named `"posteriors"` as sufficient. The existing regression method reads raw `expectation$tau`, whereas the new table-based drawers use eligible rows and aligned weights. For hierarchical regressions this also bypasses the groups altogether. This violates both the zero-weight exclusion and the two-level return contract in `man/wildcard-plots.Rd:103-110`.

**Failing input and reproduction:** a weighted hierarchical regression with five zero-weight groups.

```sh
VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 NOT_CRAN=true Rscript -e '
devtools::load_all("."); grDevices::pdf(file=NULL)
set.seed(5)
d <- data.frame(g=rep(1:20,each=5),t=rep(0:4,20))
d$y <- c(-2,2)[1+d$g%%2]+0.3*d$t+rnorm(100)
d$w <- rep(c(0,1,3,1),length.out=20)[d$g]
f <- suppressWarnings(mixture_regression(y~t,d,n_classes=2,id="g",
  n_group_classes=2,weights="w",n_starts=1,seed=1,vcov_type="none"))
p <- plot_posteriors(f)
cat("posterior plot class",class(p)[1],"rows",nrow(p$data),"\n")
print(get_entropy(f)[,c("level","n_units")])
print(vapply(plot_sizes(f),function(p)sum(p$data$count),numeric(1)))
cat("plot_all contains",paste(names(plot_all(f)),collapse=", "),"\n")
'
```

**Observed:** `plot_posteriors(f)` is a single ggplot with **100 rows**. The diagnostic table classifies **75 individuals and 15 groups**; `plot_sizes(f)` correctly returns both levels with totals 75 and 15. `plot_all(f)` includes only `posteriors`, alongside `sizes_individuals`, `sizes_groups`, `entropy_individuals`, `entropy_groups`, `avepp_individuals`, and `avepp_groups`.

**Expected:** a `latents_plots` result containing `individuals` and `groups`, with 75 and 15 eligible certainty values respectively, and corresponding entries in `plot_all()`. The zero-weight leak also reproduces for weighted binomial-count, Poisson, negative-binomial, and ordinal regressions; it does not require the hierarchical model.

## 3. Trajectory posterior plots map to the full stacked classification chart

**Severity:** P2 / medium. **Location:** `R/wildcards.R:50,97-104`.

There is a second, distinct semantic mismatch in delegation by view name: trajectory and growth methods use `"posteriors"` as an alias for `"classification"`. That chart contains every class probability for every person. The new verb promises the probability of each unit's **assigned** class (`man/wildcard-plots.Rd:85`). This mismatch occurs even without sampling weights or two classification levels.

**Failing input and reproduction:** an ordinary two-class trajectory regression.

```sh
VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 NOT_CRAN=true Rscript -e '
devtools::load_all("."); grDevices::pdf(file=NULL)
f <- mixture_regression(score~wave,subset(growth_scores,student<=30),
  n_classes=2,id="student",class_level="group",n_starts=1,seed=1)
p <- plot_posteriors(f)
q <- plot_classification(f)
cat("units",nrow(f$expectation$group_tau),"plotted rows",nrow(p$data),"\n")
cat("minimum probability",min(p$data$probability),"\n")
cat("same classification data",identical(p$data,q$data),"\n")
'
```

**Observed:** `units 30 plotted rows 60`; minimum plotted probability `7.773641e-32`; `same classification data TRUE`. The plot is the full stacked two-class membership chart. The growth-mixture method has the same alias.

**Expected:** the certainty view containing 30 assigned-class probabilities, each at least `1/2`, as supplied by the new certainty drawer. The full membership chart belongs to `plot_classification()`. Fixing the new dispatch does not require changing that existing method or alias.

## 4. One-class AvePP is incorrectly rejected as having no classified units

**Severity:** P2 / medium. **Location:** `R/wildcards.R:120-128`.

The eligibility filter permits one-class levels only for sizes. It therefore discards a perfectly defined one-by-one average-posterior matrix. This differs from the profile-family behavior and from the new getter's own result.

**Failing input and reproduction:** a one-class regression.

```sh
VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 NOT_CRAN=true Rscript -e '
devtools::load_all("."); grDevices::pdf(file=NULL)
f <- mixture_regression(score~hours,study_hours,n_classes=1,
  n_starts=1,seed=1,vcov_type="none")
print(get_avepp(f))
print(tryCatch(plot_avepp(f),error=function(e)
  list(class=class(e),message=conditionMessage(e))))
print(names(suppressMessages(plot_all(f))))
p <- lpa(faithful,"waiting",n_profiles=1,n_starts=1,seed=1)
print(plot_avepp(p)$data)
'
```

**Observed:** `get_avepp(f)` returns one cell with `n_assigned=900` and `average_posterior=1`. `plot_avepp(f)` raises `latents_nothing_to_plot`, saying the fit `has no classified units to draw`. `plot_all(f)` omits AvePP. The one-class LPA correctly draws a cell labelled `1.00`.

**Expected:** render the defined one-by-one matrix for the regression too. The same filter affects other families using `.latents_wild_classes()`. Refusing relative entropy when `K=1` is a separate, legitimate limitation; it does not justify rejecting AvePP.

## 5. `plot_all()` converts bootstrap integrity failures into “nothing to plot”

**Severity:** P2 / medium. **Location:** `R/wildcards.R:514-521,527-529`.

The catch for all `latents_bad_argument` conditions includes genuine data-integrity failures. It discards their explanation and correction. When other views remain, this branch can return a partial success; for the single-view bootstrap object it replaces the original error entirely.

**Failing input and reproduction:** a naturally reordered copy of a bootstrap table. The reordering intentionally triggers the existing integrity guard; the defect is the new wrapper masking that guard.

```sh
VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 NOT_CRAN=true Rscript -e '
devtools::load_all("."); grDevices::pdf(file=NULL)
f <- lpa(faithful,c("eruptions","waiting"),2,
  model="EEI",n_starts=1,seed=1)
b <- parameter_inference(f,method="bootstrap",iter=3,n_starts=1,seed=1)
bad <- b[rev(seq_len(nrow(b))),,drop=FALSE]
e1 <- tryCatch(plot(bad),error=identity)
e2 <- tryCatch(plot_all(bad),error=identity)
cat("direct:",class(e1)[1],conditionMessage(e1),"\n")
cat("wildcard:",class(e2)[1],conditionMessage(e2),"\n")
'
```

**Observed:** direct `plot()` raises `latents_bad_argument`: `This table no longer matches its bootstrap replicates (a subset or a reordered or modified copy loses them). Use the table parameter_inference(method = "bootstrap") returned.` The wildcard instead messages `Not drawn for this object: plot.` and raises `latents_nothing_to_plot`: `plot_all() has no view this object can draw.`

**Expected:** preserve the integrity condition and its corrective instructions. Only conditions establishing that a particular view is inapplicable should be skipped. The generic `latents_bad_argument` superclass does not establish that.

## 6. Refusal messages recommend the very function they just refused

**Severity:** P3 / low. **Location:** `R/wildcards.R:308-311,316-320`.

Recommendations depend only on the coarse object kind. Homogeneous transitions share `profile_model` with ordinary profiles, but their actual offered views differ. The resulting “For this object” list includes unsupported functions, including the failed call itself.

**Failing input and reproduction:** request a raincloud from a homogeneous transition fit.

```sh
VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 NOT_CRAN=true Rscript -e '
devtools::load_all("."); grDevices::pdf(file=NULL)
f <- suppressWarnings(lta(subset(course_engagement,student<=10),
  c("browse","lectures"),"student",n_profiles=2,time="sequence",
  n_starts=1,seed=1))
e <- tryCatch(plot_raincloud(f),error=identity)
cat(class(e)[1],"\n",conditionMessage(e),"\n")
cat("supported:",paste(.latents_plot_offered_views(f),collapse=", "),"\n")
'
```

**Observed:** `latents_nothing_to_plot` correctly refuses the request, then recommends `plot_raincloud()`, `plot_parallel()`, and `plot_pairs()` in its `For this object:` list. None is offered by this transition fit. Its offered views are `transitions, profiles, bars, heatmap, responses, sequences, sizes, posteriors, entropy, avepp`.

**Expected:** recommend functions mapped to views this particular object actually offers. The refusal class is correct; the recovery advice is misleading.

## Checked and found correct

- **New tests:** `tests/testthat/test-wildcards.R` passes **160 expectations**, with no failures, warnings, or skips. Its calibration fixture is an unweighted two-class regression; its grid fixture has no failed regression cell. The family smoke test checks that plots exist, so it does not detect the semantic mismatch in finding 3. These coverage gaps explain why the findings pass the current suite.
- **Released API:** `test-api-compat-0921.R` passes **44 expectations** and `test-compat-0921.R` passes **10**, both without failures, warnings, or skips. The reviewed NAMESPACE diff only adds 24 exports. No existing function definitions are modified by this feature's source diff.
- **Exports and documentation:** all 24 new verbs are exported, none has a `what` formal, and all three new Rd files parse successfully. The documentation accurately specifies the intended numbers and views; the discrepancies above are in their implementation.
- **Table-driven drawers:** with nonuniform positive sampling weights and zero-weight groups, checked regression, trajectory regression, additive, full cross-level, and general-transition fits. `.latents_wild_unit()` effective counts match `get_classification()` exactly; AvePP cells match `get_avepp()` exactly; weighted certainty means match the classification table's diagonal AvePP exactly; entropy subtitles match the formatted `get_entropy()` value. Maximum observed numeric differences in these checks were **zero**. Eligible row counts and modal-count labels are consistent. The drawers' ggplots build successfully. The posterior problems above arise where dispatch bypasses these drawers.
- **Profile delegation:** weighted LPA sizes match the classification table exactly, zero-weight rows are omitted from the certainty plot (40 plotted/eligible from 60 original rows), and direct diagnostics-result table access preserves the stored table. Profile plots in the new tests agree with their existing drawing methods. A two-level `multilca()` produces both classification levels and its plots build. Noise LPA retains the existing sizes view, whose caption explicitly excludes noise from the bars; `plot_all()` builds successfully and omits the unsupported separation views.
- **Grid names and stacking:** tested single-level profile, multilevel profile, additive-family, general-transition, and regression grids. Candidate names match stored fit positions; each stacked candidate subtable equals the corresponding getter called on that fit. In particular, unnamed regression fits correctly become `classes1_groups1` and `classes2_groups1`. Profile-grid `NULL` failures are skipped correctly. The regression condition-object case is finding 1.
- **Grid drawing:** criteria, entropy, AvePP, smallest-class, heatmap, delta, trust, and all-view wrappers produce buildable plots where supported. Additive-family grids correctly refuse their unavailable selection views with `latents_unsupported_view`; `plot_all()` respects that grid's advertised supported views. All five tested grid configurations support print-only overview without crashing.
- **Other fit kinds:** built every plot returned by `plot_all()` across growth mixture, trajectory regression, single-level and hierarchical regression, additive, dispersion, cross-level, pooled, general and homogeneous transitions, LPA, LCA, and covariate profiles. No additional incidental crashes were found. Weighted binomial-count, Poisson, negative-binomial, and ordinal regressions also exercise the new classification drawers; their sizes and AvePP match the getters, subject to finding 2 for their delegated certainty plots.
- **Other result kinds:** valid bootstrap inference, `bootstrap_lrt()`, sensitivity, and comparison results produce their expected `plot_all()` entry. The bootstrap-LRT plot builds and its print-only overview returns the original object invisibly. Diagnostics results unwrap to their original fit for plotting. Unsupported ordinary objects produce `latents_not_a_fit`; pooled classification views refuse with `latents_nothing_to_plot` because pooled estimates carry no classification posteriors.
- **Bare-list summaries:** cross-level `summary()` returns a bare list, but `overview(..., plots=FALSE)` correctly selects the object's own printer, produces no list-component dump, and returns the original object invisibly. Pooled objects' data-frame summaries also work. No bug found in this fallback.
- **Boundary of conclusions:** this was a targeted correctness review, not a full package test-suite run or a new statistical-engine validation. Each numbered finding has a runnable reproduction above. Repairs can be confined to the newly added dispatch/helpers, tests, and documentation; changing existing methods is unnecessary.

## Review artifacts

No package source, test, documentation, or namespace file was edited during this review. This report is the only intentional new file. An early R probe invoked an existing method that prints on call and updated the pre-existing, ignored `Rplots.pdf` through R's default graphics device; later probes used a null device. Existing `HANDOFF.md`, `LEARNINGS.md`, and other logs were left unchanged as requested.
