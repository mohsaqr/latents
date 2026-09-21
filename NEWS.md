# multilpa 0.11.8

## The introductory vignette, rewritten

`vignette("multilpa")` now opens with what the model found. The order is the
profile means, then how big each profile and class is, then the profile
probabilities, then everything else -- estimation diagnostics, the plots, the
sequences. Printing the fit gives the means one row per profile with each
profile's size, so the guide starts from a result rather than from plumbing.

Two passages are gone. One introduced the notation `pi_{k|m}` and `omega_m`,
said the probabilities lie between zero and one and sum to one, and never used
the symbols again; what it was reaching for is now said where the table is read.
The other listed the model's assumptions in a block, of which the clause that
does any work -- the default diagonal covariance treats the indicators as
independent within a profile, which a shared association between activity
measures can violate -- has moved to Limitations, where a reader can act on it.
A sentence restating the ICC formula's symbols in words went with them.

Every figure in the text was recomputed against the standardized indicators.

# multilpa 0.11.7

## Transition networks, with `get_tna()` and `get_group_tna()`

A fitted transition model hands itself to the tna package:

```r
moves <- fit_transitions(course_engagement, vars = activity, id = "student",
                         time = "sequence", n_profiles = 3, n_group_classes = 3)

get_tna(moves)        # one network for the whole sample
get_group_tna(moves)  # one network per latent group class
```

`get_tna()` returns a `tna` model and `get_group_tna()` a `group_tna`, so every
verb of that package applies: `centralities()`, `communities()`, `cliques()`,
`compare()`, `plot()`. `tna` is in Suggests and is required only by these two.

What is handed over is the **estimates**. tna's own method for a fitted mixture
model, `group_model.mhmm()`, takes the model's grouping and counts transitions
between modal assignments, discarding the estimated matrices; these verbs pass
the matrices themselves, maximised over the full posteriors, so an observation
split 0.6/0.4 between two profiles contributes to both rather than entirely to
one. The two answer different questions and will differ wherever classification
is uncertain.

The whole-sample network is the **marginal** transition matrix: the expected
counts summed across classes and then normalised, not the class-probability
average of the class matrices. A class holding a tenth of the units but a fifth
of the transitions counts for its transitions, which is what a marginal
probability means. A state no observation ever leaves is reported as a
self-transition rather than as `NaN`.

# multilpa 0.11.6

## Single-level fits, with `id = NULL`

`multilpa()` accepts `id = NULL`, which treats the observations as independent
-- each row is its own unit, `n_group_classes` becomes one -- and fits the
ordinary Gaussian or latent-class mixture the two-level model reduces to:

```r
multilpa(faithful, vars = c("eruptions", "waiting"), id = NULL, n_profiles = 3)
```

`id` has **no default**, and omitting it still raises, now with
`multilpa_bad_argument` naming both ways forward rather than R's bare
"argument \"id\" is missing". A single-level model is not the default in a
package for multilevel models, and a forgotten grouping must not quietly become
a different model: on data with 106 students in 1,422 enrolments it would report
1,422 independent observations and change every standard error. Every
`id = NULL` fit raises a `multilpa_single_level` warning saying what it fitted.

Nothing about the likelihood changes, and the fit is identical to the one you
get by numbering the rows yourself and asking for one group class; that
equivalence is asserted in the tests. Every verb of the package works on it,
including `parameter_inference()`, whose cluster bootstrap degenerates to the
ordinary nonparametric bootstrap when each unit holds one row -- which is the
right bootstrap for independent observations -- and `vcov_type = "robust"`,
whose per-group scores become per-row scores.

The unit column is fabricated internally as `.observation`. It is not handed
back by `get_data(x, "data")` or `get_data(x, "assignments")`, and data that
already carry a column of that name raise `multilpa_bad_data` rather than having
it silently overwritten. Asking for more than one group class without an `id`
raises `multilpa_bad_argument`: one observation per unit leaves no composition
for a second-level class to differ in.

A single-level fit prints as what it is --- `Latent profile analysis: 2
profiles` over `160 observations` --- rather than reporting one group class and
160 groups of one.

`fit_staged()` and `fit_transitions()` still require `id`, and say so in their
own documentation rather than inheriting `multilpa()`'s: both have a second
level by construction.

# multilpa 0.11.5

## Every covariance structure now reports uncertainty

`parameter_inference()` gains `method = "bootstrap"`. It resamples *groups*
with replacement, refits inside the same covariance family, undoes the
relabelling each refit comes back with, and reports the percentile interval and
the standard deviation of the replicates:

```r
fit <- multilpa(course_engagement, activity, id = "student", n_profiles = 2,
                n_group_classes = 1, volume = "varying", shape = "equal",
                orientation = "axis")
parameter_inference(fit, method = "bootstrap", iter = 999, seed = 1)
```

This closes the gap the fourteen structures opened. Ten of them --- EII, VII,
VEI, EVI, VEE, EVE, VVE, EEV, VEV and EVV --- constrain the volume, the shape
or the orientation, which the Wald coordinates cannot express, so the Wald path
refused them. `enumerate_classes(structure = )` could therefore select a model
the package would not then do inference on: on this package's own
`course_engagement` data the best of the fourteen is `VEI`. Every structure the
grid can select now reports an interval.

`vcov()` and `confint()` reach the same path through `...`. `confint()` returns
the percentile interval the table reports rather than rebuilding
`estimate +/- z * se` from the bootstrap standard error, so one fit does not
give two different intervals.

Groups are the resampling unit, not rows, so the interval carries the same
independence assumption as `vcov_type = "robust"` and not the stronger one
`"observed"` makes. Against the cluster-robust sandwich on `course_engagement`,
where both are available, the bootstrap standard errors agree to within 11%
across every estimated parameter. `statistic`, `p_value` and `p_adjusted` are
`NA` on this path: the interval is the inference, and putting a normal
approximation back on top of the replicates it was read from is the assumption
the path exists to avoid.

A fit that holds a measurement block is refused with
`multilpa_unsupported_inference`, because the held values came from another fit
and resampling these data does not resample them.

## Fixes

* `fit_staged()` accepted a measurement fitted with a constrained covariance
  structure and silently fitted the second stage as the nearest structure
  `variance_model` and `covariance_model` can name: a `VEI` measurement came
  back recorded as `VVI`. The held means and variances were the right numbers,
  but under the wrong model, and `n_parameters_with_measurement` counted the
  spread as `VVI`'s. `multilpa()` already refuses to hold the variances of such
  a structure, so a staged fit of one was never supported; it is now refused up
  front with `multilpa_bad_stage`, naming the structure. Refit the measurement
  as EEI, VVI, EEE or VVV, or fit both levels together with `multilpa()`, which
  estimates all fourteen directly.

* `bootstrap_lrt()` refitted each replicate from `variance_model` and
  `covariance_model` alone, which cannot express a constrained structure: a
  `VEI` model came back as `VVI`, two parameters wider, so the reference
  distribution belonged to a different pair of models than the statistic
  compared against it. Replicates are now refitted from the structure the model
  records.

## Corrections to the 0.11.3 notes below

The 0.11.3 section was written in two passes and the earlier one was left
standing. Two of its statements are wrong, and are corrected here rather than
edited away:

* "Six of `mclust`'s fourteen models --- the ones with a non-axis-parallel
  orientation --- remain unavailable" was superseded within that same release
  by "All fourteen mclust covariance structures". All fourteen ship, and their
  parameter counts match `mclust`'s own.
* The refusal list "`parameter_inference()`, `vcov()` and `confint()` refuse
  EII, VII, VEI and EVI" named four of the ten that were actually refused. The
  rule, rather than the list: Wald inference was available exactly for the four
  structures reachable by `variance_model` and `covariance_model` alone --- EEI,
  VVI, EEE and VVV. As of this release the bootstrap covers all fourteen.

# multilpa 0.11.4

## Equivalence tests no longer ship with the package

Every test that compares multilpa against other software or published output
now lives in `equivalence/` in the source repository, and that folder is not
part of the built package. This covers Mplus, mclust, tidySEM, depmixS4 and
the in-house enumeration. The shipped `tests/testthat/` keeps unit and
invariant tests only. Run the comparisons with `Rscript equivalence/run.R`.

`mclust` and `tidySEM` are no longer in `Suggests`; nothing shipped uses them.

## Latent transition references

`fit_transitions()` is now checked against the external references pinned by
the JStats project: depmixS4 Gaussian and binary panel models, two Mplus LTA
runs, and the plain-LTA row of Table 5 of Muthen & Asparouhov (2022,
*Psychological Methods*). Log-likelihoods, criteria and parameters agree to
the reference's precision. On one binary depmixS4 fixture (K = 3), the
recorded depmixS4 solution is a local maximum: multilpa reaches a likelihood
0.016 higher. A brute-force path-sum likelihood confirms both values.

# multilpa 0.11.3

## Within-person profiles, without doing the transform yourself

`multilpa(centering = "person")` subtracts each group's own mean from its rows
before the measurement model sees them, so the profiles become profiles of
*change* rather than of level. This is the person-mean-centred, within-person
design (Quintana, 2021; Voelkle, Brose, Schmiedek, & Lindenberger, 2014), and
until now it meant transforming the data frame first --- which left the fit
holding numbers that no longer matched the frame it was given, so every verb
that checks row alignment was comparing centred values against raw ones.

```r
multilpa(course_engagement,
         c("browse", "lectures", "forum_read", "forum_post", "attendance"),
         id = "student", n_profiles = 3, n_group_classes = 1,
         centering = "person")
```

The offsets stay on the fit, so nothing downstream has to guess which scale it
is looking at: `get_data(x, "data")` returns the columns you supplied,
alignment is still checked against them, and the bivariate residuals are
computed on the scale the model was estimated on. `"grand"` subtracts one mean
per indicator instead. `"none"` is the default and changes nothing.

Centring removes exactly the between-unit variation, so `"person"` refuses with
`multilpa_bad_data` when it leaves an indicator constant --- which is what
happens when a unit has one observation of it. With `"person"` the group
classes become types of *change pattern*, not types of unit.

## All fourteen mclust covariance structures

`volume`, `shape` and `orientation` are the three pieces a covariance
decomposes into --- `Sigma_k = lambda_k * D_k * A_k * D_k'` --- and naming
them reaches every model `mclust` fits, where `variance_model` and
`covariance_model` reached four.

| `volume` | `shape` | `orientation` | model |
|---|---|---|---|
| equal / varying | spherical | --- | EII, VII |
| equal / varying | equal / varying | axis | EEI, VEI, EVI, VVI |
| equal / varying | equal / varying | equal | EEE, VEE, EVE, VVE |
| equal / varying | equal / varying | varying | EEV, VEV, EVV, VVV |

All three default to `NULL`, which follows `variance_model` and
`covariance_model`, so a call that names none of them fits exactly what it did
before. An ellipsoidal structure carries a full covariance array whatever
`covariance_model` said, because an orientation cannot live in a diagonal
block.

The estimates are those of Celeux and Govaert (1995). Ten have a closed form
or a scalar fixed point; EVE and VVE share one orientation across profiles
while letting the shapes differ, which has no closed form, and use the
minorize-maximize step of Browne and McNicholas (2014) --- majorize each trace
term linearly, then solve the resulting orthogonal Procrustes problem by
singular value decomposition.

**Verified against `mclust`.** Parameter counts match its own `df` for all
fourteen. Against `mclust::mstep()` on the same responsibilities, twelve agree
to 1e-11 or better. EVE and VVE differ by at most 5.7e-4 --- and attain a
*lower* value of the objective the M-step minimizes, by 2.1e-7 and 2.8e-5, so
the difference is `mclust` stopping short rather than a disagreement about the
estimate.

## Selecting over structures as well as class counts

`enumerate_classes(structure = )` crosses the covariance structures with the
class counts, which is the grid model-based clustering is usually selected
over, and adds a `structure` column to the candidate table:

```r
enumerate_classes(course_engagement, activity, id = "student",
                  n_profiles = 2:4, n_group_classes = 1,
                  structure = c("EEI", "VEI", "EVI", "VVI"),
                  centering = "person")
```

`candidate_fit()` gains `structure` to name one candidate out of a grid where
the class counts alone name several; asking without it raises
`multilpa_unknown_candidate` listing the structures it could have meant.

## Also new

The two structures whose M-step iterates over orientations now warm start from
the previous iteration's answer, which is worth about 1.2 to 1.3 times on them.
An M-step that stops at its iteration cap is **counted and reported once per
fit, with the count**, rather than warning from inside every M-step: one step
in three hundred is a different thing from all of them, and a warning raised
hundreds of times says neither.

`get_data(x, "assignments")` gains an `uncertainty` column: one minus the
posterior of the profile each row was assigned to, which is exactly what modal
assignment discards.

## The constrained diagonal covariance family

`variance_model` ties a profile's covariance volume to its shape: `"equal"`
constrains both and `"varying"` frees both. `volume` and `shape` free them
separately, which makes the six axis-parallel `mclust` models reachable
instead of two:

| `volume` | `shape` | model | spread parameters |
|---|---|---|---|
| equal | spherical | EII | 1 |
| varying | spherical | VII | K |
| equal | equal | EEI | d |
| varying | equal | VEI | K + d − 1 |
| equal | varying | **EVI** | 1 + K(d − 1) |
| varying | varying | VVI | Kd |

Each profile's covariance decomposes as `Sigma_k = lambda_k * A_k`, a volume
`lambda_k = |Sigma_k|^(1/d)` and a shape `A_k` with determinant one; the
estimates are those of Celeux and Govaert (1995). `volume`/`shape` default to
`NULL`, which follows `variance_model`, so every existing fit is unchanged.

Verified against `mclust`: the parameter counts match its own for all six
models, and at `mclust`'s maximum-likelihood estimate the log likelihood agrees
with one written from the mixture definition to 1e-10. The maximized
likelihoods are monotone along the whole nesting lattice.

Three consequences worth knowing:

* `parameter_inference()`, `vcov()` and `confint()` refuse EII, VII, VEI and
  EVI with `multilpa_unsupported_inference`. The free coordinates this package
  differentiates are log variances, one per profile and indicator, which is the
  wrong chart for a constrained volume or shape --- it has more coordinates than
  the model has parameters. EEI, VVI, EEE and VVV are unaffected.
* A constrained structure is maximized across every profile at once, so it
  cannot be combined with a held `variances` block; `fixed = "means"` is fine.
* A `start` taken from a wider structure is projected onto the requested family
  before the first iteration, because EM only promises a likelihood that does
  not decrease *within* the family it is maximizing over.

`covariance_model = "full"` constrains neither volume nor shape, so naming
either argument with it raises `multilpa_bad_argument`. Six of `mclust`'s
fourteen models --- the ones with a non-axis-parallel orientation --- remain
unavailable.

# multilpa 0.11.2

## One verb for every table

Eleven table verbs and the `what =` argument of seven `as.data.frame()`
methods are **replaced by one generic**, `get_data(x, what = )`. A reader no
longer has to know which verb owns which table before they can ask for it:

```r
get_data(fit)                            # the primary table, the measurement model
get_data(fit, "entropy")                 # was entropy_table(fit)
get_data(fit, "transitions")             # was transitions(fit)
get_data(fit, "profile_probabilities")   # was as.data.frame(fit, what = "...")
names(get_data(fit, "all"))              # every table this fit can produce
```

`what = "all"` returns a named list of every table, built from the same
definitions a single `what` uses, so the two cannot disagree. A table this
particular fit cannot produce -- sequences for a fit made without `time`,
bivariate residuals for a family with no discrete group classes -- is left out
rather than erroring, so the names of the list say what was available.

`summary()` now carries every table and prints all of them, truncated to
`rows = 10` each with the `get_data()` call that returns the rest.
`as.data.frame()` is plain coercion to the primary table; it no longer takes
`what`, and an argument it cannot use raises `multilpa_bad_argument` naming it
rather than quietly returning a different table.

### Migration

| Was | Now |
|---|---|
| `assignments(fit, data)` | `get_data(fit, "assignments", data = data)` |
| `average_posteriors(fit)` | `get_data(fit, "average_posteriors")` |
| `bch_weights(fit)` | `get_data(fit, "bch_weights")` |
| `bivariate_residuals(fit, data)` | `get_data(fit, "residuals", data = data)` |
| `classification_errors(fit)` | `get_data(fit, "classification_errors")` |
| `classification_table(fit)` | `get_data(fit, "classification")` |
| `entropy_table(fit)` | `get_data(fit, "entropy")` |
| `information_criteria(fit)` | `get_data(fit, "information_criteria")` |
| `sequences(fit)` | `get_data(fit, "sequences")` |
| `sequence_summary(fit)` | `get_data(fit, "sequence_summary")` |
| `transitions(fit)` | `get_data(fit, "transitions")` |
| `as.data.frame(fit, what = "x")` | `get_data(fit, "x")` |
| `as.data.frame(summary(fit), what = "fit")` | `get_data(fit, "model")` |
| `as.data.frame(diagnostics(fit), what = "posteriors")` | `get_data(diagnostics(fit), "average_posteriors")` |
| `fit_covariates(...)` | `multilpa(..., profile_covariates = )` |

Three tables that could only be reached through a `summary()` object are now
on the fit itself: `"counts"` (effective class memberships), `"covariances"`
(the within-profile residual covariance matrices) and `"model"` (the one-row
fit summary, previously `what = "fit"`). An enumeration grid gains
`"criteria"`, the per-criterion minima its summary already printed.

### Two defaults changed

`"classification"`, `"average_posteriors"`, `"classification_errors"` and
`"bch_weights"` default to **both levels** on a fit that has discrete group
classes, where the old verbs defaulted to `"individuals"`. This is the rule
`diagnostics()` already followed. Pass `level = "individuals"` for the old
behaviour. `"classification_errors"` and `"bch_weights"` accept
`level = "both"` for the first time; weight a regression with one level, not
with both.

`as.data.frame()` on a fitted transition model returns the **measurement
model**, not the transition matrix, because the primary table is now the same
kind of thing for every fitted family. `get_data(fit, "transitions")` is the
transition matrix.

## Covariates are arguments to `multilpa()`, not a separate verb

`fit_covariates()` is **removed**. Membership covariates are part of the model
`multilpa()` fits, so they are arguments to it:

```r
multilpa(course_engagement, activity, id = "student", n_profiles = 2,
         n_group_classes = 2, profile_covariates = "previous_grade")
```

The result is a `multilpa_covariates` object exactly as before, and it records
the `multilpa()` call the caller wrote. The covariate likelihood has no
observed-data form, no starting-value contract of the shape the covariate-free
EM uses and no held-measurement machinery, so `start`, `missing = "fiml"` and
`fixed` are refused by name with `multilpa_bad_argument` rather than accepted
and ignored.

## Every plot a fit supports, in one call

`plot(x, what = "all")` draws every view the fit has the ingredients for,
re-issuing the caller's own call once per view so style arguments carry into
each panel. A view that refuses is named at the end rather than stopping the
sequence.

## Recovery is one call

`get_data(x, "assignments", truth = )` cross-tabulates the model's labels
against known ones, so checking a fit against a generating truth no longer
needs `xtabs()` on the assignment frame:

```r
get_data(fit, "assignments", data = course_engagement,
         truth = c("engagement", "student_type"))
```

Each truth column is compared against the level it describes, read off the
data rather than guessed: a column taking one value within every group is a
property of the group and goes against `group_class`, and one that varies
inside any group goes against `profile`. The `assignment` column records which
comparison was made.

## A fitted model prints its own estimates

`print()` on any fitted model now shows the measurement model under the header,
so `fit` alone is the whole result. `fit` followed by `as.data.frame(fit)` was
two calls showing one thing.

```r
fit <- multilpa(course_engagement,
                c("browse", "lectures", "forum_read", "forum_post", "attendance"),
                id = "student", n_profiles = 2, n_group_classes = 2)
fit
#> Two-level latent profile analysis: 2 profiles, 2 group classes
#> ...
#>  profile  indicator  mean variance standard_deviation
#>        1     browse 2.854   0.3521             0.5934
#>        ...
#> Every other table: get_data(x, what = ), or get_data(x, "all").
```

An all-categorical fit prints its response probabilities instead, rather than
an empty block under a heading promising means. `rows` bounds the printout;
`get_data()` returns any table whole.

## Also

* `get_data(x, "model")` renames the plain fit's `bic` column to `bic_groups`,
  which is what the covariate and random-intercept families already called it,
  so one verb no longer spells the same convention two ways. The table is
  family-shaped rather than a padded union: a random-intercept fit reports its
  integration diagnostics and a covariate fit its predictor counts.
  `"information_criteria"` is the table with one column set for every family,
  and is the one to compare fits across families with.
* `report()` gains `rows`, forwarded to `print(summary(x))`, because a first
  look at a fit with thousands of observations would otherwise be mostly
  posteriors.
* A covariate model with **no** covariates can no longer be requested: naming
  no covariate is a request for the covariate-free model. The two differ only
  in parameterisation, the intercept-only multinomial logits standing in for
  the mixing probabilities.
* `report()` and `diagnostics()` are unchanged in purpose. `diagnostics()`
  reports its tables through `get_data()` and renames its average-posterior
  table from `"posteriors"` to `"average_posteriors"`, which no longer
  collides with the individual posteriors.
* `multilpa_removed_argument` is **removed** from the condition catalogue. It
  could only be raised by `classification_table(detail = )`, and that verb no
  longer exists, so nothing can raise it. `?"multilpa-conditions"` no longer
  claims otherwise.
* A three-step method called on a random-intercept fit now refuses with
  `multilpa_no_group_classes` instead of an unclassed error, so
  `get_data(x, "all")` can tell a table that does not apply from a defect.

# multilpa 0.11.1

## One bundled dataset, covering every model family

`school_engagement` and `engagement_panel` are **removed** and replaced by
`course_engagement`: 1422 enrolments of 106 students across 32 courses, ordered
within each student. The same rows now carry every model the package fits.

```r
activity <- c("browse", "lectures", "forum_read", "forum_post", "attendance")

multilpa(course_engagement, activity, id = "student",
         n_profiles = 2, n_group_classes = 2)          # two-level
fit_transitions(course_engagement, activity, id = "student",
                time = "sequence", n_profiles = 2)      # the same students, in order
fit_covariates(course_engagement, activity, id = "student",
               n_profiles = 2, n_group_classes = 1,
               profile_covariates = "previous_grade")   # and what predicts it
```

The old pair could not do this. `school_engagement`'s `term` column was not
time -- it indexed twelve different students within a school -- so anything
about order needed the second dataset and a different set of people, and
neither dataset carried a covariate at all.

Two new columns make that possible. `previous_grade` is the standardised grade
the student earned in the course *before* this one, so it predicts engagement
rather than summarising it. `attendance` is the count of days the student was
active, following the same engagement as the click measures but recorded one
tick per day rather than one per click, so it carries more noise than they do.

The generating truth still ships beside the observations -- `engagement` for the
enrolment pattern and `student_type` for the kind of student -- so a fitted
model can still be checked against what produced it:

```r
xtabs(~ profile + engagement, data = assignments(fit, data = course_engagement))
```

The truth column is called `engagement` rather than `profile` because
`assignments()` adds a column named `profile`, and a dataset already carrying
that name would make the verb refuse to join its own output onto it.

The data are simulated. The simulation is calibrated to the shape of a real
learning-analytics export, using only aggregate constants from it: marginal
means and spreads on the log scale, the separation between an engaged and a
disengaged pattern, and how persistent that pattern is across a student's
courses. No row, identifier or value of any real student is present, and the
source is not distributed. See `data-raw/course-engagement.R`.

The activity indicators are natural-log counts, `log1p(events)`. Raw
learning-analytics counts are strongly right-skewed and on an unlogged draft the
skew alone looked like an extra class, exactly as a floor effect does.

This release is a correctness release. An independent review of 0.10.0 found
that several advertised paths returned plausible numbers that were wrong, rather
than failing. Every one of them is fixed below, each with a regression test that
reproduces the original defect. Nothing here adds a model family.

## Silently wrong results, now correct

### A covariate's units no longer change the fit

`fit_covariates()` ran its membership M-step on the coefficients as supplied.
BFGS stops on a *relative* tolerance measured in the coordinates it is given, so
a predictor in small units stalled the search at a point where the score was
still large, and `optim()` reported success. Rescaling one covariate by `1e-7`
moved the log likelihood from `-620.1281` to `-730.5748` and drove the slope to
zero, with no warning. The search now runs on root-mean-square-scaled
coordinates and convergence is judged by a unit-free score, not by the
optimiser's return code. Estimates are returned on your scale, unchanged:

```r
# identical to machine precision at every scale
fit_covariates(data, "y", "g", 2, profile_covariates = "z")
```

### Staged categorical measurement attaches to the right item

A held measurement was passed between stages as an unlabelled, positional list,
while the compatibility check compared only the *sorted set* of indicator names
-- which a permutation passes. Reordering `categorical = c("a", "b")` to
`c("b", "a")` swapped whole response distributions (`-587.5909` to `-671.5937`);
reversing a factor's levels swapped categories within an item (`-1047.9689`).
`starting_values()` now keeps the indicator names and category labels, and
`multilpa()` aligns on them. An encoding that cannot be aligned raises
`multilpa_bad_start` or `multilpa_bad_stage` instead of returning a number.
An unlabelled start is still read positionally, so stored starts keep working.

### Bootstrap comparisons keep the constraints they were given

`bootstrap_lrt()` accepted fixed-measurement models but dropped `fixed` and the
held values from every refit, so the observed statistic compared constrained
fits while the null distribution came from unconstrained ones -- a 1-versus-3
parameter statistic referred to a 9-versus-11 parameter null. Refits now carry
the constraint, and a pair whose constrained nesting cannot be established is
refused with `multilpa_bad_nesting`.

### Transition counts no longer overflow

Expected transition counts exponentiated a scaled forward factor and an unscaled
backward factor separately, so on a long sequence one overflowed to `Inf` where
the other underflowed to `0`. A converged 600-occasion fit reported `NaN` for
every count, and fitting 1,800-occasion sequences failed outright. The complete
pair log probability is now formed before exponentiating, which cannot overflow.
Counts agree with explicit path enumeration to `1.8e-15`.

### `descriptives()` no longer invents or replaces values

An `NA` stratifier produced `NA` row selectors, and `x[NA, ]` fabricates a row
rather than dropping one: every stratum gained a phantom missing observation
while the real row vanished. Missing strata are now a labelled stratum, and the
per-stratum counts are asserted to account for every input row. Separately, the
profile assignment was written into the same frame as the indicators, so an
indicator named `profile` or `group_class` was overwritten by its own class
labels. The assignment is now kept beside the frame, never in it.

### Cluster-robust standard errors are refused when they cannot exist

`three_step()` and `r3step()` summed estimating-equation contributions within
groups without requiring enough independent groups. Those contributions sum to
zero at the estimate -- that sum *is* the stationarity condition -- so a single
group gave a standard error of `5.9e-17` and a confidence interval of zero
width, and `r3step()`'s robust path reported p-values of exactly zero. Both now
require more independent groups than the quantities reported, raising
`multilpa_too_few_groups`. `three_step(vcov_type = "independent")` offers an
explicit, labelled unclustered alternative.

### `diagnostics()` no longer ignores the argument you passed

`diagnostics()` and `report()` documented that `...` reached the verbs they
gather, but discarded it. Because `by` is a real argument of
`bivariate_residuals()`, `diagnostics(fit, by = "overall")` silently returned
per-profile residuals under the name you asked to pool. Both verbs now take `by`
explicitly and refuse any argument they cannot forward.

### A bootstrap comparison stops withholding its p-value

`bootstrap_lrt()` judged a replicate invalid if its statistic fell below a fixed
`-1e-5`. But the statistic is `2 * (ll_alt - ll_null)` on log likelihoods of
order 10^3 to 10^4, while EM stops on a *relative* tolerance, and under the null
the alternative converges to the null solution -- so every replicate sits at the
boundary with noise of about `2 * tol * |log likelihood|`. On
`course_engagement` that is 1.6e-04, sixteen times the window. Replicates that
had all converged were reported as "Nonconvergence or reversed likelihood" and
the p-value withheld. The window now scales with `tol` and the likelihood, so
the same comparison returns every replicate valid and a p-value.

### An ordinary converged fit stops being told it has not converged

`parameter_inference()` cut a fixed `0.01` on the scaled score. On 720
observations at the default `tol = 1e-8` that fires on a fit which converged
with all 20 starts replicating. The criterion is now the displacement the score
implies, expressed in the estimate's own standard errors: a score `g` against
information `I` moves the estimate by `g / I`, and the standard error is
`sqrt(1 / I)`, so `g * SE` is that displacement and is dimensionless. The
default-tolerance fit sits 0.17% of a standard error from stationarity and is
now quiet; a fit at `tol = 1e-4` sits 5% away and still warns, now saying so in
those words.

### Misaligned data is refused, not scored

`assignments()` documented that it checked row alignment but compared only the
row count, and `bivariate_residuals()` trusted row order. A reversed frame was
paired with the wrong posteriors and returned as a result. Both now compare
every column the fit recognises and raise `multilpa_bad_inference_data` on a
mismatch; a frame sharing no column with the fit warns
`multilpa_unverified_alignment` rather than passing in silence.

## Now supported

* `parameter_inference()`, `vcov()` and `confint()` work on fixed and staged
  fits, reporting standard errors conditional on the held measurement -- which
  `?fit_staged` already promised. They previously failed with
  `length(theta) == object$n_parameters is not TRUE`. Held blocks enter the
  likelihood at their held values and contribute no row to the information
  matrix; the conditional likelihood matches an independent enumeration to
  `5.7e-14` and the standard errors match an independent numerical information
  matrix to `1.3e-07` relative.
* `assignments()` and `descriptives(by = "profile")` work on a random-intercept
  fit. Both were documented to, but the fit stored no modal assignment, so the
  first failed with an unclassed `cbind()` error.
* `multilpa_plot_types()` lists `"random_intercepts"`, which
  `plot.multilpa_random_intercept()` has always accepted.
* `fit_transitions(max_iter = 0)` returns a usable fit instead of failing with
  "attempt to set an attribute on NULL".
* `plot()` draws `what = "entropy"` and `what = "posteriors"` for covariate and
  random-intercept fits, so `plot(diagnostics(fit))` works for every family.

## Behaviour changes

* `max_iter = 0` with a supplied `start` now evaluates that start alone and
  ignores `n_starts`, which is what `?multilpa` promised. It previously scored
  all random starts and returned the best of them.
* Restart selection now breaks ties deterministically rather than by
  `which.max()`, which was choosing between equally valid optima on
  floating-point noise.
* A separated covariate fit is now reported as unconverged rather than
  converged, so `parameter_inference()` refuses it instead of returning Wald
  intervals around an unbounded coefficient.
* Every warning the package raises now carries a documented condition class, so
  a simulation loop can muffle the qualification it expects and let the rest
  through. New classes: `multilpa_extreme_coefficients`,
  `multilpa_quadrature_check`, `multilpa_empty_transition_row`,
  `multilpa_unverified_alignment`, `multilpa_no_free_parameters`,
  `multilpa_held_parameter`.
* `three_step()` gains `vcov_type`; `diagnostics()` and `report()` gain `by`;
  `plot.multilpa_covariates()` and `plot.multilpa_random_intercept()` gain
  `"entropy"` and `"posteriors"`. All are additive and keep positional calls.
* `multilpa_bootstrap_lrt` objects carry `fixed`, and
  `as.data.frame(what = "test")` gains a `fixed` column.
* `parameter_inference()` on a fixed fit reports no row for a held block and
  carries a `fixed` attribute. The table is on the natural scale, where every
  probability in a set is reported, so it has `n_parameters` rows plus one for
  each reference category the estimation scale drops. `vcov()` is
  `n_parameters` square on the unconstrained scale; on the natural scale it is
  correspondingly larger and singular.
* `descriptives(by = )` can return one additional row: the labelled missing
  stratum.
* `profile_prevalence` for a transition fit is read from the final expectation
  rather than the last M-step. These agree at convergence; the new value is the
  correct one for an unconverged or zero-iteration fit.
* `starting_values()` keeps names and column labels on categorical response
  blocks, and `as.data.frame(what = "responses")` now reports those labels
  rather than positions, matching what the fitted object's own response table
  reports for the same quantity. A start carrying no labels is still read
  positionally and its tidy view still reports positions.
* `assignments()` raises `multilpa_bad_inference_data` rather than
  `multilpa_bad_nesting` for a row-count mismatch, which is a data contract
  failure, not a model-nesting one.
* Thirteen further user-reachable guards in `bootstrap_lrt()` and
  `fit_transitions()` now raise their catalogued class instead of a bare
  `simpleError`, so `?"multilpa-conditions"` is true of them.
* A `multilpa(fixed = )` fit records `n_parameters_with_measurement`, which only
  `fit_staged()` set before. `print()` consequently reported the same count
  twice on such a fit; it now reports the two counts it names.
* `parameter_inference()` results carry `score_displacement` beside
  `scaled_score`, the quantity the convergence warning is judged on.
* `bivariate_residuals()` gains `adjust`, and its table a `p_adjusted` column.
  Every pair of indicators is a separate test, so a five-indicator fit asks ten
  questions at once; the family size is now the verb's to declare rather than
  the reader's to reconstruct. The default is `"none"`, so `p_adjusted` equals
  `p_value` unless asked otherwise.

## Infrastructure

* GitHub Actions check the package on macOS, Windows and Linux across release,
  devel and oldrel, plus an explicit R 4.1 job for the declared `Depends` floor;
  a second workflow runs `lintr`.
* `DESCRIPTION` gains `URL` and `BugReports`.

# multilpa 0.10.0

## New features

### A fit carries the data it was built from

Post-fit verbs no longer ask for data the fit is already holding. The model
frame is rebuilt from the uncentred indicators, the categorical codes and their
levels, the identifier and the occasion — all of which were stored already — so
this costs no extra memory, and the round trip is exact for every model family.

```r
parameter_inference(fit)          # was parameter_inference(fit, students)
bivariate_residuals(fit)
vcov(fit); confint(fit)
bootstrap_lrt(smaller, larger)
as.data.frame(fit, what = "data") # the columns the model was fitted to
```

Supplying `data` still works and is still the stronger check, since a frame that
does not reproduce the fitted likelihood is rejected as before.
`three_step()` and `r3step()` keep requiring `data`, because an outcome or a
covariate is a column the model never saw. The `multilpa_data_required`
condition is retired: there is now something sensible to default to.

### Three convenience verbs

* `descriptives(x, vars, id)` — one row per variable with n, missing, mean, sd,
  range, and the **intraclass correlation**: the share of variance lying between
  groups. An ICC near zero says the groups do not differ, so a model built to
  tell them apart has nothing to find. It accepts a data frame before fitting or
  a fit afterwards, and `by = "profile"` splits it by assigned profile.
* `diagnostics(x)` — entropy, classification quality, average posteriors and
  bivariate residuals in one call, as a classed object with `print()`,
  `as.data.frame(x, what = )` and `plot()`. `plots = TRUE` draws them too. A
  model family that cannot supply a table says so rather than failing.
* `report(x)` — summary, descriptives, diagnostics and every plot the fit
  supports, in one call. It skips views a given fit has no ingredients for
  rather than stopping at the first refusal.

### `assignments()`

Every observation with the profile it was assigned to, the group class, and its
posteriors -- optionally with a frame of your own columns kept alongside:

```r
labelled <- assignments(fit, data = course_engagement)
xtabs(~ profile + engagement, data = labelled)
```

Comparing an assignment with anything else means putting them in the same row,
and doing that with two separate objects assumes they share an order. This verb
owns the alignment: a frame with the wrong number of rows is refused, and a
column the assignments would overwrite is an error rather than a replacement.

It is a verb rather than another `as.data.frame(what = )` value because it
computes a join with a frame the caller supplies. The line: `as.data.frame()`
represents what a fit already holds, and anything bringing in columns the model
never saw gets its own verb with its own documented arguments.

### `plot(fit, what = "bars")` draws its intervals without being handed data

The fit carries the columns it was built from, so the 95% intervals are drawn by
default; `data` is now only an override. A model family whose standard errors
are not implemented gets bars without whiskers rather than an error.

`report()` likewise asks each `plot()` method what views it accepts instead of
assuming the full catalogue, and skips a view the method refuses -- naming what
it could not draw rather than dropping it silently.

## Bug fixes

* The intraclass correlation is now `1` rather than `NA` when a variable has no
  within-group variance. That is a well-defined case — every observation equals
  its group's own value — not a degenerate one.

# multilpa 0.9.0

## Breaking changes

### `information_criteria()` returns the shape people report

The default is now one row, with one column per criterion — the shape a
model-comparison table is published in, and the same column names
`enumerate_classes()` uses, so a chosen model's row sits under the grid it was
chosen from without renaming anything.

```r
information_criteria(fit)
#>   log_likelihood n_parameters      aic      kic bic_groups bic_individual ...
#> 1      -3778.259           15 7586.518 7604.518   7617.933       7655.207 ...

information_criteria(fit, format = "long")   # the previous shape
```

* The `penalty` column has been removed. No convention reports it — not Mplus,
  not `tidyLPA::get_fit()`, not mclust — and it is the difference between two
  numbers the table already carries.
* `definitions = TRUE` describes one criterion per row, so it now requires
  `format = "long"` and raises `multilpa_bad_argument` otherwise.
* `as.data.frame(fit, what = "information_criteria")` follows the verb, and its
  `format` argument now actually reaches it. Previously `format` was captured by
  the accessor's own signature and never forwarded.
* `convention` remains `NA_character_` for `deviance`, `aic` and `kic`, now
  documented as meaning *the question does not arise* rather than *a value is
  missing*.

## New features

* Two example datasets ship with the package, so the first example runs as
  written: `school_engagement` (720 students in 60 schools) and
  `engagement_panel` (120 students over four waves). Both are simulated and both
  carry the kind each row was generated from — `engaged` and `state` — so a
  fitted model can be checked against what produced it.
* `vignette("multilpa")` is rebuilt on those datasets and gains a section
  cross-tabulating the fitted profiles against the truth column.

# multilpa 0.8.0

## Breaking changes

### Argument names now follow the sibling package `Nestimate`

No deprecation shim: the old spellings are gone and raise
`unused argument`.

| was | is | where |
|---|---|---|
| `object` | `x` | first formal of this package's own verbs |
| `indicators` | `vars` | `multilpa()`, `fit_covariates()`, `fit_staged()`, `fit_transitions()`, `fit_random_intercept()`, `enumerate_classes()` |
| `group` | `id` | the same six |
| `level_ci` | `ci_level` | `three_step()`, `r3step()` |
| `p_adjust` | `adjust` | `parameter_inference()`, `three_step()`, `r3step()` |
| `n_boot` | `iter` | `bootstrap_lrt()` |
| `profiles`, `group_classes` | `n_profiles`, `n_group_classes` | `enumerate_classes()` |

```r
# before
multilpa(students, indicators = vars, group = "school", n_profiles = 2)
# after
multilpa(students, vars = vars, id = "school", n_profiles = 2)
```

The S3 methods for the base generics `summary()`, `coef()`, `vcov()`,
`confint()`, `logLik()` and `nobs()` keep `object`, because those generics own
that name.

**Tidy output columns did not change.** Every posterior, sequence and
classification table keeps its `group` column, `level` keeps the value
`"group"`, and fitted objects keep `group_classes` and `n_group_classes`.

## New features

* `multilpa_plot_types()` lists every view the `plot()` methods accept, with the
  group each belongs to and what it answers. It carried `@export` in 0.7.0 but
  was missing from `NAMESPACE`, so it could not be called.
* `plot()` gains four views: `"bars"` (profile means as grouped bars, with 95%
  intervals when `data` is supplied), `"heatmap"` (means as deviations from each
  indicator's grand mean), `"entropy"` and `"posteriors"` (per-case
  classification quality as ridges).

## Bug fixes

* A refused plot no longer leaves the graphics device unusable. Restoring `mfg`
  as part of `par(no.readonly = TRUE)` switches `new` on as a documented side
  effect, so a view that validated and raised before drawing left every
  subsequent plot ready to overlay onto a blank panel.
* Direct labels are no longer clipped at the right edge on narrow panels. The
  gap was measured in user units while the margin reserved for it was measured
  in inches.
* Point area in `plot(fit)` encodes prevalence absolutely rather than being
  min-max normalised within the plot, which had mapped any spread onto the full
  size range and so carried no information.

# multilpa 0.7.0

* Every public verb returns a tidy `data.frame`. See `CHANGES.md` in the source
  repository for the full sweep.
