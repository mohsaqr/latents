# multilpa

A dependency-light R package for **two-level latent profile and latent class
analysis**. All estimation runs in R; Carm, Mplus, and a server are not required.

Indicators may be continuous, categorical, or a mixture of the two. The model
identifies individual profiles from the indicators and latent group classes from
differences in profile proportions. For example, course enrolments
have engagement profiles, while students differ in the prevalence of those
profiles. The implementation uses nested expectation-maximization (EM) with
multiple starts.

## Install and run

From this project directory:

```sh
R CMD INSTALL .
```

The package ships one example dataset, so the first example runs as written.
`course_engagement` has 1,422 enrolments: 106 students, 32 courses, up to
fifteen courses per student. It replaces the two datasets earlier versions
shipped, because one frame now carries what needed two. `student` is the
nesting unit, `sequence` orders each student's own courses, and
`previous_grade` is a covariate recorded before each enrolment, so the same
1,422 rows fit the two-level model (`id = "student"`), the latent transition
model (`time = "sequence"`) and the membership-covariate model
(`profile_covariates = "previous_grade"`) without changing dataset.

The data are simulated, and they carry the truth each row was generated from
beside the indicators — `engagement` for the pattern the row came from,
`student_type` for the kind of student — so a fit can be checked against what
produced it. The five indicators are `log1p` counts of learning-analytics
events.

```r
library(multilpa)

activity <- c("browse", "lectures", "forum_read", "forum_post", "attendance")

# Before fitting: does the nesting carry any signal at all? An ICC near zero
# says the students do not differ, so a model built to tell them apart has
# nothing to find.
descriptives(course_engagement, vars = activity, id = "student")

fit <- multilpa(
  data = course_engagement,
  vars = activity,
  id = "student",
  n_profiles = 2,
  n_group_classes = 2,
  variance_model = "varying",
  n_starts = 10,
  seed = 1
)

summary(fit)
as.data.frame(fit)                                   # one row per profile and indicator
as.data.frame(fit, what = "profile_probabilities")   # one row per group class and profile
as.data.frame(fit, what = "posteriors")              # one row per individual and profile
as.data.frame(fit, what = "posteriors", format = "wide")  # one row per individual
as.data.frame(fit, what = "group_posteriors")        # one row per group and group class
as.data.frame(fit, what = "starts")                  # one row per EM start
as.data.frame(fit, what = "data")                    # the columns the model was fitted to
assignments(fit, data = course_engagement)           # every row with the class it was given

# The generating truth ships beside the indicators, so recovery is one table at
# each level. Profile and group-class numbers are arbitrary, so read the table,
# not the diagonal.
xtabs(~ profile + engagement, data = assignments(fit, data = course_engagement))
xtabs(~ group_class + student_type, data = assignments(fit, data = course_engagement))

descriptives(fit)                                    # one row per variable, with the ICC
descriptives(fit, by = "profile")                    # the same, split by assigned profile
diagnostics(fit)                                     # every classification diagnostic
diagnostics(fit, by = "overall")                     # residuals pooled rather than per profile
report(fit)                                          # summary + diagnostics + every plot

information_criteria(fit)                            # one row, the reportable shape
information_criteria(fit, format = "long")           # one row per criterion and convention
classification_table(fit, level = "both")            # classification quality
entropy_table(fit)                                   # entropy at each level
average_posteriors(fit, level = "individuals")       # mean posterior by assigned class

multilpa_plot_types()                                # every view, with what it answers
plot(fit)                                            # profile means by indicator
plot(fit, scale = "standardized")                    # comparable profile shapes
plot(fit, what = "bars")                             # the same means, with 95% Wald intervals
plot(fit, what = "heatmap")                          # which indicators define the profiles
plot(fit, what = "probabilities")                    # prevalence by group class
plot(fit, what = "entropy")                          # where the uncertainty actually sits
plot(fit, what = "posteriors")                       # modal assignment probabilities
```

`descriptives()` answers the question it was asked first: the five indicators
have ICCs of 0.09 to 0.22 across students, so the nesting carries signal and a
two-level model has something to find. That fit then converges, all ten starts
reach the same likelihood, and the two profiles separate cleanly: relative
entropy is 0.946 at the enrolment level and 0.938 at the student level. The
`xtabs()` above puts 1,398 of 1,422 enrolments in the profile that generated
them, 98.3%, with the fitted profile 1 standing
for the generated `engaged` pattern — label switching is ordinary and is not
worth reordering anything to hide. The two group classes recover the two kinds
of student: one is 83% disengaged enrolments, the other 79% engaged. Pooled
over the whole sample 58.9% of enrolments are engaged, which describes neither
kind of student. Recovering that split is what the two-level model does and a
pooled one cannot.

Every result table is a base `data.frame`, so nothing needs to be pulled out of
the fitted object by hand.

`diagnostics()` and `report()` take `by` and nothing else: an argument they
cannot forward is refused by `multilpa_bad_argument` naming it, rather than
being dropped, so a `by` that was meant to reach the residual table cannot
silently answer a different question.

A post-fit verb reads the data the fit already carries, so `data =` is optional
for `assignments()`, `diagnostics()`, `report()`, `bivariate_residuals()`,
`parameter_inference()`, `confint()` and `vcov()`. It is still required by
`three_step()` and `r3step()`, which need the outcome or covariate columns the
fit never saw. When it is supplied it is checked column by column against the
fit: `assignments()` and `bivariate_residuals()` raise
`multilpa_bad_inference_data` on a frame that is not in fitting row order, and
warn `multilpa_unverified_alignment` on a frame that shares no column with the
fit and so cannot be checked at all.

`descriptives()` also takes `by` on a plain data frame. A stratifier with
missing values gets its own labelled `NA` stratum rather than being dropped, so
the per-stratum counts always add up to the number of rows you passed in.

Everything below runs on `course_engagement` as written, with one exception:
the bundled indicators are all continuous, so the categorical examples in the
next section use `survey`, a placeholder for data with categorical indicators of
your own. For a complete reproducible example that simulates its own data and
checks recovery against the values that generated it, run:

```sh
Rscript validation/synthetic-demo.R
```

The script checks recovery, prints the data and fitted structure, and saves
the data, true classes, and fitted model to `tmp/synthetic-demo.rds`. It loads
the local sources with `pkgload::load_all()`, so package installation is
optional for this example.

## Model

For person `i` in group `j`, the vector of continuous indicators has a Gaussian
profile `k`. By default, the indicators are conditionally independent within a profile.
The latent group class `h` determines the prior profile probabilities.

```text
P(Y_j) = sum_h eta_h * product_i [sum_k pi_hk * product_d
                                Normal(y_ijd | mu_kd, variance_kd)]
```

- `means[k, d]` and `variances[k, d]` are shared across group classes.
- `profile_probabilities[h, k]` describes profile prevalence within group class h.
- `group_probabilities[h]` describes prevalence across groups, giving each
  observed group one contribution, regardless of its size.
- Individual posteriors average over uncertainty in the group's latent class.
- `variance_model = "equal"` shares each indicator's variance across profiles;
  different indicators can still have different variances.
- `n_group_classes = 1` reduces to ordinary pooled LPA.

The formulation extends the nonparametric multilevel mixture framework of
[Vermunt (2003)](https://jeroenvermunt.nl/sm2003.pdf) to Gaussian indicators.
It assumes measurement invariance across group classes. Full residual covariance,
membership covariates, latent transitions over ordered occasions, and an
alternative continuous random-intercept model are available through the
extensions below. Random slopes are not implemented.

## Categorical indicators

Name any indicators that should be treated as categorical, and they are modelled
by unrestricted, profile-specific response probabilities instead of Gaussian
means and variances. That is the latent class measurement model, so the same
verb fits two-level LPA, two-level LCA, and mixtures of the two.

```r
# Two-level latent class analysis: every indicator categorical.
lca <- multilpa(survey, c("u1", "u2", "u3", "u4", "u5"), "school_id",
                  n_profiles = 3, n_group_classes = 2,
                  categorical = c("u1", "u2", "u3", "u4", "u5"),
                  n_starts = 40, seed = 42)

as.data.frame(lca, what = "responses")   # probabilities and logit thresholds
plot(lca, what = "responses")            # response curves by profile

# Mixed measurement: name only the categorical ones.
mixed <- multilpa(survey, c("reading", "maths", "u1", "u2"), "school_id",
                    n_profiles = 3, n_group_classes = 2,
                    categorical = c("u1", "u2"), n_starts = 40, seed = 42)

# Membership covariates carry over to a categorical or mixed measurement model.
fit_covariates(survey, c("reading", "u1", "u2"), "school_id", 3, 2,
               profile_covariates = "age", categorical = c("u1", "u2"))
```

`survey` here is a placeholder: these are the only examples in this README that
do not run on the bundled data, because `course_engagement` has no categorical
indicator to name.

Binary, ordinal and unordered indicators all use the same unrestricted
parameterization, which is what mixture software estimates by default when every
threshold is free. Numeric, integer, logical, character and factor columns are
accepted; factor levels keep their declared order, and everything else is
ordered by `sort()`, so the coding never depends on row order. Category counts
come from the observed values.

`missing = "fiml"` extends to categorical indicators: an unobserved cell simply
drops out of that individual's likelihood, which is the observed-data likelihood
under an ignorable mechanism, with no imputation.

Observed-information and robust standard errors cover categorical and mixed
indicators, and the parametric bootstrap supports complete categorical and mixed
data. `min_probability` bounds every response probability
away from zero, in the same way `min_variance` bounds variances, and the bound
is applied as the exact constrained solution rather than by rescaling.

## Additional methods

Below, `course_engagement` is the complete frame and `incomplete` is the same
study with one indicator unobserved, because `missing = "fiml"` is what the
first fit is there to show and several of these model families are
complete-data only. One `transform()` makes it: suppose `attendance` was only
logged for a student's first twelve courses, which leaves 161 of the 1,422 rows
without it.

```r
incomplete <- transform(course_engagement,
                        attendance = ifelse(sequence > 12, NA, attendance))

# Observed-data ML for incomplete indicators; correlated residuals within profiles.
full_fit <- multilpa(incomplete, activity, "student",
                     n_profiles = 2, n_group_classes = 2,
                     covariance_model = "full", missing = "fiml",
                     n_starts = 10, seed = 1)

# Observed-information standard errors and Wald intervals. The fit carries the
# data it was built from, so `data =` is optional everywhere below.
parameter_inference(full_fit)
confint(full_fit)

# MLR robust sandwich errors, accumulating the score over independent groups.
parameter_inference(full_fit, vcov_type = "robust")
vcov(full_fit, vcov_type = "robust")

# Likelihood-ratio statistic with the Lo-Mendell-Rubin adjustment. The
# `p_value` column is returned as NA on purpose, because the VLMR reference
# distribution is not reproduced; use the bootstrap below for a calibrated one.
smaller_fit <- multilpa(incomplete, activity, "student",
                        n_profiles = 2, n_group_classes = 2, missing = "fiml", seed = 1)
larger_fit <- multilpa(incomplete, activity, "student",
                       n_profiles = 3, n_group_classes = 2, missing = "fiml", seed = 1)
lmr_lrt(smaller_fit, larger_fit)

# One-step membership regressions. The measurement model is the one multilpa()
# fits: Gaussian, categorical or mixed, diagonal or full covariance. Covariate
# fits are complete-data only, so they take `course_engagement`, not
# `incomplete`; an incomplete frame raises `multilpa_bad_data`. Standard
# errors cover the Gaussian and full-covariance cases; a categorical fit
# refuses them with `multilpa_unsupported_inference` rather than understating
# its parameter count. `group_covariates =` takes the same kind of predictor at
# the student level; `previous_grade` varies course to course, so it belongs at
# the enrolment level and the bundled data carry no student-constant covariate.
with_predictors <- fit_covariates(
  course_engagement, activity, "student", 2, 2,
  profile_covariates = "previous_grade", seed = 1)
as.data.frame(with_predictors, what = "coefficients")   # both levels, one row per term
parameter_inference(with_predictors)                    # with standard errors and intervals

fit_covariates(course_engagement, activity, "student", 2, 2,
               profile_covariates = "previous_grade",
               covariance_model = "full", seed = 1)

# Alternative: one continuous group intercept, loading 1 on every indicator.
# Complete data only, like the covariate fits above.
random_intercept <- fit_random_intercept(
  course_engagement, activity, "student", 2, seed = 1)
summary(random_intercept)                            # includes the quadrature check
as.data.frame(summary(random_intercept), what = "model")

# Compare profile/group-class counts; inspect diagnostics in every row. Anything
# enumerate_classes() does not name itself reaches multilpa(), so
# `missing = "fiml"` would go here to enumerate on `incomplete` instead.
candidates <- enumerate_classes(
  course_engagement, activity, "student",
  n_profiles = 2:4, n_group_classes = 1:3, n_starts = 10, seed = 1)
as.data.frame(candidates)                            # one row per candidate model
summary(candidates)                                  # the best model on each criterion
candidate_fit(candidates, n_profiles = 2, n_group_classes = 2)
plot(candidates, criterion = "sabic_individual")

# Complete-data nested models differing by one class at one level. A held
# measurement is carried into every refit, and a pair whose constrained nesting
# cannot be established is refused with `multilpa_bad_nesting`.
# bootstrap_lrt(smaller_fit, larger_fit, iter = 199, seed = 42)

# Where conditional independence fails: residual association within profiles.
# `fit` is the diagonal model from the first example, which is the one that has
# something to confess.
bivariate_residuals(fit)
bivariate_residuals(fit, by = "overall")

# What those residuals ask for: the same two profiles and two group classes,
# with the association estimated rather than assumed away.
dependent <- multilpa(course_engagement, activity, "student",
                      n_profiles = 2, n_group_classes = 2,
                      covariance_model = "full", n_starts = 10, seed = 1)
information_criteria(dependent)
bivariate_residuals(dependent, by = "overall")

# Resume from a fitted solution, or score parameters produced elsewhere. A start
# describes one measurement specification, so the refit must be given the same
# `covariance_model` and `missing` as the fit it came from; use
# `starting_values(full_fit, covariance = "drop")` to warm-start a diagonal
# model from a full-covariance solution instead.
refit <- multilpa(incomplete, activity, "student",
                  n_profiles = 2, n_group_classes = 2, n_starts = 1,
                  covariance_model = "full", missing = "fiml",
                  start = starting_values(full_fit))

# max_iter = 0 performs no update, so logLik() evaluates the supplied start
# and nothing else: `n_starts` is ignored rather than scored and beaten.
evaluated <- multilpa(incomplete, activity, "student",
                      n_profiles = 2, n_group_classes = 2,
                      covariance_model = "full", missing = "fiml",
                      start = starting_values(full_fit), max_iter = 0)
logLik(evaluated)   # identical to logLik(full_fit)
```

Those residuals are not decoration on these data. `attendance` counts the days
a student was active in a course, and a day counts as active *because*
something was clicked, so `attendance` shares variance with the click measures
beyond what the profile explains. The diagonal model has no way to say that,
and `bivariate_residuals(fit, by = "overall")` reports it: the largest residual
correlation is 0.208 (`forum_read` with `attendance`, p = 1.7e-15), and five of
the ten pairs are significant at 0.05 — the four that pair `attendance` with a
click measure, plus `browse` with `forum_read` at p = 0.035. `dependent`, the
same two profiles and two group classes with `covariance_model = "full"`,
raises the log likelihood from -7618.51 to -7494.69 for twenty more parameters,
which AIC (15283.03 to 15075.39) and group-count BIC (15344.29 to 15189.92)
both pay for, and leaves no residual above 1.5e-05. Local dependence shows up
in the enumeration above too, where the diagonal criteria keep improving out to
four profiles and three group classes: unmodelled residual association is one
of the things extra classes get recruited to absorb. Read the residual table
before reading the grid.

| Capability | Supported scope |
|---|---|
| Measurement model | Gaussian, categorical, or mixed, via `categorical`; categorical indicators use unrestricted profile-specific response probabilities |
| Missing indicators | `missing="fiml"`: observed Gaussian marginals and observed categorical responses, ignorable missingness assumption; no missing covariates |
| Residual covariance | `covariance_model="diagonal"` or `"full"`, shared or profile-specific via `variance_model`. These four combinations are mclust's EEI, VVI, EEE and VVV; the other ten mclust parameterizations, which constrain volume, shape and orientation separately, are not available |
| SEs/CIs | Observed-Hessian ML inference for Gaussian, categorical and mixed discrete models, including full covariance and FIML; membership-covariate inference for complete Gaussian models, with or without full residual covariance; and fixed or staged fits, whose standard errors are conditional on the held measurement. Not available for transition fits (`multilpa_no_inference`), random-intercept fits (`multilpa_no_inference`), covariate fits with categorical indicators (`multilpa_unsupported_inference`), or a fit that holds every parameter it has (`multilpa_no_free_parameters`), each of which refuses by condition class rather than returning a number |
| Robust SEs | `vcov_type="robust"`: Huber-White sandwich over independent groups, plus the MLR scaling correction factor. This is the same estimator Mplus `ESTIMATOR=MLR` defines, but the two have not been compared numerically; the retained Mplus comparison covers likelihoods, parameters and criteria only |
| Membership covariates | Numeric predictors at both levels; shared individual-profile slopes across group classes; Gaussian, categorical or mixed indicators with diagonal or full residual covariance; complete data only |
| Continuous random effect | One shared Gaussian group intercept with unit indicator loadings; complete diagonal model, no discrete group classes |
| Class enumeration | Grid of discrete models, every information criterion under every sample-size convention that applies to it, entropy, failed-fit and convergence diagnostics |
| Information criteria | AIC, BIC, SABIC, CAIC, AWE, ICL, CLC and KIC. The five that depend on a sample size are reported under both the group-count and individual-count conventions; AIC and KIC depend on none and are reported once; CLC is reported once per convention because its convention selects which level's classification uncertainty it penalizes |
| Classification quality | Modal and model-estimated class sizes, average posterior probabilities, odds of correct classification, relative entropy |
| Measurement estimates | `as.data.frame(fit)` returns the measurement model; passing `data =` as well is what asks for it with a standard error beside every estimate, in one call |
| LMR statistic | Likelihood-ratio statistic and the Lo-Mendell-Rubin adjustment; the `p_value` column is **always `NA`**, because the VLMR reference distribution is not reproduced |
| Bootstrap LRT | Parametric bootstrap preserving group sizes; complete discrete models differing by one class; held measurement blocks are carried into every refit and reported in a `fixed` column, and a pair whose constrained nesting cannot be established is refused with `multilpa_bad_nesting`; not Mplus TECH14 |
| Local dependence | Posterior-weighted bivariate residuals within profile or overall, with approximate unadjusted p-values; apply a multiplicity correction when comparing pairs |
| Three-step | Classification error matrix and BCH weights at either level; distal outcomes by BCH, proportional or modal assignment, with `three_step(vcov_type = "cluster")` or `"independent"`; R3STEP membership covariates with `r3step(vcov_type = "observed")` or `"robust"`. A cluster-robust request is refused with `multilpa_too_few_groups` when there are not more independent groups than reported quantities, because the contributions sum to zero at the estimate and the covariance would be singular |
| Sequences | Profile assignments in long or wide form, with group counts, sequence lengths and completeness by group class; observed transitions between assignments are not computed |
| Latent transitions | `fit_transitions()`: first-order homogeneous transition probabilities between profiles, measurement invariant across occasions, per-group-class initial distributions and transition matrices; Gaussian, categorical or mixed indicators, FIML, diagonal or full covariance, ragged sequences; no standard errors, enumeration or bootstrap |
| Staged estimation | `fit_staged()` and `multilpa(fixed=)`: hold means, variances or response probabilities at supplied values while membership is estimated; `fit_staged()` reports both parameter counts; first-stage uncertainty is not propagated |
| Warm starts | `starting_values()` round-trips any fitted solution, keeping indicator names and category labels on categorical response blocks, so a stage cannot attach a response distribution to the wrong item; an encoding that cannot be aligned raises `multilpa_bad_start` or `multilpa_bad_stage`. `max_iter = 0` evaluates a supplied parameter set without moving, and `fit_transitions()` accepts it too |
| Plots | Profile means (raw or standardized), grouped bars with Wald intervals, a standardized heatmap, categorical response curves, prevalence by group class, profile sequences, per-case entropy, modal posteriors, and any enumeration criterion; `multilpa_plot_types()` lists them. Entropy and posterior views also work for covariate and random-intercept fits. Base graphics only; a transition fit has no plot method and refuses with `multilpa_no_plot` |
| Conditions | The errors and warnings listed in `?"multilpa-conditions"` carry stable classes, so they can be caught by what went wrong. A few internal guards in `bootstrap_lrt()` and `fit_transitions()` still raise unclassed errors; match on class only where the catalogue documents one |

These are explicit model families, not every combination of Mplus options.
Random-intercept fits currently do not provide standard errors.
Random-intercept mixtures use Gaussian quadrature; increase node counts and
refit if the higher-order likelihood check fails, which warns
`multilpa_quadrature_check`. Inference rejects variance
boundaries and nonpositive information matrices. A fit whose likelihood still
carries a non-negligible score is qualified with a `multilpa_unconverged`
warning before any Wald quantity is reported, so a warning on
`parameter_inference()`, `confint()` or `plot(what = "bars")` is the package
telling you to refit with a tighter `tol` rather than an error.
Natural-scale Wald intervals
can extend beyond parameter bounds; log-variance/logit coordinates are also
available. A bootstrap p-value is withheld if any replicate fails convergence
or has a reversed likelihood. A covariate fit whose membership logits are close
to separated is reported as unconverged, so `parameter_inference()` refuses it
rather than putting a Wald interval around an unbounded coefficient.
Enumeration never automatically declares a winner.

## Profiles over time

Passing `time` records where each individual sits in an ordered sequence, so a
person measured repeatedly can be followed across profiles rather than reduced
to one assignment.

```r
over_time <- multilpa(course_engagement, activity, "student",
                      n_profiles = 2, n_group_classes = 2,
                      time = "sequence", seed = 1)

sequences(over_time)                  # one row per individual and time point
sequences(over_time, format = "wide") # one row per individual, carrying its identifier
sequence_summary(over_time)           # group counts, sequence lengths and completeness
plot(over_time, what = "sequences")
```

`sequence` is the position of a course in that student's own order, so no
argument beyond `time =` is needed and no separate panel dataset is: the rows
the cross-sectional fit used are already ordered within student. They are
ragged — the shortest student has ten courses and the longest fifteen, which
`sequence_summary()` reports per group class.

## Staged estimation: deciding the measurement model first

By default every parameter is estimated jointly, so adding group classes can
move the profiles those classes are meant to describe. `fit_staged()` estimates
the measurement model on its own and then estimates the group-class structure
with that measurement held fixed, so the profiles mean the same thing before
and after.

```r
staged <- fit_staged(course_engagement, activity, "student",
                     n_profiles = 2, n_group_classes = 2, seed = 1)

as.data.frame(staged, what = "stages")     # what each stage estimated and held
as.data.frame(staged, what = "profile_probabilities")
```

On these data the two routes nearly agree, which is the reassuring case rather
than the guaranteed one: the staged fit reaches -7618.67 against the joint
fit's -7618.51 with the same 23 parameters, so holding the measurement costs
almost nothing here and the group classes are describing the same two profiles
either way.

The result is an ordinary `multilpa` object, so every accessor, diagnostic and
method works on it unchanged. A measurement solution already fitted and
inspected can be carried in rather than refitted:

```r
measurement <- multilpa(course_engagement, activity, "student",
                        n_profiles = 2, n_group_classes = 1, seed = 1)
staged <- fit_staged(course_engagement, activity, "student",
                     n_profiles = 2, n_group_classes = 2, measurement = measurement)
```

For finer control, `multilpa()` takes `fixed` directly, naming any of `"means"`,
`"variances"` and `"response_probabilities"`, or `"measurement"` for every block
the model has. A held block stays exactly as supplied in every restart and stops
counting towards `n_parameters`, so this is a different model rather than a
different starting point for the same one.

```r
multilpa(course_engagement, activity, "student",
         n_profiles = 2, n_group_classes = 2,
         start = starting_values(measurement, what = "measurement"),
         fixed = "variances")
```

A fixed or staged fit is no longer refused by the inference verbs:
`parameter_inference()`, `vcov()` and `confint()` report the estimated
parameters only, conditional on the held measurement, which contributes no row
to the information matrix. A fit that holds *every* parameter it has raises
`multilpa_no_free_parameters`, and naming a held coordinate in `confint(parm =)`
raises `multilpa_held_parameter`, because a held value has no sampling
distribution.

```r
as.data.frame(staged, what = "stages")    # parameters and parameters_with_measurement
parameter_inference(staged)               # the estimated parameters, with their errors
```

First-stage uncertainty is **not** propagated: the second stage treats the
measurement solution as known, so its standard errors and information criteria
are conditional on that solution and are narrower than a joint fit's. Both
parameter counts are therefore reported by `fit_staged()`: compare staged fits
with one another on `parameters`, and a staged fit against a jointly estimated
one on `parameters_with_measurement`, remembering that the staged likelihood is
not the joint maximum. A fit built with `multilpa(fixed =)` directly reports
only the count it estimated here, so add the held blocks yourself when comparing
one of those against a joint fit.

## Latent transitions

`sequences()` reports where the model put each observation. It does not
estimate how observations move. `fit_transitions()` does: it fits the same
measurement model and, on top of it, the probability of moving from each
profile to each profile between consecutive occasions.

```r
moves <- fit_transitions(course_engagement, activity, "student",
                         n_profiles = 2, time = "sequence", seed = 1)

transitions(moves)                          # one row per ordered pair of profiles
as.data.frame(moves, what = "initial")      # where sequences start
as.data.frame(moves, what = "sequence_lengths")
```

Engagement is sticky in both directions here: an enrolment in the engaged
profile is followed by another with probability 0.86, and a disengaged one by
another disengaged with probability 0.81.

Profiles keep the same meaning at every occasion, because the measurement
parameters are shared across occasions, so a change of profile is a change of
state rather than a change of definition. Transitions are first order and
homogeneous: the probability of moving does not depend on the occasion or on
earlier profiles.

`n_group_classes` above one gives every group class its own initial
distribution and its own transition matrix, which separates groups that differ
in how they move from groups that differ only in where they start.

```r
mixture <- fit_transitions(course_engagement, activity, "student",
                           n_profiles = 2, n_group_classes = 2,
                           time = "sequence", seed = 1)
transitions(mixture)
```

Groups need not be observed at every occasion. `occasions = "observed"`, the
default, links each group's own consecutive observations. `occasions = "grid"`
places them on the grid of every position seen in the data, so a group that
skips a wave spends a transition crossing the gap and contributes no
measurement information at it. The two agree whenever every group is observed
at every position, which `sequence_summary(moves)` reports as a `gaps` count of
zero for every group class.

Standard errors, likelihood-ratio tests, class enumeration and plots are not
available for this model family — `parameter_inference()` and `vcov()` refuse
with `multilpa_no_inference`, and `plot()` with `multilpa_no_plot`. `logLik()`,
`information_criteria()`, `classification_table()`, `entropy_table()`,
`sequences()` and `sequence_summary()` are. `max_iter = 0` is accepted here too,
so a supplied parameter set can be scored without being moved.

## Relating classes to variables that did not define them

A covariate or outcome added to the measurement model can change the classes it
was meant to describe. The three-step approach fits the measurement model first,
then carries the classification and its error into a second stage, so the
classes stay fixed.

`previous_grade` is what makes this section runnable on bundled data:
`multilpa()` never sees it, so it is genuinely a variable the classes did not
define.

```r
# The classification error matrix, P(assigned | true), at either level.
classification_errors(fit)
classification_errors(fit, level = "groups")

# Bolck-Croon-Hagenaars weights, the inverse of that matrix by modal class.
bch_weights(fit)                                     # one row per unit and class

# A variable the classes did not define, corrected for misclassification.
three_step(fit, data = course_engagement, outcome = "previous_grade")
three_step(fit, data = course_engagement, outcome = "previous_grade",
           contrast = "pairs")
three_step(fit, data = course_engagement, outcome = "previous_grade",
           method = "modal")
three_step(fit, data = course_engagement, outcome = "previous_grade",
           vcov_type = "independent")               # unclustered, and labelled as such

# Covariates predicting class membership (Vermunt 2010 R3STEP), errors fixed.
r3step(fit, data = course_engagement, covariates = "previous_grade")
r3step(fit, data = course_engagement, covariates = "previous_grade",
       vcov_type = "robust")
```

`three_step()` defaults to `method = "bch"`, which is robust to the outcome
being unrelated to class membership. `"modal"` ignores classification error and
is biased toward the null; `"proportional"` weights by the posterior. `r3step()`
returns the multinomial logits with standard errors.

The two verbs answer different questions about the same column, which is why
both are shown on it. `three_step()` describes the classes: BCH puts mean
`previous_grade` at 0.221 in the engaged profile and -0.315 in the disengaged
one, a gap of 0.536 (SE 0.050). `r3step()` runs the regression the design
actually supports, because the grade was earned in the course *before* the
enrolment being classified: a one standard deviation higher previous grade
raises the log odds of the engaged profile by 0.574 (SE 0.062, 95% CI 0.453 to
0.695), and 0.574 (SE 0.057) with `vcov_type = "robust"`. The same covariate
estimated jointly with the measurement model by `fit_covariates()` gives 0.440
(SE 0.074) — a smaller coefficient, because there the covariate is also allowed
to move the profiles it is predicting.

`level = "groups"` is available to both verbs and needs a covariate that is
constant within a group; `previous_grade` varies from course to course, so
asking for it at the student level is refused with `multilpa_bad_outcome`
rather than silently averaged.

Both verbs default to a cluster-robust variance over the fit's groups, and both
**refuse it** with `multilpa_too_few_groups` when there are not more independent
groups than the quantities being reported. The cluster contributions sum to zero
at the estimate — that sum is the stationarity condition — so too few groups
leave the covariance singular and would hand back a standard error of
essentially zero. `three_step(vcov_type = "independent")` and
`r3step(vcov_type = "observed")` are the explicit, labelled alternatives; each
ignores the nesting, and must be reported as having done so.

## Plots

`plot()` methods are provided for fitted models, covariate fits,
random-intercept fits and enumeration grids. They use base graphics only; no
plotting package is added as a dependency. `multilpa_plot_types()` returns every
view with the question it answers; a latent transition fit has no plot method
and refuses with `multilpa_no_plot`.

Every series is distinguished by **colour, point symbol and line type together**
and labelled directly at the line end, so the plots stay readable in greyscale,
for colour-blind readers, and without a legend. Colours are Okabe-Ito.

```r
plot(fit)                                   # means in input units
plot(fit, scale = "standardized")           # shapes comparable across indicators
plot(fit, what = "probabilities")           # profile prevalence per group class
plot(fit, what = "entropy")                 # where the classification uncertainty sits
plot(fit, what = "posteriors")              # modal assignment probabilities
plot(with_predictors, what = "entropy")     # the same two views on a covariate fit
plot(random_intercept, what = "posteriors") # and on a random-intercept fit
plot(candidates, criterion = "icl_groups")  # any column of the enumeration grid
```

Use `scale = "standardized"` whenever indicators are on different scales,
otherwise the largest-scale indicator dictates the apparent profile shape. The
divisor is each indicator's observed standard deviation, not its within-profile
residual standard deviation, so the values are comparable but are not effect
sizes.

Plots show point estimates only and profile order is arbitrary. The one
exception is `what = "bars"`, which draws 95% Wald intervals and therefore
inherits the inference qualifications: on a fit whose likelihood still carries a
non-negligible score it warns `multilpa_unconverged` before drawing. In an
enumeration plot the ringed point marks the lowest criterion
value, which identifies an extremum and does not select a model; candidates that
failed to converge appear as crosses on the baseline rather than being dropped.

No visual constant is hard-coded. Pass any of them inline:

```r
plot(fit, palette = c("#0072B2", "#D55E00"),
     panel_fill = "#FFFFFF", point_size = 1.8, line_width = 3,
     main = "My title", subtitle = "My subtitle", labels = FALSE)
```

## Diagnostics and interpretation

The result retains input row order. Group order follows first occurrence in
the data; `group_values` preserves the original identifiers, and `group_ids`
provides unique display labels. Profile and group-class numbers are arbitrary:
align labels before comparing separate fits.

Row order is part of the effective seed. Initialization runs k-means on the
rows as supplied, so the same data in a different row order can converge to a
different local optimum with the same `seed`. The likelihood itself is exactly
row-order invariant given identical starting values; it is the search that
differs. Fix the row order alongside the seed when a fit must be reproducible,
and raise `n_starts` until `n_best_replicated` is comfortably above one.

Two tidy tables carry the search diagnostics, so none of this has to be pulled
out of the fitted object by hand:

```r
as.data.frame(summary(fit), what = "fit")     # converged, best start, n_best_replicated
as.data.frame(fit, what = "starts")           # one row per start, with its error if any
```

The best finite likelihood is returned even when
that start did not converge, with a `multilpa_unconverged` warning. Failed
starts are recorded with their error messages and warn `multilpa_failed_starts`.
More starts and replicated best likelihoods improve
confidence but do not establish a global maximum. Convergence is based on
relative likelihood change, not parameter stability or a Hessian test.

`min_variance` is an explicit lower bound in **squared input units**, defaulting
to `1e-6`. Gaussian mixtures with freely varying variances otherwise allow
singular solutions. The returned estimates solve a constrained likelihood
problem; the `boundary` column of `as.data.frame(summary(fit), what = "fit")`
and a `multilpa_boundary` warning identify a variance at the bound.
Choose a sensible bound for the indicator scales and examine sensitivity. The
`small_classes` column of the same table flags effective memberships below one,
warned as `multilpa_small_classes`; effective class counts are also printed by
`summary()`.

Missing indicators are rejected with `multilpa_bad_data` by default; use `missing="fiml"` to integrate
them out of the likelihood. Completely unobserved rows contribute no indicator
information, but are retained for posterior prediction. Entirely unobserved
indicators, infinite values, constant indicators, and missing group IDs are
rejected. Mean filling is used only to initialize incomplete-data fits.
Numeric, character,
and factor group identifiers are supported. Group-class identification requires
more than basic count checks: indistinguishable profile distributions or group
prevalence vectors can produce unidentified models. The package rejects some
obviously unidentified cases but does not certify general identification.

### Information criteria

Let `J` be the number of observed groups, `N` the number of individuals with at
least one observed indicator, and
`q` the number of free parameters:

| Field | Definition |
|---|---|
| `aic` / `AIC(fit)` | `-2 * logLik + 2 * q` |
| `bic` / `bic_groups` / `BIC(fit)` | `-2 * logLik + log(J) * q` |
| `bic_individual` | `-2 * logLik + log(N) * q` |

The BIC sample-size convention matters in multilevel mixtures.
[Lukočienė, Varriale, and Vermunt (2010)](https://doi.org/10.1111/j.1467-9531.2010.01231.x)
studied this issue. Both conventions are exposed; use one consistently when
comparing multilevel candidates fitted to the same observations and indicators.
For ordinary pooled-LPA enumeration (`n_group_classes = 1`), use
`bic_individual`; group-count BIC has zero penalty if there is only one observed
group. `nobs(fit)` and the `logLik` object's `nobs` attribute return `J`.

For full covariance, the variance constraint applies to every covariance
eigenvalue. New discrete fits retain `indicator_data` to verify the original
observations before inference and bootstrap comparisons.

## Vignette

`vignette("multilpa")` fits one simulated dataset end to end and prints every
diagnostic exactly as the package returns it: measurement model with standard
errors, prevalence by group class, all eight information criteria under every
convention that applies to them, classification quality, entropy, enumeration, Wald and robust
inference, bivariate residuals, classification error and BCH weights, mixed
indicators, membership covariates, staged estimation, sequences and latent
transitions. Because the data are simulated, every estimate can be read against
the value that generated it.

```r
vignette("multilpa")
```

## Verification

```sh
Rscript -e 'pkgload::load_all("."); testthat::test_dir("tests/testthat")'
R CMD build .
R CMD check --no-manual multilpa_0.11.1.tar.gz
```

`pkgload` is only a development convenience. Installed-package testing via
`R CMD check` does not require it. The package imports nothing outside base R:
`stats`, `graphics`, `grDevices` and `utils`. Everything in `Suggests` —
`testthat`, `mclust`, `tidySEM`, `knitr` and `rmarkdown` — is for testing and
the vignette, and every use of one is guarded.

Verification includes:

| Check | Observed agreement on retained synthetic examples |
|---|---|
| Exhaustive latent-assignment likelihood and posteriors | posterior differences below `1e-15` |
| Independent first EM update, both variance constraints | within `1e-9` |
| Direct `stats::optim()` BFGS likelihood optimization | likelihood difference `7.82e-11`; parameter differences below `6e-7` |
| `mclust` single-level VVI and EEI limits | likelihood differences below `6e-14` |
| Simulated parameter/class recovery | tested after aligning arbitrary labels |
| Held measurement parameters against the constrained M-step formula | exact; drift across a fit below `1e-12` |
| Latent transition likelihood against enumeration of every state path | differences below `3e-14` |
| Latent transitions against `depmixS4` | likelihood function `6.24e-16`; independently maximized likelihood `8.30e-15` |
| Input, dimension, numerical, restart, and RNG edge cases | covered by regression tests |

The `mclust` limit checks initialize at the reference solution to test that it
is a fixed point of the new engine; they do not compare global searches.
Independent likelihood enumeration and direct optimization validate the
multilevel calculations on tested examples.

The [mathematical audit](validation/MATH_AUDIT.md) documents numerical fixes,
independent regression checks, and the results of the latest full validation.

### Direct Mplus comparisons

The original six comparisons cover three datasets:

- Published Mplus examples **7.9 and 7.10**, using the same 500-person,
  four-indicator dataset: the two single-level variance models reproduce the
  published estimates within their displayed precision.
- New Mplus 9 Demo runs on **1,200 simulated individuals in 60 groups**:
  both two-level variance models agree on parameter estimates within `6e-8`
  and on marginal posterior probabilities within `2e-10`.
- New Mplus 9 Demo runs on the **public Example 10.4 dataset**, with 1,000
  individuals in 110 groups and five indicators: both matched two-level LPA
  models agree on parameter estimates within `1.1e-5` and posteriors within
  `3e-5`. These are new LPA fits; the original published CFA mixture is a
  different model and is not used as the numerical target.

A genuine `CATEGORICAL =` run validates the two-level latent class model on
1,200 individuals in 60 groups with five binary indicators: thresholds agree
within `5.8e-6`, profile probabilities within `1.1e-6`, group-class
probabilities within `8.3e-7`, and the likelihood within `2.8e-6`. The
single-level limit additionally agrees with `poLCA` to `2.7e-11` in likelihood
and `7.9e-8` in response probabilities, with identical free-parameter counts.

Additional genuine Mplus runs validate full covariance with missing indicators
(both equal and varying covariance), one-step membership covariates at both
levels, and the single-profile continuous random-intercept limit. Observed
information standard errors are also compared with Mplus.

A further genuine `ESTIMATOR = MLR` run validates the robust sandwich errors and
the new information criteria: robust standard errors agree within `1.2e-7`, the
MLR scaling correction factor within `7.1e-7`, and AIC, BIC and SABIC within
`4.1e-5`. Two genuine `TECH11` runs reproduce the Lo-Mendell-Rubin adjusted
statistic to the three decimals Mplus prints, end to end from native fits on the
published Example 7.9 data (`747.0633` against `747.063`, adjusted `723.7708`
against `723.771`). The Vuong-Lo-Mendell-Rubin reference distribution, that is
the mean, standard deviation and p-value Mplus reports beside it, is **not**
reproduced, and no p-value is returned in its place. See
[the extension report](validation/mplus/EXTENSIONS.md) for exact scope and
tolerances. The random-intercept mixture likelihood additionally agrees with
independent adaptive numerical integration; its Mplus comparison currently
covers one profile only.

All new two-level comparisons use independent starts in R and Mplus. Class labels and saved subject
IDs are aligned before comparison. Mplus's BIC is compared to
`bic_individual`, not this package's default group-count BIC. Saved Mplus
likelihoods and information criteria have finite output precision, accounted
for explicitly in the assertions. Tests run offline from retained fixtures;
neither Mplus nor internet access is needed to run them.

See [the comparison report](validation/mplus/COMPARISON.md) for numerical
tables, original URLs, raw Mplus artifacts, reproduction commands, and scope.
These checks establish agreement for the tested specifications, not full
Mplus feature parity or a guarantee for every dataset.

The independent reference implementation lives in
`tests/testthat/helper-independent-likelihood.R`; it uses explicit assignment
enumeration and direct probability arithmetic, independently of the fitting
engine's log-domain EM calculation.
