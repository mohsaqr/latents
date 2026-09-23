# multilpa

Multilevel latent profile analysis (MLPA) extends latent profile analysis to data in which observations are clustered within higher-level units, such as repeated measurements within persons, students within schools, patients within clinics, or employees within organizations. The central motivation is that observations from the same higher-level unit are generally dependent, so treating them as independent can confound within-unit heterogeneity with between-unit heterogeneity. A two-level mixture model addresses this structure by allowing Level-1 observations to belong to latent profiles while simultaneously modeling systematic differences among the Level-2 units that contain those observations. In this formulation, the lower level describes heterogeneity among observations, whereas the higher level describes heterogeneity among the units. MLPA is therefore a natural framework when the scientific question concerns both the latent structure of individual observations and the way that structure varies across clusters.

In repeated-measure applications, this distinction is particularly important. Let \(\mathbf{Y}_{it}\) denote the vector of indicators observed for person \(i\) at occasion \(t\), and let \(C_{it}\) denote a Level-1 latent profile. MLPA can estimate probabilities such as \(P(C_{it}=k)\) while accounting for the fact that multiple observations belong to the same person. A higher-level latent variable can then represent differences among persons in the distribution of these lower-level profiles. Thus, MLPA can distinguish diversity among repeated observations from systematic differences in the composition of those observations across people. The same logic applies to other nested designs: the Level-1 profile characterizes an observation, while the Level-2 structure captures dependence and heterogeneity associated with the unit in which observations are nested.

VASSTRA represents a related but distinct case in which the separation between **state diversity** and **person heterogeneity** is made explicit as the organizing principle of the analysis. Rather than treating the lower-level latent profiles and higher-level classification primarily as components of a single joint multilevel mixture, VASSTRA first establishes a common state space from the Level-1 observations. The first stage therefore asks: *what distinct states occur across observations?* Once these states have been estimated, the repeated observations belonging to each person are represented within that common state space. The second stage consequently asks: *how does each person occupy or distribute their observations across the states?* The resulting person-level representation can retain multiple states for the same individual rather than reducing the person immediately to a single latent category.

The multilpa R package provides an implementation of multilevel latent profile and related mixture models for clustered observations. It is designed to estimate latent profiles at the observation level while representing the higher-level structure of the groups containing those observations. In a typical application, the data consist of multiple Level-1 observations nested within Level-2 units, and the package provides estimation procedures for the corresponding two-level latent structure. Model fitting is likelihood-based and uses expectation-maximization estimation with multiple starting values to address the possibility of local solutions. The package also provides tools for comparing candidate solutions through information criteria, likelihood-based comparisons, classification diagnostics, uncertainty estimates, and diagnostics for residual dependence. These components allow the researcher to evaluate the latent structure rather than relying on the estimated profile assignments alone.

The multilpa package also provides a useful computational basis for examining alternative representations of hierarchical heterogeneity. Its model-fitting functions expose the estimated profile structure, posterior classification information, and model-comparison quantities needed to examine how many latent profiles and higher-level classes are supported by the data. The package is therefore useful both for conventional MLPA and for analyses in which the multilevel mixture provides a reference model for a more explicitly state-oriented approach such as VASSTRA. In the latter setting, MLPA can establish and evaluate a joint multilevel representation, while VASSTRA uses a sequential representation when the scientific objective is to identify the diversity of lower-level states first and then characterize people according to their distribution across that state space.

## Installation

```r
install.packages("multilpa")

# Development version
# install.packages("remotes")
remotes::install_github("mohsaqr/multilpa")
```

From a local clone of the repository:

```sh
R CMD INSTALL .
```

## Example Dataset and Methodological Structure

The package includes one example dataset, `course_engagement`, which is used throughout the examples. The dataset contains 1,422 course enrolments from 106 students across 32 courses, with each student contributing observations from up to fifteen courses.

The important point is that each row represents one student enrolled in one course. This gives the data a multilevel and longitudinal structure. Enrolments are grouped within students, and the courses taken by each student are ordered over time. Because all of the required variables are stored in the same data frame, the same 1,422 rows can be used for the two-level model, the latent transition model, and the membership-covariate model.

The variable `student` identifies which enrolments belong to the same student. It therefore defines the nesting structure for the multilevel model and is specified using `id = "student"`.

Within each student, the variable `sequence` records the order in which courses were taken. This ordering is used in the latent transition model through `time = "sequence"`, allowing profile membership to be examined across a student's sequence of course enrolments.

The variable `previous_grade` is recorded before each enrolment. Because it precedes the enrolment being modeled, it can be used as a covariate of profile membership through `profile_covariates = "previous_grade"`.

Taken together, these variables allow the same dataset to support several related analyses. The observations are nested within students through `student`, ordered within students through `sequence`, and accompanied by a pre-enrolment covariate through `previous_grade`.

The `course_engagement` data are simulated, which means that the data-generating structure is known. Two additional variables store this information alongside the observed indicators. The variable `engagement` records the pattern from which each row was generated, while `student_type` records the type of student from which it came. These variables make it possible to compare the profiles estimated by a model with the known structure used to generate the data.

The dataset also contains five learning-analytics indicators. These variables are event counts that were first transformed using `log1p` and then standardized within course. As a result, each value is interpreted relative to other enrolments in the same course. A value of zero corresponds to the course average, while positive and negative values indicate how many standard deviations an enrolment lies above or below that average.

## A first model

The first model estimates two engagement profiles at the enrolment level and
two classes of students at the student level. Before the fit, `descriptives()`
reports the intraclass correlation (ICC) of each indicator across students. An
ICC near zero would mean that students do not differ on that indicator, and a
model designed to separate students would have little to separate.

The last call compares the fitted classes with the columns that generated the
data. Each known label is compared with the level it describes: `engagement`
varies within a student and is compared with the enrolment-level profile,
whereas `student_type` is constant within a student and is compared with the
student-level class. Profile and class numbers are arbitrary labels, so the
table is read by its largest cells and not by its diagonal.

```r
library(multilpa)

vars <- c("browse", "lectures", "forum_read", "forum_post", "attendance")

descriptives(course_engagement, vars = vars, id = "student")

fit <- multilpa(
  data = course_engagement,
  vars = vars,
  id = "student",
  n_profiles = 2,
  n_group_classes = 2,
  variance_model = "varying",
  n_starts = 10,
  seed = 1
)

summary(fit)

get_results(fit, "assignments", data = course_engagement,
            truth = c("engagement", "student_type"))
```

 `descriptives()` provides a first check of whether the multilevel structure is worth modelling. For the five indicators, ICCs range from 0.09 to 0.22 across students, showing that the nesting contains meaningful between-student variation. The two-level model then converges successfully, with all ten starts reaching the same likelihood. Classification is also clear, with relative entropy of 0.938 at the enrolment level and 0.940 at the student level.

Because `course_engagement` is simulated, the fitted classes can be compared with the classes that generated the data. At the enrolment level, 1,389 of 1,422 observations, or 97.7%, are assigned to the correct profile. Fitted profile 2 corresponds to the generated `engaged` pattern. The numerical labels themselves are arbitrary, so this label switching does not need to be corrected.

The same recovery table compares `student_type` with the group-level class, counting each student once. Of the 62 `committed` students, 61, or 98.4%, are assigned to group class 2. Of the 44 `wavering` students, 34, or 77.3%, are assigned to group class 1. `get_results(fit, "profile_probabilities")` shows what these group classes represent: one consists of about 83% disengaged enrolments, while the other consists of about 79% engaged enrolments. By contrast, the pooled sample contains 58.8% engaged enrolments, which describes neither group well. This separation between student-level patterns is what the two-level model captures.

## Result tables and diagnostics

Every table is retrieved with one function, `get_results()`, whose `what`
argument names the table. The calls below list the tables a two-level fit
provides, and the comment beside each call states what one row describes.

```r
names(get_results(fit, "all"))                   # every table this fit provides
get_results(fit)                                 # one row per profile and indicator
as.data.frame(fit)                               # the same table
get_results(fit, "profile_probabilities")        # one row per student class and profile
get_results(fit, "counts")                       # effective class sizes at both levels
get_results(fit, "covariances")                  # within-profile residual covariances
get_results(fit, "posteriors")                   # one row per enrolment and profile
get_results(fit, "posteriors", format = "wide")  # one row per enrolment
get_results(fit, "group_posteriors")             # one row per student and student class
get_results(fit, "model")                        # one row describing the whole fit
get_results(fit, "starts")                       # one row per EM start
get_results(fit, "data")                         # the columns used in the fit
get_results(fit, "assignments")                  # each enrolment with its assigned class
get_results(fit, "information_criteria")         # one column per criterion
get_results(fit, "information_criteria", format = "long")  # one row per criterion
get_results(fit, "classification")               # classification quality by class
get_results(fit, "entropy")                      # entropy at each level
get_results(fit, "average_posteriors")           # mean posterior by assigned class
get_results(fit, "average_posteriors", level = "individuals")  # enrolments only

print(summary(fit), rows = 3)                    # every table, three rows each
descriptives(fit)                                # one row per indicator, with the ICC
descriptives(fit, by = "profile")                # the same, by assigned profile
diagnostics(fit)                                 # every classification diagnostic
diagnostics(fit, by = "overall")                 # residuals pooled across profiles
report(fit)                                      # summary, diagnostics and plots
```

All result tables are returned as base `data.frame` objects. `get_results(x, what = )` retrieves a specific table, while `what = "all"` returns all available tables as a named list. `as.data.frame(x)` returns the primary table, which is the measurement model for every fitted family. `summary(x)` prints all available tables, truncated to `rows = 10`, together with the corresponding `get_results()` call for retrieving the full result. Requesting a table that is not available raises `multilpa_bad_argument` and reports the valid choices.

`diagnostics()` and `report()` take `by`, with `report()` also forwarding `rows` to `print(summary(x))`. Unsupported arguments are rejected with `multilpa_bad_argument` rather than silently ignored.

For `get_results()`, `diagnostics()`, `report()`, `parameter_inference()`, `confint()`, and `vcov()`, `data =` is optional because the fitted object already carries the required data. It remains required for `three_step()` and `r3step()`, which need outcome or covariate columns not used in the original fit, and for recovery analyses using unseen `truth` columns.

When `data =` is supplied, shared fitted columns are checked row by row against the original fit. `get_results(x, "assignments")`, `get_results(x, "residuals")`, `three_step()`, and `r3step()` raise `multilpa_bad_inference_data` when these columns disagree. If the available shared columns cannot establish row order, the package instead warns with `multilpa_unverified_alignment`.

`descriptives()` also accepts `by` when applied to a plain data frame. Missing values in the stratifying variable are retained as a labelled `NA` stratum, so the stratum counts still sum to the total number of rows.

All examples below run directly on `course_engagement`, except the categorical examples. Because the bundled indicators are continuous, those examples use `survey` as a placeholder for a dataset with categorical indicators. For a complete reproducible example that simulates data and checks recovery against the generating values, run:

```sh
Rscript validation/synthetic-demo.R
```

The script is in the source repository and is not installed with the package.
It prints the simulated data and the fitted structure, and saves the data, the
true classes and the fitted model to `tmp/synthetic-demo.rds`. It loads the
package sources with `pkgload::load_all()`, so it runs from a clone without
installing the package.

## Model

The model has two latent variables. Each observation `i` (an enrolment) belongs
to a group `j` (a student), each observation belongs to one of `K` latent
profiles, and each group belongs to one of `H` latent group classes. Within
profile `k`, indicator `d` follows a Gaussian distribution with mean `mu_kd`
and variance `variance_kd`, and the indicators are independent within a profile
by default. Group class `h` sets the probability `pi_hk` of profile `k` for
every observation in the group, and `eta_h` is the prevalence of group class
`h`. The likelihood of group `j` is

```text
P(Y_j) = sum_h eta_h * product_i [ sum_k pi_hk * product_d Normal(y_ijd | mu_kd, variance_kd) ]
```

- The profile means and variances are shared across group classes, so a
  profile has the same definition in every group class. This is the
  measurement-invariance assumption.
- `pi_hk` is reported in the `profile_probabilities` table, and `eta_h` as the
  group-class probability. Each observed group contributes once to `eta_h`,
  whatever its size.
- An observation's posterior profile probabilities average over the
  uncertainty in its group's class.
- `variance_model = "equal"` gives each indicator one variance shared by all
  profiles. Different indicators can still have different variances.
- `n_group_classes = 1` reduces the model to single-level LPA.

The formulation applies the nonparametric multilevel mixture model of
[Vermunt (2003)](https://jeroenvermunt.nl/sm2003.pdf) to Gaussian indicators.
Full residual covariance, membership covariates and latent transitions across
ordered occasions are described below. Random slopes are not implemented.

## Categorical indicators

Continuous indicators are described within each profile by a mean and a
variance. Categorical indicators, such as yes/no or rating-scale items, are
described instead by the probability of each answer within each profile: in
one profile 80% may answer *yes*, in another 30%. The `categorical` argument
names the columns to treat this way. Naming none gives two-level latent profile
analysis (LPA), naming all gives two-level latent class analysis (LCA), and
naming some gives a model that combines both kinds of indicator.

### Two-level latent class analysis

The bundled `student_esm` data come from an experience-sampling study of
university students ([Neubauer & Schmiedek, 2024](https://doi.org/10.1007/s11618-023-01182-8)), who answered six prompts a
day for fourteen days. At each prompt when they had not studied, they reported
which of eight leisure activities they had done since the previous prompt. The
dataset holds those prompts for 100 students, 2,582 prompts in all, together
with four ratings of affect on a 1 to 7 scale. Prompts are nested in students,
so a two-level model can separate the activity profiles of prompts from
differences between students in how often they occupy each profile. The
examples use the first week, `day <= 6`, which keeps all 100 students and 1,422
prompts.

```r
activities <- c("time_with_friends", "on_social_media", "tv_video_games",
                "listened_music", "sports", "walking", "reading", "part_time_job")
first_week <- subset(student_esm, day <= 6)

# Two-level LCA: every indicator is categorical
lca <- multilpa(first_week, activities, "student",
                n_profiles = 2, n_group_classes = 2,
                categorical = activities, n_starts = 10, seed = 1)

get_results(lca, "responses")               # probability of each answer by profile
get_results(lca, "profile_probabilities")   # profile prevalence by student class
plot(lca, what = "responses")               # response probabilities by profile
```

The two profiles differ mainly in media use. In profile 2, students reported
social media at 66.3% of prompts, TV or video games at 36.9% and music at
35.4%; in profile 1 the same activities were reported at 7.7%, 13.4% and 5.9%.
Profile 1 has the higher shares of time with friends (16.8% against 6.3%) and
walking (16.3% against 11.3%), and holds 61.4% of the prompts. The student
classes differ in their mix: 92.5% of the prompts of student class 1 fall in
profile 1, whereas 70.8% of those of student class 2 fall in profile 2. The
classes are close in size, at 51.7% and 48.3% of students. Relative entropy is
0.583 for prompts and 0.813 for students, so individual prompts are classified
with considerable uncertainty, and students more clearly.

### Mixed measurement

A mixed model adds the affect ratings as continuous indicators. `worried` is
left out: about half of its ratings are 1, and a profile of low-worry prompts
would reach the variance bound on it.

```r
# Mixed measurement: affect ratings with the activity items
mixed <- multilpa(first_week, c("happy", "relaxed", "exhausted", activities),
                  "student", n_profiles = 2, n_group_classes = 2,
                  categorical = activities, n_starts = 10, seed = 1)

get_results(mixed)                         # means and variances of the ratings
get_results(mixed, "responses")            # answer probabilities of the activities
```

With affect included, the profiles become mainly profiles of mood. Profile 1
has mean ratings of 5.98 for happy, 5.77 for relaxed and 2.33 for exhausted;
profile 2 has 3.90, 3.66 and 4.10. The activities still differ between them,
but less: time with friends is reported at 19.5% of profile-1 prompts and 6.9%
of profile-2 prompts, and social media at 23.4% and 36.2%. The profiles of a
model depend on the indicators it is given, so the choice of indicators is part
of the research question.

### Coding and estimation

Binary, ordinal and nominal indicators share one parameterization: a free
response probability for every category in every profile. Mixture software
estimates the same model when every threshold is free. Numeric, integer,
logical, character and factor columns are accepted. Factor levels keep their
declared order and other types are ordered by `sort()`, so the category coding
does not depend on row order. The number of categories is taken from the
observed values.

With `missing = "fiml"`, an unobserved categorical response is omitted from
that observation's likelihood. The result is the observed-data likelihood
under an ignorable missingness mechanism, and no values are imputed.

Observed-information and robust standard errors are available for categorical
and mixed models, and the parametric bootstrap supports them when the data are
complete. `min_probability` keeps every response probability above a lower
bound, as `min_variance` does for variances. The bound enters the M-step as a
constraint, so the estimates solve the constrained problem exactly.

## Local dependence

The default diagonal model assumes that the indicators are independent within
a profile. `get_results(fit, "residuals")` reports the association left
between each pair of indicators after the profiles are fitted.

```r
get_results(fit, "residuals")
get_results(fit, "residuals", by = "overall")
```

`attendance` counts the days a student was active in a course, and a day counts as active *because*
something was clicked, so residual association is expected. In the pooled
residuals of `fit`, the largest is 0.216, between `forum_read` and
`attendance`, and four of the five pairs significant at 0.05 involve
`attendance`. A full-covariance model estimates that association directly:

```r
dependent <- multilpa(course_engagement, vars, "student",
                      n_profiles = 2, n_group_classes = 2,
                      covariance_model = "full", n_starts = 10, seed = 1)

get_results(dependent, "information_criteria")
get_results(dependent, "residuals", by = "overall")
```

It lowers AIC from 16925.63 to 16712.22 and group-count BIC from 16986.89 to
16826.75, and no pooled residual exceeds 1.4e-05. Ignored residual dependence
matters for enumeration: extra classes can be recruited to absorb it.

## Choosing the number of classes

`enumerate_classes()` fits a grid of models and reports the information
criteria, entropy and convergence diagnostics of each. It does not choose the
model.

```r
candidates <- enumerate_classes(course_engagement, vars, "student",
                                n_profiles = 2:4, n_group_classes = 1:3,
                                n_starts = 10, seed = 1)

as.data.frame(candidates)
get_results(candidates, "criteria")
summary(candidates)
candidate_fit(candidates, n_profiles = 2, n_group_classes = 2)
```

Under the diagonal model, AIC and the group-count BIC, SABIC and CAIC keep
improving up to the largest model in the grid, four profiles and three
classes, which is the pattern residual dependence produces. Criteria that
depend on a sample size are reported with both the number of groups and the
number of observations
([Lukočienė, Varriale, & Vermunt, 2010](https://doi.org/10.1111/j.1467-9531.2010.01231.x)).
For single-level LPA (`n_group_classes = 1`), use the individual-count BIC.

`bootstrap_lrt()` compares two nested models that differ by one class at one
level. It simulates data from the smaller model and refits both to each
simulated dataset; both models must be fitted to complete data.

```r
smaller_fit <- multilpa(course_engagement, vars, "student",
                        n_profiles = 2, n_group_classes = 2, seed = 1)
larger_fit <- multilpa(course_engagement, vars, "student",
                       n_profiles = 3, n_group_classes = 2, seed = 1)

bootstrap_lrt(smaller_fit, larger_fit, iter = 199, seed = 42)
```

The p-value is withheld if any replicate fails to converge.

## Missing indicators

`missing = "fiml"` fits the model to the observed part of each row. It assumes
the missingness is ignorable. In `engagement_with_na`, `attendance` is removed
after each student's twelfth course, which affects 161 of the 1,422 enrolments.

```r
engagement_with_na <- transform(course_engagement,
                                attendance = ifelse(sequence > 12, NA, attendance))

full_fit <- multilpa(engagement_with_na, vars, "student",
                     n_profiles = 2, n_group_classes = 2,
                     covariance_model = "full", missing = "fiml",
                     n_starts = 10, seed = 1)
```

Without `missing = "fiml"`, missing values raise `multilpa_bad_data`.
Membership-covariate fits and the parametric bootstrap require complete data.

## Standard errors and confidence intervals

`parameter_inference()` returns standard errors and Wald intervals. The default
uses the observed information matrix; `vcov_type = "robust"` uses a
Huber-White sandwich clustered on groups, with the MLR scaling correction.

```r
parameter_inference(full_fit)
confint(full_fit)
parameter_inference(full_fit, vcov_type = "robust")
```

The intervals are withheld when a variance is at its lower bound or the
information matrix is not positive definite, and a fit that has not converged
produces a `multilpa_unconverged` warning. Standard errors are not available
for latent transition models.

## Membership covariates

`profile_covariates` and `group_covariates` add multinomial logistic
regressions for membership, estimated jointly with the profiles.
`previous_grade` varies across enrolments, so it predicts profile membership.

```r
with_predictors <- multilpa(course_engagement, vars, "student",
                            n_profiles = 2, n_group_classes = 2,
                            profile_covariates = "previous_grade", seed = 1)

get_results(with_predictors, "coefficients")
parameter_inference(with_predictors)
```

Because the covariate enters while the profiles are estimated, it can change
them. To keep an established solution fixed, use the three-step methods below.

## Three-step analysis

The three-step methods estimate the classes first and then relate them to an
external variable, correcting for classification error, so the classes stay
fixed. `three_step()` estimates class means of an outcome with the BCH method
(Bolck, Croon, & Hagenaars, 2004); `r3step()` regresses membership on
covariates (Vermunt, 2010).

```r
get_results(fit, "classification_errors")

three_step(fit, data = course_engagement, outcome = "previous_grade")
r3step(fit, data = course_engagement, covariates = "previous_grade")
```

The BCH mean of `previous_grade` is 0.222 in the engaged profile and -0.315 in
the disengaged profile. A previous grade one standard deviation higher raises
the log odds of the engaged profile by 0.575 (SE 0.062). Fitted jointly with
`profile_covariates`, the coefficient is 0.422, smaller because the covariate
there also moves the profiles. `three_step()` clusters its variances on
groups by default; `r3step()` uses the observed information by default and
clusters with `vcov_type = "robust"`.

## Staged estimation

`fit_staged()` estimates the profiles first and then the group classes with
the profiles held fixed, so adding classes does not change what the profiles
mean.

```r
staged <- fit_staged(course_engagement, vars, "student",
                     n_profiles = 2, n_group_classes = 2, seed = 1)

get_results(staged, "stages")
parameter_inference(staged)
```

Here the staged fit reaches a log likelihood of -8440.08 against -8439.82 for
the joint fit, so holding the profiles fixed costs little. First-stage
uncertainty is not propagated, so standard errors and information criteria are
conditional on the fixed profiles. `multilpa(fixed = )` holds chosen parameter
blocks directly.

## Profiles over time

`time` records the order of each student's observations, so the fit can report
each student's sequence of assigned profiles.

```r
over_time <- multilpa(course_engagement, vars, "student",
                      n_profiles = 2, n_group_classes = 2,
                      time = "sequence", seed = 1)

get_results(over_time, "sequences", format = "wide")
get_results(over_time, "sequence_summary")
plot(over_time, what = "sequences")
```

The sequences are modal assignments. They describe observed histories but do
not estimate transition probabilities.

## Latent transition analysis

`lta()` estimates the probability of moving between profiles from one occasion
to the next. Transitions are first order and homogeneous, and the profiles have
the same definition at every occasion. With group classes, each class has its
own initial and transition probabilities.

```r
moves <- lta(course_engagement, vars, "student",
             n_profiles = 2, n_group_classes = 2, time = "sequence", seed = 1)

get_results(moves, "transitions")
get_results(moves, "initial")
plot(moves, what = "transitions")
```

`get_tna()` and `get_group_tna()` return the estimated transitions as
[tna](https://cran.r-project.org/package=tna) models, which tna's functions and
plots accept.

```r
tna::centralities(get_tna(moves))
plot(get_group_tna(moves))
```

Standard errors, likelihood-ratio tests and class enumeration are not
available for transition models. The article *Latent transition analysis and
its plots* works through a three-profile model with every plot.

## Seed sensitivity

`sensitivity()` refits a model under several seeds and reports, for each, the
log likelihood, which maximum it reached, and the agreement of its assignments
with the reference fit after matching labels.

```r
sensitivity(fit, seeds = 1:10)
```

A replicated best likelihood supports a solution but does not prove it is the
global maximum.

## Plots

`multilpa_plot_types()` lists every view and the question it answers.

```r
plot(fit)                                   # profile means
plot(fit, scale = "standardized")           # means on a common scale
plot(fit, what = "bars")                    # means with 95% Wald intervals
plot(fit, what = "heatmap")                 # which indicators separate profiles
plot(fit, what = "probabilities")           # profile prevalence by group class
plot(fit, what = "entropy")                 # classification uncertainty
plot(candidates, criterion = "icl_groups")  # any criterion across the grid
```

Series are distinguished by colour, symbol and line type together, with
Okabe-Ito colours. Standardized values are rescaled means, not effect sizes.

## Estimation details

- **Starting values.** EM starts from k-means on the rows as supplied, so row
  order is part of the seed. `starting_values()` extracts a fitted solution to
  start another fit, and `max_iter = 0` scores supplied parameters without
  updating them.
- **Bounds.** `min_variance` (default `1e-6`, in squared input units) and
  `min_probability` keep estimates away from singular solutions; a
  `multilpa_boundary` warning flags an active bound.
- **Convergence.** The best finite likelihood is returned with a
  `multilpa_unconverged` warning if its start did not converge. Convergence is
  judged by the relative change in the likelihood.
- **Labels.** Profile and class numbers are arbitrary; match labels before
  comparing fits.
- **Conditions.** Errors and warnings listed in `?"multilpa-conditions"`
  carry stable classes.

## Vignettes and articles

| Document | Contents |
|---|---|
| `vignette("multilpa")` | The longitudinal case, descriptives and a three-profile, three-class model |
| `vignette("multilpa-evaluation")` | Classification, residual dependence, enumeration and uncertainty |
| `vignette("multilpa-covariates")` | Three-step and joint covariate estimation, and staged models |
| `vignette("multilpa-transitions")` | Initial profiles and transitions across courses |
| Article: *Latent transition analysis and its plots* | A three-profile, two-class LTA with every plot and its transition networks (`vignettes/articles/` in the source repository) |

## Verification

The package imports only base R packages. `testthat`, `knitr`, `rmarkdown` and
`tna` are suggested. Comparisons with Mplus, Latent GOLD, `mclust`, `tidySEM`,
`depmixS4`, `glca`, `poLCA` and `multilevLCA` are in the `equivalence/` folder
of the source repository. Among them, two-level Mplus runs agree on parameter
estimates within `5.1e-8`, the categorical model within `5.8e-6` on
thresholds, and MLR robust standard errors within `1.8e-7`.
[The comparison report](https://github.com/mohsaqr/multilpa/blob/main/equivalence/mplus/COMPARISON.md)
gives the tables and scope.

```sh
R CMD build .
R CMD check multilpa_0.5.1.tar.gz
```

## Authors

The package is written by
[Mohammed Saqr](https://orcid.org/0000-0001-5881-3109) and
[Sonsoles López-Pernas](https://orcid.org/0000-0002-9621-1392). Mohammed Saqr
maintains it; questions are best raised as issues on the
[repository](https://github.com/mohsaqr/multilpa/issues). Run
`citation("multilpa")` to cite the package.
