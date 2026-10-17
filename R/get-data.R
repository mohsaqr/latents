#' Every table a fitted model holds, from one verb
#'
#' One verb returns every tidy table this package computes, named by `what`.
#' The tables were once eleven separate verbs and seven `as.data.frame()`
#' methods with their own `what` catalogues; a reader had to know which verb
#' owned which table before they could ask for it. Here the object decides
#' what it can offer and `what` names it, so the same call shape reads a
#' measurement model, a classification diagnostic and a transition matrix.
#'
#' Every object also answers the usual generics: `print()` shows a compact
#' overview, `summary()` collects its main tables (taking `level` and
#' `vcov_type` where the model has standard errors), and `as.data.frame()`
#' returns its primary table, the same as `get_results(x)`.
#'
#' @param x A fitted model of this package, an enumeration grid, a bootstrap
#'   comparison, a starting-value set, a classification-diagnostics object, or
#'   the summary of any of them.
#' @param what Which table to return, or `"all"` for every table the object can
#'   produce. The tables an object offers depend on its class and are listed
#'   under *Tables* below; asking for one it does not have raises
#'   `latents_bad_argument` naming the ones it does. `NULL`, the default,
#'   gives the object's primary table: the measurement model for a fit, the
#'   candidate grid for an enumeration, the test row for a bootstrap
#'   comparison, and the entropy summary for a diagnostics object.
#' @param ... Arguments for the requested table. Each table accepts only the
#'   arguments listed for it below, and any other name raises
#'   `latents_bad_argument` saying which that table takes, rather than being
#'   dropped on the way to a table that would then answer a different question.
#'
#' @section Tables:
#' Every fitted model offers these:
#' \describe{
#'   \item{`"profiles"`}{The Gaussian measurement model, one row per profile
#'     and continuous indicator: `profile`, `indicator`, `mean`, `variance`,
#'     `standard_deviation`, `mean_standard_error`, `variance_standard_error`.
#'     Zero rows, with the same columns, when every indicator is categorical.
#'     Takes `scale` and `data`.}
#'   \item{`"responses"`}{The categorical measurement model, one row per
#'     profile, categorical indicator and category: `profile`, `indicator`,
#'     `category`, `probability`, `threshold`, `probability_standard_error`.
#'     The threshold is `qlogis(P(y <= category))` and is `NA_real_` for each
#'     indicator's final category, where the cumulative probability is one; a
#'     probability held on its bound has an `NA_real_` standard error. Zero
#'     rows, with the same columns, when no indicator is categorical. Takes
#'     `data`.}
#'   \item{`"covariances"`}{The within-profile residual covariance matrices,
#'     one row per profile and ordered pair of continuous indicators:
#'     `profile`, `indicator`, `indicator_2`, `covariance`. The diagonal
#'     parameterization carries no covariance array but does state a covariance
#'     matrix -- the variances on the diagonal and exact zeros off it -- and it
#'     is reported explicitly, so one shape answers both parameterizations.}
#'   \item{`"counts"`}{Effective class memberships, one row per class at each
#'     level the fit has: `level`, `class`, `effective_count` (the summed
#'     posteriors) and `effective_proportion`. These differ from the modal
#'     counts in `"classification"` whenever entropy is below one.}
#'   \item{`"model"`}{One row describing the whole fit: its dimensions,
#'     likelihood, parameter counts, `aic`, `bic_groups`, `bic_individual` and
#'     its convergence and boundary diagnostics. The columns are those of the
#'     family: a covariate fit reports its predictor counts, because padding
#'     those onto a two-level fit would leave the common case mostly `NA`.
#'     Comparing fits
#'     across families is what `"information_criteria"` is for, and that table
#'     does have one column set for every family.}
#'   \item{`"posteriors"`}{One row per individual and profile: `row`, `group`,
#'     `profile`, `posterior`, `modal`. Takes `format`.}
#'   \item{`"assignments"`}{One row per observation, the columns the fit was
#'     built from (or those of `data`) followed by `profile`, `group_class`
#'     where the model has one, `uncertainty`, and one `posterior_profile_*`
#'     column per profile (and `posterior_noise` for a fit with a noise
#'     component, whose cases are profile 0). `uncertainty` is one minus the posterior of the
#'     assigned profile, which is what modal assignment discards. Takes `data`
#'     and `truth`; `truth` returns a recovery cross-tabulation instead,
#'     described below.}
#'   \item{`"classification"`}{Classification quality, one row per class and
#'     level: `level`, `class`, `n_modal`, `proportion_modal`, `estimated_n`,
#'     `estimated_proportion`, `average_posterior`,
#'     `odds_correct_classification`. Takes `level`. Under sampling weights
#'     `estimated_n` and `estimated_proportion` (and so the odds of correct
#'     classification) are weighted posterior totals; the modal and average
#'     posterior columns describe the units supplied.}
#'   \item{`"average_posteriors"`}{For the units assigned to each class, their
#'     mean posterior of belonging to every class: `level`, `assigned_class`,
#'     `class`, `n_assigned`, `average_posterior`. Each assigned class's rows
#'     sum to one. Conditions on the *assigned* class, where
#'     `"classification_errors"` conditions on the *true* one; the two are
#'     different numbers, not two spellings of one table. Takes `level`.}
#'   \item{`"classification_errors"`}{The probability that a unit truly in one
#'     class is assigned to another: `level`, `true_class`, `assigned_class`,
#'     `probability`, summing to one within each `true_class`. Takes `level`.
#'     On a membership-covariate fit this is a marginal summary over the fitted
#'     units, not a classification-error adjustment conditional on covariates.
#'     A fit with a noise component raises `latents_unsupported_noise`, as
#'     [three_step()] does, and `"all"` and `summary()` leave the table out.}
#'   \item{`"bch_weights"`}{The BCH inverse-error weights, one row per unit and
#'     class: `level`, `unit`, `assigned_class`, `class`, `weight`. Takes
#'     `level`; weight a regression with one level, not with `"both"`. On a
#'     membership-covariate fit these invert the marginal error matrix and are
#'     not a conditional three-step correction. Refused for a fit with a noise
#'     component, like `"classification_errors"`.}
#'   \item{`"entropy"`}{One row per level: `level`, `n_classes`, `n_units`,
#'     `entropy_sum`, `relative_entropy`.}
#'   \item{`"residuals"`}{Bivariate residuals for every indicator pair, testing
#'     the within-profile independence the measurement model assumes: `profile`,
#'     `indicator_1`, `indicator_2`, `kind`, `observed`, `expected`, `residual`,
#'     `effective_n`, `statistic`, `df`, `p_value`, `p_adjusted`. Takes `data`,
#'     `by` and `adjust`.}
#'   \item{`"information_criteria"`}{Every likelihood-penalty criterion the
#'     package computes. Takes `format` and `definitions`.}
#'   \item{`"sequences"`}{The assignments in occasion order, one row per
#'     observation: `group`, `group_class`, `time`, `profile`. Takes `format`.
#'     Raises `latents_no_time` for a fit made without `time`.}
#'   \item{`"sequence_summary"`}{One row per group class, summarizing how much
#'     data it contributes: `group_class`, `groups`, `observations`,
#'     `mean_length`, `median_length`, `shortest`, `longest`, `complete`,
#'     `gaps`. Raises `latents_no_time` for a fit made without `time`.}
#'   \item{`"starts"`}{One row per EM start: `start`, `log_likelihood`,
#'     `converged`, `iterations`, `error` and, where the family records it,
#'     `boundary`.}
#'   \item{`"data"`}{The columns of the original data the model was fitted to,
#'     one row per observation in input order: the identifier (except for a
#'     single-level fit), the occasion where there is one, the sampling
#'     weights of a weighted fit (as rescaled for fitting), and the indicators
#'     under their original names. Columns a model never saw, membership
#'     covariates included, are not here, which is why [three_step()] and
#'     [r3step()] still take `data`.}
#' }
#' Every family with discrete group classes adds:
#' \describe{
#'   \item{`"group_posteriors"`}{One row per observed group and group class:
#'     `group`, `group_size`, `log_likelihood`, `group_class`, `posterior`,
#'     `modal`. Takes `format`.}
#' }
#' And each family adds its own:
#' \describe{
#'   \item{`"profile_probabilities"`, `"stages"`}{A `multilpa` fit.
#'     `"profile_probabilities"` has one row per group class and profile:
#'     `group_class`, `profile`, `probability`, `group_class_probability`.
#'     `"stages"` has one row per estimation stage: `stage`, `group_classes`,
#'     `fixed`, `log_likelihood`, `parameters`, `parameters_with_measurement`,
#'     `converged`. An ordinary fit has a single `"joint"` row and a
#'     [fit_staged()] result has two; `fixed` is `NA_character_` for a stage
#'     that estimated every block, never a sentinel such as `"none"`.}
#'   \item{`"coefficients"`}{A covariate fit. One row per estimated membership
#'     coefficient: `level` (`"profile"` or `"group"`), `outcome`, `term`,
#'     `parameter` (always `"logit"`), `estimate`. The standard errors are not
#'     here; [parameter_inference()] reports them.}
#'   \item{`"transitions"`, `"initial"`, `"sequence_lengths"`}{A transition
#'     fit. `"transitions"` has one row per group class and ordered pair of
#'     profiles, with `group_class`, `from`, `to`, `probability`,
#'     `expected_count`, `stable`, `estimated` and
#'     `group_class_probability`; it takes `estimated` and
#'     `stable` to restrict it. `"initial"` has one row per group class and
#'     profile: `group_class`, `profile`, `probability`, `prevalence`,
#'     `group_class_probability`. `"sequence_lengths"` has one row per group:
#'     `group`, `group_class`, `occasions`, `observations`, `complete`.}
#'   \item{`"candidates"`, `"best"`, `"classification"`}{An enumeration
#'     grid. `"best"` (formerly `"criteria"`, which still works) has one row
#'     per information criterion and sample-size convention: `criterion`, `convention`, `candidate` (the
#'     candidate its minimum falls on), `n_profiles`, `n_group_classes`,
#'     `model` and `value`. `"candidates"` has one row per candidate model:
#'     its name (`candidate`, first), its class counts, its covariance
#'     `model`, log likelihood,
#'     parameter count, every criterion
#'     under both sample-size conventions, both entropies, the average
#'     posterior probability (AvePP) summaries of each level --
#'     `avepp_min`, `avepp_max`, `avepp_mean`,
#'     `avepp_below` and the same four prefixed `group_` -- and the
#'     convergence, boundary, replication, warning and error diagnostics.
#'     Failed candidates are retained with `NA` estimates and their error text.
#'     `"candidates"` takes `cutoff` (default `0.8`): `avepp_below` and
#'     `group_avepp_below` count the classes whose AvePP is below it. A class with no modally assigned
#'     unit has no AvePP; it is left out of the smallest, largest and mean and
#'     counted as below the cutoff. A level with one class has `NA`
#'     summaries. `"classification"` has one row per fitted candidate,
#'     level and class: `candidate`, `n_profiles`, `n_group_classes`,
#'     `model`, then the
#'     columns of a fit's `"classification"` table and that level's
#'     `relative_entropy`. Failed candidates and levels with a single class
#'     contribute no rows.}
#'   \item{`"deltas"`, `"knees"`}{Model-selection evidence of an enumeration
#'     grid, each taking `using` (default `"bic_groups"`). Neither selects a
#'     model. `"deltas"` has one row per candidate: `candidate`, its class
#'     counts and `model`, `criterion`, `value`, `delta_profile` (how much the
#'     criterion fell from the converged candidate with one profile fewer and
#'     the same group classes and structure; positive means this candidate is
#'     lower) with `from_profile` naming that candidate, and
#'     `delta_group_class` with `from_group_class` for one group class fewer.
#'     `value` and both deltas are `NA` for a candidate that did not converge,
#'     and a delta is `NA` when the smaller neighbour is missing or did not
#'     converge. `"knees"` has one row per
#'     line of the grid: `direction` (`"profiles"`, along the profiles at
#'     fixed group classes; `"group_classes"`, along the group classes at
#'     fixed profiles; `"parameters"`, every candidate by its number of
#'     parameters, one line per structure), `line`, `criterion`, `from`, `to`,
#'     `n_points`,
#'     `knee_candidate` (the maximum chord-distance knee: with both axes
#'     rescaled to 0-1, the
#'     candidate furthest below the straight line from the first to the last,
#'     measured perpendicular to it), `distance`,
#'     `knee_without_largest` (the knee when the last candidate of the line is
#'     left out), `stable` and `note`. This is a simplified geometric rule
#'     inspired by Kneedle, without its smoothing or sensitivity thresholds;
#'     a very small positive distance can therefore still name a knee.
#'     A line needs three candidates for a
#'     knee and four for `knee_without_largest` and `stable`; otherwise they
#'     are `NA` and `note` says why.}
#'   \item{`"comparison"`, `"classification"`}{A [compare_models()] result.
#'     `"comparison"` is its table, one row per model. `"classification"`,
#'     for profile, class, membership-covariate and transition models, has
#'     one row per model, level and class:
#'     `model`, then the columns of a fit's `"classification"` table and that
#'     level's `relative_entropy`; levels with one class are left out. The
#'     models themselves come back with [get_candidate()].}
#'   \item{`"test"`, `"replicates"`}{A bootstrap comparison. `"test"` is the
#'     one-row result; `"replicates"` has one row per simulated dataset.}
#'   \item{`"covariances"`}{A starting-value set. One row per profile and
#'     ordered pair of continuous indicators: `profile`, `indicator`,
#'     `indicator_2`, `covariance`.}
#' }
#'
#' @section Tables of the other model families:
#' Each family's first table is its default `what`; `"all"` returns a named
#' list of every table (without `"recovery"`, which needs `data` and `truth`).
#' Tables with intervals take `level`, and `vcov_type` (`"observed"`,
#' `"robust"` or `"opg"`; `NULL` reuses the type stored with a regression
#' fit). The usual generics answer as well: `coef()`, `vcov(type = )`,
#' `confint(parm = , level = )`, `logLik()` and `nobs()`.
#' \describe{
#'   \item{Mixture regressions ([mixture_regression()])}{`"coefficients"`
#'     (one row per class and term, with Wald statistics, `p_adjusted` by
#'     Benjamini-Hochberg over non-intercepts and `exp_estimate` for the
#'     binomial, Poisson and ordinal families; ordinal thresholds lead each
#'     class as terms `"threshold:<lower>|<upper>"`), `"classes"`,
#'     `"membership"` (multinomial logits and odds ratios), `"group_classes"`
#'     (rows for two-level fits only), `"fit"` (`bic` penalised by independent units,
#'     `bic_rows` by rows as flexmix reports it), `"assignments"`,
#'     `"groups"`, `"fitted"`, `"starts"` (`stage` `"screened"` or
#'     `"completed"`), `"classification"` and `"recovery"` (`by = "class"` or
#'     `"group_class"`).}
#'   \item{Growth mixtures (`mixture_regression(random = )`)}{The tables
#'     above plus `"random"` (random-effect variances, standard deviations
#'     and correlations, with intervals on the log and Fisher-z scales),
#'     `"random_effects"` (one row per person), `"trajectories"` (each class's
#'     mean trajectory with a band) and `"individual"` (each person's own
#'     predicted curve); the last two take `time`. A multilevel growth mixture
#'     adds `"group_classes"` and `"clusters"`, and its `"fit"` counts clusters
#'     in `bic` (persons in `bic_persons`).}
#'   \item{Additive fits (`multilpa(family = "additive")`)}{`"parameters"`
#'     (class means, between- and within-group variances and class weights;
#'     variance intervals on the log scale, weights on the logit scale;
#'     `adjust` sets the p-value adjustment, `"none"` by default),
#'     `"group_classes"`, `"groups"`, `"intercepts"` and
#'     `"conditional_intercepts"` (each group's posterior intercept, overall
#'     or given each class), `"assignments"`, `"classification"`, `"fit"`,
#'     `"recovery"` and `"starts"`. `effective_groups` in `"group_classes"` is
#'     `w^2 (1 - w) / Var(w)`, the number of perfectly classified groups that
#'     would estimate the class weight as precisely; below 50 the class's
#'     intervals carry a `latents_weak_class` warning.}
#'   \item{Cross-level fits (`multilpa(family = "full_cross_level")` or
#'     `"restricted_cross_level"`)}{`"profiles"` (individual profiles, then
#'     the group classes' means of group means), `"group_classes"`,
#'     `"composition"` (profile shares within each class: model-implied for
#'     the full family, a posterior-weighted description for the restricted
#'     one), `"groups"`, `"assignments"`, `"fit"` and `"starts"`. The
#'     likelihood is a working likelihood, comparable only between
#'     cross-level fits of the same data.}
#'   \item{General transition fits ([lta()] with `transitions = "occasion"`,
#'     transition or initial covariates, `measurement = "occasion"`,
#'     `order = 2`, `mover_stayer = TRUE`, `weights`, or `ordinal`/`count`
#'     indicators)}{
#'     `"transitions"` (model-implied, averaged over the groups observed,
#'     with their sampling weights in a weighted fit; `"initial"` and
#'     `"second_order_transitions"` likewise),
#'     `"transition_coefficients"` (log odds of each move against staying,
#'     with odds ratios), `"second_order_transitions"` and
#'     `"second_order_coefficients"` (`order = 2`), `"initial"` and
#'     `"initial_coefficients"`, `"profiles"`, `"responses"` (categorical
#'     indicators: `occasion`, `profile`, `indicator`, `category`,
#'     `probability`), `"ordinal"` and
#'     `"count_means"` (ordinal and count indicators), `"group_classes"`,
#'     `"assignments"`, `"fit"` and `"starts"`.}
#'   \item{Pooled imputations ([pool_imputations()])}{`"estimates"` (the
#'     pooled table), `"imputations"` (each imputation's aligned estimates)
#'     and `"fits"` (one row per imputation).}
#'   \item{Transition and group-class family grids
#'     (`enumerate_classes(time = )`, `enumerate_classes(family = )`)}{
#'     `"candidates"` (one row per candidate, its name first, with its
#'     likelihood, criteria, entropy and diagnostics) or `"best"`, the row
#'     with the lowest `criterion` (`"bic"` by groups, the default, `"aic"`
#'     or `"bic_individual"`): the lowest value of a criterion, not a choice
#'     of model. With one group class every family is the same model, so
#'     those rows share a likelihood. Their candidates are taken out by name
#'     with [get_candidate()]; transition grids can also be compared with
#'     [compare_models()].}
#' }
#' @section What `"all"` returns: A named list of every table the object can
#'   produce, in catalogue order, built from the same definitions a single
#'   `what` uses, so the two can never disagree. It takes no further arguments:
#'   every table is built with its defaults, and supplying anything else raises
#'   `latents_bad_argument`. A table this particular object cannot produce is
#'   left out rather than erroring -- a sequence table for a fit made without
#'   `time`, bivariate residuals for a family with no discrete group classes --
#'   so the names of the list say what was available. Only the classed
#'   refusals that mean *this object has no such table* are skipped; anything
#'   else propagates.
#'
#' @section Arguments the tables take:
#'   These reach the requested table through `...`, and each table takes
#'   only the ones listed for it under *Tables* above.
#' \describe{
#'   \item{`data`}{A data frame with one row per observation of the fit, in the
#'   order it was fitted in. `NULL` uses the columns the fit carries. The
#'   standard errors of `"profiles"` and `"responses"` come from
#'   [parameter_inference()] on these data. Without `data`, a fit that cannot
#'   supply them (bound-active, unconverged, singular information, a noise
#'   component) gives `NA` with a `latents_no_standard_errors` message; with
#'   `data`, it raises the underlying condition, as asked.
#'   }
#'   \item{`scale`}{For `"profiles"`, `"raw"` reports the estimates in input units
#'   and `"standardized"` divides each indicator's deviation from its grand
#'   mean by that indicator's observed standard deviation, which is what
#'   `plot(x, scale = "standardized")` draws. Every standard error beside an
#'   estimate is divided by the same constant, so the whole row is on one
#'   scale.
#'   }
#'   \item{`format`}{For the posterior tables and `"sequences"`, `"long"` gives one
#'   row per unit and class and `"wide"` one row per unit with one column per
#'   class. For `"information_criteria"`, `"wide"` is the one-row reporting
#'   shape and `"long"` gives one row per criterion and sample-size convention.
#'   }
#'   \item{`level`}{Which level of classification to report: `"individuals"`,
#'   `"groups"` or `"both"`. `NULL`, the default, gives `"both"` for a fit that
#'   has discrete group classes and `"individuals"` for one that does not,
#'   which is the same rule [diagnostics()] follows.
#'   }
#'   \item{`truth`}{For `"assignments"`, one or more column names of `data` holding
#'   a known label to check the assignments against. See *Recovery* below.
#'   }
#'   \item{`by`}{For `"residuals"`, `"profile"` assesses each profile separately
#'   and `"overall"` pools them.
#'   }
#'   \item{`adjust`}{For `"residuals"`, the multiplicity correction applied across
#'   the indicator pairs.
#'   }
#'   \item{`definitions`}{For `"information_criteria"` in long format, `TRUE` adds
#'   the formula and the reference for each criterion.
#'   }
#'   \item{`estimated`, `stable`}{For `"transitions"`, `TRUE` or `FALSE` restricts the
#'   table to the rows that were estimated from observed occupancy, or to the
#'   diagonal. `NULL` keeps every row.
#'   
#'   }
#' }
#'
#' @section Recovery: `what = "assignments"` with `truth` returns a
#'   cross-tabulation of the model's labels against known ones instead of the
#'   per-observation table: one row per class and truth value, with columns
#'   `assignment` (which model label the column was compared against),
#'   `class`, `truth` (the column name), `value`, `n` and `proportion`, the
#'   share of the units carrying that truth value that were assigned to that
#'   class. Several truth columns may be named at once and are stacked.
#'
#'   Each truth column is compared against the level it describes, which is
#'   read off the data rather than guessed: a column that takes one value
#'   within every group is a property of the group and is compared against
#'   `group_class`; a column that varies inside any group cannot be a
#'   group-level label and is compared against `profile`. The `assignment`
#'   column records which comparison was made, so the choice is never silent.
#'   Group-level recovery counts each group once, even when groups have
#'   different numbers of observations. Unused factor levels are omitted,
#'   because their within-truth proportions are undefined.
#'   Every column is compared against `profile` for a fit with no discrete
#'   group classes.
#'
#'   The labels are the model's own and carry no order, so a recovery table
#'   whose mass sits off the diagonal is a relabelling, not a failure.
#'
#' @section Information criteria: Let `q` be the number of free parameters,
#'   `n` the chosen sample size, and `EN` the classification entropy of the
#'   level matching that convention. The criteria are `deviance = -2L`,
#'   `aic = -2L + 2q`, `bic = -2L + q log(n)`,
#'   `sabic = -2L + q log((n + 2) / 24)`, `caic = -2L + q (log(n) + 1)`,
#'   `awe = -2(L - EN) + 2q (1.5 + log(n))`, `icl = -2L + q log(n) + 2 EN`,
#'   `kic = -2L + 3(q + 1)`, and `clc = -2L + 2 EN`. Note that `clc` uses the
#'   entropy sum, as `icl` and `awe` here do; `tidyLPA` reports a `CLC` built
#'   from relative entropy instead, which is bounded by one and so penalizes
#'   almost nothing, and the two numbers are not comparable.
#'
#'   Every criterion that depends on a sample size is reported under both
#'   multilevel conventions, `"groups"` and `"individuals"`, and uses
#'   group-level posteriors under the first and individual-level posteriors
#'   under the second. That is a stated per-level choice, not a unique
#'   multilevel definition: use one convention consistently across compared
#'   candidates. `deviance`, `aic` and `kic` depend on no sample size and carry
#'   `convention = NA_character_`. `clc` depends on none either, but its
#'   convention selects which level's classification uncertainty it penalizes,
#'   so it is reported once per convention. Individual sample sizes exclude
#'   rows with no observed indicators. Where the group classification entropy
#'   is undefined, group-level `awe`, `icl` and `clc` are `NA`.
#'
#' @section Classification quality: The odds of correct classification for
#'   class `k` is `(p / (1 - p)) / (r / (1 - r))`, where `p` is the average
#'   posterior in the assigned class and `r` is the model-estimated class
#'   proportion. Values near one indicate classification no better than the
#'   class proportion alone; values of five or more are conventionally read as
#'   adequate separation. It is `NA_real_` wherever it is undefined, that is
#'   whenever `average_posterior` is missing, or it or `estimated_proportion`
#'   is zero or one, so a single-class solution always reports `NA_real_`.
#'   `average_posterior` is itself `NA_real_` for a class with no modal
#'   members, where the mean is taken over nothing and is undefined rather than
#'   zero. Modal assignment discards classification uncertainty, which is why
#'   `n_modal` and `estimated_n` differ whenever entropy is below one.
#'
#'   `relative_entropy` is `1 - EN / (n log K)`, on the zero-to-one scale, so
#'   higher is sharper; `entropy_sum` is the raw `EN` the entropy-penalized
#'   criteria use, so lower is sharper. `relative_entropy` is `NA_real_` when a
#'   level has a single class, where it is undefined rather than perfect.
#'
#'   A fit with individual profiles but no discrete group classes returns
#'   individuals only for `level = "both"`, and raises
#'   `latents_no_group_classes` for `level = "groups"`.
#'
#' @section Bivariate residuals: The measurement model assumes the indicators
#'   are independent within a profile. `"residuals"` tests that pair by pair,
#'   weighting every observation by the posterior the fit holds for it, so a
#'   frame in any other row order would pair each observation with someone
#'   else's posterior; the alignment is checked rather than assumed. A large
#'   residual is evidence of local dependence between two indicators, which the
#'   enumeration grid will otherwise absorb by asking for an extra profile. The
#'   p-values are approximate, and `p_adjusted` applies the multiplicity
#'   correction named by `adjust` across the pairs.
#'   A fit with sampling weights uses those weights in its descriptive
#'   associations and `effective_n`. Design-based residual inference is not
#'   implemented: `statistic`, `df`, `p_value` and `p_adjusted` are `NA` for
#'   weighted fits rather than unweighted tests attached to weighted moments.
#'
#' @references Schwarz, G. (1978). Estimating the dimension of a model.
#'   *Annals of Statistics*, 6, 461--464.
#'   Sclove, S. L. (1987). Application of model-selection criteria to some
#'   problems in multivariate analysis. *Psychometrika*, 52, 333--343.
#'   Bozdogan, H. (1987). Model selection and Akaike's information criterion.
#'   *Psychometrika*, 52, 345--370.
#'   Banfield, J. D., & Raftery, A. E. (1993). Model-based Gaussian and
#'   non-Gaussian clustering. *Biometrics*, 49, 803--821.
#'   Biernacki, C., & Govaert, G. (1997). Using the classification likelihood
#'   to choose the number of clusters. *Computing Science and Statistics*, 29,
#'   451--457.
#'   Cavanaugh, J. E. (1999). A large-sample model selection criterion based on
#'   Kullback's symmetric divergence. *Statistics and Probability Letters*, 42,
#'   333--343.
#'   Biernacki, C., Celeux, G., & Govaert, G. (2000). Assessing a mixture model
#'   for clustering with the integrated completed likelihood. *IEEE
#'   Transactions on Pattern Analysis and Machine Intelligence*, 22, 719--725.
#'   Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent structure
#'   models with categorical variables. *Political Analysis*, 12, 3--27.
#'   Nagin, D. S. (2005). *Group-Based Modeling of Development*. Harvard
#'   University Press.
#'   Vermunt, J. K. (2010). Latent class modeling with covariates: two improved
#'   three-step approaches. *Political Analysis*, 18, 450--469.
#' @return A base `data.frame` for a single `what`, whose columns are given
#'   under *Tables*, or a named list of such data frames for `what = "all"`.
#' @seealso [as.data.frame()], which returns the primary table alone,
#'   [summary()], which prints every table, and [diagnostics()] for the
#'   classification tables gathered with a condensed reading.
#' @examples
#' fit <- multilpa(
#'   course_engagement,
#'   vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
#'   id = "student", n_profiles = 2, n_group_classes = 2, n_starts = 4,
#'   seed = 1
#' )
#' get_results(fit)
#' get_results(fit, what = "profile_probabilities")
#' get_results(fit, what = "entropy")
#' names(get_results(fit, what = "all"))
#'
#' # `engagement` and `student_type` are the kinds each row was simulated from,
#' # which the model never saw. Both recovery checks are one call.
#' get_results(fit, what = "assignments", data = course_engagement,
#'          truth = c("engagement", "student_type"))
#' @aliases latents-as-data-frame latents-print latents-summary
#' @aliases as.data.frame.latents_family_enumeration
#' @aliases as.data.frame.latents_growth_mixture
#' @aliases as.data.frame.latents_mixture_regression
#' @aliases as.data.frame.latents_pooled
#' @aliases as.data.frame.latents_table
#' @aliases as.data.frame.latents_transition_enumeration
#' @aliases as.data.frame.multilpa
#' @aliases as.data.frame.multilpa_additive
#' @aliases as.data.frame.multilpa_bootstrap_lrt
#' @aliases as.data.frame.multilpa_covariates
#' @aliases as.data.frame.multilpa_cross_level
#' @aliases as.data.frame.multilpa_enumeration
#' @aliases as.data.frame.multilpa_lta
#' @aliases as.data.frame.multilpa_start
#' @aliases as.data.frame.multilpa_transitions
#' @aliases as.data.frame.summary_multilpa
#' @aliases as.data.frame.summary_multilpa_additive
#' @aliases as.data.frame.summary_multilpa_bootstrap_lrt
#' @aliases as.data.frame.summary_multilpa_covariates
#' @aliases as.data.frame.summary_multilpa_enumeration
#' @aliases as.data.frame.summary_multilpa_transitions
#' @aliases coef.latents_growth_mixture
#' @aliases coef.latents_mixture_regression
#' @aliases coef.multilpa_additive
#' @aliases coef.multilpa_lta
#' @aliases confint.latents_growth_mixture
#' @aliases confint.latents_mixture_regression
#' @aliases confint.multilpa_additive
#' @aliases get_results.latents_family_enumeration
#' @aliases get_results.latents_growth_mixture
#' @aliases get_results.latents_mixture_regression
#' @aliases get_results.latents_pooled
#' @aliases get_results.latents_transition_enumeration
#' @aliases get_results.multilpa_additive
#' @aliases get_results.multilpa_cross_level
#' @aliases get_results.multilpa_lta
#' @aliases logLik.latents_growth_mixture
#' @aliases logLik.latents_mixture_regression
#' @aliases logLik.multilpa_additive
#' @aliases logLik.multilpa_cross_level
#' @aliases logLik.multilpa_lta
#' @aliases nobs.latents_growth_mixture
#' @aliases nobs.latents_mixture_regression
#' @aliases nobs.multilpa_additive
#' @aliases nobs.multilpa_cross_level
#' @aliases nobs.multilpa_lta
#' @aliases print.latents_family_enumeration
#' @aliases print.latents_growth_mixture
#' @aliases print.latents_mixture_regression
#' @aliases print.latents_plots
#' @aliases print.latents_pooled
#' @aliases print.latents_table
#' @aliases print.latents_transition_enumeration
#' @aliases print.multilpa
#' @aliases print.multilpa_additive
#' @aliases print.multilpa_bootstrap_lrt
#' @aliases print.multilpa_covariates
#' @aliases print.multilpa_cross_level
#' @aliases print.multilpa_enumeration
#' @aliases print.multilpa_lta
#' @aliases print.multilpa_start
#' @aliases print.multilpa_transitions
#' @aliases print.summary_latents_growth_mixture
#' @aliases print.summary_latents_mixture_regression
#' @aliases print.summary_multilpa
#' @aliases print.summary_multilpa_additive
#' @aliases print.summary_multilpa_bootstrap_lrt
#' @aliases print.summary_multilpa_covariates
#' @aliases print.summary_multilpa_enumeration
#' @aliases print.summary_multilpa_transitions
#' @aliases summary.latents_family_enumeration
#' @aliases summary.latents_growth_mixture
#' @aliases summary.latents_mixture_regression
#' @aliases summary.latents_pooled
#' @aliases summary.latents_transition_enumeration
#' @aliases summary.multilpa
#' @aliases summary.multilpa_additive
#' @aliases summary.multilpa_bootstrap_lrt
#' @aliases summary.multilpa_covariates
#' @aliases summary.multilpa_cross_level
#' @aliases summary.multilpa_enumeration
#' @aliases summary.multilpa_lta
#' @aliases summary.multilpa_transitions
#' @aliases vcov.latents_growth_mixture
#' @aliases vcov.latents_mixture_regression
#' @aliases vcov.multilpa_additive
#' @aliases vcov.multilpa_lta
#' @export
get_results <- function(x, what = NULL, ...) {
  UseMethod("get_results")
}

#' @rdname get_results
#' @export
get_results.default <- function(x, what = NULL, ...) {
  stop(errorCondition(sprintf(
    "`get_results()` has no method for an object of class %s.",
    paste(sprintf("`%s`", class(x)), collapse = ", ")),
    class = "latents_bad_argument", call = NULL))
}

#' @rdname get_results
#' @export
get_results.multilpa <- function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @rdname get_results
#' @export
get_results.multilpa_covariates <- function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @rdname get_results
#' @export
get_results.multilpa_transitions <- function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}


#' @rdname get_results
#' @export
get_results.multilpa_enumeration <- function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @rdname get_results
#' @export
get_results.latents_comparison <- function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @rdname get_results
#' @export
get_results.multilpa_bootstrap_lrt <- function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @rdname get_results
#' @export
get_results.multilpa_start <- function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @rdname get_results
#' @export
get_results.multilpa_diagnostics <- function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @rdname get_results
#' @export
get_results.summary_multilpa <- function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @rdname get_results
#' @export
get_results.summary_multilpa_covariates <- function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}


#' @rdname get_results
#' @export
get_results.summary_multilpa_transitions <- function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @rdname get_results
#' @export
get_results.summary_multilpa_enumeration <- function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @rdname get_results
#' @export
get_results.summary_multilpa_bootstrap_lrt <- function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' Coerce an object of this package to its primary table
#'
#' The one body behind every `as.data.frame()` method here, so no class can
#' coerce to something other than the first entry of its own catalogue.
#'
#' @param x The object.
#' @param row.names The caller's `row.names`.
#' @param extra `list(...)` as the method received it, which must be empty.
#' @return A base `data.frame`.
#' @noRd
.multilpa_coerce <- function(x, row.names, extra) {
  .multilpa_reject_extra_arguments(
    extra, "as.data.frame()",
    "It coerces to the primary table; get_results(x, what = ) has the others.")
  result <- get_results(x)
  row.names(result) <- row.names
  result
}

#' Resolve one `what` against an object's catalogue and build that table
#'
#' The single place `what` is matched, arguments are checked and a table is
#' built, so every method is the same three words and no class can drift from
#' the shared contract.
#'
#' @param x The object.
#' @param what The caller's `what`, possibly `NULL`.
#' @param extra `list(...)` as the method received it.
#' @return A `data.frame`, or a named list of them for `what = "all"`.
#' @noRd
.multilpa_dispatch_table <- function(x, what, extra) {
  catalogue <- .multilpa_catalogue(x)
  what <- .multilpa_match_table(what, names(catalogue))
  if (identical(what, "all")) return(.multilpa_all_tables(x, catalogue, extra))
  entry <- catalogue[[what]]
  .multilpa_check_table_arguments(what, entry$args, extra)
  do.call(entry$build, c(list(x), extra[intersect(names(extra), entry$args)]))
}

#' Match a requested table name against what an object offers
#'
#' Exact matching rather than `match.arg()`'s partial matching: the catalogues
#' run to nineteen entries and share prefixes -- `starts`, `stages`,
#' `sequences`, `sequence_summary` -- so a partial match would resolve to
#' whichever happened to come first. The refusal lists every available name, so
#' nothing is lost by being exact.
#'
#' @param what The caller's `what`, possibly `NULL`.
#' @param choices The catalogue's names.
#' @return A single table name, or `"all"`.
#' @noRd
.multilpa_match_table <- function(what, choices) {
  if (is.null(what)) return(choices[1L])
  stopifnot("`what` must be a single table name, or NULL for the first one" =
              is.character(what) && length(what) == 1L && !is.na(what))
  if (identical(what, "all") || what %in% choices) return(what)
  # An earlier name of a table still finds it, so scripts written against a
  # released version keep working.
  aliases <- c(criteria = "best")
  if (what %in% names(aliases) && aliases[[what]] %in% choices) return(aliases[[what]])
  stop(errorCondition(sprintf(
    "`what = \"%s\"` is not a table of this object. It has %s, and \"all\".",
    what, paste(sprintf("\"%s\"", choices), collapse = ", ")),
    class = "latents_bad_argument", call = NULL))
}

#' Refuse an argument the requested table does not take
#'
#' A dropped argument returns a table that answers a different question from
#' the one asked, and does so without a word. Naming it back costs the caller
#' one correction and saves them a wrong number.
#'
#' @param what The table name, for the message.
#' @param accepted The argument names that table takes.
#' @param extra `list(...)` as the method received it.
#' @return `NULL`, invisibly, when there is nothing to refuse.
#' @noRd
.multilpa_check_table_arguments <- function(what, accepted, extra) {
  if (length(extra) == 0L) return(invisible(NULL))
  supplied <- names(extra) %||% rep("", length(extra))
  unknown <- is.na(match(supplied, accepted)) | !nzchar(supplied)
  if (!any(unknown)) return(invisible(NULL))
  labels <- ifelse(nzchar(supplied[unknown]),
                   sprintf("`%s`", supplied[unknown]), "an unnamed argument")
  takes <- if (length(accepted) == 0L) "takes no further arguments" else
    sprintf("takes %s", paste(sprintf("`%s`", accepted), collapse = ", "))
  stop(errorCondition(sprintf("`what = \"%s\"` does not use %s; it %s.",
                              what, paste(labels, collapse = ", "), takes),
                      class = "latents_bad_argument", call = NULL))
}

#' Classed refusals that mean "this object has no such table"
#'
#' `what = "all"` walks the whole catalogue, and some of its entries cannot be
#' built for every object: a fit made without `time` has no sequences, a fit
#' without discrete group classes has no bivariate residuals. Those are refusals about
#' the object, not failures, and `"all"` leaves the table out. Every other
#' condition propagates, so a genuine defect is never swallowed.
#'
#' @return A character vector of condition classes.
#' @noRd
.multilpa_absent_table_conditions <- function() {
  c("latents_no_time", "latents_no_group_classes",
    "latents_no_indicator_data", "latents_incomplete_fit",
    "latents_no_continuous", "latents_no_categorical",
    # A fit whose classes cannot be separated has no three-step correction to
    # report, which is a fact about the fit and not a reason for `summary()`
    # to fail.
    "latents_inseparable_classes",
    # A noise component has no classification-error correction yet.
    "latents_unsupported_noise")
}

#' Build every table an object can produce
#'
#' @param x The object.
#' @param catalogue Its catalogue.
#' @param extra `list(...)`, which must be empty.
#' @return A named list of `data.frame`s, without the ones this object has not.
#' @noRd
.multilpa_all_tables <- function(x, catalogue, extra) {
  .multilpa_check_table_arguments("all", character(), extra)
  absent <- .multilpa_absent_table_conditions()
  tables <- lapply(catalogue, function(entry) {
    # Errors only: a table that warns is still a table, and the warning is the
    # fit's own, so it is left to reach the caller who asked for every table.
    tryCatch(entry$build(x), error = function(condition) {
      if (inherits(condition, absent)) return(NULL)
      stop(condition)
    })
  })
  Filter(Negate(is.null), tables)
}

#' The tables an object offers, and the arguments each of them takes
#'
#' One definition per table, shared by `what = <name>` and by `what = "all"`,
#' so the two cannot return differently built versions of the same table. The
#' first entry is the object's primary table, which is what `as.data.frame()`
#' and a `NULL` `what` return.
#'
#' A table the family defines but this particular object cannot produce stays
#' in the catalogue and raises its classed condition, so the refusal names the
#' reason. A table another family defines is absent, and the refusal lists what
#' this one has.
#'
#' @param x The object.
#' @return A named list of `list(build =, args =)` entries, in reporting order:
#'   the measurement model, then the mixing structure, then the unit-level
#'   tables, then the diagnostics, then the sequence tables, then provenance.
#' @noRd
.multilpa_catalogue <- function(x) {
  if (inherits(x, "multilpa_enumeration")) {
    return(list(
      candidates = .multilpa_table(
        function(x, cutoff = 0.8) .multilpa_enumeration_grid(x, cutoff),
        "cutoff"),
      best = .multilpa_table(
        function(x) .multilpa_enumeration_minima(.multilpa_enumeration_grid(x))),
      classification = .multilpa_table(.multilpa_enumeration_classification),
      deltas = .multilpa_table(
        function(x, using = "bic_groups") .latents_deltas_table(x, using), "using"),
      knees = .multilpa_table(
        function(x, using = "bic_groups") .latents_knees_table(x, using), "using")))
  }
  if (inherits(x, "latents_comparison")) {
    catalogue <- list(comparison = .multilpa_table(function(x) {
      attr(x, "models") <- NULL
      x
    }))
    models <- attr(x, "models")
    if (is.null(models) || all(vapply(models, .latents_comparison_measurement_fit,
                                    logical(1)))) {
      catalogue$classification <- .multilpa_table(.latents_comparison_classification)
    }
    return(catalogue)
  }
  if (inherits(x, "multilpa_bootstrap_lrt")) {
    return(list(
      test = .multilpa_table(.multilpa_bootstrap_test_frame),
      replicates = .multilpa_table(function(x) x$replicates)))
  }
  if (inherits(x, "multilpa_start")) {
    return(list(
      profiles = .multilpa_table(.multilpa_start_profile_frame),
      responses = .multilpa_table(.multilpa_start_response_frame),
      profile_probabilities = .multilpa_table(.multilpa_start_probability_frame),
      covariances = .multilpa_table(.multilpa_start_covariance_frame)))
  }
  if (inherits(x, "multilpa_diagnostics")) {
    return(list(
      entropy = .multilpa_table(function(x) x$entropy),
      classification = .multilpa_table(function(x) x$classification),
      average_posteriors = .multilpa_table(function(x) x$average_posteriors),
      residuals = .multilpa_table(.multilpa_diagnostic_residuals)))
  }
  if (.multilpa_any_summary(x)) return(.multilpa_stored_catalogue(x))
  stopifnot("`x` must be an object of this package" = .multilpa_any_fit(x))
  c(.multilpa_measurement_catalogue(),
    .latents_extra_catalogue(x),
    .multilpa_mixing_catalogue(x),
    .multilpa_unit_catalogue(x),
    .multilpa_diagnostic_catalogue(),
    .multilpa_sequence_catalogue(x),
    .multilpa_provenance_catalogue(x))
}

#' One catalogue entry
#' @param build A function of the object and that table's own arguments.
#' @param args The argument names it accepts.
#' @return A list of the two.
#' @noRd
.multilpa_table <- function(build, args = character()) {
  list(build = build, args = args)
}

#' Is this the summary of a fitted model of this package?
#' @param x Any object.
#' @return A single logical.
#' @noRd
.multilpa_any_summary <- function(x) {
  inherits(x, c("summary_multilpa", "summary_multilpa_covariates",
                "summary_multilpa_transitions",
                "summary_multilpa_enumeration",
                "summary_multilpa_bootstrap_lrt"))
}

#' A catalogue over the tables a summary already carries
#'
#' A summary is built by asking the fit for every table, so its own catalogue
#' serves them back rather than recomputing them. That is what makes
#' `get_results(summary(fit), what)` and `get_results(fit, what)` the same table
#' rather than two computations that could diverge.
#'
#' @param x A summary object.
#' @return A named list of catalogue entries.
#' @noRd
.multilpa_stored_catalogue <- function(x) {
  tables <- x$tables
  if (!is.list(tables) || length(tables) == 0L) {
    stop(errorCondition(paste(
      "This summary carries no tables, so it was made by an older version of",
      "the package. Take the tables from the fit instead."),
      class = "latents_incomplete_fit", call = NULL))
  }
  lapply(tables, function(table) {
    force(table)
    .multilpa_table(function(x) table)
  })
}

#' The measurement tables every fitted family offers
#' @return A named list of catalogue entries.
#' @noRd
.multilpa_measurement_catalogue <- function() {
  list(
    profiles = .multilpa_table(
      function(x, data = NULL, scale = c("raw", "standardized"))
        .multilpa_profile_frame(x, data = data, scale = match.arg(scale),
                                standard_errors = TRUE),
      c("data", "scale")),
    responses = .multilpa_table(
      function(x, data = NULL) .multilpa_response_frame(x, data = data,
                                                        standard_errors = TRUE),
      "data"),
    covariances = .multilpa_table(.multilpa_covariance_frame))
}

#' The mixing-structure tables, which differ by model family
#' @param x A fitted model of this package.
#' @return A named list of catalogue entries.
#' @noRd
.multilpa_mixing_catalogue <- function(x) {
  counts <- list(counts = .multilpa_table(.multilpa_count_frame))
  if (inherits(x, "multilpa_transitions")) {
    return(c(list(
      transitions = .multilpa_table(
        function(x, estimated = NULL, stable = NULL)
          .multilpa_transitions_table(x, estimated = estimated,
                                      stable = stable),
        c("estimated", "stable")),
      initial = .multilpa_table(.multilpa_initial_frame)), counts))
  }
  if (inherits(x, "multilpa_covariates")) {
    return(c(list(coefficients = .multilpa_table(.multilpa_coefficient_frame)),
             counts))
  }
  c(list(profile_probabilities = .multilpa_table(.multilpa_probability_frame)),
    counts)
}

#' The one-row-per-unit tables
#' @param x A fitted model of this package.
#' @return A named list of catalogue entries.
#' @noRd
.multilpa_unit_catalogue <- function(x) {
  posteriors <- list(posteriors = .multilpa_table(
    function(x, format = c("long", "wide"))
      .multilpa_posterior_frame(x, match.arg(format)),
    "format"))
  assignments <- list(assignments = .multilpa_table(
    function(x, data = NULL, truth = NULL)
      .multilpa_hide_fabricated_id(
        x, .multilpa_assignments(x, data = data, truth = truth)),
    c("data", "truth")))
  c(posteriors,
    list(group_posteriors = .multilpa_table(
      function(x, format = c("long", "wide"))
        .multilpa_group_posterior_frame(x, match.arg(format)),
      "format")),
    assignments)
}

#' The classification and model-fit tables every fitted family offers
#' @return A named list of catalogue entries.
#' @noRd
.multilpa_diagnostic_catalogue <- function() {
  # `level` follows the fit rather than a fixed default: a two-level model is
  # asked about at both levels, and a family with one level has only one to
  # report. This is the rule `diagnostics()` already applies.
  by_level <- function(verb) {
    .multilpa_table(function(x, level = NULL)
      verb(x, level = .multilpa_default_level(x, level)), "level")
  }
  list(
    classification = by_level(.multilpa_classification_table),
    average_posteriors = by_level(.multilpa_average_posteriors),
    classification_errors = by_level(.multilpa_classification_errors),
    bch_weights = by_level(.multilpa_bch_weights),
    entropy = .multilpa_table(.multilpa_entropy_table),
    residuals = .multilpa_table(
      function(x, data = NULL, by = c("profile", "overall"),
               adjust = .multilpa_p_adjust_methods)
        .multilpa_bivariate_residuals(x, data = data, by = match.arg(by),
                                      adjust = match.arg(adjust)),
      c("data", "by", "adjust")),
    information_criteria = .multilpa_table(
      function(x, definitions = FALSE, format = c("wide", "long"))
        .multilpa_information_criteria(x, definitions = definitions,
                                       format = match.arg(format)),
      c("definitions", "format")),
    model = .multilpa_table(.multilpa_fit_frame))
}

#' The occasion-ordered tables
#' @param x A fitted model of this package.
#' @return A named list of catalogue entries.
#' @noRd
.multilpa_sequence_catalogue <- function(x) {
  shared <- list(
    sequences = .multilpa_table(
      function(x, format = c("long", "wide"))
        .multilpa_sequences(x, format = match.arg(format)),
      "format"),
    sequence_summary = .multilpa_table(.multilpa_sequence_summary))
  if (!inherits(x, "multilpa_transitions")) return(shared)
  c(shared,
    list(sequence_lengths = .multilpa_table(.multilpa_sequence_length_frame)))
}

#' The tables that record how the fit was produced
#' @param x A fitted model of this package.
#' @return A named list of catalogue entries.
#' @noRd
.multilpa_provenance_catalogue <- function(x) {
  stages <- if (inherits(x, "multilpa"))
    list(stages = .multilpa_table(.multilpa_stage_frame)) else list()
  c(stages,
    list(starts = .multilpa_table(function(x) x$starts),
         data = .multilpa_table(function(x)
           .multilpa_hide_fabricated_id(x, .multilpa_model_frame(x)))))
}

#' Hide the unit column a single-level fit fabricated for itself
#'
#' `id = NULL` numbers the rows so the two-level machinery has a unit to group
#' by, and every internal rebuild of the fitting data needs that column --
#' `parameter_inference()` checks the frame against it. What the caller never
#' wrote, though, they should not be handed back, so it is dropped here, at the
#' table the caller reads, rather than at the rebuild the package relies on.
#'
#' @param x A fitted model of this package.
#' @param frame A table that may carry the fabricated column.
#' @return `frame`, without that column when the fit fabricated one.
#' @noRd
.multilpa_hide_fabricated_id <- function(x, frame) {
  if (!isTRUE(x$single_level) || !x$id %in% names(frame)) return(frame)
  frame[setdiff(names(frame), x$id)]
}

#' Which classification levels to report when the caller named none
#' @param x A fitted model of this package.
#' @param level The caller's `level`, possibly `NULL`.
#' @return A single level name.
#' @noRd
.multilpa_default_level <- function(x, level) {
  if (!is.null(level)) return(level)
  if (is.null(x$group_posteriors)) "individuals" else "both"
}

#' The bivariate residuals a diagnostics object holds, or the reason it has none
#' @param x A `multilpa_diagnostics` object.
#' @return The residual table.
#' @noRd
.multilpa_diagnostic_residuals <- function(x) {
  if (!is.null(x$residuals)) return(x$residuals)
  if (isTRUE(x$fit$noise)) {
    stop(errorCondition(paste(
      "This fit has no bivariate residuals: they do not account for a noise",
      "component yet."),
      class = "latents_unsupported_noise", call = NULL))
  }
  stop(errorCondition(
    "This fit has no bivariate residuals: they need a discrete group-class model.",
    class = "latents_no_group_classes", call = NULL))
}

#' How many occasions each group of a transition fit contributes
#' @param x A fitted `multilpa_transitions` model.
#' @return One row per group.
#' @noRd
.multilpa_sequence_length_frame <- function(x) {
  data.frame(group = x$group_values, group_class = x$group_classes,
             occasions = unname(x$sequence_lengths),
             observations = unname(tabulate(x$group_index, nbins = x$n_groups)),
             complete = unname(x$sequence_lengths) == x$n_occasions,
             row.names = NULL)
}

#' Print every table an object carries
#'
#' A summary is the whole of what an object holds, so its `print()` shows every
#' table rather than a chosen few. A table longer than `rows` is shown to that
#' depth with its remaining row count and the call that returns it whole: the
#' alternative is either a console dump or a table the reader cannot tell was
#' truncated.
#'
#' @param tables A named list of `data.frame`s, as `get_results(x, "all")` gives.
#' @param rows How many rows of each table to show.
#' @param digits Printed significant digits.
#' @return `NULL`, invisibly. Called for the printing.
#' @noRd
.multilpa_print_tables <- function(tables, rows = 10L, digits = 4L) {
  stopifnot(
    "`rows` must be a single non-negative whole number" =
      is.numeric(rows) && length(rows) == 1L && is.finite(rows) && rows >= 0 &&
      rows == as.integer(rows),
    "the tables must be a named list of data frames" =
      is.list(tables) && !is.null(names(tables)))
  invisible(lapply(names(tables), function(name) {
    frame <- tables[[name]]
    cat(sprintf("\n-- %s %s\n", name,
                strrep("-", max(3L, 64L - nchar(name)))))
    if (nrow(frame) == 0L) {
      cat("   (no rows)\n")
      return(invisible(NULL))
    }
    # `rows = 0` asks for the headings alone; printing a zero-row head would
    # dump the column names as if the table were empty.
    if (rows > 0) print(utils::head(frame, rows), digits = digits, row.names = FALSE)
    if (nrow(frame) > rows) {
      cat(sprintf("   ... %d more rows.  get_results(x, what = \"%s\")\n",
                  nrow(frame) - rows, name))
    }
    invisible(NULL)
  }))
}

#' The closing line of a summary that has printed its tables
#' @param tables The tables that were printed.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_print_table_footer <- function(tables) {
  cat(sprintf(paste0(
    "\n%d tables above, truncated to fit. get_results(x, what = ) returns any\n",
    "of them whole, and get_results(x, what = \"all\") returns every one.\n"),
    length(tables)))
  invisible(NULL)
}

#' Validate the two printing controls every summary print method takes
#' @param digits The caller's `digits`.
#' @param rows The caller's `rows`.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_check_print_arguments <- function(digits, rows) {
  stopifnot(
    "`digits` must be a single number between 1 and 22" =
      is.numeric(digits) && length(digits) == 1L && is.finite(digits) &&
      digits >= 1 && digits <= 22,
    "`rows` must be a single non-negative whole number" =
      is.numeric(rows) && length(rows) == 1L && is.finite(rows) && rows >= 0 &&
      rows == as.integer(rows))
  invisible(NULL)
}

#' The profile means as one row per profile
#'
#' The long measurement table is one row per profile *and* indicator, which is
#' the right shape to compute on and the wrong one to read: three profiles over
#' five indicators is fifteen lines of console for a table that has three rows.
#' This is the shape the means are read in. It is built by reshaping the long
#' table rather than from `x$means`, so the two cannot report different numbers.
#'
#' @param x A fitted model of this package.
#' @return A `data.frame` with one row per profile and one column per
#'   continuous indicator.
#' @noRd
.multilpa_wide_means <- function(x) {
  long <- .multilpa_profile_frame(x)
  wide <- stats::reshape(long[c("profile", "indicator", "mean")],
                         idvar = "profile", timevar = "indicator",
                         direction = "wide")
  names(wide) <- sub("^mean\\.", "", names(wide))
  row.names(wide) <- NULL
  wide
}

#' Print a fit's primary table under its header
#'
#' A fitted model auto-prints its own estimates. The alternative is a header
#' that describes a result the reader then has to fetch with a second call,
#' which is the ritual this package exists to remove: `fit` and
#' `as.data.frame(fit)` would otherwise be two calls showing one thing.
#'
#' The measurement model is the primary table, and an all-categorical fit has
#' none, so that fit prints its response probabilities instead of an empty
#' block under a heading promising means.
#'
#' @param x A fitted model of this package.
#' @param rows How many rows to print before truncating.
#' @return `NULL`, invisibly. Called for the printing.
#' @noRd
.multilpa_print_primary <- function(x, rows = 20L) {
  stopifnot("`rows` must be a single non-negative whole number" =
              is.numeric(rows) && length(rows) == 1L && is.finite(rows) &&
              rows >= 0 && rows == as.integer(rows))
  name <- names(.multilpa_catalogue(x))[1L]
  # Printing shows the estimates; the standard errors belong to the tables
  # asked for by name, and computing them on every print would be slow.
  table <- if (identical(name, "profiles")) .multilpa_profile_frame(x) else
    get_results(x)
  ## The Gaussian measurement is printed one row per profile. Everything else
  ## keeps the catalogue's own shape, which is already one row per thing.
  if (identical(name, "profiles") && nrow(table) > 0L) {
    wide <- .multilpa_wide_means(x)
    ## How big each profile is belongs beside what it looks like: a mean is read
    ## differently when it describes two percent of the sample. Matched on the
    ## class label rather than on row order, which the two tables need not share.
    sizes <- .multilpa_count_frame(x)
    individuals <- sizes[sizes$level == "individuals", , drop = FALSE]
    at <- match(wide$profile, individuals$class)
    wide$count <- individuals$effective_count[at]
    wide$proportion <- individuals$effective_proportion[at]
    cat("\n")
    if (rows > 0) print(utils::head(wide, rows), row.names = FALSE)
    if (nrow(wide) > rows) {
      cat(sprintf("   ... %d more profiles.\n", nrow(wide) - rows))
    }
    cat("\nVariances and standard errors: get_results(x, \"profiles\").",
        "\nEvery other table: get_results(x, what = ), or get_results(x, \"all\").\n")
    return(invisible(NULL))
  }
  if (nrow(table) == 0L && length(x$categorical %||% character()) > 0L) {
    table <- .multilpa_response_frame(x)
    name <- "responses"
  }
  if (nrow(table) == 0L) return(invisible(NULL))
  cat("\n")
  if (rows > 0) print(utils::head(table, rows), row.names = FALSE)
  if (nrow(table) > rows) {
    cat(sprintf("   ... %d more rows.  get_results(x, \"%s\")\n",
                nrow(table) - rows, name))
  }
  cat("\nEvery other table: get_results(x, what = ), or get_results(x, \"all\").\n")
  invisible(NULL)
}
