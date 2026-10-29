#' Extract model estimates and analysis results
#'
#' @description
#' `get_results()` provides access to the estimates and summaries produced by
#' models and analyses in `latents`. It supports examination of the estimated
#' latent classes, their membership probabilities, and model-fit statistics,
#' as well as the preparation of results for reporting and further analysis.
#'
#' The results available depend on the supplied object. Fitted models provide
#' measurement estimates, class proportions, posterior membership
#' probabilities, and other model-specific quantities, such as transition
#' probabilities. Enumeration and comparison objects provide summaries of
#' candidate models, while diagnostic and bootstrap objects provide their
#' corresponding checks and test results.
#'
#' The `what` argument selects the requested table. Each table is returned as
#' a standard R data frame, allowing its contents to be inspected, filtered,
#' combined, or exported. Omitting `what` returns the object's primary table;
#' `what = "all"` returns a named collection of the tables included by default
#' for that object.
#'
#' @details
#' Some tables calculate additional quantities when requested. For example,
#' measurement tables may obtain standard errors through
#' [parameter_inference()]. Their availability therefore depends on the fitted
#' model and the requirements of the corresponding calculation.
#'
#' The table catalogue below specifies the contents and accepted arguments of
#' each table, including results that require an explicit request.
#'
#' Every object also has `print()`, which shows a compact overview,
#' `summary()`, which collects its main tables (taking `level` and
#' `vcov_type` where the model has standard errors), and `as.data.frame()`,
#' which returns its primary table, the same as `get_results(x)`.
#'
#' Before combining supplied `data` with model results, the function checks
#' row alignment using columns retained by the fitted object, such as
#' identifiers, occasions, and indicators. Numeric values serialized by
#' [utils::write.csv()] are also accepted when serialization preserves
#' distinct values. Columns absent from the fitted object cannot establish
#' alignment. Inference and model comparison apply additional checks to
#' verify that supplied data reproduce the estimation data.
#'
#' @param x A fitted model of this package, an enumeration grid, a model
#'   comparison, a bootstrap comparison, a starting-value set, a
#'   classification-diagnostics object, or the summary of any of them.
#' @param what Name of the table to return, or `"all"` for the collection of
#'   tables included by default for the object. `NULL` (default) returns the primary table: the
#'   measurement model for a fit, the candidate grid for an enumeration, the
#'   test row for a bootstrap comparison and the entropy table for a
#'   diagnostics object. A name the object does not offer raises
#'   `latents_bad_argument`, listing the names it does offer.
#' @param ... Arguments for the requested table. Each table accepts only the
#'   arguments listed for it under *Tables*; any other argument raises
#'   `latents_bad_argument`.
#'
#' @section Tables:
#' The following tables apply to profile and class models and, where
#' supported, membership-covariate and transition models. Tables for other
#' model families are described in *Tables of the other model families*.
#' \describe{
#'   \item{`"profiles"`}{The Gaussian measurement model, one row per profile
#'     and continuous indicator: `profile`, `indicator`, `mean`, `variance`,
#'     `standard_deviation`, `mean_standard_error`, `variance_standard_error`.
#'     Zero rows, with the same columns, when every indicator is categorical.
#'     Takes `scale` and `data`.}
#'   \item{`"responses"`}{The categorical measurement model, one row per
#'     profile, categorical indicator and category: `profile`, `indicator`,
#'     `category`, `probability`, `threshold`, `probability_standard_error`.
#'     `threshold` is `qlogis(P(y <= category))` and is `NA_real_` for each
#'     indicator's final category; a probability on its bound has an
#'     `NA_real_` standard error. Zero rows, with the same columns, when no
#'     indicator is categorical. Takes `data`.}
#'   \item{`"covariances"`}{The within-profile residual covariance matrices,
#'     one row per profile and ordered pair of continuous indicators:
#'     `profile`, `indicator`, `indicator_2`, `covariance`. Under the diagonal
#'     parameterization the off-diagonal covariances are reported as exact
#'     zeros.}
#'   \item{`"ordinal"`, `"count_means"`}{Fits with ordinal or count
#'     indicators only: the ordinal and count measurement models. Take
#'     `data`.}
#'   \item{`"counts"`}{Effective class sizes, one row per class at each level
#'     of the fit: `level`, `class`, `effective_count` (the summed posteriors)
#'     and `effective_proportion`. Effective counts use posterior probabilities;
#'     modal counts in `"classification"` use discrete class assignments.}
#'   \item{`"model"`}{One row describing the fit: its dimensions, likelihood,
#'     parameter counts, `aic`, `bic_groups`, `bic_individual` and its
#'     convergence and boundary diagnostics. The columns differ by family (a
#'     covariate fit adds its predictor counts); `"information_criteria"` has
#'     the same columns for every family.}
#'   \item{`"posteriors"`}{One row per individual and profile: `row`, `group`,
#'     `profile`, `posterior`, `modal`. Takes `format`. Excluded observations
#'     retain their rows, with `posterior` and `modal` set to `NA`. In wide
#'     format their assigned `profile` and posterior columns are `NA`.}
#'   \item{`"assignments"`}{One row per observation: the columns the fit was
#'     built from (or those of `data`), then `profile`, `group_class` where the
#'     model has one, `uncertainty` (one minus the posterior of the assigned
#'     profile) and one `posterior_profile_*` column per profile, plus
#'     `posterior_noise` for a fit with a noise component (whose cases are
#'     profile 0). Takes `data` and `truth`; with `truth` it returns a recovery
#'     table instead (see *Recovery*). Row order is preserved. Under the default
#'     `drop_empty = TRUE`, observations with no observed indicators have `NA`
#'     for `profile`, `uncertainty` and the posterior columns. Group assignments
#'     are missing only when the entire group is excluded; an empty observation
#'     can therefore retain a group assignment supported by other observations.
#'     Zero-weight units are unclassified under either `drop_empty` setting.}
#'   \item{`"classification"`}{Classification quality, one row per level and
#'     class: `level`, `class`, `n_modal`, `proportion_modal`, `estimated_n`,
#'     `estimated_proportion`, `average_posterior`,
#'     `odds_correct_classification`. Takes `level`. Under sampling weights
#'     every column that estimates the population is weighted:
#'     `estimated_n` and `estimated_proportion` are weighted posterior totals
#'     and `average_posterior` a weighted mean, so the odds of correct
#'     classification compare the two on the same population; `n_modal` and
#'     `proportion_modal` stay counts of the units supplied. The odds are `Inf` for a
#'     class classified perfectly (average posterior one to within rounding)
#'     and `NA` where they are undefined: a class with no modal members, or a
#'     single class, whose proportion is one.}
#'   \item{`"average_posteriors"`}{For the units modally assigned to each
#'     class, their mean posterior for every class: `level`, `assigned_class`,
#'     `class`, `n_assigned`, `average_posterior`. Rows sum to one within each
#'     assigned class. Under sampling weights `average_posterior` is the
#'     weighted mean and `n_assigned` stays a count. This conditions on the assigned class, whereas
#'     `"classification_errors"` conditions on the true class. Takes `level`.}
#'   \item{`"classification_errors"`}{The probability that a unit in one true
#'     class is assigned to another: `level`, `true_class`, `assigned_class`,
#'     `probability`, summing to one within each `true_class`. Takes `level`.
#'     It is `sum_i w_i p_it 1(W_i = s) / sum_i w_i p_it` for true class `t`
#'     and assigned class `s` (Vermunt 2010), with `w_i` the unit's sampling
#'     weight under sampling weights and 1 otherwise, so a weighted fit reports
#'     the population's error rates.
#'     On a membership-covariate fit this is a marginal summary over the fitted
#'     units, not conditional on the covariates. A fit with a noise component
#'     raises `latents_unsupported_noise`, as [three_step()] does, and `"all"`
#'     and `summary()` omit the table.}
#'   \item{`"bch_weights"`}{The BCH inverse-error weights, one row per unit and
#'     class: `level`, `unit`, `assigned_class`, `class`, `weight`. Takes
#'     `level`; a weighted regression should use one level, not `"both"`.
#'     They invert the `"classification_errors"` matrix, weighted as it is. On a
#'     membership-covariate fit the weights invert the marginal error matrix
#'     and are not a conditional three-step correction. Not available for a fit
#'     with a noise component, as for `"classification_errors"`.}
#'   \item{`"entropy"`}{One row per level: `level`, `n_classes`, `n_units`,
#'     `entropy_sum`, `relative_entropy`. Under sampling weights
#'     `entropy_sum` is `sum_i w_i e_i` and `relative_entropy`
#'     `1 - sum_i w_i e_i / (sum_i w_i log K)`; `n_units` stays a count.}
#'   \item{`"replication"`}{One row: `n_starts` (starts run to the end),
#'     `n_best_replicated` (those whose log likelihood is within `tolerance`
#'     of the best), `share`, `tolerance`, `best_start` (the selected start's
#'     index) and `log_likelihood`. A best solution reached by one start of
#'     several may not be the global maximum: refit with more starts.}
#'   \item{`"residuals"`}{Bivariate residuals for every pair of indicators,
#'     checking within-profile independence: `profile`, `indicator_1`,
#'     `indicator_2`, `kind` (`"gaussian"`, `"categorical"`, `"mixed"` or
#'     `"count"`), `observed`, `expected`, `residual`, `effective_n`,
#'     `statistic`, `df`, `p_value`, `p_adjusted`. Takes `data`, `by` and
#'     `adjust`. See *Bivariate residuals* for what each kind computes.}
#'   \item{`"information_criteria"`}{Every information criterion the package
#'     computes (see *Information criteria*). Takes `format` and
#'     `definitions`.}
#'   \item{`"sequences"`}{The assignments in occasion order, one row per
#'     classified observation (a row with nothing observed is left out, or
#'     kept under `drop_empty = FALSE`): `group`, `group_class`, `time`,
#'     `profile`. Takes `format`.
#'     Raises `latents_no_time` for a fit made without `time`.}
#'   \item{`"sequence_summary"`}{One row per group class: `group_class`,
#'     `groups`, `observations`, `mean_length`, `median_length`, `shortest`,
#'     `longest`, `complete`, `gaps`. Summarises classified groups only;
#'     lengths describe their recorded occasions, including occasions with
#'     missing indicators. Raises `latents_no_time` for a fit made
#'     without `time`.}
#'   \item{`"starts"`}{One row per EM start: `start`, `log_likelihood`,
#'     `converged`, `iterations`, `error` and, where the family records it,
#'     `boundary`. Under `search = "emEM"`, one row per start of the pool,
#'     with its `kind`, its `stage` (`"final"` when run to convergence,
#'     `"screened"` otherwise) and its `screened_log_likelihood`.}
#'   \item{`"data"`}{The data the model was fitted to, one row per observation
#'     in input order: the identifier (except for a single-level fit), the
#'     occasion where there is one, the sampling weights of a weighted fit (as
#'     rescaled for fitting) and the indicators under their original names.
#'     Columns the model did not use, including membership covariates, are not
#'     included; [three_step()] and [r3step()] therefore take `data`.}
#' }
#' Fits with discrete group classes also offer:
#' \describe{
#'   \item{`"group_posteriors"`}{One row per observed group and group class:
#'     `group`, `group_size`, `log_likelihood`, `group_class`, `posterior`,
#'     `modal`. Takes `format`. Excluded groups retain their rows, with missing
#'     posterior probabilities and modal assignments, as for `"posteriors"`.}
#' }
#' Tables specific to one family or object:
#' \describe{
#'   \item{`"gaussian_diagnostics"`}{Gaussian profile, membership-covariate,
#'     transition and cross-level fits that record these diagnostics. One row
#'     per Gaussian component and measurement context: `context`, `profile`,
#'     `gaussian_mass`, `positive_rows`, `weight_effective_n`,
#'     `min_variance_ratio`, `min_standardized_eigenvalue`, `boundary`.
#'     `gaussian_mass` sums posterior times normalized sampling weight over
#'     rows with at least one observed Gaussian indicator; `positive_rows`
#'     counts those rows (the same for every component, since a Gaussian
#'     density is positive everywhere). `weight_effective_n` is the
#'     concentration measure `(sum(a)^2 / sum(a^2))` for those contributions,
#'     not an identification guarantee; Gaussian mass is not this effective
#'     sample size. `min_variance_ratio` is the smallest component variance
#'     divided by its observed sampling-weighted marginal variance on the
#'     fitted scale, including between-component variation.
#'     `min_standardized_eigenvalue` is the smallest eigenvalue of the
#'     component covariance after standardizing every indicator by its
#'     observed sampling-weighted marginal standard deviation. Missing reference
#'     variation gives `NA` ratios. `boundary` concerns the absolute variance
#'     or eigenvalue floor. These descriptive measures use no automatic spike
#'     cutoff and do not alter start selection.}
#'   \item{`"profile_probabilities"`, `"stages"`}{A `multilpa` fit.
#'     `"profile_probabilities"` has one row per group class and profile:
#'     `group_class`, `profile`, `probability`, `group_class_probability`.
#'     `"stages"` has one row per estimation stage: `stage`, `group_classes`,
#'     `fixed`, `log_likelihood`, `parameters`, `parameters_with_measurement`,
#'     `converged`. An ordinary fit has one `"joint"` row and a [fit_staged()]
#'     result has two; `fixed` is `NA_character_` for a stage that estimated
#'     every block.}
#'   \item{`"coefficients"`}{A covariate fit. One row per estimated membership
#'     coefficient: `level` (`"profile"` or `"group"`), `outcome`, `term`,
#'     `parameter` (always `"logit"`), `estimate`. Standard errors are given by
#'     [parameter_inference()].}
#'   \item{`"transitions"`, `"initial"`, `"sequence_lengths"`}{A transition
#'     fit. `"transitions"` has one row per group class and ordered pair of
#'     profiles: `group_class`, `from`, `to`, `probability`, `expected_count`,
#'     `stable`, `estimated`, `group_class_probability`; it takes `estimated`
#'     and `stable`. `"initial"` has one row per group class and profile:
#'     `group_class`, `profile`, `probability`, `prevalence`,
#'     `group_class_probability`. `"sequence_lengths"` has one row per group:
#'     `group`, `group_class`, `occasions`, `observations`, `complete`.}
#'   \item{`"candidates"`, `"best"`, `"classification"`}{An enumeration
#'     grid. `"candidates"` has one row per candidate model: its name
#'     (`candidate`, first), its class counts, covariance `model`, log
#'     likelihood, parameter count, every criterion under both sample-size
#'     conventions, both entropies, the average posterior probability (AvePP)
#'     summaries of each level (`avepp_min`, `avepp_max`, `avepp_mean`,
#'     `avepp_below` and the same four prefixed `group_`) and the convergence,
#'     boundary, replication, warning and error diagnostics, including
#'     `min_standardized_eigenvalue` and `min_gaussian_mass` where Gaussian
#'     diagnostics are available. These are descriptive summaries, without
#'     an automatic cutoff or change to candidate selection. Failed candidates
#'     are kept with `NA` estimates and their error text. It takes `cutoff`
#'     (default `0.8`): `avepp_below` and `group_avepp_below` count the classes
#'     whose AvePP is below it. A class with no modally assigned unit has no
#'     AvePP; it is excluded from the minimum, maximum and mean and counted as
#'     below the cutoff. A level with one class has `NA` summaries. `"best"`
#'     has one row per information criterion and sample-size convention:
#'     `criterion`, `convention`, `candidate` (the converged candidate with the
#'     lowest value), `n_profiles`, `n_group_classes`, `model` and `value`;
#'     `what = "criteria"` is accepted as a deprecated name for it.
#'     `"classification"` has one row per fitted candidate, level and class:
#'     `candidate`, `n_profiles`, `n_group_classes`, `model`, then the columns
#'     of a fit's `"classification"` table and that level's
#'     `relative_entropy`. Failed candidates and single-class levels
#'     contribute no rows.}
#'   \item{`"deltas"`, `"knees"`}{Model-selection evidence of an enumeration
#'     grid, each taking `using` (default `"bic_groups"`). Neither selects a
#'     model. `"deltas"` has one row per candidate: `candidate`, its class
#'     counts and `model`, `criterion`, `value`, `delta_profile` (the decrease
#'     in the criterion from the converged candidate with one profile fewer,
#'     the same group classes and the same structure; positive means this
#'     candidate is lower), `from_profile` naming that candidate, and
#'     `delta_group_class` with `from_group_class` for one group class fewer.
#'     `value` and both deltas are `NA` for a candidate that did not converge,
#'     and a delta is `NA` when the smaller neighbour is missing or did not
#'     converge. `"knees"` has one row per line of the grid: `direction`
#'     (`"profiles"`, along the profiles at fixed group classes;
#'     `"group_classes"`, along the group classes at fixed profiles;
#'     `"parameters"`, every candidate by its number of parameters, one line
#'     per structure), `line`, `criterion`, `from`, `to`, `n_points`,
#'     `knee_candidate`, `distance`, `knee_without_largest` (the knee when the
#'     last candidate of the line is left out), `stable` and `note`. The knee
#'     is the candidate furthest below the straight line from the first to the
#'     last candidate, measured perpendicular to it with both axes rescaled to
#'     0-1. This is a simplified rule inspired by Kneedle (Satopaa et al.,
#'     2011), without its smoothing or sensitivity thresholds, so a very small positive distance
#'     can still name a knee. A line needs three candidates for a knee and four
#'     for `knee_without_largest` and `stable`; otherwise they are `NA` and
#'     `note` gives the reason.}
#'   \item{`"local_dependence"`}{For a fit with `local_dependence`: `pair`,
#'     `kind` (`"categorical"` or `"continuous"`), `profile`, `item_1`,
#'     `category_1`, `item_2`, `category_2`, `probability`, `independent`,
#'     `ratio`, `covariance` and `correlation`. A categorical pair has one row
#'     per profile and category combination, including unobserved combinations
#'     subject to the `min_probability` floor and included in the parameter
#'     count: the fitted
#'     joint `probability`, the product of the two marginal probabilities
#'     (`independent`, the joint probability under independence within the
#'     profile) and their `ratio`; its covariance columns are `NA`. A
#'     continuous pair has one row per profile, with the pair's `covariance`
#'     and the `correlation` it implies; its category and probability
#'     columns are `NA`. A fit without pairs raises
#'     `latents_no_local_dependence`, and `what = "all"` leaves the table out.}
#'   \item{`"direct_effects"`}{For a fit with `direct_effects`, one row per
#'     effect: `indicator`, `covariate`, `estimate` (the change in the
#'     indicator's mean per unit of the covariate, the same in every profile),
#'     `standard_error`, `statistic`, `p_value`, `conf_low` and `conf_high`
#'     (Wald, from [parameter_inference()]). A fit without direct effects
#'     raises `latents_no_direct_effects`, and `what = "all"` leaves the table
#'     out.}
#'   \item{`"lmr"`}{The Vuong-Lo-Mendell-Rubin (VLMR) and adjusted
#'     Lo-Mendell-Rubin (aLMR) tests of an enumeration grid (Vuong 1989; Lo,
#'     Mendell and Rubin 2001; Vermunt 2024), one row per pair of converged
#'     candidates that differ by one profile (`direction = "profiles"`) or one
#'     group class (`"group_classes"`) at the same structure: `direction`,
#'     `null_candidate`, `alternative_candidate`, `n_profiles`,
#'     `n_group_classes` and `model` of the alternative, `n_units` (the
#'     independent units scored: groups, or observations without groups),
#'     `parameters_added`, `lr` (twice the log-likelihood difference),
#'     `vlmr_p`, `adjusted_lr` (`lr / (1 + 1 / (parameters_added * log(n_units)))`),
#'     `almr_p`, `reference_mean` and `reference_sd` of the weighted
#'     chi-square reference distribution, whose weights are the eigenvalues of
#'     Vuong's matrix of the two models' scores and Hessians, and `note`. The
#'     p-values are computed by Imhof's method to about `1e-9`. A pair whose
#'     information matrix is singular, or whose candidate is on a boundary,
#'     keeps its row with `NA` statistics and the reason in `note`; warnings
#'     raised while scoring a candidate are recorded there too. Every
#'     candidate's information is computed, so the table is built only when
#'     asked for, and `what = "all"` leaves it out. A grid with no such pair
#'     raises `latents_no_lmr`.}
#'   \item{`"tree"`}{The flows of [plot_tree()], on every grid, taking its
#'     arguments `n_group_classes`, `n_profiles` (`n_classes` on a regression
#'     grid), `model`, `candidates`, `from`, `to` and `level`
#'     (default `"individuals"`). One row per pair of classes in consecutive
#'     rows of each tree: `tree` (the tree's name, as in the plot list),
#'     `level`, `from_candidate`, `from_class`, `from_share`, `to_candidate`,
#'     `to_class`, `to_share`, `mass` (the share of all units the two classes
#'     share, the sum over units of the product of their posterior
#'     probabilities, weighted by sampling weights at the individual level),
#'     `share_of_from` and `share_of_to` (`mass` as a part of each class).
#'     Every pair is listed, including those that share no mass. Units are
#'     the rows and groups the fit counts (see `drop_empty`), weighted by
#'     their sampling weights at both levels. A grid with no tree raises
#'     `latents_no_tree`, a fit with a noise component
#'     `latents_unsupported_noise`, and `what = "all"` leaves the table out.
#'     The summary of a grid carries the default table only, without these
#'     arguments.}
#'   \item{`"comparison"`, `"classification"`}{A [compare_models()] result.
#'     `"comparison"` is its table, one row per model. `"classification"`,
#'     offered when every model is a profile, class, membership-covariate or
#'     transition model, has one row per model, level and class: `model`, then
#'     the columns of a fit's `"classification"` table and that level's
#'     `relative_entropy`; single-class levels are left out. The models
#'     themselves are returned by [get_candidate()].}
#'   \item{`"test"`, `"replicates"`}{A bootstrap comparison. `"test"` is the
#'     one-row result; `"replicates"` has one row per simulated dataset.}
#'   \item{`"profiles"`, `"responses"`, `"profile_probabilities"`,
#'     `"covariances"`}{A starting-value set: its Gaussian means and
#'     variances, response probabilities, profile probabilities and, one row
#'     per profile and ordered pair of continuous indicators, its covariances
#'     (`profile`, `indicator`, `indicator_2`, `covariance`).}
#'   \item{`"entropy"`, `"classification"`, `"average_posteriors"`,
#'     `"residuals"`, `"replication"`}{A [diagnostics()] object: the tables of the same names
#'     described above.}
#' }
#'
#' @section Tables of the other model families:
#' Each family's first table is its default `what`; `"all"` returns a named
#' list of every table (without `"recovery"`, which needs `data` and `truth`).
#' Tables with intervals take `level`, and `vcov_type` (`"observed"`,
#' `"robust"` or `"opg"`; `NULL` reuses the type stored with a regression
#' fit). These objects also have `coef()`, `vcov(type = )`,
#' `confint(parm = , level = )`, `logLik()` and `nobs()` methods. `confint()`
#' of a mixture regression or growth mixture uses the covariance stored with
#' the fit and refuses any other argument with `latents_bad_argument` (an
#' interval on another covariance is `get_results(x, "coefficients",
#' vcov_type = )`); `vcov()` warns `latents_unused_argument` for an argument
#' it does not take, such as `vcov_type =` in place of `type =`.
#' \describe{
#'   \item{Mixture regressions ([mixture_regression()])}{`"coefficients"`
#'     (one row per class and term, with Wald statistics, `p_adjusted` by
#'     Benjamini-Hochberg over non-intercepts and `exp_estimate` for the
#'     binomial, Poisson and ordinal families; ordinal thresholds lead each
#'     class as terms `"threshold:<lower>|<upper>"`), `"classes"`,
#'     `"membership"` (multinomial logits and odds ratios), `"group_classes"`
#'     (rows for two-level fits only), `"fit"` (`bic` penalised by the number
#'     of independent units, `bic_rows` by the number of rows, as flexmix
#'     reports it), `"assignments"`, `"groups"`, `"fitted"`, `"starts"`
#'     (`stage` `"screened"` or `"completed"`), `"classification"` and
#'     `"recovery"` (`by = "class"` or `"group_class"`).}
#'   \item{Growth mixtures (`mixture_regression(random = )`)}{
#'     `"coefficients"`, `"classes"` (each class's share with its interval,
#'     `expected_persons`, `assigned_persons`, `average_posterior` and
#'     `residual_sd`), `"membership"`, `"fit"`, `"assignments"` (one row per
#'     person: `class`, `probability` and one `posterior_*` column per
#'     class), `"starts"` and `"recovery"`, as for a mixture regression, plus
#'     `"random"` (random-effect variances, standard deviations
#'     and correlations, with intervals on the log and Fisher-z scales),
#'     `"random_effects"` (one row per person), `"trajectories"` (each class's
#'     mean trajectory with a band) and `"individual"` (each person's
#'     predicted curve); the last two take `time`. They have no `"groups"`,
#'     `"fitted"` or `"classification"` table. A multilevel growth mixture
#'     adds `"group_classes"` and `"clusters"`, and its `"fit"` uses the number
#'     of clusters in `bic` (persons in `bic_persons`).}
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
#'     `"transitions"` (model-implied, averaged over the observed groups, with
#'     their sampling weights in a weighted fit; `"initial"` and
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
#'     or `"bic_individual"`). `"best"` reports the lowest value of a
#'     criterion; it does not select a model. With one group class every
#'     family is the same model, so those rows share a likelihood. Candidates
#'     are extracted by name with [get_candidate()]; transition grids can also
#'     be passed to [compare_models()].}
#' }
#' @section What `"all"` returns: A named list of tables included by default
#'   for the object, in catalogue order. Some results require an explicit
#'   request, including the `"lmr"` table for an enumeration object.
#'   Each included table is built with its default arguments;
#'   supplying any further argument raises `latents_bad_argument`. A table
#'   this object cannot produce (for example a sequence table for a fit made
#'   without `time`, or bivariate residuals for a family without discrete
#'   group classes) is omitted, so the names of the list show what was
#'   available. Only the conditions that mean the table does not exist for
#'   this object are skipped; any other error is raised.
#'
#' @section Table arguments:
#'   These are passed through `...`; each table takes only those listed for it
#'   under *Tables*.
#' \describe{
#'   \item{`data`}{A data frame with one row per observation of the fit, in
#'   the order it was fitted in. `NULL` uses the columns the fit carries. The
#'   standard errors of `"profiles"` and `"responses"` come from
#'   [parameter_inference()] on these data. Without `data`, a fit that cannot
#'   supply them (bound-active, unconverged, singular information, a noise
#'   component) gives `NA` with a `latents_no_standard_errors` message; with
#'   `data`, the underlying condition is raised.
#'   }
#'   \item{`scale`}{For `"profiles"`: `"raw"` (default) reports estimates in
#'   input units; `"standardized"` divides each indicator's deviation from its
#'   grand mean by that indicator's observed standard deviation, as
#'   `plot(x, scale = "standardized")` draws. The standard errors are divided
#'   by the same constant.
#'   }
#'   \item{`format`}{For the posterior tables and `"sequences"`, `"long"`
#'   (default) gives one row per unit and class and `"wide"` one row per unit
#'   with one column per class. For `"information_criteria"`, `"wide"`
#'   (default) gives one row and `"long"` one row per criterion and
#'   sample-size convention.
#'   }
#'   \item{`level`}{The classification level to report: `"individuals"`,
#'   `"groups"` or `"both"`. `NULL` (default) gives `"both"` for a fit with
#'   discrete group classes and `"individuals"` otherwise, as in
#'   [diagnostics()].
#'   }
#'   \item{`truth`}{For `"assignments"`, one or more column names of `data`
#'   holding known labels to compare the assignments with. See *Recovery*.
#'   }
#'   \item{`by`}{For `"residuals"`, `"profile"` (default) assesses each
#'   profile separately and `"overall"` pools them.
#'   }
#'   \item{`adjust`}{For `"residuals"`, the multiplicity correction applied
#'   across the indicator pairs: `"none"` (default) or one of the
#'   [stats::p.adjust()] methods.
#'   }
#'   \item{`definitions`}{For `"information_criteria"` in long format, `TRUE`
#'   adds each criterion's formula and reference. Default `FALSE`.
#'   }
#'   \item{`estimated`, `stable`}{For `"transitions"`, `TRUE` or `FALSE`
#'   restricts the table to the rows estimated (or not) from observed
#'   occupancy, or to the diagonal (or off-diagonal) cells. `NULL` (default)
#'   keeps every row.
#'   }
#' }
#'
#' @section Recovery: `what = "assignments"` with `truth` returns a
#'   cross-tabulation of the model's labels against known ones instead of the
#'   per-observation table: one row per class and truth value, with columns
#'   `assignment` (the model label compared against), `class`, `truth` (the
#'   column name), `value`, `n` and `proportion` (the share of units with that
#'   truth value assigned to that class). Several truth columns are stacked.
#'   Excluded observations or groups appear under `class = NA`, so the table
#'   distinguishes unclassified units from assignments to an estimated class.
#'
#'   A truth column that is constant within every group is compared against
#'   `group_class`; one that varies within any group is compared against
#'   `profile`. The `assignment` column records which comparison was made.
#'   Group-level recovery counts each group once. Unused factor levels are
#'   omitted. For a fit without discrete group classes every column is
#'   compared against `profile`.
#'
#'   Class labels are arbitrary, so a truth value assigned mostly to a class
#'   with a different label indicates relabelling, not misclassification.
#'
#' @section Information criteria: Let `q` be the number of free parameters,
#'   `n` the chosen sample size, and `EN` the entropy of the complete-data
#'   classification (defined below). The criteria are `deviance = -2L`,
#'   `aic = -2L + 2q`, `bic = -2L + q log(n)`,
#'   `sabic = -2L + q log((n + 2) / 24)`, `caic = -2L + q (log(n) + 1)`,
#'   `awe = -2(L - EN) + 2q (1.5 + log(n))`, `icl = -2L + q log(n) + 2 EN`,
#'   `kic = -2L + 3(q + 1)`, and `clc = -2L + 2 EN`. `clc`, like `icl` and
#'   `awe`, uses the entropy sum; `tidyLPA` reports a `CLC` based on relative
#'   entropy, and the two values are not comparable. The `kic` form is
#'   Cavanaugh's (1999) criterion as `tidyLPA` computes it (Rosenberg et al.
#'   2018), with `q` the number of free parameters; Cavanaugh's general form
#'   `-2L + 3q` is 3 smaller for every model, so differences between models
#'   and the model with the smallest value are the same under both.
#'
#'   `EN` is the soft entropy of the complete-data classification, the
#'   entropy term of the classification log likelihood `log L_c = log L - EN`
#'   of ICL (Biernacki, Celeux and Govaert 2000), AWE (Banfield and Raftery
#'   1993) and CLC (Biernacki and Govaert 1997). The latent classification of
#'   a two-level model is the group class of every group and the profile of
#'   every row, so `EN = sum_j H(C_j | y_j) + sum_j sum_h P(C_j = h | y_j)
#'   sum_i H(X_ij | C_j = h, y_ij)`: the group-class entropy plus each row's
#'   profile entropy within its group class. With one group class, and in
#'   every single-level model, it is the entropy of the profile posteriors.
#'   It is the soft form, the posterior expectation of `-log P(classes | y)`
#'   (McLachlan and Peel 2000, section 6.10.3), not the entropy of the modal
#'   (hard) assignment that Biernacki, Celeux and Govaert (2000) also discuss.
#'   A fit with several group classes saved by latents 0.97.1 or earlier
#'   lacks the rows' within-class entropies, and its `awe`, `icl` and `clc`
#'   are `NA`; refit it to obtain them.
#'
#'   Every criterion that depends on a sample size is reported under both
#'   multilevel conventions: `"groups"` with `n` the number of groups and
#'   `"individuals"` with `n` the number of rows. `EN` is the same under both,
#'   so the conventions differ only in `n`, and in a single-level model they
#'   coincide. Use one convention consistently across compared candidates.
#'   `deviance`, `aic` and `kic` carry `convention = NA_character_`; `clc`,
#'   which uses no sample size, is the same in both rows. A transition model
#'   ([lta()]) keeps one entropy per level instead: its profiles form a Markov
#'   chain within a person, whose path entropy is not the sum of the
#'   occasions' entropies. Its `"groups"` criteria use group-class entropy,
#'   whereas its `"individuals"` criteria use occasion-level profile entropy.
#'   These definitions are reported by
#'   `get_results(x, "information_criteria", format = "long", definitions = TRUE)`.
#'   Individual sample sizes exclude rows with
#'   no observed indicators, and group sample sizes exclude groups in which
#'   every row has none (`n_informative_groups`); the entropy and the
#'   classification statistics exclude the same rows and groups, as in Mplus.
#'   Units of sampling weight zero are in no sample size and no entropy, under
#'   either `drop_empty`. For a fit with sampling weights, `L` is the pseudo
#'   log likelihood and `EN` the sampling-weighted entropy, each unit's terms
#'   multiplied by its weight, so that integer weights give the entropy of the
#'   data with each unit repeated that many times. The `"entropy"` table,
#'   `relative_entropy`, the average posteriors, the odds of correct
#'   classification, the classification errors and the BCH weights are
#'   weighted the same way; modal counts stay counts.
#'
#' @section Classification quality: The odds of correct classification for
#'   class `k` is `(p / (1 - p)) / (r / (1 - r))`, where `p` is the average
#'   posterior in the assigned class and `r` is the model-estimated class
#'   proportion. Values near one indicate classification no better than the
#'   class proportions alone. A value of five is a conventional descriptive
#'   benchmark, rather than a formal test of class separation. It is `NA_real_` when `average_posterior` is
#'   missing or zero, or when `estimated_proportion` is zero or one. A perfect
#'   `average_posterior` of one gives `Inf` at an interior estimated proportion; a
#'   single-class solution always reports `NA_real_`. `average_posterior` is
#'   `NA_real_` for a class with no modal members. `n_modal` and
#'   `estimated_n` use modal assignments and posterior probabilities,
#'   respectively, and need not coincide.
#'
#'   `relative_entropy` is `1 - EN_level / (n log K)` (Ramaswamy et al.
#'   1993), between zero and one, higher meaning sharper classification, with
#'   `EN_level` the level's own entropy sum, reported as `entropy_sum`, lower
#'   meaning sharper. With one group class the individual level's
#'   `entropy_sum` is the `EN` of the information criteria; with several, that
#'   `EN` adds the group-class entropy to the profile entropy within each
#'   group class (see *Information criteria*). `relative_entropy` is
#'   `NA_real_` for a level with a single class.
#'
#'   A fit with individual profiles but no discrete group classes reports
#'   individuals only for `level = "both"` and raises
#'   `latents_no_group_classes` for `level = "groups"`.
#'
#' @section Bivariate residuals: `"residuals"` checks, pair by pair, the
#'   assumption that indicators are independent within a profile, weighting
#'   each observation by its posterior (times its sampling weight). Every pair
#'   of indicators has a row (per profile under `by = "profile"`, once under
#'   `by = "overall"`). `data` must be in the order the model was fitted in;
#'   the alignment is checked. Large residuals suggest associations that the
#'   measurement model does not adequately represent. Calculations for each
#'   `kind` are described below:
#'
#'   * `"gaussian"`, two Gaussian indicators: `observed` is the
#'     posterior-weighted correlation of the residuals from the profile's
#'     means (pooled over the profiles under `"overall"`), `expected` the
#'     correlation the fitted covariance implies (zero for a diagonal one).
#'     `statistic` is the difference of their Fisher z transforms times
#'     `sqrt(effective_n - 3)`, referred to the standard normal (Fisher,
#'     1921), with a two-sided p-value; `df` is `NA`. This is an approximate
#'     reference distribution for a posterior-weighted diagnostic.
#'   * `"categorical"`, two categorical or ordinal indicators (an ordinal
#'     indicator's categories are its levels): `statistic` is the Pearson
#'     X^2 of the posterior-weighted two-way table against the table the
#'     model implies, which is Latent GOLD's bivariate residual times its
#'     `df`, and `p_value` refers it to chi-square with `df = (R - 1)(C - 1)`.
#'     `observed` and `expected` are the two tables' Cramer's V (Cramer,
#'     1946), zero for the model's table within a profile. The chi-square
#'     reference ignores that the model's table is estimated from the same
#'     data and can therefore be conservative (Oberski, van Kollenburg
#'     & Vermunt, 2013). Interpret `p_value` together with `residual` and
#'     `statistic`; a non-significant result does not establish conditional
#'     independence.
#'   * `"mixed"` (a Gaussian and a categorical or ordinal indicator) and
#'     `"count"` (any pair with a count indicator): `observed` is the
#'     posterior-weighted correlation of the two indicators' values within the
#'     profile, centred on the profile's weighted means (pooled over the
#'     profiles under `"overall"`); a categorical or ordinal indicator enters
#'     as its category number in category order, which for a nominal
#'     indicator with more than two categories depends on that order.
#'     `expected` is zero, the correlation of independent indicators, and
#'     `statistic`, `df` and `p_value` are `NA`: no reference distribution is
#'     computed for these pairs.
#'
#'   A half-observed row of a joined `local_dependence` pair enters that
#'   pair's tables only where both items are observed. `p_adjusted` applies
#'   the correction named by `adjust` across the pairs that have a p-value.
#'   A fit with sampling weights uses them in the descriptive associations
#'   and `effective_n`; design-based inference is not implemented, so
#'   `statistic`, `df`, `p_value` and `p_adjusted` are `NA` for weighted fits.
#'
#' @references Fisher, R. A. (1921). On the "probable error" of a coefficient
#'   of correlation deduced from a small sample. *Metron*, 1, 3--32.
#'
#'   Cramer, H. (1946). *Mathematical Methods of Statistics*. Princeton
#'   University Press.
#'
#'   Schwarz, G. (1978). Estimating the dimension of a model.
#'   *Annals of Statistics*, 6, 461--464.
#'
#'   Sclove, S. L. (1987). Application of model-selection criteria to some
#'   problems in multivariate analysis. *Psychometrika*, 52, 333--343.
#'
#'   Bozdogan, H. (1987). Model selection and Akaike's information criterion.
#'   *Psychometrika*, 52, 345--370.
#'
#'   Banfield, J. D., & Raftery, A. E. (1993). Model-based Gaussian and
#'   non-Gaussian clustering. *Biometrics*, 49, 803--821.
#'
#'   Biernacki, C., & Govaert, G. (1997). Using the classification likelihood
#'   to choose the number of clusters. *Computing Science and Statistics*, 29,
#'   451--457.
#'
#'   Cavanaugh, J. E. (1999). A large-sample model selection criterion based on
#'   Kullback's symmetric divergence. *Statistics and Probability Letters*, 42,
#'   333--343.
#'
#'   Rosenberg, J. M., Beymer, P. N., Anderson, D. J., Van Lissa, C. J., &
#'   Schmidt, J. A. (2018). tidyLPA: An R package to easily carry out latent
#'   profile analysis (LPA) using open-source or commercial software.
#'   *Journal of Open Source Software*, 3(30), 978.
#'
#'   McLachlan, G., & Peel, D. (2000). *Finite Mixture Models*. Wiley.
#'
#'   Ramaswamy, V., DeSarbo, W. S., Reibstein, D. J., & Robinson, W. T.
#'   (1993). An empirical pooling approach for estimating marketing mix
#'   elasticities with PIMS data. *Marketing Science*, 12, 103--124.
#'
#'   Biernacki, C., Celeux, G., & Govaert, G. (2000). Assessing a mixture model
#'   for clustering with the integrated completed likelihood. *IEEE
#'   Transactions on Pattern Analysis and Machine Intelligence*, 22, 719--725.
#'
#'   Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent structure
#'   models with categorical variables. *Political Analysis*, 12, 3--27.
#'
#'   Nagin, D. S. (2005). *Group-Based Modeling of Development*. Harvard
#'   University Press.
#'
#'   Vermunt, J. K. (2010). Latent class modeling with covariates: two improved
#'   three-step approaches. *Political Analysis*, 18, 450--469.
#'
#'   Satopaa, V., Albrecht, J., Irwin, D., & Raghavan, B. (2011). Finding a
#'   "kneedle" in a haystack: Detecting knee points in system behavior.
#'   *Proceedings of the 31st International Conference on Distributed
#'   Computing Systems Workshops*, 166--171.
#'
#'   Oberski, D. L., van Kollenburg, G. H., & Vermunt, J. K. (2013). A Monte
#'   Carlo evaluation of three methods to detect local dependence in binary
#'   data latent class models. *Advances in Data Analysis and
#'   Classification*, 7, 267--279.
#' @return A base `data.frame` for a single `what`, whose columns are given
#'   under *Tables*, or a named list of such data frames for `what = "all"`.
#' @seealso [as.data.frame()], which returns the primary table alone,
#'   [summary()], which prints every table, and [diagnostics()] for the
#'   classification tables with a diagnostic summary.
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
#' # Compare estimated assignments with the generating classes at both levels.
#' get_results(fit, what = "assignments", data = course_engagement,
#'          truth = c("engagement", "student_type"))
#' @export
get_results <- function(x, what = NULL, ...) {
  UseMethod("get_results")
}

#' @noRd
#' @export
get_results.default <-function(x, what = NULL, ...) {
  stop(errorCondition(sprintf(
    "`get_results()` has no method for an object of class %s.",
    paste(sprintf("`%s`", class(x)), collapse = ", ")),
    class = "latents_bad_argument", call = NULL))
}

#' @noRd
#' @export
get_results.multilpa <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @noRd
#' @export
get_results.multilpa_covariates <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @noRd
#' @export
get_results.multilpa_transitions <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}


#' @noRd
#' @export
get_results.multilpa_enumeration <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @noRd
#' @export
get_results.latents_comparison <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @noRd
#' @export
get_results.multilpa_bootstrap_lrt <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @noRd
#' @export
get_results.multilpa_start <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @noRd
#' @export
get_results.multilpa_diagnostics <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @noRd
#' @export
get_results.summary_multilpa <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @noRd
#' @export
get_results.summary_multilpa_covariates <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}


#' @noRd
#' @export
get_results.summary_multilpa_transitions <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @noRd
#' @export
get_results.summary_multilpa_enumeration <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @noRd
#' @export
get_results.summary_multilpa_bootstrap_lrt <-function(x, what = NULL, ...) {
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
  .latents_require("`what` must be a single table name, or NULL for the first one" =
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
    "latents_unsupported_noise",
    # A grid with one solution per model has no tree.
    "latents_no_tree",
    # A fit without locally dependent pairs has no table of them.
    "latents_no_local_dependence",
    # Nor does a fit without direct effects.
    "latents_no_direct_effects")
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
  catalogue <- Filter(function(entry) !isFALSE(entry$in_all), catalogue)
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
      criteria = .multilpa_table(
        function(x) .multilpa_enumeration_minima(.multilpa_enumeration_grid(x))),
      classification = .multilpa_table(.multilpa_enumeration_classification),
      deltas = .multilpa_table(
        function(x, using = "bic_groups") .latents_deltas_table(x, using), "using"),
      knees = .multilpa_table(
        function(x, using = "bic_groups") .latents_knees_table(x, using), "using"),
      tree = .multilpa_table(.latents_tree_table, .latents_tree_arguments()),
      # Every candidate's information matrix: computed when asked for, not by
      # get_results(x, "all") or summary().
      lmr = .multilpa_table(.latents_lmr_table, in_all = FALSE)))
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
    # The entropy table stays the default (as in 0.9.21); the checks follow.
    tables <- list(
      entropy = .multilpa_table(function(x) x$entropy),
      classification = .multilpa_table(function(x) x$classification),
      average_posteriors = .multilpa_table(function(x) x$average_posteriors),
      residuals = .multilpa_table(.multilpa_diagnostic_residuals),
      # A diagnostics object saved before 0.97.5 has no stored table.
      replication = .multilpa_table(function(x) {
        x$replication %||% .multilpa_replication_table(x$fit)
      }),
      checks = .multilpa_table(function(x) x$checks))
    # The family tables a fit did not produce (a family without the profile
    # classification tables, or an object saved before the checks) are absent.
    present <- c(checks = !is.null(x$checks), entropy = !is.null(x$entropy),
                 classification = !is.null(x$classification),
                 average_posteriors = !is.null(x$average_posteriors),
                 residuals = !is.null(x$fit) && .multilpa_any_fit(x$fit),
                 replication = TRUE)
    return(tables[present[names(tables)]])
  }
  if (.multilpa_any_summary(x)) return(.multilpa_stored_catalogue(x))
  stopifnot("`x` must be an object of this package" = .multilpa_any_fit(x))
  c(.multilpa_measurement_catalogue(),
    if (!is.null(x$gaussian_diagnostics) && nrow(x$gaussian_diagnostics) > 0L)
      list(gaussian_diagnostics = .multilpa_table(function(x) x$gaussian_diagnostics)),
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
#' @param in_all `FALSE` leaves the table out of `what = "all"` (and so of
#'   `summary()`), for a table too costly to build unasked.
#' @return A list of the three.
#' @noRd
.multilpa_table <- function(build, args = character(), in_all = TRUE) {
  list(build = build, args = args, in_all = in_all)
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
        .multilpa_profile_frame(x, data = data,
                                scale = .latents_match_arg(scale, c("raw", "standardized"),
                                                           "scale"),
                                standard_errors = TRUE),
      c("data", "scale")),
    responses = .multilpa_table(
      function(x, data = NULL) .multilpa_response_frame(x, data = data,
                                                        standard_errors = TRUE),
      "data"),
    covariances = .multilpa_table(.multilpa_covariance_frame),
    local_dependence = .multilpa_table(.latents_local_dependence_table),
    direct_effects = .multilpa_table(
      function(x, data = NULL) .latents_direct_frame(x, data = data), "data"))
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
      .multilpa_posterior_frame(x, .latents_match_arg(format, c("long", "wide"), "format")),
    "format"))
  assignments <- list(assignments = .multilpa_table(
    function(x, data = NULL, truth = NULL)
      .multilpa_hide_fabricated_id(
        x, .multilpa_assignments(x, data = data, truth = truth)),
    c("data", "truth")))
  # A single-level fit, or one with one group class, classifies nothing at
  # the group level: its group posteriors stay available by name but are left
  # out of "all" and summary() (API-9).
  c(posteriors,
    list(group_posteriors = .multilpa_table(
      function(x, format = c("long", "wide"))
        .multilpa_group_posterior_frame(
          x, .latents_match_arg(format, c("long", "wide"), "format")),
      "format", in_all = .multilpa_has_group_classes(x))),
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
    replication = .multilpa_table(.multilpa_replication_table),
    residuals = .multilpa_table(
      function(x, data = NULL, by = c("profile", "overall"),
               adjust = .multilpa_p_adjust_methods)
        .multilpa_bivariate_residuals(
          x, data = data, by = .latents_match_arg(by, c("profile", "overall"), "by"),
          adjust = .latents_match_arg(adjust, .multilpa_p_adjust_methods, "adjust")),
      c("data", "by", "adjust")),
    information_criteria = .multilpa_table(
      function(x, definitions = FALSE, format = c("wide", "long"))
        .multilpa_information_criteria(
          x, definitions = definitions,
          format = .latents_match_arg(format, c("wide", "long"), "format")),
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
        .multilpa_sequences(x, format = .latents_match_arg(format, c("long", "wide"),
                                                           "format")),
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
  # The rule `diagnostics()` applies: a single-level fit, or one with one
  # group class, has no group-level classification to report (API-9).
  if (.multilpa_has_group_classes(x)) "both" else "individuals"
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
  group_classes <- x$group_classes
  group_classes[!.multilpa_informative_groups(x)] <- NA_integer_
  data.frame(group = x$group_values, group_class = group_classes,
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

#' Response probabilities one row per indicator category, one column per class
#' @param x A fitted model with categorical indicators.
#' @return A data frame: `indicator`, `category`, then `profile_1`, ...
#' @noRd
.multilpa_wide_responses <- function(x) {
  long <- .multilpa_response_frame(x)
  key <- paste(long$indicator, long$category, sep = "\r")
  keys <- unique(key)
  first <- match(keys, key)
  profiles <- sort(unique(long$profile))
  columns <- lapply(profiles, function(profile) {
    mine <- long$profile == profile
    long$probability[mine][match(keys, key[mine])]
  })
  names(columns) <- paste0("profile_", profiles)
  data.frame(indicator = long$indicator[first], category = long$category[first],
             columns, check.names = FALSE, stringsAsFactors = FALSE)
}

#' Print the estimated class sizes, and the group-class sizes of a two-level fit
#' @param x A fitted model.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_print_class_sizes <- function(x) {
  sizes <- .multilpa_count_frame(x)
  line <- function(level, label, prefix) {
    rows <- sizes[sizes$level == level, , drop = FALSE]
    if (nrow(rows) == 0L) return(invisible(NULL))
    cat(sprintf("%-14s%s\n", label, paste(sprintf("%s%s %.1f%% (%.0f)", prefix, rows$class,
                                                   100 * rows$effective_proportion,
                                                   rows$effective_count),
                                           collapse = "   ")))
  }
  cat("\n")
  line("individuals", "Class sizes", "profile_")
  if (!isTRUE(x$single_level) && !identical(as.integer(x$n_group_classes), 1L)) {
    line("groups", "Group classes", "group_class_")
  }
  invisible(NULL)
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
    ## How big each profile is comes first: a mean is read differently when it
    ## describes two percent of the sample. The means follow one row per
    ## indicator and one column per profile, as a class model's probabilities
    ## do, so many indicators lengthen the table instead of wrapping it.
    .multilpa_print_class_sizes(x)
    # Read by position from the wide table, never by name: an indicator may
    # itself be called "profile".
    wide <- .multilpa_wide_means(x)
    indicators <- names(wide)[-1L]
    means <- data.frame(indicator = indicators, stringsAsFactors = FALSE)
    invisible(lapply(seq_len(nrow(wide)), function(index) {
      means[[paste0("profile_", wide[[1L]][[index]])]] <<-
        round(vapply(wide[index, -1L], as.numeric, numeric(1)), 3L)
    }))
    cat("\nMeans\n")
    if (rows > 0) print(utils::head(means, rows), row.names = FALSE)
    if (nrow(means) > rows) {
      cat(sprintf("   ... %d more indicators.  get_results(x, \"profiles\")\n",
                  nrow(means) - rows))
    }
    cat("\nVariances and standard errors: get_results(x, \"profiles\").\n",
        "Every other table: get_results(x, what = ), or get_results(x, \"all\").\n", sep = "")
    return(invisible(NULL))
  }
  if (nrow(table) == 0L && length(x$categorical %||% character()) > 0L) {
    # A class model reads as a profile of response probabilities: one row per
    # indicator category, one column per class, under the class sizes.
    .multilpa_print_class_sizes(x)
    wide <- .multilpa_wide_responses(x)
    numeric_columns <- vapply(wide, is.numeric, logical(1))
    wide[numeric_columns] <- lapply(wide[numeric_columns], round, digits = 3L)
    cat("\nResponse probabilities\n")
    if (rows > 0) print(utils::head(wide, rows), row.names = FALSE)
    if (nrow(wide) > rows) {
      cat(sprintf("   ... %d more rows.  get_results(x, \"responses\")\n", nrow(wide) - rows))
    }
    cat("\nThresholds and standard errors: get_results(x, \"responses\").\n",
        "Every other table: get_results(x, what = ), or get_results(x, \"all\").\n", sep = "")
    return(invisible(NULL))
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
