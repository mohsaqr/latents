#' Plot latents results
#'
#' `plot()` draws any result of this package and returns a ggplot object,
#' which can be printed, saved with `ggplot2::ggsave()` or restyled with
#' `+ ggplot2::theme()`. The plots of mixture regressions (for fits that
#' classify persons, the `"fitted"` view only), pooled imputations and
#' bootstrap tests are drawn when called, including inside loops, and their
#' ggplot is returned invisibly.
#'
#' The main views also have their own functions, each grouping related views
#' under `what`, with `"all"` drawing every one: `plot_profiles()` (profile
#' means as lines, bars, a heatmap, rainclouds, parallel lines or pairs, or
#' categorical means, medians and modes), `plot_responses()` (categorical
#' response probabilities), `plot_sequences()` (each group's profile at each
#' occasion), `plot_prevalence()` (profile shares within each group class),
#' `plot_diagnostics()` (classification quality), `plot_transitions()`
#' (transition matrix), `plot_enumeration()` (information criteria,
#' classification quality and model-selection evidence across a grid),
#' `plot_selection()` (the model-selection evidence of a grid) and
#' `plot_tree()` (how profiles split as more are added, one plot per tree, or
#' how chosen solutions correspond). Each draws the same
#' plot as the matching `plot(x, what = )`. `plot_views()` lists the views of
#' a fitted model, and `"enumeration"` and `"tree"` for a grid, with a
#' one-line description of each; the other grid views are listed under
#' `what` below.
#'
#' Each profile keeps the same order (largest first), Okabe-Ito colour and
#' point shape in every view. The plots require ggplot2, which is a suggested
#' package.
#'
#' @param x The object to plot: a profile, class or transition model for the
#'   named fit views and `plot_diagnostics()`; an enumeration from
#'   [enumerate_classes()], [enumerate_lpa()] or [enumerate_lca()] for
#'   `plot_enumeration()`, `plot_selection()` and `plot_tree()`; a transition
#'   fit from [lta()] for `plot_transitions()`. A general [lta()] fit
#'   (occasion-varying transitions, covariates, mover-stayer, second order,
#'   weights, ordinal or count indicators) draws only its transition matrix,
#'   with `plot_transitions()` or `plot()`. Transition and group-class family
#'   grids plot with `plot(x, what = )`.
#' @param what The view to draw. `"all"` draws every view of that function
#'   that the model supports, as a list, and names in a message any it cannot
#'   draw. For `plot_profiles()`: `"lines"` (default), `"bars"` (means from
#'   zero with 95% intervals), `"heatmap"` (means in observed standard
#'   deviations, or the response probabilities of categorical indicators),
#'   `"raincloud"` (each indicator's distribution by assigned profile),
#'   `"parallel"` (every case as a line) or `"pairs"` (scatter-plot matrix
#'   with each profile's 95% covariance ellipse).
#'
#'   For `plot_diagnostics()`: `"all"` (default), `"sizes"` (effective profile
#'   counts), `"posteriors"` (posterior probability of each case's assigned
#'   profile), `"entropy"` (each case's entropy relative to a flat posterior)
#'   or `"avepp"` (average posterior probability matrix). `"group_sizes"`,
#'   `"group_posteriors"`, `"group_entropy"` and `"group_avepp"` draw the same
#'   four views for the groups and their group classes; `"all"` includes them
#'   when the fit has two or more group classes, and with fewer they are
#'   refused with `latents_no_group_classes`. A fit with a noise component
#'   supports `"sizes"` and the group views only: `"all"` returns `"sizes"`,
#'   and the other individual-level views are refused with
#'   `latents_unsupported_noise`.
#'
#'   For `plot_selection()`: `"all"` (default), `"deltas"`, `"heatmap"`,
#'   `"smallest"`, `"trust"` or `"diagnostics"`, as described for
#'   `plot_enumeration()`.
#'
#'   For `plot_enumeration()`, three kinds of view; information criteria
#'   cannot be requested together with the other two:
#'   * Information criteria: criterion columns of `as.data.frame(x)`, such as
#'     `"bic_groups"`, or `"all"`. The default draws AIC, both BICs and ICL
#'     over individuals.
#'   * Classification quality across the candidates, at each level with two
#'     or more classes: `"entropy"` (relative entropy), `"avepp"` (the lowest
#'     AvePP, with the cutoff dashed) or `"diagnostics"` (both).
#'   * Model-selection evidence (none of which selects a model): `"deltas"`
#'     (the decrease in the criterion in `using` with each added profile and
#'     each added group class; see `get_results(x, "deltas")` and
#'     `"knees"`), `"heatmap"` (the criterion over profiles and group classes,
#'     relative to the lowest), `"smallest"` (the smallest class share, with
#'     5% dashed) or `"trust"` (convergence, boundary and replication of every
#'     candidate). `"selection"` draws all four and `"diagnostics"`, as a
#'     list.
#' @param statistic `"mean"` (default), `"median"` or `"mode"`. They coincide
#'   for continuous indicators; for categorical ones they summarise the fitted
#'   probabilities on the numeric category scale, without intervals.
#' @param scale `"raw"` (default) keeps each indicator in its own units;
#'   `"standardized"` divides its deviation from the observed mean by its
#'   observed standard deviation.
#' @param intervals For profile, bar and response views, `TRUE` (default)
#'   draws 95% intervals when the fit has standard errors; otherwise the
#'   subtitle states why there are none. They are the `conf_low` and
#'   `conf_high` of [parameter_inference()]; a response probability's interval
#'   is formed on the logit scale.
#' @param labels `TRUE` (default) labels each series directly; `FALSE` uses a
#'   legend. On sequence views, `TRUE` groups the rows by group class under
#'   labelled strips.
#' @param data The data the model was fitted to, used only to compute the
#'   intervals. `NULL` (default) uses the data stored in the fit.
#' @param category Which category's probability to draw: `"last"` (default),
#'   `"first"`, or a single label or index.
#' @param combine `TRUE` (default) draws several criteria as panels of one
#'   plot; `FALSE` returns one plot per criterion. For `plot_tree()` (and
#'   `what = "tree"`), `FALSE` (default) returns one plot per tree and `TRUE`
#'   draws every tree as a panel of one plot.
#' @param n_group_classes,n_profiles,model For `plot_tree()`, the trees and
#'   rows to draw. At the individual level a tree runs down the numbers of
#'   profiles (classes) at one number of group classes and one structure:
#'   `n_group_classes` and `model` choose trees, `n_profiles` the rows. At the
#'   group level a tree runs down the numbers of group classes at one number
#'   of profiles: `n_profiles` and `model` choose trees, `n_group_classes` the
#'   rows. `NULL` (default) keeps every converged candidate. A single count
#'   for the rows (`n_profiles` at the individual level, `n_group_classes` at
#'   the group level) is one level of every tree, and the plot compares that
#'   level's solutions instead: `plot_tree(x, n_profiles = 4,
#'   n_group_classes = c(2, 5))` draws the four-profile solution at two and at
#'   five group classes as two rows, each class in the colour of the class it
#'   matches in the first row.
#' @param candidates For `plot_tree()`, two or more candidate names drawn as
#'   the rows of one tree, in the order given, to compare those solutions.
#'   Not combined with `n_group_classes`, `n_profiles` or `model`.
#' @param from,to For `plot_tree()`, two solutions to compare as a two-row
#'   tree, each named by its class counts `c(n_profiles, n_group_classes)`
#'   (one count means one group class) or by its name:
#'   `plot_tree(x, from = c(3, 3), to = c(4, 3))`. Given together, and not
#'   with `candidates`. On a grid crossed with covariance structures, `model`
#'   names the structure of both.
#' @param mark_minimum `TRUE` (default) circles each criterion's lowest value
#'   among converged candidates, without selecting a model. Information
#'   criteria only.
#' @param cutoff For `plot_enumeration(what = "avepp")`, `"diagnostics"` or
#'   `"selection"`, the AvePP cutoff drawn as a dashed line (default 0.8).
#' @param using For `plot_enumeration(what = "deltas")`, `"heatmap"` or
#'   `"selection"`, the information criterion used (default
#'   `"bic_groups"`); `"smallest"` and `"trust"` use none.
#' @param main,subtitle Title and subtitle. `NULL` (default) uses the view's
#'   own.
#' @param cell_labels For sequence views, `TRUE` (default) prints each
#'   assigned profile in its cell when the cells are large enough for the
#'   text.
#' @param predictor For the fitted view of a mixture regression, a numeric
#'   predictor in the formula. `NULL` (default) uses the first numeric
#'   predictor.
#' @param time For growth and trajectory views, the numeric time column.
#'   `NULL` (default) uses the first numeric variable in the fitted random or
#'   regression design.
#' @param facet For trajectories, `FALSE` (default) draws all classes in one
#'   panel; `TRUE` draws one panel per class.
#' @param persons For individual trajectories, the number of persons to show
#'   (default 12), or a vector of their identifiers.
#' @param max_persons The maximum number of observed trajectories drawn
#'   behind the class means (default 120); zero hides them.
#' @param spread For growth trajectories, the band of random-effect spread of
#'   persons around the class curve, from 0 (hidden) to less than 1; default
#'   0.8.
#' @param level For `plot_tree()`, `"individuals"` (default; the cases
#'   divided into profiles or classes) or `"groups"` (the groups divided into
#'   group classes). Otherwise the confidence level for additive,
#'   mixture-regression and growth
#'   views (default 0.95). On pooled plots, the former name of `what`:
#'   `"profile"`, `"group"` or `"measurement"`.
#' @param criterion The former name of `what` for enumeration criteria. On a
#'   profile grid it is used alone or with `what = "enumeration"`; transition
#'   and group-class family grids take `"bic"`, `"aic"` or
#'   `"bic_individual"`.
#' @param palette For bootstrap plots, a colour vector for the histogram and
#'   the observed statistic; `NULL` (default) uses the standard colours.
#' @param style The former bootstrap drawing style; ignored with a
#'   `latents_unused_argument` warning.
#' @param ... With `what = "tree"` (on `plot()` of an enumeration or on
#'   `plot_enumeration()`), the arguments of `plot_tree()`. Otherwise not
#'   used. Profile, class and transition models and enumerations refuse
#'   unknown arguments with `latents_bad_argument`; the other methods ignore
#'   them with a `latents_unused_argument` warning (see *Errors*).
#' @return A ggplot object. `what = "all"` (on `plot()` of a fitted model or
#'   on a named plot function), `plot_enumeration(combine = FALSE)` with
#'   several criteria or classification-quality measures, and enumeration
#'   views when more than one is requested (including `"selection"`) return
#'   a `latents_plots` list of ggplot objects, which draws each one when
#'   printed, pausing between plots on an interactive screen device.
#'   `plot_tree()` and `what = "tree"` return one ggplot for a single tree or
#'   `combine = TRUE`, and otherwise a `latents_plots` list named by tree. `plot_views()` returns a `data.frame` with one row per view:
#'   `type` (the value for `what`), `group` and `description`.
#'
#' @section What `plot()` draws:
#' Each object takes `what =` for its views, `main`, `subtitle` and the
#' arguments named here.
#' \describe{
#'   \item{Fitted profile and class models}{From [multilpa()], [multilca()],
#'     [lpa()], [lca()] and covariate models: `"profiles"`
#'     (default), `"bars"`, `"heatmap"`, `"raincloud"`, `"parallel"`,
#'     `"pairs"` (each profile's 95% covariance ellipse), `"responses"`,
#'     `"probabilities"` (not for covariate models), `"sequences"` (fits with
#'     `time`; `cell_labels = TRUE` prints the profile in each cell),
#'     `"sizes"`, `"entropy"`, `"posteriors"`, `"avepp"` and `"all"`. The
#'     profile views name in a caption any ordinal, count or (for the means)
#'     categorical indicator they omit.}
#'   \item{Transition fits}{From [lta()]: `"transitions"` (default), then
#'     `"profiles"`, `"bars"`, `"heatmap"`, `"responses"`, `"sequences"`,
#'     `"sizes"`, `"entropy"`, `"posteriors"`, `"avepp"` and `"all"`, without
#'     `statistic`; the profile, bar and response views draw the 95%
#'     intervals of the transition fit's [parameter_inference()] (`intervals`).
#'     A general [lta()] fit takes no `what` and draws its transition matrix.}
#'   \item{Enumerations}{Criterion columns or `"all"` (with `combine`,
#'     `labels` and `mark_minimum`), `"tree"` (with the arguments of
#'     `plot_tree()`), the classification-quality
#'     views `"entropy"`, `"avepp"` and `"diagnostics"` (with `cutoff`), and
#'     the model-selection views `"deltas"` and `"heatmap"` (which use the
#'     criterion in `using`), `"smallest"` and `"trust"`; `"selection"` draws
#'     those four and `"diagnostics"`. Transition and group-class family
#'     grids take `"bic"` (default), `"aic"` or `"bic_individual"`.}
#'   \item{Additive fits}{`"means"` (default; each group class's mean with
#'     its Wald interval), `"variances"` (within- and between-group
#'     variances) or `"intercepts"` (each group's posterior intercept);
#'     `level` sets the interval.}
#'   \item{Cross-level fits}{`"profiles"` (default), `"group_means"` or
#'     `"composition"` (profile shares within each group class).}
#'   \item{Mixture regressions}{`"fitted"` (default: the outcome against a
#'     numeric `predictor`, one line per class), `"coefficients"` or
#'     `"posteriors"`; `level` sets the intervals. Fits that classify persons
#'     (`class_level = "group"`) default to `"trajectories"` (unless
#'     `predictor` is given) and also take `"individuals"` and
#'     `"classification"`, as growth mixtures do.}
#'   \item{Growth mixtures}{`"trajectories"` (default: class mean
#'     trajectories with bands; `facet`, `max_persons`, `spread`),
#'     `"individuals"` (`persons`: a number or identifiers),
#'     `"coefficients"`, `"random"` (predicted random effects with each
#'     class's 95% ellipse) or `"classification"` (also `"posteriors"`);
#'     `time` names the time variable, `level` the band.}
#'   \item{Pooled imputations}{`"profile"` (default), `"group"` or
#'     `"measurement"` parameters: each imputation's estimate beside the
#'     pooled one with its interval.}
#'   \item{Model comparisons, bootstrap tests, diagnostics}{One view each:
#'     the criteria of each compared model, the simulated null distribution
#'     with the observed statistic, and the views of `plot_diagnostics()`.}
#' }
#' @section Former argument names: Enumeration plots accept `criterion`
#'   (now `what`; on a profile grid `what = "enumeration"` is the default
#'   criteria view), pooled plots accept `level` (now `what`), and bootstrap
#'   plots accept `palette` (the histogram and observed-line colours). On an
#'   enumeration, `criterion` with any `what` other than `"enumeration"` is an
#'   error of class `latents_bad_argument` (with `"tree"`, a warning, as the
#'   tree uses no criterion); `level` together with `what` on a pooled plot is
#'   an error of the same class.
#' @section Errors: `latents_missing_package` without ggplot2;
#'   `latents_bad_argument` for an argument the object's plot does not take,
#'   except that the plots of additive, cross-level, general transition
#'   (`lta()`), growth-mixture, mixture-regression, pooled, comparison and
#'   bootstrap results and of transition or group-class family grids ignore
#'   it with a `latents_unused_argument` warning (including a bootstrap
#'   plot's `style`); `latents_bad_argument` also when a named plot function
#'   is given an object it does not draw; `latents_nothing_to_plot` when the
#'   fit lacks what a view needs (for example posteriors or entropy of a
#'   one-profile fit, prevalence of a covariate model, or the case-level
#'   views of a centred fit); `latents_unknown_category` when `category`
#'   names no category of an indicator; `latents_no_continuous`,
#'   `latents_no_categorical` and `latents_no_time` when a view needs
#'   indicators or occasions the fit does not have;
#'   `latents_unknown_criterion` for a criterion not in the grid;
#'   `latents_no_tree` (also of class `latents_nothing_to_plot`) when a tree
#'   would have fewer than two rows, or at the group level of a grid with one
#'   group class; `latents_unknown_candidate` when `candidates`,
#'   `n_group_classes`, `n_profiles` or `model` names something no converged
#'   candidate has; `latents_failed_candidate` when a named candidate was not
#'   fitted; `latents_unconverged` when a named candidate did not converge;
#'   `latents_unsupported_noise` for a tree of a fit with a noise component.
#' @references Satopaa, V., Albrecht, J., Irwin, D. and Raghavan, B. (2011). Finding a
#'   "kneedle" in a haystack: detecting knee points in system behavior.
#'   *31st International Conference on Distributed Computing Systems
#'   Workshops*, 166-171.
#'
#'   Zappia, L. and Oshlack, A. (2018). Clustering trees: a
#'   visualization for evaluating clusterings at multiple resolutions.
#'   *GigaScience*, 7(7), giy083.
#' @examples
#' plot_views()
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   fit <- multilpa(course_engagement, c("browse", "lectures", "forum_read"),
#'                   "student", n_profiles = 2, n_group_classes = 2,
#'                   n_starts = 2, seed = 1)
#'   plot_profiles(fit)
#'   plot_diagnostics(fit, what = "avepp")
#' }
#' @name latents-plots
#' @aliases plot.multilpa plot.multilpa_enumeration plot.multilpa_covariates
#'   plot.multilpa_transitions plot.multilpa_lta plot.multilpa_additive
#'   plot.multilpa_cross_level plot.latents_growth_mixture
#'   plot.latents_mixture_regression plot.latents_pooled
#'   plot.latents_comparison plot.multilpa_bootstrap_lrt
#'   plot.latents_transition_enumeration plot.latents_family_enumeration
NULL
