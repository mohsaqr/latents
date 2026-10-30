# Plot model estimates, classification summaries, and comparisons

The plotting functions support interpretation and comparison of fitted
latent variable models. Depending on the model, they display differences
between classes in indicator patterns, outcome relationships, prevalence
or transition probabilities. Classification plots describe how strongly
observations are associated with their assigned classes.

For candidate models, plots compare information criteria and
classification quality across alternative solutions. They also display
convergence diagnostics and correspondence between classes, helping
assess how the solution changes with the number of classes or the model
specification.

Available views depend on the fitted model. Plots are returned as ggplot
objects that can be customised and saved for reporting.

## Usage

``` r
plot_profiles(
  x,
  statistic = c("mean", "median", "mode"),
  scale = c("raw", "standardized"),
  intervals = TRUE,
  labels = TRUE,
  data = NULL,
  main = NULL,
  subtitle = NULL,
  what = c("lines", "bars", "heatmap", "raincloud", "parallel", "pairs", "all")
)

plot_responses(
  x,
  category = "last",
  intervals = TRUE,
  labels = TRUE,
  data = NULL,
  main = NULL,
  subtitle = NULL
)

plot_sequences(
  x,
  labels = TRUE,
  cell_labels = TRUE,
  main = NULL,
  subtitle = NULL
)

plot_prevalence(x, labels = TRUE, main = NULL, subtitle = NULL)

plot_diagnostics(
  x,
  what = c("all", "sizes", "posteriors", "entropy", "avepp", "group_sizes",
    "group_posteriors", "group_entropy", "group_avepp"),
  main = NULL,
  subtitle = NULL
)

plot_transitions(x, main = NULL, subtitle = NULL)

plot_tree(
  x,
  main = NULL,
  subtitle = NULL,
  n_group_classes = NULL,
  n_profiles = NULL,
  model = NULL,
  candidates = NULL,
  from = NULL,
  to = NULL,
  level = c("individuals", "groups"),
  combine = FALSE
)

plot_selection(
  x,
  what = c("all", "deltas", "heatmap", "smallest", "trust", "diagnostics"),
  using = "bic_groups",
  cutoff = 0.8,
  main = NULL,
  subtitle = NULL
)

plot_enumeration(
  x,
  what = c("aic", "bic_groups", "bic_individual", "icl_individual"),
  combine = TRUE,
  labels = TRUE,
  mark_minimum = TRUE,
  main = NULL,
  subtitle = NULL,
  cutoff = 0.8,
  using = "bic_groups",
  ...
)

plot_views(x = NULL)
```

## Arguments

- x:

  The object to plot. For `plot_profiles()`, `plot_responses()`,
  `plot_sequences()`, `plot_prevalence()`, `plot_diagnostics()` and
  `plot_transitions()`, any fitted model of the package: profile, class
  and transition models from
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md),
  [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
  [`lca()`](https://pak.dynasite.org/latents/reference/lca.md) and
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) (with or
  without covariates) take the views described under `what` and *Views
  by object class*; growth mixtures and mixture regressions from
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md),
  additive, dispersion and cross-level fits from
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  pooled imputations from
  [`pool_imputations()`](https://pak.dynasite.org/latents/reference/pool_imputations.md)
  and general
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fits
  (occasion-varying transitions, covariates, mover-stayer, second order,
  weights, ordinal or count indicators) take the views listed under
  *Named plot functions on other fitted models*. For
  `plot_enumeration()`, `plot_selection()` and `plot_tree()`, a grid
  from
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md),
  [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md),
  [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md)
  or
  [`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md);
  each refuses the views that grid does not offer (see above).
  Optionally a grid for `plot_views()`.

- statistic:

  `"mean"` (default), `"median"` or `"mode"`. They coincide for
  continuous indicators; for categorical ones they summarise the fitted
  probabilities on the numeric category scale, without intervals.

- scale:

  `"raw"` (default) keeps each indicator in its own units;
  `"standardized"` divides its deviation from the observed mean by its
  observed standard deviation.

- intervals:

  For profile, bar and response views, `TRUE` (default) draws 95%
  intervals when the fit has standard errors; otherwise the subtitle
  states why there are none. They are the `conf_low` and `conf_high` of
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md);
  a response probability's interval is formed on the logit scale.

- labels:

  `TRUE` (default) labels each series directly; `FALSE` uses a legend.
  On sequence views, `TRUE` groups the rows by group class under
  labelled strips.

- data:

  The data the model was fitted to, used only to compute the intervals.
  `NULL` (default) uses the data stored in the fit.

- main, subtitle:

  Title and subtitle. `NULL` (default) uses the view's default title and
  subtitle.

- what:

  The view to draw. `"all"` draws every view of that function that the
  model supports, as a list, and names in a message any it cannot draw.
  On fitted models other than profile, class and transition models,
  `plot_profiles()` draws only `"lines"` and `plot_diagnostics()` only
  `"posteriors"`, where that model has the view (see *Named plot
  functions on other fitted models*). For `plot_profiles()`: `"lines"`
  (default), `"bars"` (means from zero with 95% intervals), `"heatmap"`
  (means in observed standard deviations, or the response probabilities
  of categorical indicators), `"raincloud"` (each indicator's
  distribution by assigned profile), `"parallel"` (every case as a line)
  or `"pairs"` (scatter-plot matrix with each profile's 95% covariance
  ellipse).

  For `plot_diagnostics()`: `"all"` (default), `"sizes"` (effective
  profile counts), `"posteriors"` (posterior probability of each case's
  assigned profile), `"entropy"` (each case's entropy relative to a flat
  posterior) or `"avepp"` (average posterior probability matrix).
  `"group_sizes"`, `"group_posteriors"`, `"group_entropy"` and
  `"group_avepp"` draw the same four views for the groups and their
  group classes; `"all"` includes them when the fit has two or more
  group classes, and with fewer they are refused with
  `latents_no_group_classes`. A fit with a noise component supports
  `"sizes"` and the group views only: `"all"` returns `"sizes"`, and the
  other individual-level views are refused with
  `latents_unsupported_noise`.

  For `plot_selection()`: `"all"` (default), `"deltas"`, `"heatmap"`,
  `"smallest"`, `"trust"` or `"diagnostics"`, as described for
  `plot_enumeration()`.

  For `plot_enumeration()`, three kinds of view; information criteria
  cannot be requested together with the other two:

  - Information criteria: criterion columns of `as.data.frame(x)`, such
    as `"bic_groups"`, or `"all"`. On a profile or class grid the
    default draws AIC, both BICs and ICL over individuals; on a
    transition or family grid AIC, BIC and BIC over individuals; on a
    regression grid AIC, BIC, SABIC and ICL.

  - Classification quality across the candidates, at each level with two
    or more classes: `"entropy"` (relative entropy), `"avepp"` (the
    lowest AvePP, with the cutoff dashed) or `"diagnostics"` (both).

  - Model-selection evidence (none of which selects a model): `"deltas"`
    (the decrease in the criterion in `using` with each added profile
    and each added group class; see `get_results(x, "deltas")` and
    `"knees"`), `"heatmap"` (the criterion over profiles and group
    classes, relative to the lowest), `"smallest"` (the smallest class
    share, with 5% dashed) or `"trust"` (convergence, boundary and
    replication of every candidate). `"selection"` draws all four and
    `"diagnostics"`, as a list.

- category:

  Which category's probability to draw: `"last"` (default), `"first"`,
  or a single label or index.

- cell_labels:

  For sequence views, `TRUE` (default) prints each assigned profile in
  its cell when the cells are large enough for the text.

- n_group_classes, n_profiles, model:

  For `plot_tree()`, the trees and rows to draw. At the individual level
  a tree runs down the numbers of profiles (classes) at one number of
  group classes and one structure: `n_group_classes` and `model` choose
  trees, `n_profiles` the rows. At the group level a tree runs down the
  numbers of group classes at one number of profiles: `n_profiles` and
  `model` choose trees, `n_group_classes` the rows. `NULL` (default)
  keeps every converged candidate. A single count for the rows
  (`n_profiles` at the individual level, `n_group_classes` at the group
  level) is one level of every tree, and the plot compares that level's
  solutions instead:
  `plot_tree(x, n_profiles = 4, n_group_classes = c(2, 5))` draws the
  four-profile solution at two and at five group classes as two rows,
  each class in the colour of the class it matches in the first row.

- candidates:

  For `plot_tree()`, two or more candidate names drawn as the rows of
  one tree, in the order given, to compare those solutions. Not combined
  with `n_group_classes`, `n_profiles` or `model`.

- from, to:

  For `plot_tree()`, two solutions to compare as a two-row tree, each
  named by its class counts `c(n_profiles, n_group_classes)` (one count
  means one group class) or by its name:
  `plot_tree(x, from = c(3, 3), to = c(4, 3))`. Given together, and not
  with `candidates`. On a grid crossed with covariance structures,
  `model` names the structure of both.

- level:

  For `plot_tree()`, `"individuals"` (default; the cases divided into
  profiles or classes) or `"groups"` (the groups divided into group
  classes). For
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), see
  *Arguments of
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html)*.

- combine:

  `TRUE` (default) draws several criteria as panels of one plot; `FALSE`
  returns one plot per criterion. For `plot_tree()` (and
  `what = "tree"`), `FALSE` (default) returns one plot per tree and
  `TRUE` draws every tree as a panel of one plot.

- using:

  For `plot_enumeration(what = "deltas")`, `"heatmap"` or `"selection"`,
  the information criterion used. When not given, a profile or class
  grid uses `"bic_groups"` and the other grids `"bic"`; `"smallest"` and
  `"trust"` use none.

- cutoff:

  For `plot_enumeration(what = "avepp")`, `"diagnostics"` or
  `"selection"`, the AvePP cutoff drawn as a dashed line (default 0.8).

- mark_minimum:

  `TRUE` (default) circles each criterion's lowest value among converged
  candidates, without selecting a model. Information criteria only.

- ...:

  With `what = "tree"` (on
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of an
  enumeration or on `plot_enumeration()`), the arguments of
  `plot_tree()`. On
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a
  transition or group-class family grid, also `combine`, `labels`,
  `mark_minimum`, `cutoff` and `using`, passed to `plot_enumeration()`.
  Otherwise not used. Profile, class and transition models and
  enumerations refuse unknown arguments with `latents_bad_argument`; the
  other methods ignore them with a `latents_unused_argument` warning
  (see *Errors*).

## Value

A ggplot object. `what = "all"` (on
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a fitted
model or on a named plot function), `plot_enumeration(combine = FALSE)`
with several criteria or classification-quality measures, and
enumeration views when more than one is requested (including
`"selection"`) return a `latents_plots` list of ggplot objects, which
draws each one when printed, pausing between plots on an interactive
screen device. `plot_tree()` and `what = "tree"` return one ggplot for a
single tree or `combine = TRUE`, and otherwise a `latents_plots` list
named by tree. `plot_views()` returns a `data.frame` with one row per
view: `type` (the value for `what`), `group` and `description`.
`plot_views(x)` of a grid returns one row per `what` value of that grid:
`type`, `group` (`"criteria"`, `"quality"`, `"selection"` or `"tree"`),
`supported` (`TRUE` when the grid offers it), `description` and `reason`
(why it is not offered; `NA` when it is).

## Details

Returned plots can be printed, saved with
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html),
or customized with
[`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html).
Mixture-regression plots (only the `"fitted"` view for models that
classify persons), pooled imputation plots, and bootstrap-test plots are
displayed immediately, including inside loops, and return their ggplot
objects invisibly.

Named functions provide access to related plot types through `what`:
`plot_profiles()` displays profile summaries and indicator
distributions; `plot_responses()` displays categorical response
probabilities; `plot_sequences()` displays assigned profiles over
occasions; `plot_prevalence()` displays profile proportions within group
classes; `plot_diagnostics()` displays classification summaries; and
`plot_transitions()` displays transition probabilities. For enumeration
results, `plot_enumeration()` displays information criteria and
classification summaries, `plot_selection()` displays model-comparison
diagnostics, and `plot_tree()` displays correspondence between class
solutions.

These functions use the corresponding `plot(x, what = )` method.
Available views depend on the model family, as described below;
`what = "all"` requests all supported views for the selected function.
`plot_views()` lists plot types and their descriptions. For an
enumeration object, `plot_views(grid)` also reports whether each view is
supported and gives the reason when it is unavailable.

The enumeration functions support profile and class grids
([`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md),
[`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md),
[`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md)),
transition grids (`enumerate_classes(time = )`), group-class family
grids (`enumerate_classes(family = )`) and regression grids
([`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md)).
Available views are summarized below:

|  |  |  |  |  |
|----|----|----|----|----|
| View | Profiles, classes | Transitions | Families | Regressions |
| criteria | 14 columns | aic, bic, bic_individual | aic, bic, bic_individual | aic, bic, bic_rows, sabic, icl |
| entropy, avepp, diagnostics | yes | yes | group classes | yes |
| deltas, heatmap, smallest, trust, selection | yes | yes | no | yes |
| tree | yes | yes | group classes | yes |

Group-class family grids contain group classes without individual-level
profiles. Their criteria and classification summaries are plotted
against the number of group classes, and their trees compare group
assignments (`level = "groups"`, the default), with one tree per family.
Views requiring both profiles and group classes are unavailable.
Regression-grid trees use case-level posterior probabilities; the
smallest-class view uses `smallest_share`, and the `"trust"` view omits
boundary status because that column is unavailable. Unsupported views
raise `latents_unsupported_view`.

Profile order (largest first), Okabe-Ito colours, and point shapes are
consistent across views. Plotting requires the suggested package
ggplot2.

## Arguments of [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

Besides `x`, `what`, `main`, `subtitle` and the arguments above,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) takes these on
the objects named:

- `level`:

  The confidence level for additive, mixture-regression and growth views
  (default 0.95). On pooled plots, a compatibility alias for `what`:
  `"profile"`, `"group"` or `"measurement"`.

- `predictor`:

  For the fitted view of a mixture regression, a numeric predictor in
  the formula. `NULL` (default) uses the first numeric predictor.

- `time`:

  For growth and trajectory views, the numeric time column. `NULL`
  (default) uses the first numeric variable in the fitted random or
  regression design.

- `facet`:

  For trajectories, `FALSE` (default) draws all classes in one panel;
  `TRUE` draws one panel per class.

- `persons`:

  For individual trajectories, the number of persons to show (default
  12), or a vector of their identifiers.

- `max_persons`:

  The maximum number of observed trajectories drawn behind the class
  means (default 120); zero hides them.

- `spread`:

  For growth trajectories, the band of random-effect spread of persons
  around the class curve, from 0 (hidden) to less than 1; default 0.8.

- `criterion`:

  A compatibility alias for `what` in enumeration plots. On a profile
  grid it is used alone or with `what = "enumeration"`; transition and
  group-class family grids take `"bic"`, `"aic"` or `"bic_individual"`
  (or a unique abbreviation, such as `"bic_i"`), or any other view of
  the grid.

- `palette`:

  For bootstrap plots, a colour vector for the histogram and the
  observed statistic; `NULL` (default) uses the standard colours.

- `style`:

  An obsolete bootstrap plotting argument; ignored with a
  `latents_unused_argument` warning.

## Views by object class

Each object takes `what =` for its views, `main`, `subtitle` and the
arguments named here.

- Fitted profile and class models:

  From
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md),
  [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
  [`lca()`](https://pak.dynasite.org/latents/reference/lca.md) and
  covariate models: `"profiles"` (the default; a fit whose indicators
  are all categorical with non-numeric category labels, which the
  profile summaries cannot score, draws `"responses"` by default),
  `"bars"`, `"heatmap"`, `"raincloud"`, `"parallel"`, `"pairs"` (each
  profile's 95% covariance ellipse), `"responses"`, `"probabilities"`
  (not for covariate models), `"sequences"` (fits with `time`;
  `cell_labels = TRUE` prints the profile in each cell), `"sizes"`,
  `"entropy"`, `"posteriors"`, `"avepp"` and `"all"`. The profile views
  name in a caption any ordinal, count or (for the means) categorical
  indicator they omit.

- Transition fits:

  From [`lta()`](https://pak.dynasite.org/latents/reference/lta.md):
  `"transitions"` (default), then `"profiles"`, `"bars"`, `"heatmap"`,
  `"responses"`, `"sequences"`, `"sizes"`, `"entropy"`, `"posteriors"`,
  `"avepp"` and `"all"`, without `statistic`; the profile, bar and
  response views draw the 95% intervals of the transition fit's
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  (`intervals`). A general
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit takes
  no `what` and draws its transition matrix, which `plot_transitions()`
  draws too.

- Enumerations:

  Criterion columns or `"all"` (with `combine`, `labels` and
  `mark_minimum`), `"tree"` (with the arguments of `plot_tree()`), the
  classification-quality views `"entropy"`, `"avepp"` and
  `"diagnostics"` (with `cutoff`), and the model-selection views
  `"deltas"` and `"heatmap"` (which use the criterion in `using`),
  `"smallest"` and `"trust"`; `"selection"` draws those four and
  `"diagnostics"`. Transition and group-class family grids draw `"bic"`
  by default and take every view they offer (the table under
  *Description*); a regression grid draws its default criteria and takes
  the same arguments as `plot_enumeration()`.

- Additive fits:

  `"means"` (default; each group class's mean with its Wald interval),
  `"variances"` (within- and between-group variances) or `"intercepts"`
  (each group's posterior intercept); `level` sets the interval.
  `plot_profiles()` draws `"means"`.

- Cross-level fits:

  `"profiles"` (default), `"group_means"` or `"composition"` (profile
  shares within each group class). `plot_profiles()` draws `"profiles"`
  and `plot_prevalence()` `"composition"`.

- Mixture regressions:

  `"fitted"` (default: the outcome against a numeric `predictor`, one
  line per class), `"coefficients"` or `"posteriors"`; `level` sets the
  intervals. Fits that classify persons (`class_level = "group"`)
  default to `"trajectories"` (unless `predictor` is given) and also
  take `"individuals"` and `"classification"`, as growth mixtures do.
  `plot_profiles()` draws `"fitted"`, or `"trajectories"` for fits that
  classify persons, and `plot_diagnostics(what = "posteriors")` draws
  `"posteriors"`.

- Growth mixtures:

  `"trajectories"` (default: class mean trajectories with bands;
  `facet`, `max_persons`, `spread`), `"individuals"` (`persons`: a
  number or identifiers), `"coefficients"`, `"random"` (predicted random
  effects with each class's 95% ellipse) or `"classification"` (also
  `"posteriors"`); `time` names the time variable, `level` the band.
  `plot_profiles()` draws `"trajectories"` and
  `plot_diagnostics(what = "posteriors")` draws `"posteriors"`.

- Pooled imputations:

  `"profile"` (default; the probability of each profile within each
  group class), `"group"` or `"measurement"` (profile means and
  variances) parameters: each imputation's estimate beside the pooled
  one with its interval. `plot_profiles()` draws `"measurement"` and
  `plot_prevalence()` draws `"profile"`.

- Model comparisons, bootstrap tests, diagnostics:

  One view each: the criteria of each compared model, the simulated null
  distribution with the observed statistic, and the views of
  `plot_diagnostics()`.

## Named plot functions on other fitted models

Named plotting functions dispatch to the corresponding view for each
model family. For example, `plot_profiles(fit)` and
`plot(fit, what = "trajectories")` produce the same growth-mixture plot:

- Growth mixtures, and mixture regressions that classify persons
  (`class_level = "group"`):

  `plot_profiles()` draws `"trajectories"` (each class's curve over
  time) and `plot_diagnostics(what = "posteriors")` draws `"posteriors"`
  (each person's class probabilities).

- Other mixture regressions:

  `plot_profiles()` draws `"fitted"` (each class's line against the
  first numeric predictor) and `plot_diagnostics(what = "posteriors")`
  draws `"posteriors"`.

- Additive and dispersion fits:

  `plot_profiles()` draws `"means"` (each group class's means across
  indicators).

- Cross-level fits:

  `plot_profiles()` draws `"profiles"` and `plot_prevalence()` draws
  `"composition"` (profile shares within each group class).

- Pooled imputations:

  `plot_profiles()` draws `"measurement"` (the pooled profile means and
  variances beside each imputation's) and `plot_prevalence()` draws
  `"profile"` (the pooled probability of each profile within each group
  class).

- General [`lta()`](https://pak.dynasite.org/latents/reference/lta.md)
  fits:

  `plot_transitions()` draws the transition matrix.

Every other combination, including the `what` values of
`plot_profiles()` other than `"lines"` and those of `plot_diagnostics()`
other than `"posteriors"`, is refused with `latents_nothing_to_plot`
naming the views that fit's
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) offers;
`what = "all"` draws the mapped views and names the others in a message.
These functions pass on only the arguments given, and refuse with
`latents_bad_argument` one the mapped view does not take (for example
`intervals` or `scale` for a growth trajectory); use
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) for that view's
own arguments, such as `time`, `facet` or `predictor`, and for the views
no named function draws.

## Compatibility arguments

Enumeration plots accept `criterion` as an alternative to `what` (on a
profile grid, `what = "enumeration"` requests the default criteria
view). Pooled plots accept `level` as an alternative to `what`, and
bootstrap plots accept `palette` (the histogram and observed-line
colours). On an enumeration, `criterion` with any `what` other than
`"enumeration"` is an error of class `latents_bad_argument` (with
`"tree"`, a warning, as the tree uses no criterion); `level` together
with `what` on a pooled plot is an error of the same class.

## Errors

`latents_missing_package` without ggplot2; `latents_bad_argument` for an
argument the object's plot does not take, except that the plots of
additive, cross-level, general transition
([`lta()`](https://pak.dynasite.org/latents/reference/lta.md)),
growth-mixture, mixture-regression, pooled, comparison and bootstrap
results and of transition or group-class family grids ignore it with a
`latents_unused_argument` warning (including a bootstrap plot's
`style`); `latents_bad_argument` also when a named plot function is
given an object it does not draw (a named fit function given anything
other than a fitted model, or a grid function given anything other than
an enumeration); `latents_nothing_to_plot` (also of class
`latents_bad_argument`) when a named fit function, or its `what`, has no
view of that fitted model, naming the views its
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) offers (for
example `plot_transitions()` on a model without transitions, or
`plot_diagnostics(what = "sizes")` on a growth mixture);
`latents_nothing_to_plot` when the fit lacks what a view needs (for
example posteriors or entropy of a one-profile fit, prevalence of a
covariate model, or the case-level views of a centred fit);
`latents_unknown_category` when `category` names no category of an
indicator; `latents_no_continuous`, `latents_no_categorical` and
`latents_no_time` when a view needs indicators or occasions the fit does
not have; `latents_unknown_criterion` (also `latents_bad_argument`) for
a criterion not in the grid; `latents_unsupported_view` (also
`latents_bad_argument`) for a view the kind of grid does not offer;
`latents_no_tree` (also of class `latents_nothing_to_plot`) when a tree
would have fewer than two rows, or at the group level of a grid with one
group class; `latents_unknown_candidate` when `candidates`,
`n_group_classes`, `n_profiles` or `model` names something no converged
candidate has; `latents_failed_candidate` when a named candidate was not
fitted; `latents_unconverged` when a named candidate did not converge;
`latents_unsupported_noise` for a tree of a fit with a noise component.

## References

Satopaa, V., Albrecht, J., Irwin, D. and Raghavan, B. (2011). Finding a
"kneedle" in a haystack: detecting knee points in system behavior. *31st
International Conference on Distributed Computing Systems Workshops*,
166-171.

Zappia, L. and Oshlack, A. (2018). Clustering trees: a visualization for
evaluating clusterings at multiple resolutions. *GigaScience*, 7(7),
giy083.

## Examples

``` r
plot_views()
#>             type       group
#> 1       profiles measurement
#> 2           bars measurement
#> 3        heatmap measurement
#> 4      raincloud measurement
#> 5       parallel measurement
#> 6          pairs measurement
#> 7      responses measurement
#> 8  probabilities   structure
#> 9      sequences   structure
#> 10   transitions   structure
#> 11         sizes diagnostics
#> 12       entropy diagnostics
#> 13    posteriors diagnostics
#> 14         avepp diagnostics
#> 15   enumeration   selection
#> 16          tree   selection
#> 17           all       every
#>                                                                                                                               description
#> 1                                                                 Profile means, medians or modes across indicators, one line per profile
#> 2                                                                             Profile means as grouped bars from zero, with 95% intervals
#> 3                                                                    Profile means in observed standard deviations from the observed mean
#> 4                                                           Each indicator's distribution by assigned profile: density, box, observations
#> 5                                                                           Every case as a line across indicators, one panel per profile
#> 6                                                                          Scatter-plot matrix with each profile's 95% covariance ellipse
#> 7                                                                                Categorical response probabilities, one line per profile
#> 8                                                                      Profile prevalence within each group class, the two-level quantity
#> 9                                                                                Each group's profile at each occasion, one row per group
#> 10                                                              Estimated transition matrix, one panel per group class (a transition fit)
#> 11                                                                              Effective number of cases in each profile, with its share
#> 12                                                                           Each case's entropy relative to a flat posterior, by profile
#> 13                                                                      Posterior probability of each case's assigned profile, by profile
#> 14                                                                   Average posterior probability: assigned profile by posterior profile
#> 15 Information criteria across a candidate grid (plot any enumeration; or name criteria in `what`; plot_views(grid) lists a grid's views)
#> 16                                           How classes split as more are added, one plot per tree (a profile, class or transition grid)
#> 17                                                                    Every view above that this fit has the ingredients for, in one call
if (requireNamespace("ggplot2", quietly = TRUE)) {
  fit <- multilpa(course_engagement, c("browse", "lectures", "forum_read"),
                  "student", n_profiles = 2, n_group_classes = 2,
                  n_starts = 2, seed = 1)
  plot_profiles(fit)
  plot_diagnostics(fit, what = "avepp")
}
```
