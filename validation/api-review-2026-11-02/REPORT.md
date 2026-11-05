# API review: diagnostics, plotting and reporting after a fit (2026-11-02)

Question: once a model is fitted, can a user get diagnostics, plots and a
report through easy, named arguments, the same way for every model family?

Evidence: `probe1.R` (refusals and accepted values for a multilpa fit),
`probe2.R` (the same verbs on lpa, multilca, lta, mixture regression, growth
mixture and an enumeration grid), `probe3.R` (what `report()` prints and
draws). Run from the repo root on the working tree (1.1.13). Vignette usage
counts come from every chunk in `vignettes/*.Rmd` and `vignettes/articles/*.Rmd`.

## What already works

- `get_results(x, what = )` is the single table accessor. It lists the
  available tables when a name is wrong, and individual tables take named
  arguments (`level`, `format`, `adjust`, `by`). Vignettes use it 189 times.
- The named plot functions (`plot_profiles()`, `plot_diagnostics()`,
  `plot_prevalence()` and the rest) go through one registry
  (`R/plot-named.R`), so `plot_profiles()` and `plot_diagnostics()` already
  draw mixture regressions and growth mixtures.
- `diagnostics()` returns a checks table with `status` and a reference for
  every check. The fit-family split is internal and invisible to the user.
- `plot_views(grid)` lists a grid's views with a `supported` column.

## Findings, in order of how much they cost the user

### A. One-call reporting stops at a third of the families

1. **`report()` refuses five of the families.** For a mixture regression,
   growth mixture, general LTA, additive or cross-level fit, and for any
   enumeration grid, it raises `latents_not_a_fit` (probe2). Every one of
   these already has `summary()`, `diagnostics()` and `plot()`. A user
   therefore has to know which families `report()` accepts before calling
   it.
2. **`report()` cannot be scoped.** On a two-level LPA it prints 300 console
   lines and draws 11 plots (probe3: profiles, bars, heatmap, raincloud,
   parallel, pairs, probabilities, sizes, entropy, posteriors, avepp). The
   only controls are `plots = TRUE/FALSE` and `rows`. There is no way to ask
   for only the diagnostics, or only two plots. This probably explains why
   the vignettes call `report()` once in 1,278 lines of code.
3. **`report()` cannot be written to a file.** Its output exists only in the
   console. knitr and rmarkdown are already in Suggests, so an HTML report
   adds no dependency.
4. **`descriptives()` refuses mixture regressions and growth mixtures**
   (probe2). It sends the user to `descriptives(data, vars = )`, so they have
   to re-specify the columns the fit already holds.

### B. Diagnostics: the right tables, but no arguments to narrow them

5. **`diagnostics(fit)` draws four plots by default**, including when the
   result is assigned (`R/report.R:184`). In a script or a knitted document,
   `quality <- diagnostics(fit)` emits figures nobody asked for. `plots`
   accepts only TRUE or FALSE, so the user cannot pick one plot, such as the
   AvePP matrix.
6. **`plot(diagnostics(fit))` refuses `what`** (probe1). The same views are
   available as `plot_diagnostics(fit, what = "avepp")`, but from the fit,
   not from the diagnostics object. That makes two paths for the same picture.
7. **No argument narrows the checks or residuals to the ones that need
   attention.** The LPA workflow article filters by hand:
   `head(subset(get_results(independent, "residuals"), statistic > 5))`. A
   `get_results(d, "checks", status = )` and a
   `get_results(d, "residuals", flagged = TRUE)` would replace that.
8. **`parameter_inference()` accepts a growth mixture but refuses a mixture
   regression** (probe2). It sends the user to `confint()` and `vcov()`. The
   two are sibling families with the same coefficient table.

### C. Plotting: inconsistent vocabulary across the methods

9. **`plot_views(fit)` refuses a fit** (probe1). It lists views for a grid or
   for nothing. For a fit, the list of views it can draw already exists
   internally (`.multilpa_supported_views()`, which `report()` uses) but is
   not exposed.
10. **Three plot methods return an unclassed `match.arg` error for a bad
    `what`**: `plot.multilpa_transitions`, `plot.latents_mixture_regression`
    and `plot.latents_growth_mixture` raise `'arg' should be one of ...`.
    `plot.multilpa` raises `latents_bad_argument` with a full sentence.
    Tests and user code cannot catch these errors by class.
11. **The same view has two names.** It is `plot_profiles(what = "lines")`
    but `plot(fit, what = "profiles")`.
12. **The criterion argument has two names and two vocabularies.**
    - Name: `plot.multilpa_enumeration(criterion = )` and
      `plot.latents_family_enumeration(criterion = )`, but
      `plot_enumeration(what = )` and
      `plot.latents_regression_enumeration(what = )`.
    - Values: `"bic"` is valid for regression, family and transition grids
      but refused for a profile grid (`latents_unknown_criterion`, probe2),
      where it must be spelled `"bic_groups"` or `"bic_individual"`.
13. **`level` means three different things.** It is a confidence level
    (0.95) in the growth, mixture-regression and additive plots, a model
    level (`"profile"`, `"group"`, `"measurement"`) in the pooled plot, and a
    classification level (`"individuals"`, `"groups"`) in `plot_tree()` and
    `get_results()`.
14. **A refusal points to an argument the function does not have.**
    `plot_profiles(lca_fit)` says `Use what = "responses" instead`, but
    `plot_profiles()` has no such view. The correct call is
    `plot_responses(fit)` (`R/plot-gg.R:363`).
15. **`plot(fit, what = )` takes one view or `"all"`**, never a chosen few.
    It already returns a `latents_plots` list for `"all"`, so a vector of
    views needs no new return type.

## Proposed surface (additive only; 0.9.21 code keeps working)

```r
report(fit)                                        # every family and every grid
report(fit, sections = c("summary", "diagnostics"))
report(fit, plots = c("profiles", "avepp"))        # TRUE = all, FALSE = none, or view names
report(fit, file = "engagement-report.html")       # self-contained HTML via rmarkdown (Suggests)

diagnostics(fit, plots = "avepp")                  # TRUE / FALSE / view names
plot(diagnostics(fit), what = "avepp")
get_results(diagnostics(fit), "checks", status = c("check", "fail"))
get_results(diagnostics(fit), "residuals", flagged = TRUE)

plot_views(fit)                                    # views this fit can draw, `supported` column
plot(fit, what = c("profiles", "avepp"))           # several views -> latents_plots
plot_profiles(fit, what = "profiles")              # accepted beside "lines"
plot_enumeration(grid, what = "bic")               # the grid's primary BIC, documented per grid
descriptives(growth_fit)                           # outcome and predictors by class
```

Compatibility: every item widens what an argument accepts or turns a refusal
into a result. No formal is removed or reordered, and new formals go before
`...`. Items 10 and 14 change only the class and the wording of an error.
Item 5's default (`plots = TRUE`) stays, because changing it would change
0.9.21 output. `level` (item 13) stays as it is, because renaming it would
break released code. The fix there is a sentence in each Rd stating which
meaning applies.

## Suggested order

1. Quick fixes that need no design decision: 10, 14, 9 (`plot_views(fit)`),
   11 (alias), 6 (`plot(diagnostics, what = )`), 15 (vector `what`).
2. Reporting: 1 (`report()` for every family and grid), 2
   (`sections =`, `plots =` as view names), 5 (`diagnostics(plots =)` view
   names), 7 (`status =`, `flagged =`), 4 (`descriptives()` for the
   regression families).
3. Needs owner decisions: 3 (`report(file = )`: HTML only, or also .md/.docx;
   figure size), 12 (what `"bic"` means on a two-level profile grid:
   `bic_groups`, which is the package default `using`), 8 (Wald
   `parameter_inference()` for mixture regressions or not).
