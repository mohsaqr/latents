# Changelog

## latents 0.8.1

- A documentation website, built with pkgdown, is published at
  <https://pak.dynasite.org/latents/>.

## latents 0.8.0

### Package renamed

- The package is renamed from multilpa to **latents**, a name that
  covers latent profile, latent class and latent transition models and
  the planned extensions. The model functions keep their names
  ([`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md),
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md)), as do
  the result classes. Condition classes now use the `latents_` prefix
  (for example `latents_bad_argument`, previously
  `multilpa_bad_argument`), `multilpa_plot_types()` is
  [`plot_views()`](https://pak.dynasite.org/latents/reference/plot_views.md),
  and the conditions catalogue is
  [`?"latents-conditions"`](https://pak.dynasite.org/latents/reference/latents-conditions.md).
  The vignettes are
  [`vignette("lpa", package = "latents")`](https://pak.dynasite.org/latents/articles/lpa.md),
  `"evaluation"`, `"covariates"`, `"lca"` and `"lta"`.

### New features

- [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
  fits a two-level latent class model:
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  with every indicator categorical, so the items are named once. Mixed
  models remain `multilpa(categorical = )`.
