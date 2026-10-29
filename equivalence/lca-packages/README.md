# lca() against five independent R implementations of latent class analysis

Until this folder `lca()` had one software reference with current results
(poLCA, in the harness). Here every case is fitted by `lca()` and by up to five
independent maximum-likelihood implementations, each from many starts, and
compared on the maximized log likelihood, the class shares, every item response
probability and, where the package returns them per row, every posterior.

```
Rscript equivalence/lca-packages/compare.R      # comparison.csv, RUN.txt
```

`equivalence/test-lca-packages.R` runs the same comparison under testthat.

**Result (2026-10-29): 211 quantities compared, 203 agree, 8 known
differences in multilevLCA (below), 0 unexplained.** With the HELP data:
glca, BayesLCA and depmixS4 agree to 7.3e-7 in the log likelihood and 3.4e-6
in every probability, and all 11 published values agree.

| reference | algorithm and settings | compared on |
|---|---|---|
| glca 1.4.2 | EM, 50 starts, `eps = 1e-12` | likelihood, shares, probabilities, posteriors |
| randomLCA | EM with `penalty = 0` (its default 0.01 penalizes the likelihood), 50 trials, `EMtol = 1e-12`, one core | likelihood, shares, probabilities |
| BayesLCA | `blca.em` with uniform priors (`alpha = beta = delta = 1`), so the posterior mode is the ML estimate; 50 restarts, `conv = 1e-12` | likelihood, shares, probabilities |
| multilevLCA | `multiLCA()` single level, `tol = 1e-12`, sequential (one core) | likelihood, shares, probabilities |
| depmixS4 | `mix()` with multinomial responses, 20 random starts, EM to 1e-15 | likelihood, shares, probabilities, posteriors |

Datasets: poLCA's `carcinoma` (2 and 3 classes), `values`, `cheating` (the four
cheating items) and `gss82` (polytomous: 2 and 3 classes; glca and depmixS4
only); BayesLCA's `Alzheimer` (2 and 3 classes); randomLCA's `myocardial`,
`dentistry` (2 and 3 classes) and `hivtests`, expanded from their frequency
tables to one row per respondent; and the HELP study's four binary indicators
of the SAS and R blog's Example 8.21 (`help-lca.csv`, 431 rows, three classes;
`PROVENANCE.txt`). randomLCA 1.1.4 refuses the HELP model as not identifiable
(one residual degree of freedom, a probability at 1), so it is left out there.

**Published (HELP, three classes):** Kleinman and Horton (2011), *SAS and R*
blog, Example 8.21 prints poLCA's fit (log likelihood -1032.484, 14
parameters, AIC 2092.967, BIC 2149.893, shares .0785 / .2162 / .7053) and SAS
PROC LCA's (log likelihood -1032.48, P(housed) 1.0000 / .2703 / .5625). All
agree to their printed digits. SAS's prevalences (.2163 / .0785 / .7052)
differ from poLCA's in the fourth decimal and are not compared.

Tolerances, fixed before any comparison was run: the log likelihood to 1e-6;
shares, probabilities and posteriors to 1e-4. Classes are matched once per
package by the ordering closest in shares and probabilities.

## Reading each package on the same scale

- **BayesLCA** reports `logpost`, the log likelihood plus the log density of
  its priors. With uniform priors that density is constant: 0 for each
  Beta(1, 1) and lgamma(G) for the Dirichlet(1, ..., 1) on G shares, so its
  three-class values were exactly log 2 = 0.693 above every other program's
  until the constant was subtracted.
- **depmixS4** stops on a relative criterion. At 1e-12 it stopped 3e-4 short in
  the posteriors of the flat three-class dentistry model; at 1e-15 every gap is
  within 1e-4. There it ends when the likelihood falls by floating-point noise
  (1e-11), which it reports as a warning: those notices are counted in the log,
  not hidden. Its shares and probabilities are read from its EM fixed point:
  the mean posteriors and the posterior-weighted category frequencies.
- **randomLCA** penalizes the likelihood by default (`penalty = 0.01`), which
  moves the optimum; it is run with `penalty = 0`.

## Known differences (`known-differences.csv`)

- **multilevLCA bounds every response probability to [1e-5, 1 - 1e-5].** Where
  the maximum has probabilities at 0 or 1 (carcinoma at 2 and 3 classes,
  myocardial, hivtests, HELP) its log likelihood is lower by 4e-5 to 0.0058;
  its shares and probabilities agree to within the bound.
- **multilevLCA runs one k-means start.** On Alzheimer with three classes it
  stops at a lower maximum (-744.968 against -743.484). glca, randomLCA,
  BayesLCA and depmixS4 all reach latents' maximum with the same parameters.

Nothing in latents disagreed with any package.
