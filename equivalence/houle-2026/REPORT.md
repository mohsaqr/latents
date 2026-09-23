# Houle et al. (2026): package comparison and numerical checks

Reviewed 19 September 2026 against multilpa 0.5.2 and the corrected source tree
described in the [whole-package audit](../MATH_AUDIT.md).

**The core `multilpa()` model corresponds structurally to the paper's
dispersion-heterogeneity family. It does not implement all six families, and its
joint estimation does not reproduce the supplement's fixed-measurement stages.**
The supplied fit table passes independent information-criterion checks. The
participant data and fitted output needed to reproduce the published estimates
were not present in the supplied materials or the listed supplements checked.

## Sources and scope

- Houle, S. A., Morin, A. J. S., and Harvey, J.-F. (2026), *Multilevel Latent
  Profile Analyses: A Comprehensive Guide*, Organizational Research Methods,
  [DOI 10.1177/10944281261469432](https://doi.org/10.1177/10944281261469432).
  The user supplied the complete 35-page article PDF.
- The user supplied `sj-docx-1-orm-10.1177_10944281261469432.docx`: Table S1
  and the complete illustrative Mplus syntax. Both documents were read in full.
- The matching local `sj-csv-2-orm-10.1177_10944281261469432.csv` contains six
  variable descriptions in two text columns, **not participant records**. The
  [publisher's supplement listing](https://journals.sagepub.com/doi/10.1177/10944281261469432)
  lists the DOCX and this 0.24 KB CSV. The syntax references `Data.dat`,
  `LPA.csv`, and `Fscores.dat`; those data are not embedded in the DOCX.

Page references below use the article's printed page numbers. Syntax findings
were checked against the original DOCX XML, not just converted plain text.
No Mplus executable was run for this comparison. The source materials are not
copied into the repository; the small numerical fixtures record provenance.

## Which model is implemented?

For continuous indicators, the package's core observed-data likelihood is

\[
L=\prod_j\sum_{h=1}^{H}\omega_h
       \prod_{i\in j}\sum_{k=1}^{K}\pi_{hk}
       \phi(y_{ij};\mu_k,\Sigma_k).
\]

Here \(h\) is a group class and \(k\) an individual profile. Individual
measurement parameters \(\mu_k,\Sigma_k\) are shared across group classes;
the conditional profile proportions \(\pi_{hk}\) vary. This is the key
dispersion-heterogeneity structure on pp. 15–16. Diagonal covariance and the
appropriate equal/varying variance restrictions must also match before a
particular specification is comparable.

| Article family | Defining feature | Current package scope |
|---|---|---|
| Additive | Classes of groups differ in between-level means; within-level variances remain common across those classes. | No dedicated between-level measurement model. |
| Dispersion | Classes of groups differ in within-level variances, with between-level means/variances constrained. | No group-class-specific residual dispersion model. |
| Additive-dispersion | Group classes can differ in both location and dispersion. | Not implemented as this two-level measurement family. |
| Dispersion-heterogeneity | Group classes differ in their mix of individual profiles. | Structural match to `multilpa()`, with shared individual measurement. |
| Restricted cross-level | Separate individual and group profiles; individual profile composition does not define group profiles. | No separate group-indicator density. |
| Full cross-level | Group indicators and individual profile composition both inform group profiles. | No separate group-indicator density. |

The first three families are described on pp. 11–14; the cross-level families
on pp. 15–16. Means or covariances derived from a mixture's composition are
not extra between-level measurement parameters. `variance_model = "varying"`
varies measurement variances across **individual profiles**, not across the
paper's group-level classes of additive/dispersion models.

`fit_random_intercept()` adds one scalar normally distributed group shift with
unit loadings on all continuous indicators. It has no discrete group classes
and does not supply the missing general between-level measurement families.

## Estimation and preprocessing are separate compatibility requirements

The paper recommends retaining a primary-level solution and **fixing its
measurement estimates** while estimating the other level (p. 16). Its examples
use Mplus `@` constraints. Contextual applications start with the individual
level; climate applications can start with the group level.

In `multilpa()`, `start` supplies initial values; subsequent EM iterations update
the measurement parameters. `max_iter = 0` scores a complete supplied parameter
set but does not optimize mixing probabilities while holding measurement fixed.
Consequently, even a structurally compatible dispersion-heterogeneity fit is
not a replication of that staged procedure. Fixed-measurement estimation would
need an explicit API and corresponding parameter-count and inference changes.

The predictor/outcome extensions on pp. 22–23 also fix measurement and free
selected membership regressions or outcomes. `fit_covariates()` jointly
re-estimates measurement and uses individual predictor slopes shared across
group classes. The supplement also illustrates slopes that vary across group
classes. The package's three-step/BCH procedures are separate methods and must
not be labeled as that same fixed-measurement likelihood.

The paper distinguishes contextual and climate constructs, grand-mean and
group-mean centering, manifest aggregation, and latent/factor-score approaches
(pp. 4–8). The package's internal centering is only a reversible numerical
translation. It does not perform the paper's group-mean centering, factor-score
estimation, measurement-error correction, or separate group-level aggregation.

The examples' MLR and `STARTS = 10000 1000 100` settings are also not identical
to simply choosing a package `n_starts` value. Entropy remains descriptive, and
agreement of information criteria does not establish agreement of estimators,
standard errors, likelihood-ratio reference distributions, or fitted classes.

## Numerical checks

### Table S1: all 24 rows pass

The table contains eight additive-dispersion models, eight individual-level
LPA models, and eight cross-level models. With \(q\) the printed number of
free parameters and \(N=10{,}000\), every row is jointly consistent with

\[
\begin{aligned}
\mathrm{AIC}&=-2\ell+2q,\\
\mathrm{BIC}&=-2\ell+q\log N,\\
\mathrm{CAIC}&=\mathrm{BIC}+q,\\
\mathrm{ABIC}&=-2\ell+q\log((N+2)/24).
\end{aligned}
\]

The checker intersects intervals for the unknown unrounded deviance implied
by all four criteria and the printed likelihood. This accounts for three-decimal
rounding rather than demanding spurious exact equality. The sample-size
interval inferred from BIC minus AIC uniquely identifies integer N = 10,000.
The package calls the paper's ABIC `sabic`; use the `individuals` convention
when comparing this table. The package's default group-count BIC remains a
different, explicitly documented convention.

This validates **penalty arithmetic**, not the published likelihoods or their
maximization. The validation script compares the public `information_criteria()`
output with independently computed formulas using the table's likelihoods and
parameter counts as inputs: **96/96 criterion values agree**, with maximum
absolute difference zero. All 24 fixture rows were also verified directly
against the original DOCX table.

### Independent synthetic likelihood check

`check-likelihood.R` uses the five printed individual-profile means/variances
from the dispersion-heterogeneity example, with explicitly invented mixing
weights and nine observations in three groups. It enumerates every latent
assignment (two group classes and five possible profiles for each member),
then compares likelihood and both levels of posterior probabilities with
`multilpa(..., max_iter = 0)`.

- Log-likelihood difference: **0** at printed machine precision.
- Maximum group posterior difference: **3.33e-16**.
- Maximum individual posterior difference: **5.55e-16**.
- One EM iteration changes a supplied measurement mean by **0.22807735**,
  confirming that `start` does not fix measurement. The deliberately truncated
  run issues the expected nonconvergence warning.

These are synthetic checks of the model structure and API behavior; they are
not a reproduction of the study or evidence that this nine-row fit is identified.

## Source issues that affect reproducibility

Several illustrative syntax blocks cannot be executed unchanged:

| Supplement section | Confirmed issue |
|---|---|
| Dispersion-heterogeneity, second stage | Declares `C(4)` but specifies `C#5`; declares X1/X2 but uses undeclared `burn`/`eng`. |
| Fixed-L1 cross-level examples | Declare `WC(2)` but specify five within classes; some model variable names differ from the declarations. |
| Contextual additive-dispersion/dispersion examples | Use aggregate MX1/MX2 without including them in `USEV`. |
| Doubly latent factor-score cross-level example | Heading says L2 measurement is freely estimated, while numeric `@` constraints fix it. |
| Manifest climate cross-level example | Contains malformed `MX2;@.3` syntax. |
| Doubly latent factor-score predictor/outcome example | Includes `X1_W@1 X2_W@-1` and `X1_W@-1 X2_W@1` outside mean brackets: these literally specify negative Gaussian variances. Variable names also disagree with declarations. |

These observations concern the distributed illustrative code. They do not
establish that the authors used those erroneous lines in their fitted analyses.
In particular, the intended correction to a negative variance cannot be
inferred safely without clarification or the actual fitted input/output.

There is also a **parameter-count ambiguity** in Table S1's final block. For
H group classes, its counts are \(5H+3\). Given two group indicators with free
diagonal means/variances and five fixed individual measurement profiles, that
count corresponds to \(4H\) group measurement parameters, \(H-1\) group
weights, and four *shared* individual-profile weights. If individual-profile
weights vary freely across group classes as the heading states, the count
would instead be \(4H+(H-1)+4H=9H-1\). Additional constraints or different
fitted syntax could resolve this ambiguity; the table alone cannot.

## Reproduce the checks

Run from the package root, with the existing development dependencies installed:

```sh
Rscript equivalence/houle-2026/check-fit-indices.R
Rscript equivalence/houle-2026/check-likelihood.R
```

No estimator source was changed in this paper comparison. The earlier full
audit's 1,381 passing assertions, 974 external comparisons, and successful
package check apply to the implemented package models; they do not certify all
six families in this article. Exact study replication remains contingent on
participant data, fitted input/output, preprocessing details, and matching
measurement constraints.
