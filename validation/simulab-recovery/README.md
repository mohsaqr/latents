# Recovery of a known two-level latent profile model

`recovery-study.R` generates data from the model `multilpa()` estimates,
refits it, and asks how close the estimates get to the values that produced the
data. The generator is `simulab::simulate_multilpa()`; the comparison is
`simulab::validate_recovery()`. Both are sibling packages, so the study needs no
machinery of its own beyond matching the arbitrary mixture labels.

Run from the project root:

```
Rscript validation/simulab-recovery/recovery-study.R
```

It needs `simulab` >= 0.4.3 installed, which is where `simulate_multilpa()` was
added, and it writes `tmp/simulab-recovery-terms.csv`,
`tmp/simulab-recovery-diagnostics.csv`, `tmp/simulab-recovery-enumeration.csv`
and `tmp/simulab-recovery-study.rds`.

## The generating model

Four continuous indicators, three individual profiles, two latent cluster
classes. Profile means are `-1.3 / 0.0 / 1.3` on reading and correspondingly
spaced on maths, science and writing, all with unit standard deviations and no
within-profile correlation.

The two cluster classes share that measurement model entirely. They differ only
in profile prevalence:

| | Profile 1 | Profile 2 | Profile 3 |
|---|---|---|---|
| Class 1 | 0.60 | 0.30 | 0.10 |
| Class 2 | 0.10 | 0.30 | 0.60 |

The middle profile is equally common in both classes, so the cluster classes are
not a relabelled copy of the profiles and cannot be recovered from the
measurement model alone. Cluster classes are equally likely.

Two conditions, 100 datasets each, 20 EM starts per fit:

- **60 clusters x 20** individuals (n = 1200)
- **30 clusters x 10** individuals (n = 300)

Mixture labels are arbitrary, so every fit is matched to the generating model
before anything is compared: profiles by the permutation of estimated mean
vectors closest to the true ones, then cluster classes by the permutation of
prevalence rows closest to the true ones, in the profile order just fixed. With
three profiles and two classes the search is exhaustive and therefore exact.

## What was recovered

Every one of the 200 fits converged and none hit the variance bound.

### Fit diagnostics, averaged over 100 replications

| condition | individual entropy | group entropy | individual accuracy | group accuracy |
|---|---|---|---|---|
| 60 clusters x 20 | 0.692 | 0.989 | 0.844 | 0.997 |
| 30 clusters x 10 | 0.744 | 0.922 | 0.790 | 0.960 |

Accuracy is the share of individuals assigned to their generating profile, and
of clusters assigned to their generating class, after label matching. The
cluster classes are recovered almost perfectly at the larger size: 81 of 100
replications placed every one of the 60 clusters correctly, and the worst
replication still reached 0.967. Individual assignment sits near 0.84 because
the profiles genuinely overlap -- a relative entropy of 0.69 is what a
four-indicator design with unit-variance profiles 1.3 apart looks like, and no
estimator can do better than the posterior allows.

### Parameter recovery

Mean bias over replications, with its Monte Carlo standard error, and the root
mean squared error of the bias:

| condition | quantity | bias | MC SE | RMSE | within 0.1 |
|---|---|---|---|---|---|
| 60 clusters x 20 | profile mean | 0.0013 | 0.0025 | 0.088 | 0.79 |
| 60 clusters x 20 | profile sd | -0.0046 | 0.0014 | 0.048 | 0.96 |
| 60 clusters x 20 | prevalence | -2e-17 | 0.0016 | 0.038 | 0.98 |
| 60 clusters x 20 | class proportion | -3e-18 | 0.0049 | 0.070 | 0.86 |
| 30 clusters x 10 | profile mean | -0.0223 | 0.0076 | 0.265 | 0.40 |
| 30 clusters x 10 | profile sd | -0.0226 | 0.0037 | 0.094 | 0.65 |
| 30 clusters x 10 | prevalence | -2e-17 | 0.0045 | 0.110 | 0.75 |
| 30 clusters x 10 | class proportion | -4e-18 | 0.0069 | 0.098 | 0.64 |

At 60 clusters of 20 the estimator is unbiased for every quantity: no bias
exceeds its own Monte Carlo standard error by a meaningful margin, the largest
absolute bias in any of the 32 individual terms is 0.016 (on
`sd[Profile 2, reading]`), and no cell of the prevalence matrix is off by more
than 0.005. The prevalence and class
proportion biases are exactly zero to machine precision because both truth and
estimate sum to one within a cluster class, so the errors must cancel; the RMSE
column, not the bias column, is what carries information for those two rows.

At 30 clusters of 10 the picture is the familiar one for mixtures at small n:
still centred, but much more variable, and standard deviations are pulled down
by about 2%. The variability is concentrated in the middle profile, whose mean
terms have RMSEs of 0.30 to 0.41 against 0.17 to 0.26 for the outer profiles --
it is the profile with no distinctive indicator pattern, defined only by sitting
between the other two.

Local optima are the mechanism. At 60 clusters of 20, all 20 starts reached the
same best likelihood in all 100 replications. At 30 clusters of 10 the mean was
16.7 of 20, and in at least one replication only a single start found the best
solution. The worst replication in the study (30 clusters of 10, replication 20)
converged, with an unremarkable individual entropy of 0.819, yet classified only
48.3% of individuals into their generating profile -- convergence and entropy
are not evidence that the right solution was found.

### Structure selection

`enumerate_classes()` was run on the first 25 datasets of each condition over
2 to 4 profiles and 1 to 3 cluster classes. The table is the proportion of
datasets on which the criterion's minimum fell on the true (3 profiles,
2 cluster classes).

| criterion | 60 clusters x 20 | 30 clusters x 10 |
|---|---|---|
| BIC (groups) | 1.00 | 0.28 |
| BIC (individuals) | 1.00 | 0.00 |
| SABIC (individuals) | 1.00 | 0.60 |
| CAIC (individuals) | 0.96 | 0.00 |
| AIC | 0.24 | 0.24 |
| ICL (individuals) | 0.00 | 0.00 |

Every criterion identified two cluster classes in essentially every dataset;
the disagreement is entirely about the number of profiles. At the larger size
the three BIC-family criteria are exact or near-exact. AIC over-extracts, taking
four profiles on 11 of 25 datasets in both conditions, which is the expected
behaviour of a criterion whose penalty does not grow with n.

ICL never selected the true model, in either condition, always settling on two
profiles. That is not a defect: ICL adds the estimated classification entropy to
the BIC penalty, and with individual entropy at 0.69 the three-profile solution
is penalised for exactly the overlap that the design builds in. ICL answers
"how many well-separated profiles are there", which is a different question from
"how many profiles generated the data".

The two sample-size conventions come apart where it matters. At 30 clusters of
10 the group-count convention selected the true model on 7 of 25 datasets while
the individual-count convention selected it on none, because the individual
count of 300 imposes a penalty heavy enough to collapse the middle profile. The
package documents the two as alternative conventions rather than interchangeable
criteria, and this is what that difference costs at small n.

## Reading the output

`validate_recovery()` returns one row per parameter per replication, which the
script writes to `tmp/simulab-recovery-terms.csv` with a `condition` and
`replication` column. Terms are named `mean[Profile 1, reading]`,
`sd[Profile 2, maths]`, `probability[Class 1, Profile 3]` and
`class proportion[Class 2]`, so the file can be aggregated any way needed
without re-running the study.

Eight warnings were raised over the whole run, all during enumeration of
misspecified structures: effective memberships below one when four profiles were
fitted to three-profile data, and k-means start failures on the same models. No
warning arose from a fit at the true structure.
