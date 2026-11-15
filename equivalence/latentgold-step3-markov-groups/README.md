# Three-step latent Markov model with group classes against Latent GOLD 6.1 (pending)

`r3step(fit, data, initial_covariates =, transition_covariates =)` on an
`lta(n_group_classes = 2)` fit is the bias-adjusted three-step mixed latent
Markov model: each occasion is classified from its own indicators and
assigned to its modal profile, and step three is a mixture of hidden Markov
chains for the assignments, one per group class, with the classification
errors held fixed and each group class's initial and transition logits
regressed on the covariates. Latent GOLD's reference is `step3 modal ml`
with a static latent `GClass` beside the dynamic `State` whose posteriors it
reads (`make-kit.R` writes the syntax):

```
latent GClass nominal 2 coding=last,
   State nominal dynamic coding=last posterior=(post1 post2);
equations
   GClass <- 1;
   State[=0] <- 1 | GClass + w | GClass;
   State <- (~tra) 1 | State[-1] GClass + (~tra) z | State[-1] GClass;
```

**Status: pending, not run.** Latent GOLD was not available on the machine
where this kit was prepared (2026-11-15), so no listing is kept and no
Latent GOLD number is claimed anywhere. `returned/` is to hold `g01.lst` and
`g02.lst`; until then `equivalence/test-latentgold-step3-markov-groups.R`
skips, and a skipped test earns no credit in the verification matrix. Two
things are to be confirmed on the first run: that `step3` accepts a latent
variable without `posterior=` (GClass), and the layout of the Regression
Parameters table with a second conditioning variable, for which the test
then needs a coefficient parser (it compares the log likelihood only).

To run it:

```
Rscript equivalence/latentgold-step3-markov-groups/make-kit.R   # writes kit/
# copy kit/*.dat and kit/*.lgs to the Latent GOLD machine, run each .lgs in
# batch, and copy the listings back as returned/<case>.lst
```

| case | design | persons |
|---|---|---|
| g01 | 2 profiles, 5 occasions, 2 group classes (.5/.5); initial on w, transitions on z, every coefficient per group class | 500 |
| g02 | 2 profiles, 6 occasions, 2 group classes (.6/.4); transitions on z | 500 |

Until the listings return, the group-class chain is verified by an
independent likelihood (a sum over every group class and latent path) and by
depmixS4 on the expanded state space (`tests/equivalence/test-lta-three-step-groups.R`),
and its recovery of known parameters by
`validation/lta-three-step-groups-simulation.R`.
