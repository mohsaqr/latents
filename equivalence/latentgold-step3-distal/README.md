# Distal outcome of a latent Markov chain against Latent GOLD 6.1 (pending)

**Status: kit prepared 2026-11-15, not yet run.** Latent GOLD was not
installed on the machine where `three_step(occasion =)` was written, so this
directory holds the inputs only: no listing is in `returned/`, no test reads
it, and no claim of agreement with Latent GOLD is made for this feature. The
evidence that does back it is listed at the end.

`three_step(fit, data, outcome, occasion = t, method = "ml")` on an `lta()`
fit relates a person-level outcome to the profile the chain occupies at
occasion `t`. Each occasion is classified from its own indicators (the step
two of the three-step latent Markov model, as in `../latentgold-step3-markov/`),
the assignments at every occasion are indicators of the chain with the pooled
classification errors held fixed, the chain's initial and transition
probabilities are held at the `lta()` estimates, and the outcome is normal
with a mean and variance per profile. Latent GOLD's reference is
`step3 modal ml` on a `dynamic` latent variable with the outcome as a
`dependent` variable observed on the record of occasion `t` only.

| case | design | persons | occasion | outcome by profile |
|---|---|---|---|---|
| d01 | 2 profiles, 4 occasions | 500 | 3 | means 0, 1.5; sds 1, 0.8 |
| d02 | 3 profiles, 4 occasions | 600 | 2 | means -1, 0, 1.2; sds 1, 1.2, 0.9 |

Files in `kit/`: `<case>.dat` (tab-separated; each occasion's step-two
posteriors `post1...`, the outcome `distal` on the chosen occasion's record and
`.` elsewhere), `<case>.lgs` (chain fixed at the `lta()` logits with `~fix`,
the comparison of record), `<case>-free.lgs` (chain re-estimated by Latent
GOLD, expected to differ slightly), `<case>-fit.rds`, `<case>-data.rds`,
`cases.rds`, and `latents-results.csv`, latents' answer when the kit was
written.

## To run

1. `Rscript equivalence/latentgold-step3-distal/make-kit.R` (only if the kit
   must be rebuilt; the inputs above are retained).
2. Run each `.lgs` in batch with Latent GOLD 6.1 and copy the listings to
   `returned/<case>.lst` and `returned/<case>-free.lst`.
3. First check the `~fix` syntax took effect: the listing must print the
   initial and transition logits equal to the values in the `.lgs` braces,
   without standard errors. The label syntax and the order of the braced
   values were written from the manual, not executed.
4. Write `equivalence/test-latentgold-step3-distal.R` on the pattern of
   `test-latentgold-step3-markov.R`: means and variances of `distal` by state
   within half a unit of Latent GOLD's last printed digit, standard errors with
   `vcov_type`-equivalent information. Add its rule to
   `verification-matrix/test-attribution.csv`, its datasets to
   `verification-matrix/test-datasets.R`, and remove the `pending` row for it
   from `verification-matrix/candidates.csv`.

## Evidence in the meantime

- `tests/equivalence/test-lta-distal.R`: the chain posterior by a sum over
  every path, the ML step three maximized by `optim()` on that sum with a
  numDeriv sandwich, BCH means and their sandwich from the inverse error rows,
  the uncorrected methods and a categorical outcome; all agree to 1e-4 or
  better, and the test fails when the transition matrix is transposed or the
  error matrix inverted the wrong way round.
- `validation/lta-distal-simulation.R`: predeclared recovery of the profile
  and pair means by BCH and ML (bias within 3.5 Monte Carlo standard errors,
  95% coverage within its 99.9% binomial range).
