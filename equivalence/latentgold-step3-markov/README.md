# Three-step latent Markov model against Latent GOLD 6.1

`r3step(fit, data, initial_covariates =, transition_covariates =)` on an
`lta()` fit is the bias-adjusted three-step latent Markov model (Vermunt
2010; Di Mari, Oberski and Vermunt 2016): each occasion is classified from
its own indicators, assigned to its modal profile, and the step-three chain
regresses the initial profile and the transitions on covariates with the
classification errors held fixed. Latent GOLD's reference is `step3 modal ml`
on a `dynamic` latent variable, which reads the posteriors from the data
file, makes the modal assignments and the error matrix itself, and fits the
chain. With `coding=last` (initial logits against the last profile) and
`(~tra)` (transition logits against staying), its coefficients are ours.

`make-kit.R` simulates three cases, fits `lta()` and writes each occasion's
step-two posteriors with the covariates; Latent GOLD was run in batch under
Wine on 2026-10-17 and its listings are kept in `returned/`, so the test
(`equivalence/test-latentgold-step3-markov.R`) runs offline.

| case | design | persons | latents LL | Latent GOLD LL | coefficients | max estimate gap | max SE gap |
|---|---|---|---|---|---|---|---|
| m01 | 2 profiles, 5 occasions, initial on w, transitions on z | 400 | -1216.789662 | -1216.7897 | 6 | 4.5e-05 | 4.2e-05 |
| m02 | 3 profiles, 4 occasions, transitions on z | 500 | -1943.345779 | -1943.3458 | 14 | 5.0e-05 | 4.6e-05 |
| m03 | 2 profiles, up to 6 occasions, sequences ending early | 450 | -1122.407980 | -1122.4080 | 6 | 3.9e-05 | 4.8e-05 |

Every gap is inside Latent GOLD's printed rounding (four decimals). Standard
errors are compared with `vcov_type = "observed"`, the information Latent
GOLD's `standarderrors` reports; `r3step()`'s default is the sandwich over
persons.

The step-two posteriors are each occasion's own: p(k | y_it) proportional to
w_k f_k(y_it), with the transition fit's measurement and its profile shares
pooled over occasions. The transition fit's smoothed posteriors would bias
the transitions: in a first draft that used them, a move-out intercept came
out -1.51 where the one-step `lta()` gave -0.98, because the chain had
already removed apparent moves that the step-three correction then removed
again.
