# Independent covariance follow-up review

Disposition: **block pending repair of active-constraint identification**. No tracked source or tests were edited by this reviewer.

The log-eigenvalue design, gradient/Hessian, null-space search, multiplier sign, and `expm1` objective-difference identity are consistent with the stated convex optimization problem. The relative descent acceptance protects the true objective. Positive multipliers alone, however, do not establish the KKT conditions if a supposedly active constraint has nonzero slack.

## Reproduced blocker

```r
scatter <- rbind(c(1e-5, 1), c(1e-4, 100))
result <- .multilpa_bounded_variances(scatter, c(1, 1), "EVI", 1)
```

The exact optimum has both rows `c(1, 50.5)`. Both first variances are at the floor. Given those floors, equal volume forces a common second variance, whose maximum-likelihood estimate is `(1 + 100) / 2 = 50.5`. The two active-floor multipliers are strictly positive; these feasible KKT conditions establish the unique global optimum of the strictly convex log-coordinate objective.

The new solver returns approximately:

```
     [,1]        [,2]
[1,] 1.00002996  50.49850206
[2,] 1.00000000  50.50001498
```

Maximum variance error is `0.00149794`; objective excess is `5.9296e-7`. It reports `settled = TRUE` after seven iterations. Its retained active rows are 1 and 3, but row 1 still has log slack `2.9959204659e-5` and multiplier `0.00989628`. Complementary slackness fails.

The active-boundary selection compares step-length differences with `1e-10 * max(1, alpha)`. For a very small feasible Newton step, this admits a constraint whose step limit is near zero in absolute terms while its actual remaining slack is appreciable. Adding it to the active set freezes that slack, so the subsequent projected-gradient/multiplier checks certify the wrong equality-constrained problem.

Evidence: `active-kkt-probe.R/.log` records the public helper, exact reference, objective excess, final active set, slacks, multipliers, iteration count, and settlement status. `tiny-scatter-sweep.log` shows the same systematic issue over nearby positive scatter scales. Source inspection at review time places the boundary selector near line 860 of `R/covariance-structure.R`.

The parent and covariance owner received the concrete reproducer. A safe correction should identify genuinely binding constraints from appropriately scaled final slack/step error and check complementary slackness before treating an active-set stationary point as a constrained optimum. Tests should include this exact KKT oracle without loosening existing tolerances.
