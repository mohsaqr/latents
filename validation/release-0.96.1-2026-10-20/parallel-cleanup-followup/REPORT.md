# Socket worker cleanup follow-up

The final Windows CI failure identifies a production cleanup defect: `parLapply()` detects a dead worker, but `on.exit(parallel::stopCluster(cluster))` then tries to send DONE to the broken socket. That cleanup error replaces `latents_parallel_failure`, and aborts shutdown before all owned connections are closed. The two remaining connections later cause GC warnings.

Added `.latents_stop_cluster()` in `R/parallel.R`. It stops nodes independently through the public `parallel::stopCluster()` API, preserving graceful termination for live workers. Each node's connection is explicitly closed in `finally`, including when the graceful stop errors. Errors from this best-effort cleanup cannot replace the original condition. The helper does not suppress warnings or close unrelated connections.

`tests/testthat/test-parallel.R` now includes a deterministic mock of simultaneous worker-read and cleanup-write failures: before the fix, four assertions fail (condition class/message and two leaked connections). It also exercises a real cluster containing an already-closed socket and a live worker, and repeated cleanup of the resulting closed sockets. The original worker-death test remains enabled.

Evidence:

- `pre-fix.log`: four reproduced deterministic failures against original cleanup (CRAN-skipped socket cases omitted only in this initial mock reproduction).
- `lint.log`: configured source lint clean.
- `targeted-final.log`, `results.rds`: authoritative full parallel-file run using local socket access and an absolute R_LIBS path to the matching installed source package.
- `cross-review/probe.R`, `cross-review/probe.log`: independent review accepts normal results, original remote errors, real SIGTERM worker death, per-node stop failures, absence of warnings, restoration of baseline connection counts, preservation of an unrelated open connection and repeated cleanup.

The initial `targeted.log` is a restricted-sandbox run with socket-creation failures. `targeted-socket.log` uses a mismatched worker namespace and is superseded by the matching-source run. Neither is the final validation result.

Final validation: **37 passing expectations**, zero failures, errors, warnings or skips. Source is frozen for release validation.
