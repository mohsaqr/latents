## Resubmission

This addresses the review of 0.9.14 ("reset to user's options() ... in your
inst folder"). Every vignette now saves the options it changes in its first
chunk and restores them in its last, so the installed `inst/doc` scripts
leave the user's options unchanged. No example, vignette or demo changes
`par()` or the working directory.

## Test environments

* local macOS (aarch64), R 4.5.2
* win-builder, R-release (4.6.1)

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.
* "Vermunt" (possibly misspelled) is the surname of the cited author.

`student_esm` is a subset of an openESM dataset under CC-BY 4.0; the source
and licence are in `?student_esm`.
