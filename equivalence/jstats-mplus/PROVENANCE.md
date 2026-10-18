# Mplus mixture runs from JStats

Copied on 2026-10-18 from the sibling JStats repository (`../JStats`, HEAD
`8dfac8d1ecc9481c2bde4f6b7cd90f4c061e67ed`), directory
`validation/r-reference/`. Every file is a byte-for-byte copy; the md5 of
each copy equals the md5 of its source at that HEAD. The outputs are Mplus 9
runs (`TYPE = MIXTURE`); each `.dat` is the file named in its output's
`DATA: FILE IS` line, and `mplus-demo-savedata-cprob-out.dat` is the
`SAVEDATA: SAVE = CPROBABILITIES` file of `mplus-demo-savedata-cprob`.
They are compared with this package by `equivalence/test-jstats-mplus.R`.

| File | JStats path | md5 |
|---|---|---|
| `mplus-demo-lca-4cat-tech10.dat` | `validation/r-reference/mplus-demo-lca-4cat-tech10.dat` | `305b3ca784110117512a25abf0851b4f` |
| `mplus-demo-lca-4cat-tech10.out.txt` | `validation/r-reference/mplus-demo-lca-4cat-tech10.out.txt` | `455aec5a4a40b4c7013833579a2ba63b` |
| `mplus-demo-lca-count-tech11.dat` | `validation/r-reference/mplus-demo-lca-count-tech11.dat` | `db4b64b55d96cd74a47e45515fe55005` |
| `mplus-demo-lca-count-tech11.out.txt` | `validation/r-reference/mplus-demo-lca-count-tech11.out.txt` | `f2ff588e260bd5d2ec2e5ade5ca122f3` |
| `mplus-demo-lca-cov-entropy.dat` | `validation/r-reference/mplus-demo-lca-cov-entropy.dat` | `9c0521e0981f2d2ea9849a6e1c7aff11` |
| `mplus-demo-lca-cov-entropy.out.txt` | `validation/r-reference/mplus-demo-lca-cov-entropy.out.txt` | `f9041af6bdf67f4074d25dd95ae3d8e1` |
| `mplus-demo-lca-cov-k3-tech14.dat` | `validation/r-reference/mplus-demo-lca-cov-k3-tech14.dat` | `9c0521e0981f2d2ea9849a6e1c7aff11` |
| `mplus-demo-lca-cov-k3-tech14.out.txt` | `validation/r-reference/mplus-demo-lca-cov-k3-tech14.out.txt` | `9d2da7e5beac1d64d83aa756c6cd95ad` |
| `mplus-demo-lca-cov-tech10.dat` | `validation/r-reference/mplus-demo-lca-cov-tech10.dat` | `9c0521e0981f2d2ea9849a6e1c7aff11` |
| `mplus-demo-lca-cov-tech10.out.txt` | `validation/r-reference/mplus-demo-lca-cov-tech10.out.txt` | `f707e55008de9bd75b9bdd431989c61e` |
| `mplus-demo-lca-covariate-tech11.dat` | `validation/r-reference/mplus-demo-lca-covariate-tech11.dat` | `4d97e644ed09e3efa444cbce0d4539b3` |
| `mplus-demo-lca-covariate-tech11.out.txt` | `validation/r-reference/mplus-demo-lca-covariate-tech11.out.txt` | `7e68715368cf5ca39b859ad4e1dddcee` |
| `mplus-demo-lca-de3step-multi.dat` | `validation/r-reference/mplus-demo-lca-de3step-multi.dat` | `312cc9b6a263cc8c62e310409eb6e19b` |
| `mplus-demo-lca-de3step-multi.out.txt` | `validation/r-reference/mplus-demo-lca-de3step-multi.out.txt` | `91b937b3ff2fbcb864f8be1c13d69258` |
| `mplus-demo-lca-du3step-multi.dat` | `validation/r-reference/mplus-demo-lca-du3step-multi.dat` | `312cc9b6a263cc8c62e310409eb6e19b` |
| `mplus-demo-lca-du3step-multi.out.txt` | `validation/r-reference/mplus-demo-lca-du3step-multi.out.txt` | `fe216c8be5167ebd2410dfaf0738c2ba` |
| `mplus-demo-lca-entropy.dat` | `validation/r-reference/mplus-demo-lca-entropy.dat` | `b50a787cb65c76aa1aae1d9c9c300f78` |
| `mplus-demo-lca-entropy.out.txt` | `validation/r-reference/mplus-demo-lca-entropy.out.txt` | `afc3dedbc95faecda702e00de236c7f8` |
| `mplus-demo-lca-k1starts-tech11.dat` | `validation/r-reference/mplus-demo-lca-k1starts-tech11.dat` | `c59cf5183edcae9a6f8a60da0df4014f` |
| `mplus-demo-lca-k1starts-tech11.out.txt` | `validation/r-reference/mplus-demo-lca-k1starts-tech11.out.txt` | `812de63c12da1ef6cfb9509c8c624627` |
| `mplus-demo-lca-k3-tech10.dat` | `validation/r-reference/mplus-demo-lca-k3-tech10.dat` | `c59cf5183edcae9a6f8a60da0df4014f` |
| `mplus-demo-lca-k3-tech10.out.txt` | `validation/r-reference/mplus-demo-lca-k3-tech10.out.txt` | `055f02b1dd4dee79823d8da4d50cc5fe` |
| `mplus-demo-lca-k3-tech11.dat` | `validation/r-reference/mplus-demo-lca-k3-tech11.dat` | `c59cf5183edcae9a6f8a60da0df4014f` |
| `mplus-demo-lca-k3-tech11.out.txt` | `validation/r-reference/mplus-demo-lca-k3-tech11.out.txt` | `927a1a7ed36a15a62b8977b19405cd6c` |
| `mplus-demo-lca-k3-tech14.dat` | `validation/r-reference/mplus-demo-lca-k3-tech14.dat` | `c59cf5183edcae9a6f8a60da0df4014f` |
| `mplus-demo-lca-k3-tech14.out.txt` | `validation/r-reference/mplus-demo-lca-k3-tech14.out.txt` | `d039082f8e6c9ab04a9ec0ab60c02bc9` |
| `mplus-demo-lca-k4-tech11.dat` | `validation/r-reference/mplus-demo-lca-k4-tech11.dat` | `a4809527bf207e1a5b686f3de77fff5f` |
| `mplus-demo-lca-k4-tech11.out.txt` | `validation/r-reference/mplus-demo-lca-k4-tech11.out.txt` | `beee9e8542108c6589e9166255d5f29d` |
| `mplus-demo-lca-k5-tech11.dat` | `validation/r-reference/mplus-demo-lca-k5-tech11.dat` | `8998f4a93968660fa2e4e933f4a90782` |
| `mplus-demo-lca-k5-tech11.out.txt` | `validation/r-reference/mplus-demo-lca-k5-tech11.out.txt` | `0f373340aa175464b1dd0b273a4c6df2` |
| `mplus-demo-lca-mixed-tech10.dat` | `validation/r-reference/mplus-demo-lca-mixed-tech10.dat` | `ca46541926e0ae6a553081415366c867` |
| `mplus-demo-lca-mixed-tech10.out.txt` | `validation/r-reference/mplus-demo-lca-mixed-tech10.out.txt` | `c864cbcc9a07a05dd1d96ac7fa990209` |
| `mplus-demo-lca-mixed-tech14.dat` | `validation/r-reference/mplus-demo-lca-mixed-tech14.dat` | `ca46541926e0ae6a553081415366c867` |
| `mplus-demo-lca-mixed-tech14.out.txt` | `validation/r-reference/mplus-demo-lca-mixed-tech14.out.txt` | `326ab84c0421c56e7dd9f7b0ece88322` |
| `mplus-demo-lca-nominal-k3-tech11.dat` | `validation/r-reference/mplus-demo-lca-nominal-k3-tech11.dat` | `32b7485948988bd9e375b085495161f0` |
| `mplus-demo-lca-nominal-k3-tech11.out.txt` | `validation/r-reference/mplus-demo-lca-nominal-k3-tech11.out.txt` | `4a1dba58f8c1458bf2fe7bad5bfd7a41` |
| `mplus-demo-lca-nominal-tech11.dat` | `validation/r-reference/mplus-demo-lca-nominal-tech11.dat` | `1d5409c53dc8a098ed5c315fc5b134fa` |
| `mplus-demo-lca-nominal-tech11.out.txt` | `validation/r-reference/mplus-demo-lca-nominal-tech11.out.txt` | `d700aff2146f34ee9eeb5699371d12f1` |
| `mplus-demo-lca-poly-tech10.dat` | `validation/r-reference/mplus-demo-lca-poly-tech10.dat` | `360fa7b23fcdf0ade5deb91bd7244e87` |
| `mplus-demo-lca-poly-tech10.out.txt` | `validation/r-reference/mplus-demo-lca-poly-tech10.out.txt` | `352c7ba4b6a452ef51ea9ac5335a89ad` |
| `mplus-demo-lca-r3step-joint.dat` | `validation/r-reference/mplus-demo-lca-r3step-joint.dat` | `9c0521e0981f2d2ea9849a6e1c7aff11` |
| `mplus-demo-lca-r3step-joint.out.txt` | `validation/r-reference/mplus-demo-lca-r3step-joint.out.txt` | `efb32d32e48ca25d677f1028fbcd8b5b` |
| `mplus-demo-lca-tech10.dat` | `validation/r-reference/mplus-demo-lca-tech10.dat` | `ce20c0b21ea2c90bfc8aa8305afa9f93` |
| `mplus-demo-lca-tech10.out.txt` | `validation/r-reference/mplus-demo-lca-tech10.out.txt` | `ca3e856c67bdeb856530da6935bb9e76` |
| `mplus-demo-lca-tech14-lrtstarts.dat` | `validation/r-reference/mplus-demo-lca-tech14-lrtstarts.dat` | `ce20c0b21ea2c90bfc8aa8305afa9f93` |
| `mplus-demo-lca-tech14-lrtstarts.out.txt` | `validation/r-reference/mplus-demo-lca-tech14-lrtstarts.out.txt` | `6b88c277086545d94ca30ced8cc3d23d` |
| `mplus-demo-lca-tech14-tech11.dat` | `validation/r-reference/mplus-demo-lca-tech14-tech11.dat` | `ce20c0b21ea2c90bfc8aa8305afa9f93` |
| `mplus-demo-lca-tech14-tech11.out.txt` | `validation/r-reference/mplus-demo-lca-tech14-tech11.out.txt` | `74e54907803dc20eb233b3189f46677b` |
| `mplus-demo-lca-tech14.dat` | `validation/r-reference/mplus-demo-lca-tech14.dat` | `ce20c0b21ea2c90bfc8aa8305afa9f93` |
| `mplus-demo-lca-tech14.out.txt` | `validation/r-reference/mplus-demo-lca-tech14.out.txt` | `b0858c5e443e7a49004a4d896a815e6c` |
| `mplus-demo-lpa-constraint-new.dat` | `validation/r-reference/mplus-demo-lpa-constraint-new.dat` | `9192bc71f45a5f0570df971a3dc8bff8` |
| `mplus-demo-lpa-constraint-new.out.txt` | `validation/r-reference/mplus-demo-lpa-constraint-new.out.txt` | `1dec28bd44ffe61de1af351eb5414051` |
| `mplus-demo-lpa-cov-tech14.dat` | `validation/r-reference/mplus-demo-lpa-cov-tech14.dat` | `fe25d208ceb8004fd160beb7fbd724b1` |
| `mplus-demo-lpa-cov-tech14.out.txt` | `validation/r-reference/mplus-demo-lpa-cov-tech14.out.txt` | `739971bcdf34c3fadd297c0188cdb103` |
| `mplus-demo-lpa-eei-k3-tech11.dat` | `validation/r-reference/mplus-demo-lpa-eei-k3-tech11.dat` | `3dcd10d774927112a3f2047f9e7dec65` |
| `mplus-demo-lpa-eei-k3-tech11.out.txt` | `validation/r-reference/mplus-demo-lpa-eei-k3-tech11.out.txt` | `759eed2b6d2b9097227ae5d12e6e78f9` |
| `mplus-demo-lpa-ex710-tech14.dat` | `validation/r-reference/mplus-demo-lpa-ex710-tech14.dat` | `aeef8cb1c4cd476048fe328792b2b1fa` |
| `mplus-demo-lpa-ex710-tech14.out.txt` | `validation/r-reference/mplus-demo-lpa-ex710-tech14.out.txt` | `a71e63c57fc5ca67cf250bd32937f549` |
| `mplus-demo-lpa-k3-tech11.dat` | `validation/r-reference/mplus-demo-lpa-k3-tech11.dat` | `6b77f20ef0c9c1c7210f1aeac7a9b84f` |
| `mplus-demo-lpa-k3-tech11.out.txt` | `validation/r-reference/mplus-demo-lpa-k3-tech11.out.txt` | `91e5ea987953ba23236940e582d0a47e` |
| `mplus-demo-lpa-k3-vvi-tech14.dat` | `validation/r-reference/mplus-demo-lpa-k3-vvi-tech14.dat` | `6b77f20ef0c9c1c7210f1aeac7a9b84f` |
| `mplus-demo-lpa-k3-vvi-tech14.out.txt` | `validation/r-reference/mplus-demo-lpa-k3-vvi-tech14.out.txt` | `124bc0b068bffa81af678204d75596a8` |
| `mplus-demo-lpa-tech14.dat` | `validation/r-reference/mplus-demo-lpa-tech14.dat` | `aeef8cb1c4cd476048fe328792b2b1fa` |
| `mplus-demo-lpa-tech14.out.txt` | `validation/r-reference/mplus-demo-lpa-tech14.out.txt` | `77a1e56f35354158883a29b4e8829aa4` |
| `mplus-demo-lta-classvarying-measurement.dat` | `validation/r-reference/mplus-demo-lta-classvarying-measurement.dat` | `c764b60decc5f08e92f5e47885addaa2` |
| `mplus-demo-lta-classvarying-measurement.out.txt` | `validation/r-reference/mplus-demo-lta-classvarying-measurement.out.txt` | `02de5587b44d5cc2ce4a881f983c7d1d` |
| `mplus-demo-miterations-uncapped.dat` | `validation/r-reference/mplus-demo-miterations-uncapped.dat` | `cd8bb519e1d84a8c7bbbf0984d3bb81e` |
| `mplus-demo-miterations-uncapped.out.txt` | `validation/r-reference/mplus-demo-miterations-uncapped.out.txt` | `a14a703df88e0702bd3651f19f58cd77` |
| `mplus-demo-mixreg-classspecific.dat` | `validation/r-reference/mplus-demo-mixreg-classspecific.dat` | `ec2b135e55d7de985bc52f7458dc2418` |
| `mplus-demo-mixreg-classspecific.out.txt` | `validation/r-reference/mplus-demo-mixreg-classspecific.out.txt` | `4a00f4ee5489e1e08fbd903713405f17` |
| `mplus-demo-savedata-cprob-out.dat` | `validation/r-reference/mplus-demo-savedata-cprob-out.dat` | `455c71cf25c4840356fe98a51727240f` |
| `mplus-demo-savedata-cprob.dat` | `validation/r-reference/mplus-demo-savedata-cprob.dat` | `be4bd9f7008866b534edfd71639d91c0` |
| `mplus-demo-savedata-cprob.out.txt` | `validation/r-reference/mplus-demo-savedata-cprob.out.txt` | `2c3cdbb9bf40ddfacdf57d32f212670f` |
| `mplus-demo-training-partial-off.dat` | `validation/r-reference/mplus-demo-training-partial-off.dat` | `059dc73c5ea36169484ddb88f7c77573` |
| `mplus-demo-training-partial-off.out.txt` | `validation/r-reference/mplus-demo-training-partial-off.out.txt` | `a5fa4f8bb0b8cee47a7f6b2e824fec31` |
| `mplus-ug-ex7.3-de3step.out.txt` | `validation/r-reference/mplus-ug-ex7.3-de3step.out.txt` | `00dc109a8fd2a0d82097c652f4f03113` |
| `mplus-ug-ex7.3-du3step.out.txt` | `validation/r-reference/mplus-ug-ex7.3-du3step.out.txt` | `384dd662c0d6b7d2ae771360a4f654e7` |
| `mplus-ug-ex7.3-r3step.out.txt` | `validation/r-reference/mplus-ug-ex7.3-r3step.out.txt` | `8fbe3c34b3363eba81f0386c38bf7cd7` |
| `mplus-ug-ex7.3-three-step.dat` | `validation/r-reference/mplus-ug-ex7.3-three-step.dat` | `c9f1356aa46eb5b2c3b19948f94ebf66` |
| `mplus-ug-ex7.3-two-step-du3step-both.out.txt` | `validation/r-reference/mplus-ug-ex7.3-two-step-du3step-both.out.txt` | `8bdbf1245e2f8addb1af90a660554498` |
| `mplus-ug-ex7.3-two-step-du3step.out.txt` | `validation/r-reference/mplus-ug-ex7.3-two-step-du3step.out.txt` | `aa74eee4690517dd2962e79883269722` |
| `mplus-ug-ex7.3-two-step.dat` | `validation/r-reference/mplus-ug-ex7.3-two-step.dat` | `297bc76bfa159a22dd50b198400a8305` |

## Runs not copied

Of the 95 mixture outputs in `validation/r-reference/`, 61 have the data
file their `DATA: FILE IS` line names (a 62nd, `mplus-lca-cprob`, names
`lca_cprob.dat`, which is not there; `mplus-lca-cprob.dat` is, but nothing
records that it is the same file). The 19 of those 61 not copied are outside
what this package fits as the same model: factor mixtures (`BY`: fmm-tech11, fmm-tied-k3,
fmm-classvarying-k3, fmm-classvarying-loadings, factor-irt-tech8-4item,
skewnormal-cfa, the last also skew-normal), growth mixtures whose occasion
residual variances Mplus frees by occasion while `mixture_regression()` has
one per class (gmm-covariate, gmm-k3, gmm-quadratic, ug-ex8.1-tech11),
censored indicators (lca-censored), KNOWNCLASS (knownclass-lpa), TRAINING
(training-partial-on), equality labels on thresholds or variances
(lca-poly-threshold-equality, lca-threshold-equality-tech10-tech14,
mixture-classspecific-constraints), class-varying direct effects `u ON x`
(lca-classvarying-direct-effect), two latent class variables with different
numbers of classes (lta-three-status), and a run Mplus stopped after three
iterations (miterations-capped, MITERATIONS = 3).
