# GillespieSSA 0.6.2

* MINOR CHANGE: Allow using `.t` as parameter in the propensity functions.

## Test environments
* local Fedora install, R 4.1
* ubuntu 20.04 (on travis-ci), R 4.0
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 note

## Revdepcheck

revdepcheck resulted in no new errors or warnings for reverse dependencies.

```
> revdepcheck::revdep_check(timeout = as.difftime(600, units = "mins"), num_workers = 30)
── CHECK ───────────────────────────────────────────────────────────────────────────────── 2 packages ──
✓ hybridModels 0.3.7                     ── E: 0     | W: 0     | N: 0                                                                                             
✓ GillespieSSA2 0.2.8                    ── E: 0     | W: 0     | N: 1                                                                                             
OK: 2                                                                                                                                                            
BROKEN: 0
Total time: 14 min
```

