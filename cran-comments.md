# GillespieSSA 0.6.1

This release contains a major rewrite of the internal code, to make sure
the code is readable and that the algorithm doesn't continuously update
the local environment.

* MAJOR CHANGE: Instead of passing `"D"`, `"ETL"`, `"OTL"`, or `"BTL"` to `ssa()`,
  it is expected to pass `ssa.d()`, `ssa.etl()`, `ssa.otl()`, or `ssa.btl()`. 
  This cleans up parameter setting clutter in the `ssa()` function.
  
* MAJOR CHANGE: Rewrite `ssa.*()` and `ssa.*.diag()` as
  `ssa_step.ssa_*()` and `ssa_step_diag.ssa_*()` S3 functions.
  
* MAJOR CHANGE: Do not save the current state in the function environment. 
  Instead, simply save it in a local variable. 
  
* MAJOR CHANGE: Precompile propensity functions instead of evaluating them 
  as R code at each iteration.
  
* MAJOR CHANGE: Clean up and merge `ssa.run()`, `ssa.terminate()`, `ssa.check.args()` 
  and `ssa.check.method()` into `ssa()`.

## Test environments
* local Fedora install, R 3.6.0
* ubuntu 16.04 (on travis-ci), R 3.6.0
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 note

## Revdepcheck

revdepcheck was run with the following command:
```
revdepcheck::revdep_check()
```

It resulted in no new errors or warnings for reverse dependencies.

```
✔ hybridModels 0.3.5                     ── E: 0     | W: 0     | N: 0                                                   
OK: 1                                                                                                                  
BROKEN: 0
Total time: 5 min
```
