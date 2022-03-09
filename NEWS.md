# GillespieSSA 0.6.2

* MINOR CHANGE: Allow using `.t` as parameter in the propensity functions.

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
  
# GillespieSSA 0.6.0 (2019-07-15)

* MAINTAINER: Maintainer has been changed to Robrecht Cannoodt.

* DOCUMENTATION: Documentation was roxygenised and markdownised.

* DOCUMENTATION: Port demo's to vignettes.

* DOCUMENTATION: Added NEWS and README.

* DOCUMENTATION: Remove \dontrun's from examples.

* MINOR CHANGE: Many functions were refactorised in order to clean up the code.

* MINOR CHANGE: Functions which are marked "Not intended to be invoked stand alone."
  are no longer being exported.

* BUG FIX: Fix warning and potential error in OTL.

# GillespieSSA 0.5-4 (2010-08-16)

* DOCUMENTATION: Fix typos in documentation.

# GillespieSSA 0.5-3 (2008-01-17)

No changes, it seems.

# GillespieSSA 0.5-2 (2008-01-17)

* MINOR CHANGE `ssa.run()`: Throw error when propensity values are negative.

* DOCUMENTATION: Add radioactiveDecay demo

# GillespieSSA 0.5-1 (2008-01-15)

* MINOR CHANGE `ssa()`: Move chunks of code to separate functions `ssa.check.args()`, 
  `ssa.check.method()` and `ssa.run()`.
  
* DOCUMENTATION: Expand documentation of many functions.

# GillespieSSA 0.5-0 (2007-10-19)

* MINOR CHANGE `ssa()`: Minor changes to how the states and the propensity functions are evaluated.

* DOCUMENTATION: Rework demos.

# GillespieSSA 0.4-0 (2007-10-12)

* MINOR CHANGES: Minor updates to package to bring it up to date with the JSS manuscript.

* DOCUMENTATION: Added the metapopulation SIRS model described in the JSS manuscript as a demo model to the package

# GillespieSSA 0.3-1 (2007-10-04)

* BUG FIX: bug fixes.

* DOCUMENTATION: fixed error in url of package home page.

* MINOR CHANGES: Revised how model parameters are passed to the `ssa()` (per Ben Bolker's suggestion). 
  Model parameters are now passed directly to `ssa()` as a named vector and 
  hence they do not have to be declared in the global environment.

# GillespieSSA 0.2-0 (2007-09-21)

* BUG FIX: mainly bug fixes.

# GillespieSSA 0.1-0 (2007-08-01)

* Initial release 

* Presented at useR! 2007 (as poster)
