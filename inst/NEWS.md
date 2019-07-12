# GillespieSSA 0.6.0 (2019-07-12)

* MAINTAINER: Maintainer has been changed to Robrecht Cannoodt.

* DOCUMENTATION: Documentation was roxygenised and markdownised.

* DOCUMENTATION: Added NEWS and README.

* MINOR CHANGE: Functions which are marked "Not intended to be invoked stand alone."
  are no longer being exported.
  
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

* Presented at useR! 2007 (as [poster](http://user2007.org/program/posters/pineda-krch.pdf))
