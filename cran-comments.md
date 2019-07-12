I'm working on refactoring the GillespieSSA package
such that the `ssa.run()` function does not use
its environment in order to store the state values.
This will prevent name clashes with user-defined
state names and any variables used in the `ssa.run()`
function.

This update is the first step in updating this package
to modern R standards.

## Test environments
* local Fedora install, R 3.6.0
* ubuntu 14.04 (on travis-ci), R 3.6.0
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 note
