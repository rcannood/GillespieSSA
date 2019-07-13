One of our lab's projects is closely related to this 
package. As such, I would like to provide custody to
this package in order to ensure it is well maintained
over time. 

Special care will be taken to ensure the interface of the package
remains the unchanged, as not to break the package for legacy users.

## List of changes
* DOCUMENTATION: Documentation was roxygenised and markdownised.

* DOCUMENTATION: Port demo's to vignettes.

* DOCUMENTATION: Added NEWS and README.

* DOCUMENTATION: Remove \dontrun's from examples.

* MINOR CHANGE: Many functions were refactorised in order to clean up the code.

* MINOR CHANGE: Functions which are marked "Not intended to be invoked stand alone."
  are no longer being exported.

* BUG FIX: Fix warning and potential error in OTL.

## Test environments
* local Fedora install, R 3.6.0
* ubuntu 16.04 (on travis-ci), R 3.6.0
* win-builder (devel and release)

Win-builder produced the following NOTE:
```
Possibly mis-spelled words in DESCRIPTION:
  Gillespie's (3:9, 20:10)
```

I tried multiple variations of
Gillespie's, 'Gillespie's', and 'Gillespie's, 
all resulted in the same NOTE. Feel free to 
suggest alternatives.


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
