One of our lab's project is closely related to this 
package. As such, I would like to provide custody to
this package in order to ensure it is well maintained
over time.

Current improvements are simply refactorisations in 
order to update the code to modern R standards. 

Future improvements include:

* Major speed-ups due to the way code is currently evaluated
  at run-time.
* Avoid using the environment of `ssa.run()` in order to
  store the algorithms state vector.

These improvements will be made whilst keeping the
interface of the package fixed as not to break legacy users.

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
