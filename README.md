
<!-- README.md is generated from README.Rmd. Please edit that file -->

<a href="https://travis-ci.org/rcannood/GillespieSSA"><img src="https://travis-ci.org/rcannood/GillespieSSA.svg" align="left"></a>
<a href="https://codecov.io/gh/rcannood/GillespieSSA">

# `GillespieSSA`: Gillespie’s Stochastic Simulation Algorithm (SSA)

**GillespieSSA** provides a simple to use, intuitive, and extensible
interface to several stochastic simulation algorithms for generating
simulated trajectories of finite population continuous-time model.
Currently it implements Gillespie’s exact stochastic simulation
algorithm (Direct method) and several approximate methods (Explicit
tau-leap, Binomial tau-leap, and Optimized tau-leap). The package also
contains a library of template models that can be run as demo models and
can easily be customized and extended. Currently the following models
are included, decaying-dimerization reaction set, linear chain system,
logistic growth model, Lotka predator-prey model, Rosenzweig-MacArthur
predator-prey model, Kermack-McKendrick SIR model, and a metapopulation
SIRS model.

## Install

You can install **GillespieSSA** from CRAN using

``` r
install.packages("GillespieSSA")
```

Or, alternatively, you can install the development version of
GillespieSSA from GitHub
using

``` r
devtools::install_github("rcannood/GillespieSSA", build_vignettes = TRUE)
```

## Examples

Individual demo models can be run by issuing `demo(<model name>)`,
alternatively all of the demo models can be run using
`demo(GillespieSSA)`. The following example models are available:

  - [Decaying-Dimerisation Reaction Set](demo/decayingDimer.R)
    (Gillespie, 2001)
  - [Metapopulation SIRS model](demo/epiChain.R) (Pineda-Krch, 2008)
  - [Linear Chain System](demo/linearChain.R) (Cao et al., 2004)
  - [Logistic Growth (Pearl-Verhulst model)](demo/logisticGrowth.R)
    (Kot, 2001; Pineda-Krch, 2008)
  - [Lotka Predator-Prey model](demo/lotka.R) (Gillespie, 1977; Kot,
    2001)
  - [Radioactive Decay model](demo/radioactiveDecay.R) (Gillespie, 1977)
  - [Rosenzweig-MacArthur Predator-Prey model](demo/rma.R) (Pineda-Krch
    et al., 2007, Pineda-Krch, 2008)
  - [Kermack-McKendrick SIR model](demo/sir.R) (Brown & Rothery, 1993)

## Latest changes

Check out `news(package = "GillespieSSA")` or [NEWS.md](inst/NEWS.md)
for a full list of
changes.

<!-- This section gets automatically generated from inst/NEWS.md, and also generates inst/NEWS -->

### Recent changes in GillespieSSA 0.6.0 (2019-07-12)

  - MAINTAINER: Maintainer has been changed to Robrecht Cannoodt.

  - DOCUMENTATION: Documentation was roxygenised and markdownised.

  - DOCUMENTATION: Added NEWS and README.

  - MINOR CHANGE: Functions which are marked “Not intended to be invoked
    stand alone.” are no longer being exported.

### Recent changes in GillespieSSA 0.5-4 (2010-08-16)

  - DOCUMENTATION: Fix typos in documentation.

## References

  - Brown D. and Rothery P. 1993. Models in biology: mathematics,
    statistics, and computing. John Wiley & Sons.
  - Cao Y., Li H., and Petzold L. 2004. Efficient formulation of the
    stochastic simulation algorithm for chemically reacting systems. J.
    Chem. Phys. 121:4059-4067.
    [doi:10.1063/1.1778376](http://dx.doi.org/10.1063/1.1778376)
  - Cao Y., Gillespie D.T., and Petzold L.R. 2006. Efficient step size
    selection for the tau-leaping method. J. Chem. Phys. 124:044109.
    [doi:10.1063/1.2159468](http://dx.doi.org/10.1063/1.2159468)
  - Cao Y., Gillespie D.T., and Petzold L.R. 2007. Adaptive explicit
    tau-leap method with automatic tau selection. J. Chem. Phys.
    126:224101.
    [doi:10.1063/1.2745299](http://dx.doi.org/10.1063/1.2745299)
  - Chatterjee A., Vlachos D.G., and Katsoulakis M.A. 2005. Binomial
    distribution based tau-leap accelerated stochastic simulation. J.
    Chem. Phys. 122:024112.
    [doi:10.1063/1.1833357](http://dx.doi.org/10.1063/1.1833357)
  - Gillespie D.T. 1977. Exact stochastic simulation of coupled chemical
    reactions. J. Phys. Chem. 81:2340.
    [doi:10.1021/j100540a008](http://dx.doi.org/10.1021/j100540a008)
  - Gillespie D.T. 2001. Approximate accelerated stochastic simulation
    of chemically reacting systems. J. Chem. Phys. 115:1716-1733.
    \<<http://dx.doi.org/10.1063/1.1378322> \>
  - Gillespie D.T. 2007. Stochastic simulation of chemical kinetics.
    Annu. Rev. Chem. 58:35
    [doi:10.1146/annurev.physchem.58.032806.104637](http://dx.doi.org/10.1146/annurev.physchem.58.032806.104637)
  - Kot M. 2001. Elements of mathematical ecology. Cambridge University
    Press.
    [doi:10.2277/052180213X](http://dx.doi.org/10.2277/052180213X)
  - Pineda-Krch M. 2008. Implementing the stochastic simulation
    algorithm in R. Submitted to the Journal of Statistical Software
    25(12): 1-18.
    [doi: 10.18637/jss.v025.i12](http://dx.doi.org/10.18637/jss.v025.i12)
  - Pineda-Krch M., Blok H.J., Dieckmann U., and Doebeli M. 2007. A tale
    of two cycles — distinguishing quasi-cycles and limit cycles in
    finite predator-prey populations. Oikos 116:53-64.
    [doi:10.1111/j.2006.0030-1299.14940.x](http://dx.doi.org/10.1111/j.2006.0030-1299.14940.x)
