
<!-- README.md is generated from README.Rmd. Please edit that file -->

<a href="https://travis-ci.org/rcannood/GillespieSSA"><img src="https://travis-ci.org/rcannood/GillespieSSA.svg" align="left"></a>

# `GillespieSSA`: Gillespie’s Stochastic Simulation Algorithm (SSA)

**GillespieSSA** provides a simple to use, intuitive, and extensible
interface to several stochastic simulation algorithms for generating
simulated trajectories of finite population continuous-time model.
Currently it implements Gillespie’s exact stochastic simulation
algorithm (Direct method) and several approximate methods (Explicit
tau-leap, Binomial tau-leap, and Optimized tau-leap).

The package also contains a library of template models that can be run
as demo models and can easily be customized and extended. Currently the
following models are included, decaying-dimerization reaction set,
linear chain system, logistic growth model, Lotka predator-prey model,
Rosenzweig-MacArthur predator-prey model, Kermack-McKendrick SIR model,
and a metapopulation SIRS model.

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

The following example models are available:

  - [Decaying-Dimerization Reaction Set (Gillespie,
    2001)](vignettes/decaying_dimer.md): `vignette("decaying_dimer",
    package="GillespieSSA")`
  - [SIRS metapopulation model (Pineda-Krch,
    2008)](vignettes/epi_chain.md): `vignette("epi_chain",
    package="GillespieSSA")`
  - [Linear Chain System (Cao et al., 2004)](vignettes/linear_chain.md):
    `vignette("linear_chain", package="GillespieSSA")`
  - [Pearl-Verhulst Logistic growth model (Kot,
    2001)](vignettes/logistic_growth.md): `vignette("logistic_growth",
    package="GillespieSSA")`
  - [Lotka predator-prey model (Gillespie, 1977; Kot,
    2001)](vignettes/lotka_predator_prey.md):
    `vignette("lotka_predator_prey", package="GillespieSSA")`
  - [Radioactive decay model (Gillespie,
    1977)](vignettes/radioactive_decay.md):
    `vignette("radioactive_decay", package="GillespieSSA")`
  - [Rosenzweig-MacArthur predator-prey model (Pineda-Krch et al.,
    2007)](vignettes/rm_predator_prey.md): `vignette("rm_predator_prey",
    package="GillespieSSA")`
  - [Kermack-McKendrick SIR model (Brown & Rothery,
    1993)](vignettes/sir.md): `vignette("sir", package="GillespieSSA")`

## Latest changes

Check out `news(package = "GillespieSSA")` or [NEWS.md](inst/NEWS.md)
for a full list of
changes.

<!-- This section gets automatically generated from inst/NEWS.md, and also generates inst/NEWS -->

### Recent changes in GillespieSSA 0.6.0

  - MAINTAINER: Maintainer has been changed to Robrecht Cannoodt.

  - DOCUMENTATION: Documentation was roxygenised and markdownised.

  - DOCUMENTATION: Port demo’s to vignettes.

  - DOCUMENTATION: Added NEWS and README.

  - DOCUMENTATION: Remove ’s from examples.

  - MINOR CHANGE: Many functions were refactorised in order to clean up
    the code.

  - MINOR CHANGE: Functions which are marked “Not intended to be invoked
    stand alone.” are no longer being exported.

  - BUG FIX: Fix warning and potential error in OTL.

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
    [doi:10.1063/1.1378322](http://dx.doi.org/10.1063/1.1378322)
  - Gillespie D.T. 2007. Stochastic simulation of chemical kinetics.
    Annu. Rev. Chem. 58:35
    [doi:10.1146/annurev.physchem.58.032806.104637](http://dx.doi.org/10.1146/annurev.physchem.58.032806.104637)
  - Kot M. 2001. Elements of mathematical ecology. Cambridge University
    Press.
    [doi:10.1017/CBO9780511608520](http://dx.doi.org/10.1017/CBO9780511608520)
  - Pineda-Krch M. 2008. Implementing the stochastic simulation
    algorithm in R. Journal of Statistical Software 25(12): 1-18.
    [doi: 10.18637/jss.v025.i12](http://dx.doi.org/10.18637/jss.v025.i12)
  - Pineda-Krch M., Blok H.J., Dieckmann U., and Doebeli M. 2007. A tale
    of two cycles — distinguishing quasi-cycles and limit cycles in
    finite predator-prey populations. Oikos 116:53-64.
    [doi:10.1111/j.2006.0030-1299.14940.x](http://dx.doi.org/10.1111/j.2006.0030-1299.14940.x)
