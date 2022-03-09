
-   [`GillespieSSA`: Gillespie’s Stochastic Simulation Algorithm
    (SSA)](#gillespiessa-gillespies-stochastic-simulation-algorithm-ssa)
    -   [Install](#install)
    -   [Examples](#examples)
    -   [Latest changes](#latest-changes)
        -   [Recent changes in GillespieSSA
            0.6.2](#recent-changes-in-gillespiessa-062)
        -   [Recent changes in GillespieSSA
            0.6.1](#recent-changes-in-gillespiessa-061)
    -   [References](#references)

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
GillespieSSA from GitHub using

``` r
devtools::install_github("rcannood/GillespieSSA", build_vignettes = TRUE)
```

## Examples

The following example models are available:

-   Decaying-Dimerization Reaction Set (Gillespie, 2001):
    `vignette("decaying_dimer", package="GillespieSSA")`
-   SIRS metapopulation model (Pineda-Krch, 2008):
    `vignette("epi_chain", package="GillespieSSA")`
-   Linear Chain System (Cao et al., 2004):
    `vignette("linear_chain", package="GillespieSSA")`
-   Pearl-Verhulst Logistic growth model (Kot, 2001):
    `vignette("logistic_growth", package="GillespieSSA")`
-   Lotka predator-prey model (Gillespie, 1977; Kot, 2001):
    `vignette("lotka_predator_prey", package="GillespieSSA")`
-   Radioactive decay model (Gillespie, 1977):
    `vignette("radioactive_decay", package="GillespieSSA")`
-   Rosenzweig-MacArthur predator-prey model (Pineda-Krch et al., 2007):
    `vignette("rm_predator_prey", package="GillespieSSA")`
-   Kermack-McKendrick SIR model (Brown & Rothery, 1993):
    `vignette("sir", package="GillespieSSA")`

## Latest changes

Check out `news(package = "GillespieSSA")` or [NEWS.md](NEWS.md) for a
full list of changes.

<!-- This section gets automatically generated from inst/NEWS.md -->

### Recent changes in GillespieSSA 0.6.2

-   MINOR CHANGE: Allow using `.t` as parameter in the propensity
    functions.

### Recent changes in GillespieSSA 0.6.1

This release contains a major rewrite of the internal code, to make sure
the code is readable and that the algorithm doesn’t continuously update
the local environment.

-   MAJOR CHANGE: Instead of passing `"D"`, `"ETL"`, `"OTL"`, or `"BTL"`
    to `ssa()`, it is expected to pass `ssa.d()`, `ssa.etl()`,
    `ssa.otl()`, or `ssa.btl()`. This cleans up parameter setting
    clutter in the `ssa()` function.

-   MAJOR CHANGE: Rewrite `ssa.*()` and `ssa.*.diag()` as
    `ssa_step.ssa_*()` and `ssa_step_diag.ssa_*()` S3 functions.

-   MAJOR CHANGE: Do not save the current state in the function
    environment. Instead, simply save it in a local variable.

-   MAJOR CHANGE: Precompile propensity functions instead of evaluating
    them as R code at each iteration.

-   MAJOR CHANGE: Clean up and merge `ssa.run()`, `ssa.terminate()`,
    `ssa.check.args()` and `ssa.check.method()` into `ssa()`.

## References

-   Brown D. and Rothery P. 1993. Models in biology: mathematics,
    statistics, and computing. John Wiley & Sons.
-   Cao Y., Li H., and Petzold L. 2004. Efficient formulation of the
    stochastic simulation algorithm for chemically reacting systems. J.
    Chem. Phys. 121:4059-4067.
    [doi:10.1063/1.1778376](https://doi.org/10.1063/1.1778376)
-   Cao Y., Gillespie D.T., and Petzold L.R. 2006. Efficient step size
    selection for the tau-leaping method. J. Chem. Phys. 124:044109.
    [doi:10.1063/1.2159468](https://doi.org/10.1063/1.2159468)
-   Cao Y., Gillespie D.T., and Petzold L.R. 2007. Adaptive explicit
    tau-leap method with automatic tau selection. J. Chem. Phys.
    126:224101.
    [doi:10.1063/1.2745299](https://doi.org/10.1063/1.2745299)
-   Chatterjee A., Vlachos D.G., and Katsoulakis M.A. 2005. Binomial
    distribution based tau-leap accelerated stochastic simulation. J.
    Chem. Phys. 122:024112.
    [doi:10.1063/1.1833357](https://doi.org/10.1063/1.1833357)
-   Gillespie D.T. 1977. Exact stochastic simulation of coupled chemical
    reactions. J. Phys. Chem. 81:2340.
    [doi:10.1021/j100540a008](https://doi.org/10.1021/j100540a008)
-   Gillespie D.T. 2001. Approximate accelerated stochastic simulation
    of chemically reacting systems. J. Chem. Phys. 115:1716-1733.
    [doi:10.1063/1.1378322](https://doi.org/10.1063/1.1378322)
-   Gillespie D.T. 2007. Stochastic simulation of chemical kinetics.
    Annu. Rev. Chem. 58:35
    [doi:10.1146/annurev.physchem.58.032806.104637](https://doi.org/10.1146/annurev.physchem.58.032806.104637)
-   Kot M. 2001. Elements of mathematical ecology. Cambridge University
    Press.
    [doi:10.1017/CBO9780511608520](https://doi.org/10.1017/CBO9780511608520)
-   Pineda-Krch M. 2008. Implementing the stochastic simulation
    algorithm in R. Journal of Statistical Software 25(12): 1-18. [doi:
    10.18637/jss.v025.i12](https://doi.org/10.18637/jss.v025.i12)
-   Pineda-Krch M., Blok H.J., Dieckmann U., and Doebeli M. 2007. A tale
    of two cycles — distinguishing quasi-cycles and limit cycles in
    finite predator-prey populations. Oikos 116:53-64.
    [doi:10.1111/j.2006.0030-1299.14940.x](https://doi.org/10.1111/j.2006.0030-1299.14940.x)
