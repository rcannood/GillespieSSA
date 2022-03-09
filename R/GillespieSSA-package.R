#' Gillespie Stochastic Simulation Algorithm package
#'
#' @name GillespieSSA-package
#' @aliases GillespieSSA-package GillespieSSA
#' @docType package
#'
#' @description Package description and overview of basic SSA theory
#'
#' \pkg{GillespieSSA} is a versatile and extensible framework for stochastic
#' simulation in R and provides a simple interface to a number of Monte Carlo
#' implementations of the stochastic simulation algorithm (\acronym{SSA}). The
#' methods currently implemented are: the Direct method, Explicit
#' tau-leaping (\acronym{ETL}), Binomial tau-leaping (\acronym{BTL}), and
#' Optimized tau-leaping (\acronym{OTL}). The package also provides a library
#' of ecological, epidemiological, and evolutionary continuous-time (demo)
#' models that can easily be customized and extended. Currently the following
#' models are included, Decaying-Dimerization Reaction Set, Linear Chain
#' System, single-species logistic growth model, Lotka predator-prey model,
#' Rosenzweig-MacArthur predator-prey model, Kermack-McKendrick \acronym{SIR}
#' model, and a metapopulation \acronym{SIRS} model.
#'
#' @concept Direct method
#' @concept Explicit tau-leap method
#' @concept Binomial tau-leap method
#' @concept Optimized tau-leap method
#' @concept Poisson
#' @concept stochastic simulation algorithm
#' @concept Gillespie
#' @concept ecology
#' @concept epidemiology
#' @concept evolution
#' @concept biology
#'
#' @section The stochastic simulation algorithm:
#' The stochastic simulation algorithm (\acronym{SSA}) is a procedure
#' for constructing simulated trajectories of finite populations in continuous time.
#' If \eqn{X_i(t)} is the number of individuals in population \eqn{i}
#' (\eqn{i=1,\ldots,N}{i=1,...,N}) at time \eqn{t} the \acronym{SSA} estimates
#' the state vector \eqn{ \mathbf{X}(t) \equiv (X_1(t),\ldots,X_N(t)) }{ X(t) =
#' (X_1(t),...,X_N(t))}, given that the system initially (at time \eqn{t_0})
#' was in state \eqn{\mathbf{X}(t_0)=\mathbf{x_0}}{X(t_0)=x_0}. Reactions,
#' single instantaneous events changing at least one of the populations (e.g.
#' birth, death, movement, collision, predation, infection, etc), cause the
#' state of the system to change over time. The \acronym{SSA} procedure samples
#' the time \eqn{\tau}{tau} to the next reaction \eqn{R_j}
#' (\eqn{j=1,\ldots,M}{j=1,...,M}) and updates the system state
#' \eqn{\mathbf{X}(t)}{X(t)} accordingly. Each reaction \eqn{R_j} is
#' characterized mathematically by two quantities; its state-change vector
#' \eqn{\bm{\nu}_j \equiv ( \nu_{1j},\ldots,\nu_{Nj} )}{nu_j =
#' (nu_1j,...,nu_Nj)}, where \eqn{ \nu_{ij} }{nu_ij} is the change in the
#' number of individuals in population \eqn{i} caused by one reaction of type
#' \eqn{j} and its propensity function \eqn{a_j(\mathbf{x})}{a_j(x)}, where
#' \eqn{a_j(\mathbf{x})dt}{a_j(x)dt} is the probability that a particular
#' reaction \eqn{j} will occur in the next infinitesimal time interval
#' \eqn{\left[t,t+dt\right]}{[t,t+dt]}.
#'
#' @section SSA implementations:
#' There are numerous exact Monte Carlo procedures implementing the \acronym{SSA}.
#' Perhaps the simplest is the Direct method of Gillespie (1977. The Direct method is
#' an exact continuous-time numerical realization of the corresponding stochastic
#' time-evolution equation. Because the Direct method simulates one reaction at a time
#' it is often, however, computationally too slow for practical applications.
#'
#' Approximate implementations of the \acronym{SSA} sacrifices exactness for large
#' improvements in computational efficiency. The most common technique used is
#' tau-leaping where reaction-bundles are attempted in coarse-grained time increments
#' \eqn{\tau}{tau}. Speed-ups of several orders of magnitude compared to the Direct
#' method are common. Tau-leaping must be used with care, however, as it is not as
#' foolproof as the Direct method.
#'
#' @section Example models:
#' Individual demo models can be run by issuing `demo(<model name>)`,
#' alternatively all of the demo models can be run using `demo(GillespieSSA)`.
#' The following example models are available:
#'
#' \tabular{l}{
#'   Decaying-Dimerization Reaction Set (Gillespie, 2001)\cr
#'   `vignette("decaying_dimer", package = "GillespieSSA")` \cr
#'   \cr
#'   SIRS metapopulation model (Pineda-Krch, 2008)\cr
#'   `vignette("epi_chain", package = "GillespieSSA")` \cr
#'   \cr
#'   Linear Chain System (Cao et al., 2004)\cr
#'   `vignette("linear_chain", package = "GillespieSSA")` \cr
#'   \cr
#'   Pearl-Verhulst Logistic growth model (Kot, 2001, Pineda-Krch, 2008) \cr
#'   `vignette("logistic_growth", package = "GillespieSSA")` \cr
#'   \cr
#'   Lotka predator-prey model (Gillespie, 1977; Kot, 2001)\cr
#'   `vignette("lotka_predator_prey", package = "GillespieSSA")` \cr
#'   \cr
#'   Radioactive decay model (Gillespie, 1977)\cr
#'   `vignette("radioactive_decay", package = "GillespieSSA")` \cr
#'   \cr
#'   Rosenzweig-MacArthur predator-prey model (Pineda-Krch et al., 2007, Pineda-Krch, 2008) \cr
#'   `vignette("rm_predator_prey", package = "GillespieSSA")` \cr
#'   \cr
#'   Kermack-McKendrick SIR model (Brown & Rothery, 1993)\cr
#'   `vignette("sir", package = "GillespieSSA")` \cr
#'   \cr
#' }
#'
#' @references
#'  * Brown D. and Rothery P. 1993. Models in biology: mathematics, statistics, and computing. John Wiley & Sons.
#'  * Cao Y., Li H., and Petzold L. 2004. Efficient formulation of the stochastic simulation algorithm for chemically reacting systems. J. Chem. Phys. 121:4059-4067. \doi{10.1063/1.1778376}
#'  * Cao Y., Gillespie D.T., and Petzold L.R. 2006. Efficient step size selection for the tau-leaping method. J. Chem. Phys. 124:044109. \doi{10.1063/1.2159468}
#'  * Cao Y., Gillespie D.T., and Petzold L.R. 2007. Adaptive explicit tau-leap method with automatic tau selection. J. Chem. Phys. 126:224101. \doi{10.1063/1.2745299}
#'  * Chatterjee A., Vlachos D.G., and Katsoulakis M.A. 2005. Binomial distribution based tau-leap accelerated stochastic simulation. J. Chem. Phys. 122:024112. \doi{10.1063/1.1833357}
#'  * Gillespie D.T. 1977. Exact stochastic simulation of coupled chemical reactions. J. Phys. Chem. 81:2340. \doi{10.1021/j100540a008}
#'  * Gillespie D.T. 2001. Approximate accelerated stochastic simulation of chemically reacting systems. J. Chem. Phys. 115:1716-1733. \doi{10.1063/1.1378322}
#'  * Gillespie D.T. 2007. Stochastic simulation of chemical kinetics. Annu. Rev. Chem. 58:35 \doi{10.1146/annurev.physchem.58.032806.104637}
#'  * Kot M. 2001. Elements of mathematical ecology. Cambridge University Press. \doi{10.1017/CBO9780511608520}
#'  * Pineda-Krch M. 2008. Implementing the stochastic simulation algorithm in R. Journal of Statistical Software 25(12): 1-18. \doi{10.18637/jss.v025.i12}
#'  * Pineda-Krch M., Blok H.J., Dieckmann U., and Doebeli M. 2007. A tale of two cycles --- distinguishing quasi-cycles and limit cycles in finite predator-prey populations. Oikos 116:53-64. \doi{10.1111/j.2006.0030-1299.14940.x}
#'
#' @section Acknowledgements:
#' * Heinrich zu Dohna for many caffein induced discussions on the package and reference manual, and for providing comments on the vignette documentation.
#' * Ben Bolker for comments on the initial release of the package and for providing a hint for how to more elegantly handle model parameters as arguments to the ssa() function.
#' * Josh Obrien for copy editing and feedback on the JSS manuscript.
#' * Thomas Petzoldt for comments on the package, the JSS manuscript and for preparing version 0.5-4.
#' * Three anonymous referees whose comments substantially improved some of the functionality.
#'
#' @seealso [ssa()], [ssa.d()], [ssa.etl()], [ssa.btl()], [ssa.otl()], [ssa.plot()]
#'
#' @keywords package distribution
NULL

