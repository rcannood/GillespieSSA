#' Explicit tau-leap method (ETL) for nu-diagonalized systems
#'
#' Explicit tau-leap method for nu-diagonalized systems.
#'
#' @details
#' Performs one time step using the Explicit tau-leap method. It is usually called from
#' within [ssa()], but can be invoked directly, see [ssa.etl()]
#' for Examples.
#'
#' @param a vector of evaluated propensity functions.
#' @param nu_tile state-change-matrix.
#' @param tau tau-leap.
#'
#' @return A list with two elements:
#' * the time leap (`tau`) and
#' * the realized state change vector (`nu_j`).
#'
#' @seealso [ssa.etl()]
#'
#' @keywords misc datagen ts
#'
#' @importFrom stats rpois
ssa.etl.diag <- function(a, nu_tile, tau) {
  MU <- length(a)          # Toto nr of reaction channels
  k  <- rpois(MU, (a*tau)) # Nr of firings per channel
  M  <- dim(nu_tile)[2]    # Nr of reaction channel per patch (nu_tile)
  U  <- MU / M             # Nr of tilings
  nu_j <- NULL

  for (f in (seq_len(U)-1)) {
    nu_j_ <- rowSums(matrix(rep(k[1:M+f*M],dim(nu_tile)[1]),byrow=TRUE,ncol=M)*nu_tile)
    nu_j <- c(nu_j, nu_j_)
  }

  list(
    tau = tau,
    nu_j = nu_j
  )
}

