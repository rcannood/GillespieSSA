#' Direct method nu-diagonalization mapping
#'
#' Auxiliary function for `ssa.d.diag` performing virtual mapping of nu-diagonalized systems.
#'
#' @param a vector of evaluated propensity functions.
#' @param nu state-change matrix.
#' @param j Reaction index to map
#'
#' @return The virtual realized state change vector (`nu_j`).
#'
#' @seealso [ssa.d.diag()], [ssa.d()]
#'
#' @keywords misc datagen ts
ssa.nutiling <- function(a,nu,j) {
  M  <- dim(nu)[2]        # Number of reaction channels in nu-tile
  N  <- dim(nu)[1]        # Number of states in nu tile
  U  <- length(a)/M       # Number of tessallations of nu tile
  f  <- ceiling((j/M)-1)  # Frameshift factor
  jp <- j-f*M             # Relative reaction channel index
  nu_jp <- nu[,jp]
  nu_j <- c(rep(0,f*N),   # Leading zeros
            nu_jp,        # Relative state-change matrix
            rep(0,(U*N-(f*N+N)))) # Lagging zeros
  return(nu_j)
}
