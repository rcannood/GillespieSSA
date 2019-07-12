#' Explicit tau-leap method (ETL)
#'
#' Explicit tau-leap method implementation of the \acronym{SSA} as described by
#' Gillespie (2001). It is usually called from within [ssa()], but
#' can be invoked directly.
#'
#' @details Performs one time step using the Explicit tau-leap method.
#' Intended to be invoked by [ssa()].
#'
#' @param a vector of evaluated propensity functions.
#' @param nu state-change matrix.
#' @param tau tau-leap.
#'
#' @return A list with two elements:
#' * the time leap (`tau`) and
#' * the realized state change vector (`nu_j`).
#'
#' @references Gillespie (2001)
#'
#' @seealso [GillespieSSA-package], [ssa()]
#'
#' @examples
#' a <- function(parms,x){
#'   b <- parms[1]
#'   d <- parms[2]
#'   K <- parms[3]
#'   N <- x[1]
#'   return(c(b*N , N*b + (b-d)*N/K))
#' }
#' parms <- c(2,1,1000,500)
#' x <- 500
#' nu <- matrix(c(+1, -1),ncol=2)
#' t <- 0
#' for (i in seq(100)) {
#'   out <- ssa.etl(a(parms,x),nu,tau=0.3)
#'   x <- x + out$nu_j
#'   t <- t + 1
#'   cat("t:",t,", x:",x,"\n")
#' }
#'
#' @keywords misc datagen ts
#'
#' @importFrom stats rpois
#'
#' @export
ssa.etl <- function(
  a = stop("missing propensity vector (a)"),
  nu = stop("missing state-change matrix (nu)"),
  tau = stop("missing step size (tau)")
) {
  M <- length(a)
  k <- rpois(M,(a*tau))
  # MPK: Strictly speaking it is not correct to call the realized state-change
  # vector nu_j here since, in leap methods, actually is not reaction specific.
  # In Pineda-Krch (JSS ms) it is refered to as the \bm{\tilde{nu}}
  nu_j <- rowSums(matrix(rep(k,dim(nu)[1]),byrow=TRUE,ncol=M)*nu)

  list(
    tau = tau,
    nu_j = nu_j
  )
}

