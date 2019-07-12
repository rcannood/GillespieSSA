#' Direct method (D)
#'
#' Direct method implementation of the \acronym{SSA} as described by Gillespie (1977).
#' It is usually called from within [ssa()], but can be invoked directly.
#'
#' @param a vector of evaluated propensity functions.
#' @param nu state-change matrix.
#'
#' @return A list with two elements:
#' * the time leap (`tau`) and
#' * the realized state change vector (`nu_j`).
#'
#' @seealso [GillespieSSA-package], [ssa()]
#' @references Gillespie (1977)
#' @keywords misc datagen ts
#'
#' @importFrom stats runif
#'
#' @examples
#' ## Logistic growth model
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
#'   out <- ssa.d(a(parms,x),nu)
#'   x <- x + out$nu_j
#'   t <- t + 1
#'   cat("t:",t,", x:",x,"\n")
#' }
#'
#' @export
ssa.d <- function(
  a = stop("missing propensity vector (a)"),
  nu = stop("missing state-change matrix (nu)")
) {
  j <- sample.int(length(a), size = 1, prob = a)
  nu_j <- nu[,j]
  tau <- -log(runif(1)) / sum(a)

  list(
    tau = tau,
    nu_j = nu_j
  )
}
