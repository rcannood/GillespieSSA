#' Binomial tau-leap method (BTL)
#'
#' Binomial tau-leap method implementation of the \acronym{SSA} as described by
#' Chatterjee et al. (2005). It is usually called from within
#' [ssa()], but can be invoked directly.
#'
#' Performs one time step using the Binomial tau-leap method. Intended to be
#' invoked by [ssa()].
#'
#' @param x state vector.
#' @param a vector of evaluated propensity functions.
#' @param nu state-change matrix.
#' @param f coarse-graining factor (see page 4 in Chatterjee et al. 2005).
#'
#' @return A list with two elements:
#' * the time leap (`tau`) and
#' * the realized state change vector (`nu_j`).
#'
#' @seealso [GillespieSSA-package], [`ssa()`][ssa]
#' @references Chatterjee et al. (2005)
#' @keywords datagen misc ts
#'
#' @importFrom stats rbinom rpois
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
#'   out <- ssa.btl(x, a(parms, x), nu, f = 10)
#'   x <- x + out$nu_j
#'   t <- t + 1
#'   cat("t:",t,", x:",x,"\n")
#' }
#'
#'
#' @export ssa.btl
ssa.btl <- function(
  x = stop("missing state vector (x)"),
  a = stop("missing propensity vector (a)"),
  nu = stop("missing state-change matrix (nu)"),
  f = stop("missing coarse-graining factor (f)")
) {
  coercing <- FALSE

  # Calculate tau
  tau <- f / sum(a)   # Calculate the current tau
  if (tau>1) tau<-1 # Tau cannot be larger than unity!

  M <- length(a)    # Number of reaction channels
  tilde_x <- x
  nu_j <- matrix(rep(0,length(x)))

  # Loop over all reaction channels having propensity fun>0
  for (j in seq(M)[a>0]) {
    if (any(nu[,j]<0)) { # do this if there are limiting reactions
      mask <- nu[,j]<0
      L <- min(floor(tilde_x[mask]/abs(nu[mask,j])))
      if (a[j]*tau>L) {
        p <- 1
        coercing <- TRUE
      }
      else p <- a[j]*tau/L
      k <- rbinom(1,L,p)
    }
    else { # do this if there are no limiting reactions
      k <- rpois(1,(a[j]*tau))
    }

    # Update tilde_x for the current reaction j
    tmp_nu_j <- matrix(rep(k,dim(nu)[1]), byrow=TRUE, ncol=1)*nu[,j]
    tilde_x <- tilde_x + tmp_nu_j

    # Record the current state change (it is returned by ssa.btl)
    nu_j <- nu_j + tmp_nu_j
  } # for()

  # Throw a warning message if p was coerced to unity. Coercing implies too
  # large steps-size due to a high coarse-graining factor (f)
  if (coercing) warning("coerced p to unity - consider lowering f")

  list(
    tau = tau,
    nu_j = nu_j
  )
}

