#' Binomial tau-leap method (BTL) for nu-diagonalized systems
#'
#' Binomial tau-leap method for nu-diagonalized systems
#'
#' Performs one time step using the Binomial tau-leap method. It is usually
#' called from within [ssa()], but can be invoked directly, see
#' [ssa.btl()] for Examples.
#'
#' @param x state vector.
#' @param a vector of evaluated propensity functions.
#' @param nu_tile state-change matrix.
#' @param f coarse-graining factor (see page 4 in Chatterjee et al. 2005).
#'
#' @return A list with two elements:
#' * the time leap (`tau`) and
#' * the realized state change vector (`nu_j`).
#'
#' @importFrom stats rbinom rpois
#'
#' @seealso [ssa.btl()],
#' @keywords misc datagen ts
ssa.btl.diag <- function(
  x,
  a,
  nu_tile,
  f
) {
  coercing <- FALSE

  # Calculate tau
  tau <- f/sum(a)   # Calculate the current tau
  if (tau>1) tau<-1 # Tau cannot be larger than unity!

  M <- dim(nu_tile)[2] # Number of reaction channels per nu-tile
  N <- dim(nu_tile)[1] # Number of states per nu-tile
  MU <- length(a)      # Toto nr of reaction channels
  U <- MU/M            # Nr of tilings
  tilde_x <- x
  nu_j <- rep(0,(N*U))

  # Identify potential limiting reactions
  #mask <- apply(nu_tile,2,function(x) any(x<0))

  # Loop over all reaction channels having a non-zero (>0) propensity fun
  for (j in seq(U*M)[a>0]) {
    f <- ceiling((j/M)-1)
    jp <- j-f*M  # Intra-patch reaction channel index (j->jp)
    x1 <- 1+f*N
    x2 <- 1+f*N+(N-1)
    if (any(nu_tile[,jp]<0)) {  # Do this if there are limiting reactions
      mask <- nu_tile[,jp]<0    # Which species has the limiting reaction
      tilde_xt <- tilde_x[x1:x2]
      L <- min(floor(tilde_xt[mask]/abs(nu_tile[mask,jp])))
      if (a[j]*tau>L) {
        p <- 1
        coercing <- TRUE
      }
      else {
        p <- a[j]*tau/L
        k <- rbinom(1,L,p)
      }
    } else { # do this if there are no limiting reactions
      k <- rpois(1,(a[j]*tau))
    }

    # Update tilde_x for the current reaction j
    tmp_nu_j <- rep(k,dim(nu_tile)[1])*nu_tile[,jp]
    tilde_x[x1:x2] <- tilde_x[x1:x2] + tmp_nu_j

    # Record the current state change
    nu_j[x1:x2] <- nu_j[x1:x2] + tmp_nu_j
    if (any(is.na(nu_j))) browser() # MPK: Just in case!
  } # for()

  # Throw a warning message if p was coerced to unity. Coercing implies too
  # large steps-size due to a high coarse-graning factor (f)
  if(coercing) warning("coerced p to unity - consider lowering f")

  list(
    tau = tau,
    nu_j = nu_j
  )
}

