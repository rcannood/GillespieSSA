#' Binomial tau-leap method (BTL)
#'
#' Binomial tau-leap method implementation of the \acronym{SSA} as described by
#' Chatterjee et al. (2005). Should be passed as `method` argument for `ssa()`.
#'
#' Performs one time step using the Binomial tau-leap method. Intended to be
#' invoked by [ssa()].
#'
#' @param f coarse-graining factor (see page 4 in Chatterjee et al. 2005).
#'
#' @seealso [GillespieSSA-package], [ssa()]
#' @references Chatterjee et al. (2005)
#' @keywords datagen misc ts
#'
#' @importFrom stats rbinom rpois
#'
#' @examples
#' ssa.btl(f = 40)
#'
#' @export
ssa.btl <- function(f = 10) {
  if (!is.numeric(f)) stop("'f' is not numeric")
  if (f <= 1) stop("f has to be >1")

  ssa_method(
    "BTL",
    list(
      f = f
    )
  )
}


ssa_step.ssa_BTL <- function (method, x, a, nu, method_state) {
  f <- method$params$f
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
    tmp_nu_j <- matrix(rep(k,nrow(nu)), byrow=TRUE, ncol=1)*nu[,j]
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

ssa_step_diag.ssa_BTL <- function (method, x, a, nu, method_state) {
  f <- method$params$f
  coercing <- FALSE

  # Calculate tau
  tau <- f/sum(a)   # Calculate the current tau
  if (tau>1) tau<-1 # Tau cannot be larger than unity!

  M <- ncol(nu) # Number of reaction channels per nu-tile
  N <- nrow(nu) # Number of states per nu-tile
  MU <- length(a)      # Toto nr of reaction channels
  U <- MU/M            # Nr of tilings
  tilde_x <- x
  nu_j <- rep(0,(N*U))

  # Identify potential limiting reactions
  #mask <- apply(nu,2,function(x) any(x<0))

  # Loop over all reaction channels having a non-zero (>0) propensity fun
  for (j in seq(U*M)[a>0]) {
    f <- ceiling((j/M)-1)
    jp <- j-f*M  # Intra-patch reaction channel index (j->jp)
    x1 <- 1+f*N
    x2 <- 1+f*N+(N-1)
    if (any(nu[,jp]<0)) {  # Do this if there are limiting reactions
      mask <- nu[,jp]<0    # Which species has the limiting reaction
      tilde_xt <- tilde_x[x1:x2]
      L <- min(floor(tilde_xt[mask]/abs(nu[mask,jp])))
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
    tmp_nu_j <- rep(k,nrow(nu))*nu[,jp]
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
