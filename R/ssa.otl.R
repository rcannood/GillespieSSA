#' Optimized tau-leap method (OTL)
#'
#' Optimized tau-leap method implementation of the \acronym{SSA} as described
#' by Cao et al. (2006). It is usually called from within [ssa()],
#' but can be invoked directly.
#'
#' Performs one time step using the Explicit tau-leap method. Intended to be
#' invoked by [ssa()].
#'
#' @param x state vector.
#' @param a vector of evaluated propensity functions.
#' @param nu state-change matrix.
#' @param hor highest order reaction vector (one entry per species in `x`)
#' @param nc number of critical reactions threshold parameter.
#' @param epsilon error control parameter.
#' @param dtf Direct method threshold factor for temporarily suspending the
#' `OTL` method.
#' @param nd number of Direct method steps to perform during an `OTL`
#' suspension.
#' @return A list with three elements:
#' * the time leap (`tau`),
#' * the realized state change vector (`nu_j`) and
#' * a boolean value (`suspendedTauLeapMethod`) indicating if the simulation should revert to the Direct method for `nd` time steps.
#'
#' @note Third order-reactions (\eqn{S_1+S_2+S_3 \rightarrow \ldots}{S_1 + S_2
#' + S_3 ---> ...}) are not supported currently since they are approximations
#' to sets of coupled first- and second-order reactions). See Cao et al. (2006)
#' for more details.
#' @seealso [GillespieSSA-package], [ssa()]
#' @references Cao et al. (2006)
#' @keywords misc datagen ts
#' @examples
#'
#' a = function(parms,x){
#'  b <- parms[1]
#'  d <- parms[2]
#'  K <- parms[3]
#'  N <- x[1]
#'  return(c(b*N , N*b + (b-d)*N/K))
#' }
#' parms <- c(2,1,1000,500)
#' x <- 500
#' nu <- matrix(c(+1, -1),ncol=2)
#' t <- 0
#' for (i in seq(100)) {
#'   out <- ssa.otl(x,a(parms,x),nu,hor=1,nc=10,epsilon=0.03,dtf=10,nd=100)
#'   x <- x + out$nu_j
#'   t <- t + 1
#'   cat("t:",t,", x:",x,"\n")
#' }
#' @export
#' @importFrom stats runif rbinom
ssa.otl <- function(
  x = stop("missing state vector (x)"),
  a = stop("missing propensity vector (a)"),
  nu = stop("missing state-change matrix (nu)"),
  hor = stop("missing highest order reaction vector (hot)"),
  nc = stop("missing critical reactions threshold parameter (nc)"),
  epsilon = stop("missing error control parameter"),
  dtf = stop("missing direct method threshold factor (dtf)"),
  nd = stop("missing OTL suspension duration parameter (nd)")
) {
  verbose <- FALSE
  if (verbose) cat("Starting OTL...\n")

  # 1. Identify the current critical reactions
  # Calculate the minimum number of firings for reaction before one of it's
  # reactants becomes extinct (L). The 'L' notation is from Eq. 10 in Cao
  # et al. (2006), J. Chem. Phys. 124:044109.
  tmp_nu <- nu
  tmp_nu[nu>=0] <- NA  # We are only interested in negative state changes

  # Turn off warning temporarily because min() throws a warning if it tries to
  # evaluate only 'NA's, which it will for reaction channels that have no
  # negative entries. Despite the warning the end result is correct, i.e. the
  # number of firings for such a channel becomes Inf.
  options(warn = -1) # warnings off
  L <- apply(floor(x/abs(tmp_nu)),2,min,na.rm=TRUE)
  options(warn = 0) # warnings on
  Jncr <- L >= nc                     # Indices of the non-critical reactions

  # 2. Compute the first candidate time leap, tau1
  if (sum(Jncr) == 0) {
    tau1 <- Inf                       # It is simple if there are no critical reactions present
  } else {                            # It is complicate if there are critical reactions present
    Irs <- apply((nu != 0),1,any)     # Subset the reactant species
    g <- rep(NA,length(x))
    g[hor==1]  <- 1                   # First-order reactions (S1->...)
    g[hor==2]  <- 2                   # Interspecific second-order reaction, first type (S1+S2->...)
    g[hor==22] <- (2+1/(x[hor==2]-1)) # Intraspecific second-order reaction (S1+S1->...)

    # Define mu ($\hat{\mu$}_i(\matnbf{x}) in Eq. 32a)
    tmp_nu <- nu[apply(nu, 1, function(x) any(x != 0)), , drop = FALSE] # Remove non-reacting species from nu
    tmp <- tmp_nu[,Jncr]*a[Jncr]
    if (is.matrix(tmp)) mu <- rowSums(tmp)
    else mu <- tmp

    # Define sigma ($\hat{\sigma}^2_i(\mathbf{x})$ in Eq. 32b)
    tmp <- tmp_nu[,Jncr]^2*a[Jncr]
    if (is.matrix(tmp)) sigmaSq <- rowSums(tmp)
    else sigmaSq <- tmp

    # Calculate tau1 (Eq. 33). If there are no noncritical reactions (Jncr only
    # has FALSE elements) tau1<-Inf (see #2 in paper)
    leftTerm  <- max(epsilon*x/g,1) / abs(mu)
    rightTerm <- max(epsilon*x/g,1)^2 / abs(sigmaSq)
    tau1      <- min(leftTerm[Irs],rightTerm[Irs])
    if (is.infinite(tau1)) cat("tau1=Inf\n") # DEBUG
    if (is.na(tau1)) browser() # DEBUG
  } # if (sum(Jncr) == 0)

  # We need to the 'while' loop with it's constructs so that we can recaulate
  # tau if we end up with negative population sizes (see step #6 in paper, page 4)
  calculateTau <- TRUE
  while (calculateTau) {
    if (verbose) cat("Calculating tau...\n")

    # 3. If tau1 is "too small" return to stochRxn() and execute a number of direct method steps.
    if (verbose) cat("tau1=",tau1,",(",dtf,"*1/sum(a))=",(dtf*1/sum(a)),", a=",a,"\n",sep=" ")
    if (tau1 < (dtf*1/sum(a))) {
      if (verbose) cat("*** Suspending tauLeap method (tau1=",tau1,", (",dtf,"/sum(a)=",(dtf*1/sum(a)),")...\n")
      return(list(tau=NA, nu_j=NA, suspendedTauLeapMethod=nd))
    }

    # 4. Compute the second candidate time leap from the set of critical reactions, tau2. If there are no critical reactions tau2=Inf
    tau2 <- -log(runif(1))/sum(a[!Jncr])

    # 5. Select the actual tau from the two candidate tau (the smaller of the two)
    # and determine the number of firings each reaction will have
    if (verbose) cat("tau1=",tau1,", tau2=",tau2," -> ")
    if (tau1 < tau2) {                                           # Step 5a
      if (verbose) cat("Selecting tau1...\n")
      tau <- tau1
      k <- as.numeric(!Jncr)                                     # Sets all critical reactions to one firings and non-critical to zero firings
      k[k==0] <- rpois(sum(Jncr),(a[Jncr]*tau))        # Sets the number of firings for non-critical reactions
    } else {                                                     # Step 5b
      if (verbose) cat("Selecting tau2...\n")
      tau <- tau2
      jc <- sample(seq(dim(nu)[2]),size=1,prob=(a/sum(a[!Jncr]))) # Pick one of the critical reactions that will fire once
      k <- rep(0,dim(nu)[2])                                     # Setting up an empty vector
      k[jc] <- 1                                                 # Add the selected critical reaction that is firing
      k[Jncr %in% TRUE] <- rpois(sum(Jncr),(a*tau))              # The number of firings of non-critical reactions is drawn from a Poisson distribution
    } # if (tau1 < tau2)

    # 6. Update the state vector and check for negative elements. If negative elements are found reduce
    # tau1 by half and return to step 3
    nu_j <- rowSums(matrix(rep(k,dim(nu)[1]),byrow=TRUE,ncol=length(a))*nu)
    if (verbose) cat("x=",x,", nu_j=",nu_j,"\n")
    if (any((x+nu_j)<0)) {
      tau1 <- tau1/2
      calculateTau <- TRUE
      if (verbose) cat("Detected negative elements in 'x'...\n")
    } else {
      calculateTau <- FALSE
    }
    if (verbose) cat("tau=",tau,"\n")
  } # while
  if (verbose) cat("Done with optimizedTauLeap()...\n")
  return(list(tau=tau, nu_j=nu_j, suspendedTauLeapMethod=FALSE))
}

