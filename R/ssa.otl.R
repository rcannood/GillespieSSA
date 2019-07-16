#' Optimized tau-leap method (OTL)
#'
#' Optimized tau-leap method implementation of the \acronym{SSA} as described
#' by Cao et al. (2006). Should be passed as `method` argument for `ssa()`.
#'
#' @param hor Highest order reaction vector. There must be one entry per species in `x`.
#'   Must be one of `1`: first-order, `2`: second-order or `22`: homo-dimer.
#'   If `hor` is `NA`, defaults are all second-order.
#' @param nc number of critical reactions threshold parameter.
#' @param epsilon error control parameter.
#' @param dtf Direct method threshold factor for temporarily suspending the
#' `OTL` method.
#' @param nd number of Direct method steps to perform during an `OTL`
#' suspension.
#'
#' @note Third order-reactions (\eqn{S_1+S_2+S_3 \rightarrow \ldots}{S_1 + S_2
#' + S_3 ---> ...}) are not supported currently since they are approximations
#' to sets of coupled first- and second-order reactions). See Cao et al. (2006)
#' for more details.
#'
#' @seealso [GillespieSSA-package], [ssa()]
#' @references Cao et al. (2006)
#' @keywords misc datagen ts
#' @examples
#' ssa.otl(
#'   hor = 1,
#'   nc = 10,
#'   epsilon = .03,
#'   dtf = 10,
#'   nd = 100
#' )
#' @export
#' @importFrom stats runif rbinom
ssa.otl <- function(
  epsilon = 0.03,
  nc = 10,
  hor = NA_real_,
  dtf = 10,
  nd = 100
) {
  if (!is.numeric(epsilon)) stop("'epsilon' is not numeric")
  if (!is.numeric(nc)) stop("'nc' is not numeric")
  if (!is.numeric(hor)) stop("'hor' is not numeric")
  if (!is.numeric(dtf)) stop("'dtf' is not numeric")
  if (!is.numeric(nd)) stop("'nd' is not numeric")

  ssa_method(
    "OTL",
    list(
      epsilon = epsilon,
      nc = nc,
      hor = hor,
      dtf = dtf,
      nd = nd
    )
  )
}

ssa_initialise_state.ssa_OTL <- function(method) {
  list(
    suspensions = 0,
    total_suspensions = 0
  )
}

ssa_validate_parameters.ssa_OTL <- function(method, x, a, nu) {
  params <- method$params

  # Check if hor vector is defined as NA, in which the conservative
  # default value of 2 is used for each species.
  hor <- params$hor

  if (any(is.na(hor))) {
    hor <- rep(2, length(x)) # Undefined hor - use default values
  }

  if (length(hor) != length(x)) {
    stop("length of hor vector is different from length of 'x0'")
  }

  if (any(hor != 1 & hor != 2 & hor != 22)) {
    stop("wrong value(s) in hor vector (can only be 1, 2, or 22)")
  }

  params$hor <- hor

  params
}

ssa_step.ssa_OTL <- function (method, x, a, nu, method_state) {
  if (method_state$suspensions > 0) {
    method_state$suspensions <- method_state$suspensions - 1
    method_state$total_suspensions <- method_state$total_suspensions + 1
    out <- ssa_step.ssa_D(method, x, a, nu, method_state)
    out$method_state <- method_state
    return(out)
  } else {
    epsilon <- method$params$epsilon
    nc <- method$params$nc
    hor <- method$params$hor
    dtf <- method$params$dtf
    nd <- method$params$nd

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
        method_state$suspensions <- nd
        return(ssa_step.ssa_OTL(method, x, a, nu, method_state))
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
        jc <- sample(seq(ncol(nu)),size=1,prob=(a/sum(a[!Jncr]))) # Pick one of the critical reactions that will fire once
        k <- rep(0,ncol(nu))                                     # Setting up an empty vector
        k[jc] <- 1                                                 # Add the selected critical reaction that is firing
        k[Jncr %in% TRUE] <- rpois(sum(Jncr),(a*tau))              # The number of firings of non-critical reactions is drawn from a Poisson distribution
      } # if (tau1 < tau2)

      # 6. Update the state vector and check for negative elements. If negative elements are found reduce
      # tau1 by half and return to step 3
      nu_j <- rowSums(matrix(rep(k,nrow(nu)),byrow=TRUE,ncol=length(a))*nu)
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

    list(
      tau = tau,
      nu_j = nu_j,
      method_state
    )
  }
}

ssa_step_diag.ssa_OTL <- function (method, x, a, nu, method_state) {
  if (method_state$suspensions > 0) {
    method_state$suspensions <- method_state$suspensions - 1
    method_state$total_suspensions <- method_state$total_suspensions + 1
    out <- ssa_step_diag.ssa_D(method, x, a, nu, method_state)
    out$method_state <- method_state
    return(out)
  } else {
    epsilon <- method$params$epsilon
    nc <- method$params$nc
    hor <- method$params$hor
    dtf <- method$params$dtf
    nd <- method$params$nd

    # 1. Identify the current critical reactions
    # Calculate the minimum number of firings for reaction before one of it's
    # reactants becomes extinct (L). The 'L' notation is from Eq. 10 in Cao
    # et al. (2006). We have to turn off warning messages temporarily because
    # 'min()' throws a warning if it tries to evaluate only 'NA's, which it will
    # for reaction channels that have no negative entries. Despite the warning the
    # end result is correct, i.e. the number of firings for such a channel becomes
    # 'Inf'.
    N <- nrow(nu)  # Nr of states per tile
    M <- ncol(nu)  # Nr of reaction channels per tile
    U <- length(a)/M      # Nr of tilings
    nu_negative <- nu
    nu_negative[nu>=0] <- NA  # We are only interested in negative state changes
    L <- NULL
    options(warn = -1)
    for(f in (seq(U)-1))
      L <- c(L, apply(floor(x[1:N+f*N]/abs(nu_negative)),2,min,na.rm=TRUE))
    options(warn = 0)
    Jncr <- L >= nc  # Indices of the non-critical reactions

    # 2. Compute the first candidate time leap, tau1
    if (sum(Jncr) == 0) {
      tau1 <- Inf                       # No critical reactions present
    } else {                            # Critical reactions are present
      Irs <- rep(apply((nu != 0),1,any),U) # Subset the reactant species
      g <- rep(NA,length(x))
      g[hor==1]  <- 1                   # First-order reactions (S1->...)
      g[hor==2]  <- 2                   # Interspecific 2nd order reaction, first type (S1+S2->...)
      g[hor==22] <- (2+1/(x[hor==2]-1)) # Intraspecific 2nd order reaction (S1+S1->...)

      # Define mu ($\hat{\mu$}_i(\matnbf{x}) in Eq. 32a)
      # Define sigma ($\hat{\sigma}^2_i(\mathbf{x})$ in Eq. 32b)
      nu_reacting <- nu[apply(nu, 1, function(x) any(x != 0)), , drop = FALSE] # Remove non-reacting species from nu
      mu <- NULL
      sigmaSq <- NULL
      for (f in (seq(U)-1)) {
        a_current_frame <- a[1:M+f*M]
        Jncr_current_frame <- Jncr[1:M+f*M]
        mu_tmp <- nu_reacting[,Jncr_current_frame]*a_current_frame[Jncr_current_frame]
        sigmaSq_tmp <- nu_reacting[,Jncr_current_frame]^2*a_current_frame[Jncr_current_frame]
        if (is.matrix(mu_tmp)) mu <- c(mu, rowSums(mu_tmp))
        else mu <- c(mu, mu_tmp)
        if (is.matrix(sigmaSq_tmp)) sigmaSq <- c(sigmaSq, rowSums(sigmaSq_tmp))
        else sigmaSq <- c(sigmaSq, sigmaSq_tmp)
      }

      # Calculate tau1 (Eq. 33). If there are no noncritical reactions (Jncr only
      # has FALSE elements) tau1<-Inf (see #2 in paper)
      leftTerm  <- max(epsilon*x/g,1) / abs(mu)
      rightTerm <- max(epsilon*x/g,1)^2 / abs(sigmaSq)
      tau1      <- min(leftTerm[Irs],rightTerm[Irs])
      if (is.infinite(tau1)) cat("tau1=Inf\n") # DEBUG
      if (is.na(tau1)) browser() # DEBUG
    } # if (sum(Jncr) == 0)

    # We need to the 'while' loop with it's constructs so that we can recaulate
    # tau if we end up with negative population sizes (see #6 in paper, page 4)
    calculateTau <- TRUE
    while (calculateTau) {

      # 3. If tau1 is "too small" return to stochRxn() and execute a number of
      # direct method steps.
      if (tau1 < (dtf*1/sum(a))) {
        method_state$suspensions <- nd
        return(ssa_step_diag.ssa_OTL(method, x, a, nu, method_state))
      }

      # 4. Compute the second candidate time leap from the set of critical
      # reactions, tau2. If there are no critical reactions tau2=Inf
      tau2 <- -log(runif(1))/sum(a[!Jncr])

      # 5. Select the actual tau from the two candidate tau (the smaller of the
      # two) and determine the number of firings each reaction will have
      if (tau1 < tau2) {                           # Step 5a
        tau <- tau1
        k <- as.numeric(!Jncr)                     # Sets all critical reactions to one firings and non-critical to zero firings
        lambda <- (a[Jncr]*tau) # Fudge for negative probabilities
        lambda[lambda<0] <- 0
        if (any(lambda<0)) {cat("1\n"); browser()}
        k[k==0] <- rpois(sum(Jncr),lambda)  # Sets the number of firings for non-critical reactions
      } else {                                     # Step 5b
        tau <- tau2
        pr <- (a/sum(a[!Jncr])) # Fudge for negative probabilities
        pr[pr<0] <- 0
        if (any(pr<0)) { cat("3\n"); browser()}
        jc <- sample(seq(M*U),size=1,prob=pr) # Pick one of the critical reactions that will fire once
        k <- rep(0,(M*U))                          # Setting up an empty vector
        k[jc] <- 1                                 # Add the selected critical reaction that is firing
        lambda <- (a*tau) # Fudge for negative probabilities
        lambda[lambda<0] <- 0
        if (any(lambda<0)) {cat("2\n"); browser()}
        k[Jncr %in% TRUE] <- rpois(sum(Jncr),(a*tau))  # The number of firings of non-critical reactions is drawn from a Poisson distribution
      }

      # 6. Update the state vector and check for negative elements. If negative
      # elements are found reduce tau1 by half and return to step 3

      # Update the state-change vector by nu-tiling
      nu_j <- NULL
      for(f in (seq(U)-1))
        nu_j <- c(nu_j, rowSums(matrix(rep(k[1:M+f*M],nrow(nu)),byrow=TRUE,ncol=M)*nu))

      if (any((x+nu_j)<0)) {
        tau1 <- tau1/2
        calculateTau <- TRUE
      } else {
        calculateTau <- FALSE
      }
    } # while

    list(
      tau = tau,
      nu_j = nu_j,
      method_state = method_state
    )
  }
}
