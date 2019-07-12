#' Higher-level interface to the method functions
#'
#' Higher-level interface to the method functions.
#'
#' Invokes a specific method function until the termination criteria are
#' fulfilled. Updates the state vector, time, and re-evaluates the propensity
#' functions in-between time steps. Also collects simulation data and returns
#' it as a list object. This function is called from within [ssa()]
#' and is not intended to be invoked stand alone.
#'
#' @param x0 initial states vector.
#' @param a vector of propensity functions.
#' @param nu state-change matrix.
#' @param parms vector of model parameters.
#' @param tf final time.
#' @param method ssa method to use.
#' @param tau step size for the `ETL` method (\eqn{>0}).
#' @param f coarse-graining factor for the `BTL` method (\eqn{>1}) where a
#' higher value results in larger step-size.
#' @param epsilon accuracy control parameter for the `OTL` method
#' (\eqn{>0}).
#' @param nc critical firing threshold for the `OTL` method (positive
#' integer).
#' @param hor numerical vector of the highest order reaction for each species
#' where \eqn{\mathtt{hor} \in \{1,2,22\}}{hor=(1,2,22)}. Only applicable in
#' the `OTL` method.
#' @param dtf `D` method threshold factor for the `OTL` method. The
#' `OTL` method is suspended if `tau` it estimates is smaller than
#' the `dtf` multiple of the `tau` that the `D` method would
#' have used (i.e. \eqn{\tau_{\mathtt{OTL}} < \mathtt{dtf} \times
#' \tau_{\mathtt{D}}}{tau_OTL<dtf*\tau_D}) (See step 3, page 3 in Cao et al.
#' 2006).
#' @param nd number of single-reaction steps performed using the Direct method
#' during `otl` suspension (See step 3, page 3, Cao et al. 2006).
#' @param ignoreNegativeState boolean object indicating if negative state
#' values should be ignored (this can occur in the `etl` method). If
#' `ignoreNegativeState=TRUE` the simulation finishes gracefully when
#' encountering a negative population size (i.e. does not throw an error). If
#' `ignoreNegativeState=FALSE` the simulation stops with an error message
#' when encountering a negative population size.
#' @param consoleInterval (approximate) interval at which `ssa` produces
#' simulation status output on the console (assumes `verbose=TRUE`). If
#' `consoleInterval=0` console output is generated each time step (or
#' tau-leap). If `consoleInterval=Inf` no console output is generated.
#' Note, `verbose=FALSE` disables all console output. **Console
#' output drastically slows down simulations.**
#' @param censusInterval (approximate) interval between recording the state of
#' the system. If `censusInterval=0` \eqn{(t,x)} is recorded at each time
#' step (or tau-leap). If `censusInterval=Inf` only
#' \eqn{(t_0,x_0)}{(t0,x0)} and \eqn{(t_f,x_t)}{(tf,xf)} is recorded. Note, the
#' size of the time step (or tau-leaps) ultimatelly limits the interval between
#' subsequent recordings of the system state since the state of the system
#' cannot be recorded at a finer time interval the size of the time steps (or
#' tau-leaps).
#' @param verbose boolean object indicating if the status of the simulation
#' simulation should be displayed on the console. If `verbose=TRUE` the
#' elapsed wall time and \eqn{(t,x)} is displayed on the console every
#' `consoleInterval` time step and a brief summary is displayed at the end
#' of the simulation. If `verbose=FALSE` the simulation runs
#' *entirely* silent (overriding `consoleInterval`). **Verbose
#' runs drastically slows down simulations.**
#' @param maxWallTime maximum wall time duration (in seconds) that the
#' simulation is allowed to run for before terminated. This option is useful,
#' in particular, for systems that can end up growing uncontrollably.
#' @return Returns a list object with the following elements,
#' \item{timeSeries}{a numerical matrix object of the simulation time series
#' where the first column is the time vector and subsequent columns are the
#' state frequencies.} \item{eval_a}{vector of the evaluated propensity
#' functions.} \item{elapsedWallTime}{elapsed wall time in seconds.}
#' \item{startWallTime}{start wall clock time (YYYY-mm-dd HH:MM:SS)}.
#' \item{endWallTime}{end wall clock time (YYYY-mm-dd HH:MM:SS).}
#' \item{stepSize}{vector of step sizes (i.e. time increments).}
#' \item{nSuspendedTauLeaps}{number of steps performed using the Direct method
#' due to `OTL` suspension (only applicable for the `OTL` method).}
#' @seealso [ssa()]
#' @keywords misc datagen ts
#'
#' @importFrom utils flush.console
ssa.run <- function(x0,a,nu,parms,tf,method,tau,f,epsilon,nc,hor,dtf,nd,
                    ignoreNegativeState,consoleInterval,censusInterval,
                    verbose,maxWallTime) {

  # In this function the user defined state variables and model parameters get
  # assigned. To reduce the risk for potential name clashed with internal
  # variables some of the function arguments and internal variables are
  # (re)named using dot notation (i.e. .*).
  .x0 <- x0
  .a <- a
  .nu <- nu
  .tau <- tau
  .f <- f
  .epsilon <- epsilon
  .nc <- nc

  # Assign the state variables defined in the state vector
  .x <- .x0
  varNames <- names(.x)
  for (.i in seq(length(varNames))) assign(varNames[.i],.x[[.i]])

  # Assign the parameters defined in the parms vector
  if (!is.null(parms)) {
    parmsNames <- names(parms)
    for (.i in seq(length(parmsNames))) assign(parmsNames[.i],parms[[.i]])
  }

  # Initialize miscellaneous counters
  .t <- 0 # Initialize the simulation time
  timeOfNextCensus <- .t + censusInterval  # Time of first data census
  timeForConsole   <- .t + consoleInterval # Time of first console output
  timeToTerminate  <- FALSE

  # Add the initial state of the system to the time series matrix and 'pre-grow'
  # with NAs
  timeSeries <- c(.t, .x) # First data point
  numCols    <- length(timeSeries)
  timeSeries <- rbind(timeSeries, matrix(nrow=1000, ncol=(numCols)))

  # Set up empty vectors for the evaluated propensity functions
  .M <- length(.a)
  eval_a <- rep(0,.M)

  # Evaluate the initial transition rates by evaluating the propensity functions
  parse_a <- parse(text=.a)
  for (.i in seq(length(parse_a))) eval_a[.i] <- eval(parse_a[.i])
  if (any(eval_a<0)) stop("negative propensity function")

  # If required (depends on the solver method) check if hor vector is defined
  # as NA in which case the user did not submitt his/her own hor vector as an
  # argument to ssa(). Fall back to the conservative default vaule of 2 for
  # each species. If the hor vector was defined by the user check it's length
  # (should have the same number of elements as the state vector) and check that
  # is only consists of 1, 2, or 22 (i.e. first-, second-order, or the homo-dimer
  # reactions).
  if (method == "OTL") {
    if (any(is.na(hor)))
      hor <- rep(2,length(.x0)) # Undefined hor - use default values
    else if (length(hor) != length(.x0))
      stop("length of hor vector is different from length of 'x0'")
    else if (any(hor!=1 & hor!=2 & hor!=22))
      stop("wrong value(s) in hor vector (can only be 1, 2, or 22)")
  }

  #############################################################################
  # We are ready to roll, start the simulation loop...
  # Continue the simulation loop as long as we have not reached the end time,
  # all of the populations have not gone extincs, as long as no population is
  # negative (occasinal by-product of the ETL method), and as long as at
  # least one propensity function is larger than zero
  #############################################################################

  # Start the timer and make an announcement if running silent
  procTimeStart   <- proc.time()
  elapsedWallTime <- 0
  startWallTime   <- format(Sys.time())
  if (verbose) {
    cat("Running ",method,
        " method with console output every ",consoleInterval,
        " time step\n",sep="")
    cat("Start wall time: ",startWallTime,"...\n",sep="")
    flush.console()
  }

  # Display the first time step on the console (not necessary if
  # consoleInterval=0 since it is displayed every time step anyway and hence is
  # already taken care of below
  if ((verbose) & (consoleInterval>0)) {
    cat("t=",.t," : ",sep="")
    cat(.x,sep=",")
    cat("\n")
    flush.console()
  }

  stepSize <- NULL
  currentRow <- 2 # First row contains (t0,x0)
  suspendedTauLeapMethod <- FALSE
  nSuspendedTauLeaps <- 0

  while( (.t<tf) & (any(.x>0)) &
         (all(.x>=0)) & (any(eval_a>0)) &
         (elapsedWallTime<=maxWallTime) ) {
    doCalc <- TRUE
    if ((verbose) & (timeForConsole<=.t)) {
      cat("(",elapsedWallTime,"s) t=",.t," : ",sep="")
      cat(.x,sep=",")
      cat("\n")
      flush.console()
      timeForConsole <- timeForConsole + consoleInterval
    }

    switch( method,
            "D" = { out <- ssa.d(eval_a, .nu)
                    if (suspendedTauLeapMethod) {
                      suspendedTauLeapMethod <- suspendedTauLeapMethod - 1
                      nSuspendedTauLeaps <- nSuspendedTauLeaps + 1
                      if (!suspendedTauLeapMethod) method <- "OTL"
                    }
                  },
          "ETL" = { out <- ssa.etl(eval_a, .nu, .tau) },
          "BTL" = { out <- ssa.btl(.x, eval_a, .nu, .f) },
          "OTL" = { out <- ssa.otl(.x, eval_a, .nu, hor, .nc, .epsilon, dtf, nd)
                    suspendedTauLeapMethod <- out$suspendedTauLeapMethod
                    if (suspendedTauLeapMethod) {
                      method <- "D"
                      doCalc <- FALSE
                    }
                  },
       "D.diag" = { out <- ssa.d.diag(eval_a, .nu)
                    if (suspendedTauLeapMethod) {
                      suspendedTauLeapMethod <- suspendedTauLeapMethod - 1
                      nSuspendedTauLeaps <- nSuspendedTauLeaps + 1
                      if (!suspendedTauLeapMethod) method <- "OTL.diag"
                    }
                  },
     "ETL.diag" = { out <- ssa.etl.diag(eval_a, .nu, .tau) },
     "BTL.diag" = { out <- ssa.btl.diag(.x, eval_a, .nu, .f) },
     "OTL.diag" = { out <- ssa.otl.diag(.x, eval_a, .nu, hor, .nc, .epsilon, dtf, nd)
                    suspendedTauLeapMethod <- out$suspendedTauLeapMethod
                    if (suspendedTauLeapMethod) {
                      method <- "D.diag"
                      doCalc <- FALSE
                    }
                  },
                  stop("unknown SSA method")
    )

    if (doCalc) {
      .t <- .t + out$tau  # Update the time
      .x <- .x + out$nu_j # Update the state vector

      # Check that no states are negative (can occur in some tau-leaping methods)
      if ((any(.x<0)) & (!ignoreNegativeState)) {
        cat("at least one population in 'x' is negative. Bailing to browser...\n")
        browser()
      }

      # We need to record the step size separatelly from the resultMatrix (below)
      # since the history may not be recorded at each step (depending on the value
      # of 'censusInterval')
      stepSize <- c(stepSize, out$tau)

      # If it's time record the current state of the system (t,x)
      if (timeOfNextCensus <= .t) {
        timeSeries[currentRow,] <- c(.t, .x)
        currentRow              <- currentRow + 1
        timeOfNextCensus        <- .t + censusInterval

        # If necessary add empty rows to the time series matrix
        if (currentRow > dim(timeSeries)[1])
          timeSeries <- rbind(timeSeries, matrix(nrow=1000, ncol=(numCols)))
      } # if()

      # Evaluate the transition rates for the next time step
      for (.i in seq(length(varNames))) assign(varNames[.i],.x[[.i]])
      for (.i in seq(length(parse_a))) eval_a[.i] <- eval(parse_a[.i])

      eval_a[is.na(eval_a)] <- 0 # Replace NA with zero (0/0 gives NA)
      if(any(eval_a<0))
        warning("negative propensity function - coersing to zero\n")
      eval_a[eval_a<0] <- 0
    } # if (!suspendedTauLeapMethod)
    procTimeEnd <- proc.time()
    elapsedWallTime <- procTimeEnd[3] - procTimeStart[3]
  } # while()

  # If applicable, display the last time step on the console
  if (verbose) {
    cat("t=",.t," : ",sep="")
    cat(.x,sep=",")
    cat("\n")
    flush.console()
  }

  # Remove all the remaining "pre-grown" NA rows
  timeSeries <- timeSeries[!is.na(timeSeries[,1]),]

  # Record the final state of the system
  timeSeries <- rbind(timeSeries, c(.t, .x))
  endWallTime <- format(Sys.time())

  return(list(timeSeries=timeSeries, eval_a=eval_a, elapsedWallTime=elapsedWallTime, startWallTime=startWallTime, endWallTime=endWallTime, stepSize=stepSize, nSuspendedTauLeaps=nSuspendedTauLeaps))
}
