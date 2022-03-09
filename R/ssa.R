#' Invoking the stochastic simulation algorithm
#'
#' Main interface function to the implemented \acronym{SSA} methods. Runs a
#' single realization of a predefined system.
#'
#' @usage
#' ssa(
#'                    x0,            # initial state vector
#'                     a,            # propensity vector
#'                    nu,            # state-change matrix
#'                 parms = NULL,     # model parameters
#'                    tf,            # final time
#'                method = ssa.d(),  # SSA method
#'               simName = "",
#'                   tau = 0.3,      # deprecated
#'                     f = 10,       # deprecated
#'               epsilon = 0.03,     # deprecated
#'                    nc = 10,       # deprecated
#'                   hor = NA_real_, # deprecated
#'                   dtf = 10,       # deprecated
#'                    nd = 100,      # deprecated
#'   ignoreNegativeState = TRUE,
#'       consoleInterval = 0,
#'        censusInterval = 0,
#'               verbose = FALSE,
#'           maxWallTime = Inf
#' )
#'
#' @details Although `ssa` can be invoked by only specifying the system
#' arguments (initial state vector `x0`, propensity vector `a`,
#' state-change matrix `nu`), the final time (`tf`), and the
#' \acronym{SSA} method to use, substantial improvements in speed and accuracy
#' can be obtained by adjusting the additional (and optional) `ssa`
#' arguments. By default `ssa` (tries to) use conservative default values
#' for the these arguments, prioritizing computational accuracy over
#' computational speed. These default values are, however, **not** fool
#' proof for the approximate methods, and occasionally one will have to hand
#' tweak them in order for a stochastic model to run appropriately.
#'
#' @param x0 numerical vector of initial states where the component elements
#'   must be named using the same notation as the corresponding state variable in
#'   the propensity vector, `a`.
#' @param a character vector of propensity functions where state variables
#'   correspond to the names of the elements in `x0`.
#' @param nu numerical matrix of change if the number of individuals in each
#'   state (rows) caused by a single reaction of any given type (columns).
#' @param parms named vector of model parameters.
#' @param tf final time.
#' @param method an SSA method,
#' the valid options are:
#' * [ssa.d()] --- Direct method (default method),
#' * [ssa.etl()] - Explicit tau-leap,
#' * [ssa.btl()] --- Binomial tau-leap, or
#' * [ssa.otl()] --- Optimized tau-leap.
#' @param simName optional text string providing an arbitrary name/label for
#' the simulation.
#' @param tau \[DEPRECATED\], see [ssa.etl()]
#' @param f \[DEPRECATED\], see [ssa.btl()]
#' @param epsilon \[DEPRECATED\], see [ssa.otl()]
#' @param nc \[DEPRECATED\], see [ssa.otl()]
#' @param hor \[DEPRECATED\], see [ssa.otl()]
#' @param dtf \[DEPRECATED\], see [ssa.otl()]
#' @param nd \[DEPRECATED\], see [ssa.otl()]
#' @param ignoreNegativeState boolean object indicating if negative state
#'  values should be ignored (this can occur in the `etl` method).
#'  If `ignoreNegativeState=TRUE` the simulation finishes gracefully when
#'  encountering a negative population size (i.e. does not throw an error).
#'  If `ignoreNegativeState=FALSE` the simulation stops with an error
#'  message when encountering a negative population size.
#' @param consoleInterval (approximate) interval at which `ssa` produces
#'   simulation status output on the console (assumes `verbose=TRUE`).
#'   If `consoleInterval=0` console output is generated each time step (or
#'   tau-leap). If `consoleInterval=Inf` no console output is generated.
#'   Note, `verbose=FALSE` disables all console output. **Console
#'   output drastically slows down simulations.**
#' @param censusInterval (approximate) interval between recording the state of the system.
#'   If `censusInterval=0` \eqn{(t,x)} is recorded at each time step (or
#'   tau-leap). If `censusInterval=Inf` only \eqn{(t_0,x_0)}{(t0,x0)}
#'   and \eqn{(t_f,x_t)}{(tf,xf)} is recorded. Note, the size of the time
#'   step (or tau-leaps) ultimately limits the interval between subsequent
#'   recordings of the system state since the state of the system cannot be
#'   recorded at a finer time interval the size of the time steps (or tau-leaps).
#' @param verbose  boolean object indicating if the status of the simulation
#'   simulation should be displayed on the console. If `verbose=TRUE`
#'   the elapsed wall time and \eqn{(t,x)} is displayed on the console every
#'   `consoleInterval` time step and a brief summary is displayed at
#'   the end of the simulation. If `verbose=FALSE` the simulation runs
#'   *entirely* silent (overriding `consoleInterval`).
#'   **Verbose runs drastically slows down simulations.**
#' @param maxWallTime maximum wall time duration (in seconds) that the
#' simulation is allowed to run for before terminated. This option is
#' useful, in particular, for systems that can end up growing uncontrolably.
#'
#' @return Returns a list object with the following elements,
#'  * `data`: a numerical matrix object of the simulation time series where the first column is the time vector and subsequent columns are the state frequencies.
#'  * `stats`: sub-list object with elements containing various simulation statistics. The of the sub-list are:
#'  * `stats$startWallTime`: start wall clock time (YYYY-mm-dd HH:MM:SS).
#'  * `stats$endWallTime`: end wall clock time (YYYY-mm-dd HH:MM:SS).
#'  * `stats$elapsedWallTime`: elapsed wall time in seconds.
#'  * `stats$terminationStatus`: string vector listing the reason(s) for the
#'    termination of the realization in 'plain words'. The possible termination statuses are:
#'    - `finalTime` = if the simulation reached the maximum simulation time `tf`,
#'    - `extinction` = if the population size of all states is zero,
#'    - `negativeState` = if one or several states have a negative population size (can occur in the ETL method),
#'    - `zeroProp` = if all the states have a zero propensity function,
#'    - `maxWallTime` = if the maximum wall time has been reached. Note the termination status may have more than one message.
#'  * `stats$nSteps`` total number of time steps (or tau-leaps) executed.
#'  * `stats$meanStepSize`: mean step (or tau-leap) size.
#'  * `stats$sdStepSize`: one standard deviation of the step (or tau-leap) size.
#'  * `stats$SuspendedTauLeaps`: number of steps performed using the Direct method due to `OTL` suspension (only applicable for the `OTL` method).
#'  * `arg$...`: sub-list with elements containing all the arguments and their values used to invoke `ssa` (see Usage and Arguments list above).
#'
#' @note Selecting the appropriate \acronym{SSA} method is a trade-off between
#'   computational speed, accuracy of the results, and which \acronym{SSA}
#'   actually works for a given scenario. This depends on the characteristics of
#'   the defined system (e.g. number of reaction channels, number of species, and
#'   the absolute and relative magnitude of the propensity functions).
#'   **Not all methods are appropriate for all models.** When selecting a
#'   \acronym{SSA} method all of these factors have to be taken into
#'   consideration. The various tau-leap methods accept a number of additional
#'   arguments. While the default values of these arguments may work for some
#'   scenarios they may have to be adjusted for others. The default values for
#'   the tau-leap methods are conservative in terms of computational speed and
#'   substantial increase in efficiency may be gained by optimizing their values
#'   for a specific system.
#' @section Preparing a run: In order to invoke \acronym{SSA} the stochastic
#'   model needs at least four components, the initial state vector (`x0`),
#'   state-change matrix (`nu`), propensity vector (`a`), and the final
#'   time of the simulation (`tf`). The initial state vector defines the
#'   population sizes in all the states at \eqn{t=0}, e.g. for a system with two
#'   species `X1` and `X2` where both have an initial population size
#'   of 1000 the initial state vector is defined as `x0 <-
#'   c(X1=1000,X2=1000)`. The elements of the vector have to be labelled using
#'   the same notation as the state variables used in the propensity functions.
#'   The state-change matrix defines the change in the number of individuals in
#'   each state (rows) as caused by one reaction of a given type (columns). For
#'   example, the state-change matrix for system with the species \eqn{S_1}{S1}
#'   and \eqn{S_2}{S2} with two reactions \deqn{S_1
#'   \stackrel{c_1}{\longrightarrow} S_2}{S1 --c1--> S2} \deqn{S_2
#'   \stackrel{c_2}{\longrightarrow} 0}{S2 --c2--> 0}
#'
#'   is defined as `nu <- matrix(c(-1,0,+1,-1),nrow=2,byrow=TRUE)` where
#'   \eqn{c_1}{c1} and \eqn{c_2}{c2} are the per capita reaction probabilities.
#'   The propensity vector, `a`, defines the probabilities that a particular
#'   reaction will occur over the next infinitesimal time interval \eqn{\left[
#'   t,t+dt \right]}{[t,t+dt]}. For example, in the previous example the
#'   propensity vector is defined as `a <- c("c1*X1","c2*X2")`. The
#'   propensity vector consists of character elements of each reaction's
#'   propensity function where each state variable requires the corresponding
#'   named element label in the initial state vector (`x0`).
#'
#' @section Example: Irreversible isomerization:
#' Perhaps the simplest model that can be formulated using the \acronym{SSA}
#' is the irreversible isomerization (or radioactive decay) model. This model
#' is often used as a first pedagogic example to illustrate the \acronym{SSA}
#' (see e.g. Gillespie 1977). The deterministic formulation of this model is
#'
#' \deqn{\frac{dX}{dt}=-cX}{dX/dt=-cX}
#'
#' where the single reaction channel is
#'
#' \deqn{S \stackrel{c}{\longrightarrow} 0}{S --c--> 0.}
#'
#' By setting \eqn{X_0=1000} and \eqn{c=0.5} it is now simple to define this model
#' and run it for 10 time steps using the Direct method,
#'
#' \preformatted{
#'   out <- ssa(x0=c(X=1000),a=c("c*X"),nu=matrix(-1),parms=c(c=0.5),tf=10)
#' }
#' The resulting time series can then be displayed by,
#' \preformatted{
#'   ssa.plot(out)
#' }
#'
#' @seealso [GillespieSSA-package], [ssa.d()], [ssa.etl()], [ssa.btl()], [ssa.otl()]
#'
#' @keywords misc datagen ts
#' @examples
#'
#' ## Irreversible isomerization
#' ## Large initial population size (X=1000)
#' parms <- c(c=0.5)
#' x0  <- c(X=10000)
#' a   <- c("c*X")
#' nu  <- matrix(-1)
#' out <- ssa(x0,a,nu,parms,tf=10,method=ssa.d(),simName="Irreversible isomerization") # Direct method
#' plot(out$data[,1],out$data[,2]/10000,col="red",cex=0.5,pch=19)
#'
#' ## Smaller initial population size (X=100)
#' x0  <- c(X=100)
#' out <- ssa(x0,a,nu,parms,tf=10,method=ssa.d()) # Direct method
#' points(out$data[,1],out$data[,2]/100,col="green",cex=0.5,pch=19)
#'
#' ## Small initial population size (X=10)
#' x0  <- c(X=10)
#' out <- ssa(x0,a,nu,parms,tf=10,method=ssa.d()) # Direct method
#' points(out$data[,1],out$data[,2]/10,col="blue",cex=0.5,pch=19)
#'
#' ## Logistic growth
#' parms <- c(b=2, d=1, K=1000)
#' x0  <- c(N=500)
#' a   <- c("b*N", "(d+(b-d)*N/K)*N")
#' nu  <- matrix(c(+1,-1),ncol=2)
#' out <- ssa(x0,a,nu,parms,tf=10,method=ssa.d(),maxWallTime=5,simName="Logistic growth")
#' ssa.plot(out)
#'
#' ## Kermack-McKendrick SIR model
#' parms <- c(beta=0.001, gamma=0.1)
#' x0  <- c(S=499,I=1,R=0)
#' a   <- c("beta*S*I","gamma*I")
#' nu  <- matrix(c(-1,0,+1,-1,0,+1),nrow=3,byrow=TRUE)
#' out <- ssa(x0,a,nu,parms,tf=100,method=ssa.d(),simName="SIR model")
#' ssa.plot(out)
#'
#' ## Lotka predator-prey model
#' parms <- c(c1=10, c2=.01, c3=10)
#' x0  <- c(Y1=1000,Y2=1000)
#' a   <- c("c1*Y1","c2*Y1*Y2","c3*Y2")
#' nu  <- matrix(c(+1,-1,0,0,+1,-1),nrow=2,byrow=TRUE)
#' out <- ssa(x0,a,nu,parms,tf=100,method=ssa.etl(),simName="Lotka predator-prey model")
#' ssa.plot(out)
#'
#' @keywords misc datagen ts
#'
#' @importFrom stats sd
#' @importFrom utils flush.console
#' @importFrom methods is
#'
#' @export
ssa <- function(
  x0 = stop("undefined 'x0'"),
  a = stop("undefined 'a'"),
  nu = stop("undefined 'nu'"),
  parms = NULL,
  tf = stop("undefined 'tf'"),
  method = ssa.d(),
  simName = "",
  tau = 0.3,
  f = 10,
  epsilon = 0.03,
  nc = 10,
  hor = NA_real_,
  dtf = 10,
  nd = 100,
  ignoreNegativeState = TRUE,
  consoleInterval = 0,
  censusInterval = 0,
  verbose = FALSE,
  maxWallTime = Inf
) {

  ##############################################
  ###              INPUT CHECKS              ###
  ##############################################
  # Do some basic check of the argument types
  if (!is.numeric(x0)) stop("'x0' is not numeric")
  if (!is.character(a)) stop("'a' is not of character type")
  if (!is.numeric(nu)) stop("'nu' is not numeric")
  if (!is.numeric(tf)) stop("'tf' is not numeric")
  if (!is.numeric(consoleInterval)) stop("'consoleInterval' is not numeric")
  if (!is.numeric(censusInterval)) stop("'censusInterval' is not numeric")
  if (!is.logical(ignoreNegativeState) || is.na(ignoreNegativeState)) stop("'ignoreNegativeState' is not boolean")
  if (!is.logical(verbose) || is.na(verbose)) stop("'verbose' is not boolean")
  if (is.null(names(x0))) stop("'x0' is missing element names")

  # preserve old interface by allowing method to be a character
  if (is.character(method)) {
    method <- match.arg(toupper(method), choices = c("D", "ETL", "BTL", "OTL"))

    method <- switch(
      method,
      "D" = ssa.d(),
      "ETL" = ssa.etl(tau = tau),
      "BTL" = ssa.btl(f = f),
      "OTL" = ssa.otl(hor = hor, nc = nc, epsilon = epsilon, dtf = dtf, nd = nd),
      stop("unknown SSA method")
    )
  }

  if (!is(method, "ssa_method")) {
    stop("'method' should be an object created by ssa.d(), ssa.etl(), ssa.btl(), or ssa.otl()")
  }

  # Check the consistency of the system dimensions, i.e. number of rows and
  # columns in the state-change matrix and the number of elements in the initial
  # state vector and the vector of propensity functions
  if ((length(a) / ncol(nu)) != (length(x0) / nrow(nu)))
    stop("inconsistent system dimensions (unequal 'nu' tessallation)")
  if (((length(a) %% ncol(nu))>0) || ((length(x0) %% nrow(nu))>0))
    stop("inconsistent system dimensions (fractional tessallation)")

  # Take a snapshot of all the options so they can be saved later
  args <- as.list(environment())

  # Is the system nu-tiled along the diagonal?
  is_nu_tiled <- length(a) > ncol(nu) && length(x0) > nrow(nu)

  ##############################################
  ###             RUN SIMULATION             ###
  ##############################################

  # Assign the state variables defined in the state vector
  x <- x0

  # Initialize miscellaneous counters
  t <- 0 # Initialize the simulation time
  timeOfNextCensus <- t + censusInterval  # Time of first data census
  timeForConsole   <- t + consoleInterval # Time of first console output

  # Add the initial state of the system to the time series matrix and 'pre-grow'
  # with NAs
  timeSeries <- vector("list", length = 1000)
  timeSeries[[1]] <- c(t, x)

  # add current time to parms
  parms <- c(parms, .t = t)

  # Set up empty vectors for the evaluated propensity functions
  a.funs <- parse.propensity.functions(a, x0, parms)

  # Evaluate the initial transition rates by evaluating the propensity functions
  eval_a <- a.funs(x, parms)
  if (any(eval_a < 0)) stop("negative propensity function")

  # validate parameters
  method$params <- ssa_validate_parameters(method, x, eval_a, nu)

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
    cat(
      "Running ", method$name, " method with console output every ", consoleInterval, " time step\n",
      "Start wall time: ", startWallTime, "...\n",
      sep = ""
    )
    flush.console()
  }

  # Display the first time step on the console (not necessary if
  # consoleInterval=0 since it is displayed every time step anyway and hence is
  # already taken care of below
  if (verbose && consoleInterval > 0) {
    cat("t=", t, " : ", paste0(x, collapse = ","), "\n", sep = "")
    flush.console()
  }

  stepSize <- integer(0)
  currentRow <- 2 # First row contains (t0,x0)

  method_state <- ssa_initialise_state(method)

  while(
    t < tf &&
    any(x > 0) &&
    all(x >= 0) &&
    any(eval_a > 0) &&
    elapsedWallTime <= maxWallTime
  ) {
    if (verbose && timeForConsole <= t) {
      cat("(", elapsedWallTime, "s) t=", t, " : ", paste0(x, collapse = ","), "\n", sep="")
      flush.console()
      timeForConsole <- timeForConsole + consoleInterval
    }

    out <-
      if (!is_nu_tiled) {
        ssa_step(method, x = x, a = eval_a, nu = nu, method_state = method_state)
      } else {
        ssa_step_diag(method, x = x, a = eval_a, nu = nu, method_state = method_state)
      }

    t <- t + out$tau
    x <- x + out$nu_j

    if ("method_state" %in% names(out)) {
      method_state <- out$method_state
    }

    # Check that no states are negative (can occur in some tau-leaping methods)
    if (any(x < 0) && !ignoreNegativeState) {
      stop("at least one population in 'x' is negative.\n")
    }

    # We need to record the step size separatelly from the resultMatrix (below)
    # since the history may not be recorded at each step (depending on the value
    # of 'censusInterval')
    stepSize <- c(stepSize, out$tau)

    # If it's time record the current state of the system (t,x)
    if (timeOfNextCensus <= t) {
      timeSeries[[currentRow]] <- c(t, x)
      currentRow <- currentRow + 1
      timeOfNextCensus <- t + censusInterval

      # If necessary add empty rows to the time series matrix
      if (currentRow > length(timeSeries)) {
        length(timeSeries) <- length(timeSeries) * 2
      }
    }

    # Evaluate the transition rates for the next time step
    parms[[".t"]] <- t
    eval_a <- a.funs(x, parms)
    eval_a[is.na(eval_a)] <- 0 # Replace NA with zero (0/0 gives NA)
    if (any(eval_a < 0))
      warning("negative propensity function - coersing to zero\n")
    eval_a[eval_a < 0] <- 0

    procTimeEnd <- proc.time()
    elapsedWallTime <- procTimeEnd[3] - procTimeStart[3]
  }

  # If applicable, display the last time step on the console
  if (verbose) {
    cat("t=", t, " : ", paste(x, collapse = ","), "\n", sep="")
    flush.console()
  }

  # Record the final state of the system
  timeSeries[[currentRow]] <- c(t, x)
  endWallTime <- format(Sys.time())

  # Remove all the remaining "pre-grown" NA rows
  timeSeries <- do.call(rbind, timeSeries)
  colnames(timeSeries) <- c("t", names(x0))

  ##############################################
  ###             RETURN OUTPUT              ###
  ##############################################

  # Figure out all the reasons why the simulation terminated
  terminationStatus <- character(0)
  if (t >= tf) {
    terminationStatus <- c(terminationStatus, "finalTime")
  }
  if (all(x == 0)){
    terminationStatus <- c(terminationStatus, "extinction")
  }
  if (any(x < 0)) {
    terminationStatus <- c(terminationStatus, "negativeState")
  }
  if (all(eval_a == 0)) {
    terminationStatus <- c(terminationStatus, "zeroProp")
  }
  if (elapsedWallTime >= maxWallTime) {
    terminationStatus <- c(terminationStatus, "maxWallTime")
  }

  # Calculate some stats for the used method
  stats <- list(
    startWallime = startWallTime,
    endWallTime = endWallTime,
    elapsedWallTime = elapsedWallTime,
    terminationStatus = terminationStatus,
    nSteps = length(stepSize),
    meanStepSize = mean(stepSize),
    sdStepSize = sd(stepSize)
  )
  if (method$name == "OTL") {
    stats$nSuspendedTauLeaps <- method_state$total_suspensions
  }

  # Print some info/stats
  if (verbose) {
    cat(
      "tf: ", t, "\n",
      "TerminationStatus: ",         paste(stats$terminationStatus, collapse = ","), "\n",
      "Duration: ",                  stats$elapsedWallTime, " seconds\n",
      "Method: ",                    method$name, "\n",
      "Nr of steps: ",               stats$nSteps, "\n",
      "Mean step size: ",            stats$meanStepSize, "+/-", stats$sdStepSize, "\n",
      ifelse(method$name != "OTL", "", paste0(
        "Nr suspended tau leaps: ",  stats$nSuspendedTauLeaps, "(", 100*(round(stats$nSuspendedTauLeaps/stats$nSteps)), "%)\n"
      )),
      "End wall time: ",             stats$endWallTime,"\n",
      "--------------------\n",
      sep = ""
    )
  }

  # Return simulation results ('chopping' off any rows in the timeSeries matrix that have no values (NA))
  list(
    data = timeSeries,
    stats = stats,
    args = args
  )
}
