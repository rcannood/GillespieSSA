#' Invoking the stochastic simulation algorithm
#'
#' Main interface function to the implemented \acronym{SSA} methods. Runs a
#' single realization of a predefined system.
#'
#' @usage
#' ssa(
#'                    x0,        # initial state vector
#'                     a,        # propensity vector
#'                    nu,        # state-change matrix
#'                 parms = NULL, # model parameters
#'                    tf,        # final time
#'                method = "D",  # SSA method
#'               simName = "",
#'                   tau = 0.3,  # only applicable for ETL
#'                     f = 10,   # only applicable for BTL
#'               epsilon = 0.03, # only applicable for OTL
#'                    nc = 10,   # only applicable for OTL
#'                   hor = NaN,  # only applicable for OTL
#'                   dtf = 10,   # only applicable for OTL
#'                    nd = 100,  # only applicable for OTL
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
#' @param method text string indicating the \acronym{SSA} method to use,
#' the valid options are:
#' * `"D"` --- Direct method (default method),
#' * `"ETL"` - Explicit tau-leap,
#' * `"BTL"` --- Binomial tau-leap, or
#' * `"OTL"` --- Optimized tau-leap.
#' @param simName optional text string providing an arbitrary name/label for
#' the simulation.
#' @param tau step size for the `ETL` method (\eqn{>0}).
#' @param f coarse-graining factor for the `BTL` method (\eqn{>1}) where a
#'   higher value results in larger step-size.
#' @param epsilon accuracy control parameter for the `OTL` method (\eqn{>0}).
#' @param nc critical firing threshold for the `OTL` method (positive integer).
#' @param hor numerical vector of the highest order reaction for each species where
#'   \eqn{\mathtt{hor} \in \{1,2,22\}}{hor=(1,2,22)}. Setting `hor=NaN` uses
#'   the default `hor=rep(22,N)` where `N` is the number of species (See
#'   page 6 in Cao et al. 2006). Unless `hor=NaN` the number of elements must
#'   equal the number of states \eqn{N}. Only applicable in the `OTL` method.
#' @param dtf `D` method threshold factor for the `OTL` method. The
#' `OTL` method is suspended if `tau` it estimates is smaller than the
#'   `dtf` multiple of the `tau` that the `D` method would have used
#'   (i.e. \eqn{\tau_{\mathtt{OTL}} < \mathtt{dtf} \times \tau_{\mathtt{D}}}{tau_OTL<dtf*\tau_D})
#'   (See step 3, page 3 in Cao et al. 2006).
#' @param nd number of single-reaction steps performed using the Direct method
#'  during `otl` suspension (See step 3, page 3, Cao et al. 2006).
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
#' out <- ssa(x0,a,nu,parms,tf=10,simName="Irreversible isomerization") # Direct method
#' plot(out$data[,1],out$data[,2]/10000,col="red",cex=0.5,pch=19)
#'
#' ## Smaller initial population size (X=100)
#' x0  <- c(X=100)
#' out <- ssa(x0,a,nu,parms,tf=10) # Direct method
#' points(out$data[,1],out$data[,2]/100,col="green",cex=0.5,pch=19)
#'
#' ## Small initial population size (X=10)
#' x0  <- c(X=10)
#' out <- ssa(x0,a,nu,parms,tf=10) # Direct method
#' points(out$data[,1],out$data[,2]/10,col="blue",cex=0.5,pch=19)
#'
#' ## Logistic growth
#' parms <- c(b=2, d=1, K=1000)
#' x0  <- c(N=500)
#' a   <- c("b*N", "(d+(b-d)*N/K)*N")
#' nu  <- matrix(c(+1,-1),ncol=2)
#' out <- ssa(x0,a,nu,parms,tf=10,method="D",maxWallTime=5,simName="Logistic growth")
#' ssa.plot(out)
#'
#' ## Kermack-McKendrick SIR model
#' parms <- c(beta=0.001, gamma=0.1)
#' x0  <- c(S=499,I=1,R=0)
#' a   <- c("beta*S*I","gamma*I")
#' nu  <- matrix(c(-1,0,+1,-1,0,+1),nrow=3,byrow=TRUE)
#' out <- ssa(x0,a,nu,parms,tf=100,simName="SIR model")
#' ssa.plot(out)
#'
#' ## Lotka predator-prey model
#' parms <- c(c1=10, c2=.01, c3=10)
#' x0  <- c(Y1=1000,Y2=1000)
#' a   <- c("c1*Y1","c2*Y1*Y2","c3*Y2")
#' nu  <- matrix(c(+1,-1,0,0,+1,-1),nrow=2,byrow=TRUE)
#' out <- ssa(x0,a,nu,parms,tf=100,method="ETL",simName="Lotka predator-prey model")
#' ssa.plot(out)
#'
#' @keywords misc datagen ts
#'
#' @export
ssa <- function(
  x0 = stop("undefined 'x0'"),
  a = stop("undefined 'a'"),
  nu = stop("undefined 'nu'"),
  parms = NULL,
  tf = stop("undefined 'tf'"),
  method = "D",
  simName = "",
  tau = 0.3,
  f = 10,
  epsilon = 0.03,
  nc = 10,
  hor = NaN,
  dtf = 10,
  nd = 100,
  ignoreNegativeState = TRUE,
  consoleInterval = 0,
  censusInterval = 0,
  verbose = FALSE,
  maxWallTime = Inf
) {
  ssa.check.args(
    x0 = x0,
    a = a,
    nu = nu,
    tf = tf,
    method = method,
    tau = tau,
    f = f,
    epsilon = epsilon,
    nc = nc,
    hor = hor,
    dtf = dtf,
    nd = nd,
    ignoreNegativeState = ignoreNegativeState,
    consoleInterval = consoleInterval,
    censusInterval = censusInterval,
    verbose = verbose
  )

  # Convert lower case method names to upper case (undocumented featurette)
  if (method=="d")   method <- "D"
  if (method=="etl") method <- "ETL"
  if (method=="btl") method <- "BTL"
  if (method=="otl") method <- "OTL"

  ssa.check.method(
    x0 = x0,
    a = a,
    nu = nu,
    method = method,
    tau = tau,
    f = f
  )

  # Is the system nu-tiled along the diagonal?
  if (length(a) > ncol(nu) && length(x0) > nrow(nu)){
    if (method=="D")   method <- "D.diag"
    if (method=="ETL") method <- "ETL.diag"
    if (method=="BTL") method <- "BTL.diag"
    if (method=="OTL") method <- "OTL.diag"
  }

  # Take a snapshot of all the options so they can be saved later
  args <- as.list(environment())

  # Run the simulation
  out.rxn <- ssa.run(
    x0 = x0,
    a = a,
    nu = nu,
    parms = parms,
    tf = tf,
    method = method,
    tau = tau,
    f = f,
    epsilon = epsilon,
    nc = nc,
    hor = hor,
    dtf = dtf,
    nd = nd,
    ignoreNegativeState = ignoreNegativeState,
    consoleInterval = consoleInterval,
    censusInterval = censusInterval,
    verbose = verbose,
    maxWallTime = maxWallTime
  )

  # Wrap up the simulation
  ssa.terminate(args, out.rxn, tf, method, maxWallTime, verbose)
}
