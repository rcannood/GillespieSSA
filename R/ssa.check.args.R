#' Validates the arguments for the ssa wrapper function
#'
#' Validates the arguments for the ssa wrapper function.
#'
#' Performs basic type checking of many of the arguments passed to the
#' [ssa()] wrapper function. Note that no logical checking is
#' currently performed, e.g. which arguments are required with which method
#' (see [ssa.check.method()]). This function is called from within
#' [ssa()] and is not intended to be invoked stand alone.
#'
#' @param x0 numerical vector of initial states where the component elements
#' must be named using the same notation as the corresponding state variable in
#' the propensity vector, `a`.
#' @param a character vector of propensity functions where state variables
#' correspond to the names of the elements in `x0`.
#' @param nu numerical matrix of change if the number of individuals in each
#' state (rows) caused by a single reaction of any given type (columns).
#' @param tf final time.
#' @param method text string indicating the \acronym{SSA} method to use, the
#' valid options are: `D` --- Direct method (default method), `ETL` -
#' Explicit tau-leap, `BTL` --- Binomial tau-leap, or `OTL` ---
#' Optimized tau-leap.
#' @param tau step size for the `ETL` method (\eqn{>0}).
#' @param f coarse-graining factor for the `BTL` method (\eqn{>1}) where a
#' higher value results in larger step-size.
#' @param epsilon accuracy control parameter for the `OTL` method
#' (\eqn{>0}).
#' @param nc critical firing threshold for the `OTL` method (positive
#' integer).
#' @param hor numerical vector of the highest order reaction for each species
#' where \eqn{\mathtt{hor} \in \{1,2,22\}}{hor=(1,2,22)}. Setting
#' `hor=NaN` uses the default `hor=rep(22,N)` where `N` is the
#' number of species (See page 6 in Cao et al. 2006). Unless `hor=NaN` the
#' number of elements must equal the number of states \eqn{N}. Only applicable
#' in the `OTL` method.
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
#' size of the time step (or tau-leaps) ultimately limits the interval between
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
#' @seealso [ssa()] [ssa.check.method()]
#' @keywords misc datagen ts
ssa.check.args <- function(
  x0,
  a,
  nu,
  tf,
  method,
  tau,
  f,
  epsilon,
  nc,
  hor,
  dtf,
  nd,
  ignoreNegativeState,
  consoleInterval,
  censusInterval,
  verbose
) {

  # Do some basic check of the argument types
  if (!is.numeric(x0))              stop("'x0' is not numeric")
  if (!is.character(a))             stop("'a' is not of character type")
  if (!is.numeric(nu))              stop("'nu' is not numeric")
  if (!is.numeric(tf))              stop("'tf' is not numeric")
  if (!is.character(method))        stop("'method' is not of character type")
  if (!is.numeric(tau))             stop("'tau' is not numeric")
  if (!is.numeric(f))               stop("'f' is not numeric")
  if (!is.numeric(epsilon))         stop("'epsilon' is not numeric")
  if (!is.numeric(nc))              stop("'nc' is not numeric")
  if (!is.numeric(hor))             stop("'hor' is not numeric")
  if (!is.numeric(dtf))             stop("'dtf' is not numeric")
  if (!is.numeric(nd))              stop("'nd' is not numeric")
  if (!is.numeric(consoleInterval)) stop("'consoleInterval' is not numeric")
  if (!is.numeric(censusInterval))  stop("'censusInterval' is not numeric")
  if ((ignoreNegativeState != TRUE) & (ignoreNegativeState != FALSE))
    stop("'ignoreNegativeState' is not boolean")
  if ((verbose != TRUE) & (verbose != FALSE)) stop("'verbose' is not boolean")
  if (is.null(names(x0))) stop("'x0' is missing element names")
}
