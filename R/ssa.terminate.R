#' Terminates a simulation that was invoked using ssa
#'
#' Terminates a simulation that was invoked using ssa.
#'
#' Terminates an invocation of the \code{link{ssa}} wrapper function. Returns
#' the same list object as [ssa()]. This function is called from
#' within [ssa()] and is not intended to be invoked stand alone.
#'
#' @param args list of arguments and their values passed to [ssa()].
#' @param out.rxn list object as returned from [ssa.run()].
#' @param tf final time.
#' @param method ssa method to use.
#' @param maxWallTime maximum simulation wall time duration (in seconds) that
#' the simulation is allowed to run for before terminated.
#' @param verbose boolean value indicating if some basic simulation statistics
#' should be displayed on the console.
#' @seealso [ssa()]
#' @keywords misc datagen ts
#' @importFrom stats sd
ssa.terminate <- function(
  args,
  out.rxn,
  tf,
  method,
  maxWallTime,
  verbose
) {

  # Get the final time and state vector
  t <- out.rxn$timeSeries[dim(out.rxn$timeSeries)[1],1]
  x <- out.rxn$timeSeries[dim(out.rxn$timeSeries)[1],2:dim(out.rxn$timeSeries)[2]]

  # Figure out all the reasons why the simulation terminated
  terminationStatus <- NULL
  if (t>=tf)          terminationStatus <- c(terminationStatus, "finalTime")
  if (all(x==0))      terminationStatus <- c(terminationStatus, "extinction")
  if (any(x<0))       terminationStatus <- c(terminationStatus, "negativeState")
  if (all(out.rxn$eval_a==0)) terminationStatus <- c(terminationStatus, "zeroProp")
  if (out.rxn$elapsedWallTime>=maxWallTime) terminationStatus <- c(terminationStatus, "maxWallTime")

  # Calculate some stats for the used method
  stats <- list(startWallime       = out.rxn$startWallTime,
                endWallTime        = out.rxn$endWallTime,
                elapsedWallTime    = out.rxn$elapsedWallTime,
                terminationStatus  = terminationStatus,
                nSteps             = length(out.rxn$stepSize),
                meanStepSize       = mean(out.rxn$stepSize),
                sdStepSize         = sd(out.rxn$stepSize),
                nSuspendedTauLeaps = out.rxn$nSuspendedTauLeaps)

  # Figure out why the simulation terminated and print some info/stats
  if (verbose) {
    cat("tf: ",t,"\n","TerminationStatus: ",sep="")
    cat(stats$terminationStatus,sep=",")
    cat("\nDuration: ",stats$elapsedWallTime," seconds\n",
        "Method: ",method,"\n",
        "Nr of steps: ",stats$nSteps,"\n",
        "Mean step size: ",stats$meanStepSize,"+/-",stats$sdStepSize,"\n",sep="")
    if (method=="OTL") {
      cat("Nr suspended tau leaps: ",stats$nSuspendedTauLeaps,
          "(",100*(round(stats$nSuspendedTauLeaps/stats$nSteps)),"%)\n",sep="")
    }
    cat("End wall time: ",stats$endWallTime,"\n",sep="")
    cat("--------------------\n")
  }

  # Return simulation results ('chopping' off any rows in the timeSeries matrix that have no values (NA))
  return(list(data  = out.rxn$timeSeries,
              stats = stats,
              args  = args))
}
