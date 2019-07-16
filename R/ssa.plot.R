#' Simple plotting of ssa output
#'
#' Provides basic functionally for simple and quick time series plot of
#' simulation output from [ssa()].
#'
#' @param out data object returned from [ssa()].
#' @param file name of the output file (only applicable if
#' `plot.device!="x11"`.
#' @param by time increment in the plotted time series
#' @param plot.from first population to plot the time series for (see note)
#' @param plot.to last population to plot the time series for (see note)
#' @param plot.by increment in the sequence of populations to plot the time
#' series for (see note)
#' @param show.title boolean object indicating if the plot should display a
#' title
#' @param show.legend boolean object indicating if the legend is displayed
#' @note The options `by`, `plot.from`, `plot.to`, and
#' `plot.by` can be used to plot a sparser sequence of data points. To
#' plot the population sizes using a larger time interval the `by` option
#' can be set, e.g. to plot only every 10th time point `by=10`. To plot
#' only specific populations the `plot.from`, `plot.to`, and
#' `plot.by` options can be set to subset the state vector. Note that the
#' indexing of the populations is based on the \eqn{(t,\mathbf{X})}{(t,X)}
#' vector, i.e. the first column is the time vector while the first population
#' is index by 2 and the last population by \eqn{N+1}. Display of a plot title
#' above the plot and legend is optional (and are set with the arguments
#' show.title and show.legend. Above the plot panel miscellaneous information
#' for the simulation are displayed, i.e. method, elapsed wall time, number of
#' time steps executed, and the number of time steps per data point.
#' @seealso [GillespieSSA-package], [ssa()]
#' @keywords misc datagen ts device utilities hplot
#'
#' @examples
#' ## Define the Kermack-McKendrick SIR model and run once using the Direct method
#' parms <- c(beta=.001, gamma=.100)
#' x0 <- c(S=500, I=1, R=0)                         # Initial state vector
#' nu <- matrix(c(-1,0,1,-1,0,1),nrow=3,byrow=TRUE) # State-change matrix
#' a  <- c("beta*S*I", "gamma*I")                   # Propensity vector
#' tf <- 100                                        # Final time
#' simName <- "Kermack-McKendrick SIR"
#' out <- ssa(x0,a,nu,parms,tf,method="D",simName,verbose=TRUE,consoleInterval=1)
#'
#' ## Basic ssa plot
#' ssa.plot(out)
#'
#' # Plot only the infectious class
#' ssa.plot(out,plot.from=3,plot.to=3)
#'
#' @importFrom grDevices rainbow
#' @importFrom graphics matplot title mtext legend
#'
#' @export
ssa.plot <- function(
  out = stop("requires simulation output object"),
  file = "ssaplot",
  by = 1,
  plot.from = 2,
  plot.to = ncol(out$data),
  plot.by = 1,
  show.title = TRUE,
  show.legend = TRUE
) {
  if (plot.from == 1 || plot.from > ncol(out$data) || plot.from > plot.to) stop("error in plot.from/plot.to arguments")

  # Render the plot(s)
  colorVector <- rainbow(ncol(out$data)-1)
  mask <- seq(1, nrow(out$data), by = by)
  matplot(
    out$data[mask,1],
    out$data[mask, seq(plot.from, plot.to, plot.by)],
    pch = 19,
    cex = 0.1,
    col = colorVector,
    bty = "n",
    xlab = "Time",
    ylab = "Frequency"
  )
  if (show.title) title(out$args$simName)
  legendTxt <- names(out$arg$x0)

  # If there are more states than 20 the legend starts to look crazy, so we don't show it...
  if (length(legendTxt) < 20 & show.legend)
    legend("topright", legend = legendTxt, bty = "y", pch = 19, col = colorVector)

  stepShowStr <- paste0("(", by, " steps/point)")
  textStr <- paste0(
    out$args$method$name, ", ",
    round(out$stats$elapsedWallTime, 2), " sec, ",
    out$stats$nSteps, " steps ", stepShowStr
  )
  mtext(textStr, line = 0, cex = 0.75)
}
