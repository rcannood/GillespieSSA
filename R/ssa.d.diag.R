#' Direct method (D) for nu-diagonalized systems
#'
#' Performs one time step using the Direct method. It is usually called
#' from within [ssa()], but can be invoked directly, see [ssa.d()] for Examples.
#'
#' @param a vector of evaluated propensity functions.
#' @param nu state-change matrix.
#'
#' @return A list with two elements, 1) the time leap (`tau`) and 2) the realized state change vector (`nu_j`).
#'
#' @seealso [ssa.d()]
#' @keywords misc datagen ts
#'
#' @importFrom stats runif
#'
#' @examples
#' ## Not intended to be invoked stand alone
#'
#' @export
ssa.d.diag <- function(a, nu) {
  j <- sample.int(length(a), size = 1, prob = a)
  nu_j <- ssa.nutiling(a, nu, j)
  tau <- -log(runif(1)) / sum(a)

  list(
    tau = tau,
    nu_j = nu_j
  )
}
