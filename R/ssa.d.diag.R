#' Direct method (D) for nu-diagonalized systems
#'
#' Direct method (D) for nu-diagonalized systems.
#'
#' @details Performs one time step using the Direct method. It is usually called
#' from within [ssa()], but can be invoked directly, see [ssa.d()] for Examples.
#'
#' @param a vector of evaluated propensity functions.
#' @param nu state-change matrix.
#'
#' @return A list with two elements:
#' * the time leap (`tau`) and
#' * the realized state change vector (`nu_j`).
#'
#' @seealso [ssa.d()]
#' @keywords misc datagen ts
#'
#' @importFrom stats runif
ssa.d.diag <- function(a, nu) {
  j <- sample.int(length(a), size = 1, prob = a)
  nu_j <- ssa.nutiling(a, nu, j)
  tau <- -log(runif(1)) / sum(a)

  list(
    tau = tau,
    nu_j = nu_j
  )
}
