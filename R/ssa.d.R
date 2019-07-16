#' Direct method (D)
#'
#' Direct method implementation of the \acronym{SSA} as described by Gillespie (1977).
#' Should be passed as `method` argument for `ssa()`.
#'
#' @seealso [GillespieSSA-package], [ssa()]
#' @references Gillespie (1977)
#' @keywords misc datagen ts
#'
#' @importFrom stats runif
#'
#' @examples
#' ssa.d()
#'
#' @export
ssa.d <- function() {
  ssa_method("D", list())
}

ssa_step.ssa_D <- function (method, x, a, nu, method_state) {
  j <- sample.int(length(a), size = 1, prob = a)
  nu_j <- nu[,j]
  tau <- -log(runif(1)) / sum(a)

  list(
    tau = tau,
    nu_j = nu_j
  )
}

ssa_step_diag.ssa_D <- function (method, x, a, nu, method_state) {
  j <- sample.int(length(a), size = 1, prob = a)
  nu_j <- ssa.nutiling(a, nu, j)
  tau <- -log(runif(1)) / sum(a)

  list(
    tau = tau,
    nu_j = nu_j
  )
}
