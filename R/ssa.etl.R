#' Explicit tau-leap method (ETL)
#'
#' Explicit tau-leap method implementation of the \acronym{SSA} as described by
#' Gillespie (2001). Should be passed as `method` argument for `ssa()`.
#'
#' @details Performs one time step using the Explicit tau-leap method.
#' Intended to be invoked by [ssa()].
#'
#' @param tau tau-leap.
#'
#' @references Gillespie (2001)
#'
#' @seealso [GillespieSSA-package], [ssa()]
#'
#' @examples
#' ssa.etl(tau = .1)
#'
#' @keywords misc datagen ts
#'
#' @importFrom stats rpois
#'
#' @export
ssa.etl <- function(tau = 0.3) {
  if (!is.numeric(tau)) stop("'tau' is not numeric")
  if (tau <= 0) stop("ETL method requires tau>0")

  ssa_method(
    "ETL",
    list(
      tau = tau
    )
  )
}

ssa_step.ssa_ETL <- function (method, x, a, nu, method_state) {
  tau <- method$params$tau
  M <- length(a)
  k <- rpois(M, a*tau)
  # MPK: Strictly speaking it is not correct to call the realized state-change
  # vector nu_j here since, in leap methods, actually is not reaction specific.
  # In Pineda-Krch (JSS ms) it is refered to as the \bm{\tilde{nu}}
  nu_j <- rowSums(matrix(rep(k,nrow(nu)),byrow=TRUE,ncol=M)*nu)

  list(
    tau = tau,
    nu_j = nu_j
  )
}

ssa_step_diag.ssa_ETL <- function (method, x, a, nu, method_state) {
  tau <- method$params$tau
  MU <- length(a)          # Nr of reaction channels
  k  <- rpois(MU, a*tau)   # Nr of firings per channel
  M  <- ncol(nu)    # Nr of reaction channel per patch (nu)
  U  <- MU / M             # Nr of tilings
  nu_j <- NULL

  for (f in (seq_len(U)-1)) {
    nu_j_ <- rowSums(matrix(rep(k[1:M+f*M],nrow(nu)),byrow=TRUE,ncol=M)*nu)
    nu_j <- c(nu_j, nu_j_)
  }

  list(
    tau = tau,
    nu_j = nu_j
  )
}
