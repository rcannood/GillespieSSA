ssa_method <- function(
  name,
  params
) {
  out <- list(
    name = name,
    params = params
  )
  class(out) <- c(paste0("ssa_", name), "ssa_method", "list")
  out
}


ssa_initialise_state <- function (method, ...) {
  UseMethod("ssa_initialise_state", method)
}

ssa_initialise_state.ssa_method <- function(method) {
  list()
}

ssa_validate_parameters <- function (method, x, a, nu) {
  UseMethod("ssa_validate_parameters", method)
}

ssa_validate_parameters.ssa_method <- function(method, x, a, nu) {
  method$params
}

ssa_step <- function (method, x, a, nu, method_state) {
  UseMethod("ssa_step", method)
}

ssa_step_diag <- function (method, x, a, nu, method_state) {
  UseMethod("ssa_step_diag", method)
}
