# Copyright 2007, 2008, 2010 Mario Pineda-Krch.
#
# This file is part of the R package GillespieSSA.
#
# GillespieSSA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# GillespieSSA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with GillespieSSA.  If not, see <http://www.gnu.org/licenses/>.



#' Validates consistency of the system definition
#'
#' Validates consistency of the system definition.
#'
#' Performs a few basic consistency checks the defined system, e.g. that the
#' number of rows and columns in the state-change matrix and the number of
#' elements in the initial state vector and the vector of propensity functions
#' are consistent. This function is called from within [ssa()] and is
#' not intended to be invoked stand alone.
#'
#' @param x0 numerical vector of initial states where the component elements
#' must be named using the same notation as the corresponding state variable in
#' the propensity vector, `a`.
#' @param a character vector of propensity functions where state variables
#' correspond to the names of the elements in `x0`.
#' @param nu numerical matrix of change if the number of individuals in each
#' state (rows) caused by a single reaction of any given type (columns).
#' @param method text string indicating the \acronym{SSA} method to use, the
#' valid options are: `D` --- Direct method (default method), `ETL` -
#' Explicit tau-leap, `BTL` --- Binomial tau-leap, or `OTL` ---
#' Optimized tau-leap.
#' @param tau step size for the `ETL` method (\eqn{>0}).
#' @param f coarse-graining factor for the `BTL` method (\eqn{>1}) where a
#' higher value results in larger step-size.
#' @seealso [ssa()] [ssa.check.args()]
#' @keywords misc datagen ts
ssa.check.method <- function(x0,a,nu,method,tau,f) {

  # Check the consistency of the system dimensions, i.e. number of rows and
  # columns in the state-change matrix and the number of elements in the initial
  # state vector and the vector of propensity functions
  if ((length(a)/dim(nu)[2]) != (length(x0)/dim(nu)[1]))
    stop("inconsistent system dimensions (unequal 'nu' tessallation)")
  if (((length(a)%%dim(nu)[2])>0) || ((length(x0)%%dim(nu)[1])>0))
  stop("inconsistent system dimensions (fractional tessallation)")

  # For the ETL method tau>0
  if ((method=="ETL") & (!(tau>0))) stop("ETL method requires tau>0")

  # Check that f (used in the BTL method) is >1
  if (method=="BTL" & f<=1) stop("f has to be >1")
}
