#' Approximate Gaussian Copula Log-Likelihood
#'
#' Computes an approximate log-likelihood for a Gaussian copula count
#' time series model from latent lower and upper truncation bounds.
#' The approximation method is selected through the argument
#' \code{method}.
#'
#' The package currently supports three likelihood approximations:
#' continuous extension (CE), Geweke--Hajivassiliou--Keane simulation
#' (GHK), and Time Series Minimax Exponential Tilting (TMET).
#'
#' @param lower Numeric vector of length \code{n} giving the lower
#'   truncation bounds of the latent variables.
#' @param upper Numeric vector of length \code{n} giving the upper
#'   truncation bounds of the latent variables.
#' @param tau Numeric vector of ARMA dependence parameters ordered as
#'   \code{c(phi_1, ..., phi_p, theta_1, ..., theta_q)}.
#' @param od Integer vector \code{c(p, q)} specifying the AR and MA
#'   orders of the latent ARMA process.
#' @param method Character string specifying the likelihood
#'   approximation method. Must be one of \code{"CE"},
#'   \code{"GHK"}, or \code{"TMET"}.
#' @param c Smoothing parameter for the CE approximation. Used only when
#'   \code{method = "CE"}. Default is \code{0.5}.
#' @param pm Integer specifying the number of past lags used to
#'   approximate an ARMA(\eqn{p,q}) process by a finite-order AR
#'   representation. Used only when \code{method = "TMET"}.
#' @param M Positive integer specifying the number of Monte Carlo or
#'   quasi-Monte Carlo samples. Used by the simulation-based methods
#'   \code{"GHK"} and \code{"TMET"}.
#' @param QMC Logical; if \code{TRUE} (default), quasi-Monte Carlo
#'   integration is used when applicable. Otherwise, standard Monte
#'   Carlo sampling is used.
#' @param ret_llk Logical; if \code{TRUE} (default), returns the
#'   approximate log-likelihood. If \code{FALSE}, method-specific
#'   internal quantities are returned for diagnostic or research use.
#'
#' @return A numeric scalar giving the approximate log-likelihood. If
#'   \code{ret_llk = FALSE}, method-specific diagnostic output is
#'   returned.
#'
#' @details
#' The function \code{pmvn()} provides a unified interface for Gaussian
#' copula likelihood approximation. The argument \code{method} selects
#' among:
#' \itemize{
#'   \item \code{"CE"}: continuous extension approximation,
#'   \item \code{"GHK"}: sequential importance sampling via the GHK simulator,
#'   \item \code{"TMET"}: minimax exponential tilting approximation.
#' }
#' The arguments \code{c}, \code{pm}, \code{M}, and \code{QMC} are used
#' only by the methods to which they apply.
#'
#' @seealso \code{\link{pmvt}}, \code{\link{sim_poisson}},
#'   \code{\link{poisson.marg}}
#'
#' @references
#'
#'Nguyen, Q. N. and De Oliveira, V. (2026). Approximating Gaussian Copula Models for Count Time Series: 
#'  Connecting the Distributional Transform and a Continuous Extension, \emph{Journal of Applied Statistics}, \strong{53}, 1--22.
#'
#'Nguyen, Q. N., and De Oliveira, V. (2026), Likelihood Inference in Gaussian Copula Models for Count Time Series via Minimax Exponential Tilting,
#'\emph{Computational Statistics and Data Analysis}, \strong{218}: 108344.

#'
#' @examples
#' mu <- 10
#' tau <- 0.2
#' arma_order <- c(1, 0)
#'
#' sim_data <- sim_poisson(mu = mu, tau = tau, arma_order = arma_order,
#'                         nsim = 500, family = "gaussian", seed = 1)
#' y <- sim_data$y
#'
#' lower <- qnorm(ppois(y - 1, lambda = mu))
#' upper <- qnorm(ppois(y, lambda = mu))
#'
#' ## Continuous extension
#' pmvn(lower, upper, tau = tau, od = arma_order, method = "CE", c = 0.5)
#'
#' ## GHK approximation
#' pmvn(lower, upper, tau = tau, od = arma_order, method = "GHK", M = 1000)
#'
#' ## TMET approximation
#' pmvn(lower, upper, tau = tau, od = arma_order, method = "TMET",
#'      pm = 30, M = 1000)
#' @export
pmvn <- function(lower, upper, tau, od, method = c("CE", "GHK", "TMET"),
                 c = 0.5, pm = 30, M = 1000, QMC = TRUE, ret_llk = TRUE) {
  
  method <- match.arg(method)
  
  switch(
    method,
    CE = sum(ce_core(lower, upper, tau, od, family = "gaussian", c = c,
                 ret_llk = ret_llk)),
    GHK = sum(ghk_core(lower, upper, tau, od, family = "gaussian",
                   M = M, QMC = QMC, ret_llk = ret_llk)),
    TMET = sum(tmet_core(lower, upper, tau, od, family = "gaussian",
                     pm = pm,M = M, QMC = QMC))
  )
}

