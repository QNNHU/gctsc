#' Simulate from Gaussian and t Copula Time Series Models
#'
#' These functions simulate time series data from Gaussian and t copula models
#' with various discrete marginals and an ARMA dependence structure.
#'
#' @param mu Mean parameter(s) for Poisson-, ZIP-, and negative
#'   binomial-type marginals. Must satisfy \eqn{\mu > 0}. May be specified
#'   as a scalar or as a numeric vector of length \code{nsim} to allow
#'   time-varying means.
#'
#' @param prob Success probability parameter(s) for binomial-type marginals.
#'   Must satisfy \eqn{0 < p < 1}. May be a scalar or a numeric vector of
#'   length \code{nsim}.
#'
#' @param tau Numeric vector of ARMA dependence coefficients, ordered as
#'   \code{c(phi_1, ..., phi_p, theta_1, ..., theta_q)}, where
#'   \eqn{\phi_i} are autoregressive (AR) coefficients and
#'   \eqn{\theta_j} are moving-average (MA) coefficients.
#'   The model \code{ARMA(0, 0)} is not supported.
#'
#' @param arma_order Integer vector \code{c(p, q)} specifying the AR and MA orders.
#'
#' @param nsim Positive integer giving the number of time points to simulate.
#'
#' @param seed Optional integer used to set the random seed.
#'
#' @param dispersion Overdispersion parameter for negative binomial marginals.
#'   Must satisfy \eqn{\kappa > 0}, where
#'   \eqn{\mathrm{Var}(Y) = \mu + \kappa \mu^2}.
#'   May be a scalar or a numeric vector of length \code{nsim}.
#'
#' @param pi0 Zero-inflation probability for ZIP, ZIB, and ZIBB marginals.
#'   Must satisfy \eqn{0 \le \pi_0 < 1}. May be a scalar or a numeric vector
#'   of length \code{nsim}.
#'
#' @param rho Intra-class correlation parameter for beta-binomial and ZIBB
#'   marginals. Must satisfy \eqn{0 < \rho < 1}, where
#'   \eqn{\mathrm{Var}(Y) = n p (1-p)\{1 + (n-1)\rho\}} and \eqn{n} is the
#'   number of trials. May be a scalar or a numeric vector of length
#'   \code{nsim}.
#'
#' @param size Number of trials for binomial-type marginals; a positive
#'   integer scalar.
#'
#' @param df Degrees of freedom for the t copula. Must be a single numeric
#'   value greater than 2. Required only when \code{family = "t"}.
#'
#' @param family Character string specifying the copula family:
#'   \code{"gaussian"} or \code{"t"}.
#'
#' @details
#'**Marginal types:**
#'\itemize{
#'  \item {Poisson}: Counts with mean  \eqn{\mu}.
#'  \item {Negative binomial (NB)}: Overdispersed counts with mean  \eqn{\mu} and dispersion parameter \eqn{\kappa}.
#'  \item {Binomial}: Number of successes in \eqn{n} trials with success probability \eqn{p}.
#'  \item {Beta–-binomial (BB)}: Binomial with success probability \eqn{p} following a beta distribution, allowing intra-cluster correlation  \eqn{\rho}.
#'  \item {Zero--inflated Poisson (ZIP)}: Poisson with extra probability \eqn{\pi_0} of an excess zero.
#'  \item {Zero--inflated binomial (ZIB)}: Binomial with extra probability \eqn{\pi_0} of an excess zero.
#'  \item {Zero--inflated beta–binomial (ZIBB)}: Beta–binomial with extra probability \eqn{\pi_0} of an excess zero.
#'}
#'
#' **Parameterization notes:**
#' \itemize{
#'   \item Negative binomial uses \code{dispersion} (\eqn{\kappa}) to model
#'         overdispersion: larger \code{dispersion} increases variance for a fixed mean.
#'   \item Beta--binomial and ZIBB use \code{rho} as the overdispersion parameter:
#'         \eqn{\rho} is the intra-class correlation, with \eqn{\rho \to 0}
#'         giving the binomial model.
#'   \item Zero--inflated marginals include a separate \code{pi0} parameter that
#'         controls the extra probability mass at zero.
#' }
#' 
#' \strong{Worked examples.}
#' Additional worked examples, including Gaussian and Student--t copula
#' models with zero-inflated marginals, are provided in the installed
#' example scripts; see \code{\link{gctsc-examples}}.
#' 
#' 
#' @return A list with components:
#' \itemize{
#'   \item \code{y}: Simulated time series data.
#'   \item \code{z}: Latent Gaussian process values.
#'   \item \code{marginal}: Marginal distribution name.
#'   \item \code{parameters}: List of parameters used.
#'   \item \code{cormat}: ARMA structure.
#' }
#' 
#' @examples
#' # Poisson example
#' sim_poisson(mu = 10, tau = c(0.2, 0.2), 
#'   arma_order = c(1, 1), nsim = 100, 
#'   family = "gaussian", seed = 42)
#'
#' # Negative Binomial example
#' sim_negbin(mu = 10, dispersion = 2, tau = c(0.5, 0.5),
#'   arma_order = c(1, 1),family = "gaussian",
#'   nsim = 100, seed =1)
#'
#' # Zero Inflated Beta-Binomial example with seasonal covariates
#' n <- 100
#' xi <- numeric(n)
#' zeta <- rnorm(n)
#' for (j in 3:n) {
#'   xi[j] <- 0.6 * xi[j - 1] - 0.4 * xi[j - 2] + zeta[j]
#' }
#' prob <- plogis(0.2 + 0.3 * sin(2 * pi * (1:n) / 12) +
#'              0.5 * cos(2 * pi * (1:n) / 12) + 0.3 * xi)
#' sim_zibb(prob, rho = 1/6, pi0 = 0.2, size = 24, tau = 0.5,
#'  arma_order = c(1, 0),family = "t", df = 10, nsim = 100)
#' @seealso \code{\link{gctsc}}, \code{\link{sim_gctsc}}, \code{\link{marginal.gctsc}},
#'   \code{\link{pmvn}}, \code{\link{pmvt}}, \code{\link{predict.gctsc}}
#' @name sim_gctsc
#' @export
sim_poisson <- function(mu, tau, arma_order, nsim,
                        family = c("gaussian","t"),
                        df = NULL, seed = NULL) {
  family <- match.arg(family)
  .check_common(nsim, tau, arma_order, seed, family , df, "sim_poisson")
  .check_mu(mu, nsim, "sim_poisson")
  
  par <- list(mu = mu, arma_order = arma_order)
  if (family == "t") par$df <- df
  
  sim(object = list(marg = "poisson", par = par, tau = tau),
      nsim = nsim, seed = seed, cop = family)
}




#' @rdname sim_gctsc
#' @export
sim_negbin <- function(mu, dispersion, tau, arma_order, nsim = 100, 
                       family = c("gaussian","t"),
                       df = NULL, seed = NULL) {
  family <- match.arg(family)
  .check_common(nsim, tau, arma_order, seed,  family , df, "sim_negbin")
  .check_mu(mu, nsim, "sim_negbin")
  .check_dispersion(dispersion, nsim, "sim_negbin")
  
  par = list(mu = mu, dispersion = dispersion, arma_order = arma_order)
  if (family == "t") par$df <- df
  
  sim(object = list(marg = "negbin",par=par,
                    tau = tau), nsim = nsim, seed = seed,  cop = family)
  
}

#' @rdname sim_gctsc
#' @export
sim_zip <- function(mu, pi0, tau, arma_order, nsim = 100, family = c("gaussian","t"), 
                    df = NULL, seed = NULL) {
  family <- match.arg(family)
  .check_common(nsim, tau, arma_order, seed,  family , df, "sim_zip")
  .check_mu(mu, nsim, "sim_zip")
  .check_pi0(pi0, nsim, "sim_zip")
  
  par = list(mu = mu, pi0 = pi0, arma_order = arma_order)
  if (family == "t") par$df <- df
  
  sim(object = list(marg = "zip",par = par,
                    tau = tau), nsim = nsim, seed = seed,  cop = family)
}

#' @rdname sim_gctsc
#' @export
sim_binom <- function(prob, size, tau, arma_order, nsim = 100, family = c("gaussian","t"),
                      df = NULL, seed = NULL) {
  family <- match.arg(family)
  .check_common(nsim, tau, arma_order, seed, family , df,  "sim_binom")
  .check_prob(prob, nsim, "sim_binom")
  .check_size(size, "sim_binom", scalar_only = TRUE)
  
  par = list(prob = prob, size = size, arma_order = arma_order)
  if (family == "t") par$df <- df
  
  sim(object = list(marg = "binom",par=par, tau = tau),
      nsim = nsim, seed = seed,  cop = family)
}

#' @rdname sim_gctsc
#' @export
sim_bbinom <- function(prob, rho, size, tau, arma_order, nsim = 100, family = c("gaussian","t"),
                  df = NULL, seed = NULL) {
  family <- match.arg(family)
  .check_common(nsim, tau, arma_order, seed, family , df,  "sim_bbinom")
  .check_prob(prob, nsim, "sim_bbinom")
  .check_rho(rho, nsim, "sim_bbinom")
  .check_size(size, "sim_bbinom", scalar_only = TRUE)
  
  par = list(prob = prob, rho = rho, size = size, arma_order = arma_order)
  if (family == "t") par$df <- df
  
  sim(object = list(marg = "bbinom", par= par ,
                    tau = tau),  nsim = nsim, seed = seed,  cop = family)
}


#' @rdname sim_gctsc
#' @export
sim_zib <- function(prob, pi0, size, tau, arma_order, nsim = 100, family = c("gaussian","t"),
                    df = NULL, seed = NULL) {
  family <- match.arg(family)
  .check_common(nsim, tau, arma_order, seed, family , df,  "sim_zib")
  .check_prob(prob, nsim, "sim_zib")
  .check_pi0(pi0, nsim, "sim_zib")
  .check_size(size, "sim_zib", scalar_only = TRUE)
  
  par = list(prob = prob, pi0 = pi0, size = size, arma_order = arma_order)
  if (family == "t") par$df <- df
  
  sim(object = list(marg = "zib",
                    par = par,
                    tau = tau),nsim = nsim, seed = seed,  cop = family)
}

#' @rdname sim_gctsc
#' @export
sim_zibb <- function(prob, rho, pi0, size, tau, arma_order, nsim = 100, 
                     family = c("gaussian", "t"),  df = NULL, seed = NULL) {
  family <- match.arg(family)
  .check_common(nsim, tau, arma_order, seed,  family , df, "sim_zibb")
  .check_prob(prob, nsim, "sim_zibb")
  .check_rho(rho, nsim, "sim_zibb")
  .check_pi0(pi0, nsim, "sim_zibb")
  .check_size(size, "sim_zibb", scalar_only = TRUE)
  
  par = list(prob = prob, rho = rho, pi0 = pi0, size = size, 
             arma_order = arma_order)
  if (family == "t") par$df <- df
  
  sim(object = list(marg = "zibb",
                    par = par,
                    tau = tau), nsim = nsim, seed = seed,  cop = family)
}


#' Simulate data from a gctsc or specification list
#'
#' @keywords internal
#' @noRd
sim <- function(object, nsim = 100, cop , seed = NULL, ...) {
  if (missing(object)) {
    stop("Argument 'object' is required.")
  }
  if (!is.null(seed)) set.seed(seed)
  
  # Validate components
  if (!is.list(object) || is.null(object$marg) || is.null(object$par) || is.null(object$tau)) {
    stop("object must be a list with components 'marg', 'par', and 'tau'.")
  }
  
  marg <- object$marg
  par  <- object$par
  tau  <- object$tau
  
  
  if (is.null(par$arma_order)) stop("'arma_order' must be provided in object$par.")
  if (marg %in% c("poisson","negbin","zip") && is.null(par$mu)) stop("'mu' must be provided in object$par.")
  if (marg %in% c("binom","zib","bbinom","zibb") && is.null(par$prob)) stop("'prob' must be provided in object$par.")
  
  mu         <- par$mu
  prob       <- par$prob
  od         <- par$arma_order
  dispersion <- par$dispersion
  rho        <- par$rho
  size   <- par$size
  pi0   <- par$pi0
  
  if (length(tau) != sum(od)) {
    stop("Length of 'tau' must equal sum of AR and MA orders: length(tau) = ",
         length(tau), ", expected = ", sum(od), ".")
  }
  
  if (all(od == 0)) {
    stop("ARMA(0,0) (white noise) is not supported.")
  }
  
  # Convert ARMA params
  p <- od[1]
  q <- od[2]
  iar <- if (p > 0) 1:p else NULL
  ima <- if (q > 0) (p + 1):(p + q) else NULL
  phi <- if (length(iar)) tau[iar] else numeric(0)
  theta <- if (length(ima)) tau[ima] else numeric(0)
  tau_list <- list(phi = phi, theta = theta)
  sigma2 <- 1 / sum(ma.inf(tau_list)^2)
  
  # Simulate latent process
  z <- switch(cop,
              "gaussian" = gau_latent(tau_list, sigma2, nsim),
              "t"        = t_latent(tau_list, sigma2, nsim, df = par$df))
  u <- if (cop == "t") pt(z, df = par$df) else pnorm(z)
  
  # Parameters for marginals
  
  lambda <- switch(marg,
                   "poisson" = list(mu = mu),
                   "zip"     = list(mu = mu,pi0 = pi0),
                   "negbin"  = list(mu = mu, dispersion = dispersion),
                   "binom"   = list(prob = prob, size = size),
                   "zib"     = list(prob = prob, size = size, pi0 = pi0),
                   "bbinom"  = list(prob = prob, rho = rho, size = size),
                   "zibb"    = list(prob = prob, rho = rho, size = size, pi0 = pi0),
                   stop("Invalid marginal type: ", marg)
  )
  
  y <- simulate_marginal(marg, u, lambda)
  
  return(list(
    y = y,
    z = z,
    marginal = marg,
    parameters = lambda,
    cormat = list(arma_order = od, tau = tau)
  ))
}



#' Simulate latent ARMA Gaussian process
#'
#' @name sim_latent
#' @keywords internal
#' @noRd
gau_latent <- function(a, sigma2, n = 100) {
  p <- length(a$phi)
  q <- length(a$theta)
  burnin <- 1000
  z <- arima.sim(
    model = list(
      ar = if (p > 0) a$phi else NULL,
      ma = if (q > 0) a$theta else NULL
    ),
    n  = burnin + n,
    sd = sqrt(sigma2)   # sd is standard deviation
  )
  z <- as.numeric(z[(burnin + 1):(burnin + n)])
  return(z)
}



#' Simulate latent ARMA t process
#'
#' @name sim_latent
#' @keywords internal
#' @noRd
t_latent <- function(a, sigma2, n = 100, df) {
  stopifnot(df > 0)
  burnin <- 1000
  p <- length(a$phi)
  q <- length(a$theta)
  
  z <- arima.sim(
    model = list(
      ar = if (p > 0) a$phi else NULL,
      ma = if (q > 0) a$theta else NULL
    ),
    n  = burnin + n,
    sd = sqrt(sigma2)   # sd is standard deviation
  )
  z <- as.numeric(z[(burnin + 1):(burnin + n)])
  
  S <- rchisq(1, df = df)
  v <- z / sqrt(S / df)
  
  v
}



#' Inverse CDF simulation for copula-based marginals
#'
#' @name simulate_marginal
#' @keywords internal
#' @noRd
simulate_marginal <- function(marg, u, lambda) {
  marg <- match.arg(marg, choices = c("poisson", "negbin", "zip", "binom", "zib", "bbinom", "zibb"))
  
  switch(marg,
         "poisson" = {
           mu <- lambda$mu
           qpois(u, lambda = mu)
         },
         "negbin" = {
           mu <- lambda$mu
           size <- 1 / lambda$dispersion
           qnbinom(u, mu = mu, size = size)
         },
         "zip" = {
           mu <- lambda$mu
           pi0 <- lambda$pi0
           v <- (u - pi0) / (1 - pi0)    #remove the effect of u <= pi0
           v <- pmax(0, pmin(1, v))
           ifelse(u <= pi0, 0, qpois(v, lambda = mu))
         }
         ,
         "binom" = {
           prob <- lambda$prob
           size <- lambda$size
           qbinom(u, size, prob)
         }
         ,
         "zib" = {
           size <- lambda$size
           prob <- lambda$prob
           pi0 <- lambda$pi0
           v <- (u - pi0) / (1 - pi0)
           v <- pmax(0, pmin(1, v))
           ifelse(u <= pi0, 0, qbinom(v, size,prob))
         }
         ,
         "bbinom" = {
           if (!requireNamespace("VGAM", quietly = TRUE)) {
             stop("Please install the 'VGAM' package for beta-binomial simulation.")
           }
           prob <- .recyclen(lambda$prob, n, "prob")
           rho <- lambda$rho
           size <- lambda$size
           qbbinom_custom(u, size, prob, rho)
         },
         "zibb" = {
           if (!requireNamespace("VGAM", quietly = TRUE)) {
             stop("Please install the 'VGAM' package for zero inflated beta-binomial simulation.")
           }
           n <- length(u)
           size <- lambda$size
           prob <- .recyclen(lambda$prob, n, "prob")
           rho <- lambda$rho
           pi0 <- lambda$pi0
           v <- (u - pi0) / (1 - pi0)
           v <- pmax(0, pmin(1, v))
           ifelse(u <= pi0, 0, qbbinom_custom(v, size,prob, rho))
         },
         stop("Unknown marginal distribution: ", marg)
  )
}


#' Quantile function for beta-binomial using inversion
#'
#' @name qbbinom_custom
#' @keywords internal
#' @noRd
qbbinom_custom <- function(p, size, prob, rho) {
  alpha <- prob * (1 - rho) / rho
  beta <- (1 - prob) * (1 - rho) / rho
  
  sapply(seq_along(p), function(i) {
    support <- 0:size
    cdf_vals <- VGAM::pbetabinom.ab(support, size = size,
                                    shape1 = alpha[i], shape2 = beta[i])
    idx <- which(cdf_vals >= p[i])[1]
    if (is.na(idx)) size else support[idx]
  })
}





