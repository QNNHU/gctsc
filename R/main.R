#' Fit a Gaussian Copula Time Series Model for Count Data
#'
#' Fits a Gaussian copula model to univariate count time series using discrete
#' marginals (Poisson, negative binomial, binomial/beta–binomial, and zero–inflated
#' variants) and latent dependence via ARMA correlation structures. The
#' multivariate normal rectangle probability is evaluated using Minimax
#' Exponential Tilting (TMET), Geweke–Hajivassiliou–Keane (GHK), or Continuous
#' Extension (CE).
#'
#' The interface mirrors \code{glm()}. Zero–inflated marginals accept a list of
#' formulas, e.g., \code{list(mu = y ~ x, pi0 = ~ z)}. Non–zero–inflated
#' marginals accept a single formula or \code{list(mu = ...)}.
#'
#' @param formula A formula (e.g., \code{y ~ x1 + x2}) or, for zero–inflated
#'   marginals, a named list of formulas \code{list(mu = ..., pi0 = ...)}.
#' @param data A data frame containing \code{y} and covariates referenced in the formula(s).
#' @param marginal A marginal model object such as \code{\link{poisson.marg}},
#'   \code{\link{negbin.marg}}, \code{\link{zib.marg}}, or \code{\link{zibb.marg}}.
#' @param cormat A correlation structure such as \code{\link{arma.cormat}}.
#' @param method One of \code{"TMET"}, \code{"GHK"}, or \code{"CE"}.
#' @param c Smoothing/tilting constant for CE (ignored otherwise). Default \code{0.5}.
#' @param QMC Logical; use quasi–Monte Carlo for simulation–based methods.
#' @param pm Integer; truncated AR order used to approximate ARMA(\eqn{p,q}) when \eqn{q>0} (TMET only).
#' @param start Optional numeric vector of starting values (marginal then dependence parameters).
#' @param options Optional list of tuning/optimization controls. If \code{NULL},
#'   defaults from \code{gctsc.opts()} are used. Any fields supplied override the defaults.
#'
#' @details
#'
#' ### Formulas
#' For zero–inflated marginals, if neither \code{mu} nor \code{pi0} is supplied,
#' both default to intercept–only: \code{mu ~ 1}, \code{pi0 ~ 1}.  
#' If \code{mu} is supplied but \code{pi0} is missing, \code{pi0 ~ 1} is used.
#'
#' ### Dependence
#' The ARMA parameters are encoded in \code{cormat}; models must be
#' stationary/invertible. ARMA(0,0) is not supported.
#'
#' ### Method-specific notes
#' CE ignores \code{QMC} and \code{M}.  
#' GHK/TMET require \code{options$M} (a positive integer).  
#' TMET may also use \code{pm} when \eqn{q>0}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{coef} — parameter estimates,
#'   \item \code{maximum} — approximate log-likelihood,
#'   \item \code{se} — standard errors (if available),
#'   \item \code{terms}, \code{model}, \code{call} — model metadata.
#' }
#'
#' @examples
#' set.seed(42)
#' n <- 200
#' y <- sim_poisson(mu = 10, tau = 0.3, arma_order = c(1,0), nsim = n)$y
#' fit <- gctsc(y ~ 1,
#'   marginal = poisson.marg(lambda.lower = 0),
#'   cormat = arma.cormat(p = 1, q = 0),
#'   method = "CE",
#'   options = gctsc.opts(M = 1000, seed = 42))
#' summary(fit)
#'
#' @seealso \code{\link{arma.cormat}}, \code{\link{poisson.marg}},
#'   \code{\link{zib.marg}}, \code{\link{zibb.marg}}, \code{\link{gctsc.opts}}
#' @export
#'
gctsc <- function(formula=NULL, data, marginal, cormat,
                  method = c("TMET", "GHK", "CE","VMET"),
                  c = 0.5, QMC = TRUE, pm = 30, start = NULL,
                 options = gctsc.opts()) {

  method <- match.arg(method)
  .validate_method(method, "gctsc")
  objs <- .validate_marg_cormat(marginal, cormat, "gctsc")
  marginal <- objs$marginal; cormat <- objs$cormat
  formula  <- .validate_formula_input(formula, marginal, "gctsc")
  .validate_options(method, QMC, options, "gctsc")
  
  des <- .build_design(formula, data, marginal, "gctsc")
  y <- des$y; x <- des$x
  validate_x_structure(x, marginal, "gctsc")
  check_x_nrow_matches_y(x, y, marginal, "gctsc")
  
  # hand off
  fit <- gctsc.fit(x = x, y = y, marginal = marginal, cormat = cormat,
                   method = method, c = c, QMC = QMC, pm = pm,
                   start = start, options = options)
  
  fit$call <- match.call(expand.dots = FALSE)
  fit$formula <- formula
  fit$terms <- des$terms
  fit$model <- des$model
  class(fit) <- "gctsc"
  fit
}


#' Fit a Gaussian Copula Time Series Model (Internal)
#'
#' Internal workhorse called by \code{\link{gctsc}}. Validates inputs, builds
#' starting values and bounds from the marginal and correlation structures, and
#' maximizes the approximate log–likelihood for the chosen method.
#'
#' @inheritParams gctsc
#' @param x Design matrix (non–ZI) or list of design matrices \code{list(mu = X_mu, pi0 = X_pi0)} (ZI).
#' @param y Numeric response vector of non–negative integer counts.
#' @return A list with estimates, log–likelihood, (optionally) Hessian, and diagnostics.
#' @keywords internal
#' @seealso \code{\link{gctsc}}
#' @noRd
gctsc.fit <- function(x = NULL, y, marginal, cormat,
                      method = "GHK", c = 0.5, QMC = TRUE,
                      start = NULL, pm = 30, options = gctsc.opts()) {

  objs <- .validate_marg_cormat(marginal, cormat, "gctsc.fit")
  marginal <- objs$marginal; cormat <- objs$cormat
  .validate_method(method, "gctsc.fit")
  .validate_options(method, QMC, options, "gctsc.fit")
  
  if (is.null(x)) {
    if (has_ZI(marginal)) {
      x <- list(mu = matrix(1, nrow = length(y), ncol = 1L),
                pi0 = matrix(1, nrow = length(y), ncol = 1L))
    } else {
      x <- matrix(1, nrow = length(y), ncol = 1L)
    }
  }
  validate_x_structure(x, marginal, "gctsc.fit")
  check_x_nrow_matches_y(x, y, marginal, "gctsc.fit")
  
  # Missing handling
  x_mat <- if (has_ZI(marginal)) do.call(cbind, x) else x
  not.na <- rowSums(is.na(cbind(y, x_mat))) == 0
  if (!any(not.na)) stop("gctsc.fit(): Have NA after combining y and x.", call. = FALSE)
  if (sum(!not.na) > 0) warning(sprintf("gctsc.fit(): dropping %d row(s) with NA.", sum(!not.na)))
  
  y <- as.matrix(y)[not.na, , drop = FALSE]
  x <- if (has_ZI(marginal)) {
    lapply(x, function(col) as.matrix(col)[not.na, , drop = FALSE])
  } else {
    as.matrix(x)[not.na, , drop = FALSE]
  }
  
  nbeta <- marginal$npar(x)
  ntau  <- cormat$npar
  
  # Starting values & bounds (always read template attrs)
  beta_tmpl <- marginal$start(y, x)
  tau_tmpl  <- cormat$start(y)
  lb <- c(attr(beta_tmpl, "lower") %||% rep(-Inf, length(beta_tmpl)),
          attr(tau_tmpl,  "lower") %||% rep(-Inf, length(tau_tmpl)))
  ub <- c(attr(beta_tmpl, "upper") %||% rep( Inf, length(beta_tmpl)),
          attr(tau_tmpl,  "upper") %||% rep( Inf, length(tau_tmpl)))
  
  if (is.null(start)) {
    init_eta <- c(beta_tmpl, tau_tmpl)
  } else {
    if (!is.numeric(start) || length(start) != (nbeta + ntau) || any(!is.finite(start)))
      stop("gctsc.fit(): 'start' must be numeric, finite, length nbeta + ntau.", call. = FALSE)
    init_eta <- start
  }
  
  # Optional AR/MA admissibility check at initial tau
  if (ntau > 0 && !is.null(cormat$p) && !is.null(cormat$q)) {
    p <- cormat$p; q <- cormat$q
    tau_init <- tail(init_eta, ntau)
    .check_arima_admissibility(tau_init, p, q, "gctsc.fit")
  }
  
  f <- structure(list(
    y = y, x = x, c = c, n = sum(not.na), method = method,
    marginal = marginal, cormat = cormat,
    ibeta = 1:nbeta, itau = (nbeta + 1):(nbeta + ntau),
    nbeta = nbeta, ntau = ntau, QMC = QMC, pm = pm,
    call = match.call(), init_eta = init_eta, coef = init_eta,
    lower = lb, upper = ub, options = options
  ), class = "gctsc")
  
  gctsc.estimate(f)
}



#' @keywords internal
#' @noRd
gctsc.estimate <- function(cf) {

  if ( cf$method != "CE" ) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      seed.keep <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
    }
  }

  start <- cf$init_eta
  low <- cf$lower
  up <- cf$upper
  penalty <- -sqrt(.Machine$double.xmax)
  M <- cf$options$M
  log.lik <- build_loglik(cf,M, penalty)

  # saving/restoring the random seed (only for methods that need it)
  ans <- if (cf$method != "CE") {
    preserve_seed({
      suppressWarnings(cf$options$opt(start, log.lik, low, up))
    })
  } else {
    suppressWarnings(cf$options$opt(start, log.lik, low, up))
  }

  eta <- ans$estimate
  names(eta) <- names(cf$coef)
  cf$coef <- eta
  cf$maximum <- ans$maximum
  cf$convergence <- ans$convergence

  # Store Hessian if available
  if (!is.null(ans$hessian) && is.matrix(ans$hessian) && all(is.finite(ans$hessian))) {
    cf$hessian <- ans$hessian
    if (all(eigen(cf$hessian, symmetric = TRUE, only.values = TRUE)$values > 0)) {
      vcov <- try(solve(cf$hessian), silent = TRUE)
      if (!inherits(vcov, "try-error")) {
        cf$se <- sqrt(diag(vcov))
      }
    } else {
      warning("Hessian is not positive definite. Standard errors not computed.")
      cf$se <- rep(NA, length(cf$coef))
    }
  } else {
    warning("Hessian not available from optimization.")
    cf$se <- rep(NA, length(cf$coef))
  }

  return(cf)
}








