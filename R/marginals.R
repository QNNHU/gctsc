#' Marginal Models for Copula Count Time Series
#'
#' @name marginal.gctsc
#' @aliases
#' poisson.marg
#' negbin.marg
#' binom.marg
#' bbinom.marg
#' zip.marg
#' zib.marg
#' zibb.marg
#'
#' @title Marginal Model Constructors for \pkg{gctsc}
#'
#' @description
#' These functions construct marginal model objects for use with
#' \code{\link{gctsc}}. Each constructor returns an object of class
#' \code{"marginal.gctsc"} containing the information needed to initialize
#' marginal parameters and compute the latent truncation bounds used in the
#' copula likelihood.
#'
#' The following marginal families are currently supported:
#' \itemize{
#'   \item Poisson: \code{poisson.marg()}
#'   \item Negative Binomial: \code{negbin.marg()}
#'   \item Binomial: \code{binom.marg()}
#'   \item Beta--Binomial: \code{bbinom.marg()}
#'   \item Zero-Inflated Poisson: \code{zip.marg()}
#'   \item Zero-Inflated Binomial: \code{zib.marg()}
#'   \item Zero-Inflated Beta--Binomial: \code{zibb.marg()}
#' }
#'
#' Supported link functions depend on the marginal family:
#' \itemize{
#'   \item \code{poisson.marg()}, \code{zip.marg()}, and
#'         \code{negbin.marg()} support \code{"identity"} and \code{"log"}.
#'   \item \code{binom.marg()}, \code{bbinom.marg()},
#'         \code{zib.marg()}, and \code{zibb.marg()} currently support
#'         \code{"logit"} only.
#' }
#'
#' @param link Link function used for the main marginal component.
#'   Supported links depend on the marginal family; see Description.
#' @param size Number of trials for Binomial, Zero-Inflated Binomial,
#'   Beta--Binomial, and Zero-Inflated Beta--Binomial marginals. For these
#'   marginals, \code{size} is currently assumed to be a single fixed positive
#'   integer. For Beta--Binomial and Zero-Inflated Beta--Binomial marginals,
#'   \code{size} should be greater than 1.
#' @param lambda.lower Optional lower bounds on the marginal parameters.
#' @param lambda.upper Optional upper bounds on the marginal parameters.
#'
#' @details
#' The marginal constructors are designed to be supplied to \code{\link{gctsc}},
#' rather than called during likelihood evaluation by the user. Each returned
#' marginal object contains internal functions for:
#' \itemize{
#'   \item computing starting values for the marginal parameters;
#'   \item determining the number of marginal parameters;
#'   \item converting observed counts into lower and upper latent truncation
#'         bounds for the copula likelihood.
#' }
#'
#' Internally, the design input \code{x} is represented as a named list of
#' design matrices. For non-zero-inflated marginals, \code{x} contains
#' \code{x$mu}. For zero-inflated marginals, \code{x} contains both
#' \code{x$mu} and \code{x$pi0}. These matrices are constructed automatically
#' by \code{\link{gctsc}} from the supplied formula and data.
#'
#' For zero-inflated marginals, the \code{mu} component controls the main
#' count distribution, while the \code{pi0} component controls the structural
#' zero probability through a logit link.
#'
#' The optional bounds \code{lambda.lower} and \code{lambda.upper} are attached
#' to the starting values and used during numerical optimization.
#'
#' @return
#' A marginal model object of class \code{"marginal.gctsc"}.
#'
#' @seealso
#' \code{\link{gctsc}}, \code{\link{arma.cormat}}
#'
#' @examples
#' poisson.marg(link = "log")
#' negbin.marg(link = "log")
#' binom.marg(link = "logit", size = 10)
#' bbinom.marg(link = "logit", size = 24)
#' zip.marg(link = "log")
#' zib.marg(link = "logit", size = 10)
#' zibb.marg(link = "logit", size = 24)
#' @rdname marginal.gctsc
#'@export
#' @usage poisson.marg(link = "log", lambda.lower = NULL, lambda.upper = NULL)
poisson.marg <- function(link = "log", lambda.lower = NULL, lambda.upper = NULL) {


  invlink <- switch(link,
                    "identity" = function(eta) eta,
                    "log" = exp,
                    stop("Unsupported link function")
  )


  obj <- list(
    name = "poisson",
    
    zero_inflated = FALSE,
    
    start = function(y, x) {
      x  <- x$mu
      fit <- glm.fit(x = x, y = y,family = poisson(link = link))
      
      lambda <- fit$coefficients
      if (anyNA(lambda)) stop("NA detected in coefficient estimates. Check design matrix.",call. = FALSE)

      names(lambda) <- prefixed_names(x, "mu_")
      lambda <- attach_bounds(lambda, lambda.lower, lambda.upper, "poisson.marg()")
      lambda
    },

    npar = function(x) NCOL(x$mu),

    bounds = function(y, x, lambda,family = "gaussian", df=NULL) {
      x  <- x$mu
      mu <- invlink(x %*% lambda)
      if (any(mu < 0)) stop("Negative mean detected. Use 'log' link or check predictors.",call. = FALSE)
      pdf <- dpois(y, mu)
      cdf <- ppois(y, mu)
      bounds <- safe_cdf_bounds(pdf, cdf, family, df)
      cbind(bounds$lower, bounds$upper)
    }
    


  )
  class(obj) <- "marginal.gctsc"
  obj
}


#' Binomial marginal model (supports y as vector or cbind(success, failure))
#'
#' @rdname marginal.gctsc
#'@export
#' @usage binom.marg(link = "logit", size = NULL, lambda.lower = NULL, lambda.upper = NULL)
binom.marg <- function(link = "logit", size = NULL, lambda.lower = NULL, lambda.upper = NULL) {
  if (missing(size) || is.null(size)) {
    stop("binom.marg(): 'size' must be provided.", call. = FALSE)
  }
  
  ok_size <- is.numeric(size) && length(size) == 1L && is.finite(size) && size > 0 &&
    size == as.integer(size)
  
  if (!ok_size) {
    stop(sprintf("%s(): 'size' must be a single positive integer.", "binom.marg"),
         call. = FALSE)
  }
  
  size <- as.integer(size)
  
  if (link != "logit") {
    stop("binom.marg(): currently only link = 'logit' is supported.",
         call. = FALSE)
  }
  
  
  invlink <- plogis
                   

  obj <- list(
    name = "binom",
    zero_inflated = FALSE,
    
    start = function(y, x) {
      x  <- x$mu
      size_vec  <- rep(size, length(y))
      
      if (any(y < 0 | y > size_vec)) {
        stop("binom.marg(): 'y' must be between 0 and 'size'.",
             call. = FALSE)
      }
      
      fit <- glm.fit( x = x, y = cbind(y, size_vec  - y),family = binomial(link = link))
      
      lambda <- coef(fit)
      
      if (anyNA(lambda)) {
        stop("NA in coefficient estimates. Check covariates for collinearity or data issues.",
             call. = FALSE)
      }
      
      names(lambda) <- prefixed_names(x, "mu_")
      lambda <- attach_bounds(lambda, lambda.lower, lambda.upper, "binom.marg()")
      lambda
    },

    npar = function(x) NCOL(x$mu),

    bounds = function(y, x, lambda,family = "gaussian", df=NULL) {
      x  <- x$mu
      size_vec <- rep(size, length(y))
      successes <- y
      mu <- invlink(x %*% lambda)
      pdf <- dbinom(successes, size_vec, mu)
      cdf <- pbinom(successes, size_vec, mu)
      bounds <- safe_cdf_bounds(pdf, cdf,family, df)
      cbind(bounds$lower, bounds$upper)
    }
  )

  class(obj) <- "marginal.gctsc"
  obj
}



#' Zero-Inflated Binomial marginal model
#'
#' @rdname marginal.gctsc
#'@export
#' @usage zib.marg(link = "logit", size = NULL, lambda.lower = NULL, lambda.upper = NULL)
zib.marg <- function(link = "logit", size = NULL, lambda.lower = NULL, lambda.upper = NULL) {
  if (missing(size) || is.null(size)) {
    stop("zib.marg(): 'size' must be provided.", call. = FALSE)
  }
  
  ok_size <- is.numeric(size) && length(size) == 1L && is.finite(size) && size > 0 &&
    size == as.integer(size)
  
  if (!ok_size) {
    stop(sprintf("%s(): 'size' must be a single positive integer.", "zib.marg"),
         call. = FALSE)
  }
  
  size <- as.integer(size)
  
  if (link != "logit") {
    stop("zib.marg(): currently only link = 'logit' is supported.",
         call. = FALSE)
  }
  
  
  invlink <- plogis
  
  obj <- list(
    name = "zib",
    zero_inflated = TRUE,
    
    start = function(y, x) {
      X_mu  <- x$mu
      X_pi0 <- x$pi0
      
      size_vec <- rep(size, length(y))
      
      if (any(y < 0 | y > size_vec)) {
        stop("zib.marg(): 'y' must be between 0 and 'size'.",
             call. = FALSE)
      }
      
      # Fit binomial only on non-zero observations
      nonzero_idx <- (y > 0)
      if (any(nonzero_idx)) {
        fit_mu <- glm.fit(x = X_mu[nonzero_idx, , drop = FALSE],
            y = cbind(y[nonzero_idx], size_vec[nonzero_idx] - y[nonzero_idx]),
            family = binomial(link = link))
        beta <- coef(fit_mu)
      } else {
        beta <- rep(0, ncol(X_mu))  # fallback
      }
      
      if (anyNA(beta)) {
        stop(sprintf(
          "%s(): NA in mean-component coefficient estimates. Check covariates for collinearity or data issues.",
          "zib.marg"
        ), call. = FALSE)
      }
      
      prob <- as.vector(invlink(X_mu %*% beta))
      
      # Fit zero-inflation model (logit link)
      f0 <- dbinom(0, size_vec, prob )
      
      # Crude estimate of structural zero probability
      pi0_est <- pmax((as.numeric(y == 0) - f0) / (1 - f0), 1e-6)
      
      # Fit alpha on this adjusted probability
      
      fit_pi0 <- glm.fit(x = X_pi0, y = pi0_est, family = quasibinomial(link = "logit"))
      alpha <- coef(fit_pi0)
      
      if (anyNA(alpha)) {
        stop(sprintf(
          "%s(): NA in zero-inflation coefficient estimates. Check covariates for collinearity or data issues.",
          "zib.marg"
        ), call. = FALSE)
      }
      
      lambda <- c(beta, alpha)
      names(lambda) <- c(prefixed_names(X_mu, "mu_"),
                         prefixed_names(X_pi0, "pi0_"))
      
      lambda <- attach_bounds(lambda, lambda.lower, lambda.upper, "zib.marg()")
      lambda
      },

    npar = function(x) {
      if (!is.list(x) || is.null(x$mu) || is.null(x$pi0)) {
        stop("x must be a list with elements 'mu' and 'pi0'")
      }
      ncol(x$mu) + ncol(x$pi0)
    },

    bounds = function(y, x, lambda, family = "gaussian", df=NULL) {
      if (!is.list(x)) stop("x must be a list with 'mu' and 'pi0'")
      X_mu <- x$mu
      X_pi0 <- x$pi0
      p_mu <- ncol(X_mu)
      beta <- lambda[1:p_mu]
      alpha <- lambda[(p_mu + 1):length(lambda)]
      
      size_vec <- rep(size, length(y))

      prob <- as.vector(invlink(X_mu %*% beta))
      pi0 <- plogis(X_pi0 %*% alpha)
      
      f0 <- dbinom(0, size_vec, prob)
      pmf_binom <- dbinom(y, size_vec, prob)
      cdf_binom <- pbinom(y, size_vec, prob)

      pdf <- ifelse(y == 0,
                    pi0 + (1 - pi0) * f0,
                    (1 - pi0) * pmf_binom)
      cdf <- pi0 + (1 - pi0) * cdf_binom

      bds <- safe_cdf_bounds(pdf, cdf, family, df)
      cbind(bds$lower, bds$upper)
    }
  )

  class(obj) <- "marginal.gctsc"
  obj
}


#' Negative binomial marginal model
#'
#' @rdname marginal.gctsc
#'@export
#' @usage negbin.marg(link = "log", lambda.lower = NULL, lambda.upper = NULL)
negbin.marg <- function(link = "log" ,lambda.lower = NULL, lambda.upper = NULL) {
  invlink <- switch(link,
                    "identity" = function(eta) eta,
                    "log" = exp,
                    stop("Unsupported link function")
  )

  obj <- list(
    name = "negbin",
    zero_inflated = FALSE,
    
    start = function(y, x) {
      x  <- x$mu
      eps <- sqrt(.Machine$double.eps)
      fit <- glm.fit(x, y, family = poisson(link=link))
      
      beta <- coef(fit)
      mu <- as.vector(fitted(fit))
      
      if (anyNA(beta)) {
        stop("negbin.marg(): NA in coefficient estimates. Check covariates for collinearity or data issues.",
             call. = FALSE)
      }
      
      if (any(mu <= 0)) {
        stop("negbin.marg(): nonpositive fitted mean detected when estimating dispersion.",
             call. = FALSE)
      }
      

      dispersion <- mean(((y - mu)^2 - mu) / mu^2, na.rm = TRUE)
      dispersion <- max(10 * eps, dispersion)
      lambda <- c(beta, dispersion = dispersion)
      
      names(lambda) <- c(prefixed_names(x, "mu_"),"dispersion")

      # Check bounds
      lambda <- attach_bounds(lambda, lambda.lower, lambda.upper, "negbin.marg()")
      lambda
    }
    ,

    npar = function(x) NCOL(x$mu) + 1 ,

    bounds = function(y, x, lambda,family = "gaussian", df=NULL) {
      x  <- x$mu
      beta <- lambda[1:NCOL(x)]
      dispersion <- lambda[length(lambda)]
      
      mu <- as.vector(invlink(x %*% beta))
      
      if (any(mu < 0)) {
        stop("negbin.marg(): mean must be positive. Consider using link = 'log'.",
             call. = FALSE)
      }
      
      size <- 1 / dispersion
      pdf <- dnbinom(y, mu = mu, size = size)
      cdf <- pnbinom(y, mu = mu, size = size)
      bounds <- safe_cdf_bounds(pdf, cdf,family = family, df = df)
      cbind(bounds$lower, bounds$upper)
    }
    

  )
  class(obj) <- "marginal.gctsc"
  obj
}

#' Zero-Inflated Poisson marginal model
#'
#' @rdname marginal.gctsc
#'@export
#' @usage zip.marg(link = "log", lambda.lower = NULL, lambda.upper = NULL)
zip.marg <- function(link = "log", lambda.lower = NULL, lambda.upper = NULL) {
  invlink <- switch(link,
                    "identity" = function(eta) eta,
                    "log" = exp,
                    stop("Unsupported link function")
  )

  obj <- list(
    name = "zip",
    zero_inflated = TRUE,
    
    
    start = function(y, x) {

      X_mu <- x$mu
      X_pi0 <- x$pi0
      if (!is.numeric(y)) stop("y must be numeric.")

      # Crude starting value: estimate Poisson mean using positive observations
      # to reduce the influence of structural zeros.
      nonzero_idx <- y > 0
      
      if (any(nonzero_idx)) {
        fit_mu <- glm.fit(
          x = X_mu[nonzero_idx, , drop = FALSE], y = y[nonzero_idx],
          family = poisson(link = link)
        )
        
        beta <- coef(fit_mu)
      } else {
        beta <- rep(0, ncol(X_mu))
      }
      
      if (anyNA(beta)) {
        stop(sprintf(
          "%s(): NA in mean-component coefficient estimates. Check covariates for collinearity or data issues.",
          "zip.marg"
        ), call. = FALSE)
      }
      
      mu <- as.vector(invlink(X_mu %*% beta))
      
      if (any(mu <= 0)) {
        stop(sprintf(
          "%s(): nonpositive fitted mean detected. Consider using link = 'log'.",
          "zip.marg"
        ), call. = FALSE)
      }
      
      # Crude estimate of structural-zero probability
      f0 <- dpois(0, lambda = mu)
      
      pi0_est <- (as.numeric(y == 0) - f0) /
        pmax(1 - f0, .Machine$double.eps)
      
      pi0_est <- pmin(pmax(pi0_est, 1e-6), 1 - 1e-6)
      
      fit_pi0 <- glm.fit( x = X_pi0,y = pi0_est, family = quasibinomial(link = "logit") )
      
      alpha <- coef(fit_pi0)
      
      if (anyNA(alpha)) {
        stop(sprintf(
          "%s(): NA in zero-inflation coefficient estimates. Check covariates for collinearity or data issues.",
          "zip.marg"
        ), call. = FALSE)
      }
      
      lambda <- c(beta, alpha)
      
      names(lambda) <- c(
        prefixed_names(X_mu, "mu_"),
        prefixed_names(X_pi0, "pi0_")
      )
      
      lambda <- attach_bounds(lambda, lambda.lower, lambda.upper, "zip.marg()")
      lambda
    },

    npar = function(x) {
      if (!is.list(x) || is.null(x$mu) || is.null(x$pi0)) {
        stop("x must be a list with elements 'mu' and 'pi0'")
      }
      ncol(x$mu) + ncol(x$pi0)
    },

    bounds = function(y, x, lambda,family = "gaussian", df=NULL) {
      X_mu <- x$mu
      X_pi0 <- x$pi0
      p_mu <- ncol(X_mu)
      
      beta  <- lambda[seq_len(p_mu)]
      alpha <- lambda[(p_mu + 1L):length(lambda)]
      
      mu  <- as.vector(invlink(X_mu %*% beta))
      pi0 <- as.vector(plogis(X_pi0 %*% alpha))

      
      if (any(mu < 0)) {
        stop(sprintf("%s(): mean must be positive. Consider using link = 'log'.", "zip.marg"),
             call. = FALSE)
      }
      
      
      f0 <- dpois(0, mu)
      pmf_pois <- dpois(y, mu)
      cdf_pois <- ppois(y, mu)

      pdf <- ifelse(y == 0,
                    pi0 + (1 - pi0) * f0,
                    (1 - pi0) * pmf_pois)
      cdf <- pi0 + (1 - pi0) * cdf_pois

      bds <- safe_cdf_bounds(pdf, cdf, family, df)
      cbind(bds$lower, bds$upper)
      }
  )

  class(obj) <- "marginal.gctsc"
  obj
}


#' Beta-Binomial marginal model
#'
#' @rdname marginal.gctsc
#'@export
#' @usage bbinom.marg(link = "logit", size, lambda.lower = NULL, lambda.upper = NULL)
bbinom.marg <- function(link = "logit", size, lambda.lower = NULL, lambda.upper = NULL) {
  
  if (missing(size) || is.null(size)) {
    stop("bbinom.marg(): 'size' must be provided.", call. = FALSE)
  }
  
  ok_size <- is.numeric(size) && length(size) == 1L && is.finite(size) && size > 0 &&
    size == as.integer(size)
  
  if (!ok_size) {
    stop(sprintf("%s(): 'size' must be a single positive integer.", "bbinom.marg"),
         call. = FALSE)
  }
  
  size <- as.integer(size)
  
  if (link != "logit") {
    stop("bbinom.marg(): currently only link = 'logit' is supported.",
         call. = FALSE)
  }
  
  
  invlink <- plogis

  obj <- list(
    name = "bbinom",
    zero_inflated = FALSE,

    
    start = function(y, x) {
      x <- x$mu

      size_vec <- rep(size, length(y))
      
      if (any(y < 0 | y > size_vec)) {
        stop(sprintf("%s(): 'y' must be between 0 and 'size'.", "bbinom.marg"),
             call. = FALSE)
      }
      
      if (any(abs(y - round(y)) > sqrt(.Machine$double.eps))) {
        stop(sprintf("%s(): 'y' must contain integer-valued counts.", "bbinom.marg"),
             call. = FALSE)
      }
      
      fit <- glm.fit( x = x, y = cbind(y, size_vec - y), family = binomial(link = link))
      
      beta <- coef(fit)
      
      if (anyNA(beta)) {
        stop(sprintf(
          "%s(): NA in mean-component coefficient estimates. Check covariates for collinearity or data issues.",
          "bbinom.marg"
        ), call. = FALSE)
      }
      
      prob_hat <- as.vector(invlink(x %*% beta))
      
      # Crude moment estimate of beta-binomial overdispersion.
      # Var(Y_t) = m p_t(1-p_t) {1 + (m - 1) rho}
      var_binom <- size_vec * prob_hat * (1 - prob_hat)
      resid2 <- (y - size_vec * prob_hat)^2
      
      rho_hat <- mean((resid2 / pmax(var_binom, .Machine$double.eps) - 1) /
                        (size_vec - 1), na.rm = TRUE)
      
      rho_hat <- pmin(pmax(rho_hat, 1e-6), 1 - 1e-6)
      rho_logit <- qlogis(rho_hat)
      
      lambda <- c(beta, rho_logit)
      
      names(lambda) <- c(prefixed_names(x, "mu_"),"logit_rho")
      
      lambda <- attach_bounds(lambda, lambda.lower, lambda.upper, "bbinom.marg()")
      lambda
    },

    npar = function(x) NCOL(x$mu) + 1L,  # slopes + rho

    bounds = function(y, x, lambda,family = "gaussian", df=NULL) {
      x <- x$mu
      p_mu <- NCOL(x)
      beta <- lambda[seq_len(p_mu)]
      rho_logit <- lambda[p_mu + 1L]
      
      prob <- as.vector(invlink(x %*% beta))
      rho <- plogis(rho_logit)

      kap1 <- prob * (1 - rho) / rho
      kap2 <- (1 - prob) * (1 - rho) / rho

      pdf <- VGAM::dbetabinom.ab(y, size = size, shape1 = kap1, shape2 = kap2, log = FALSE)
      cdf <- VGAM::pbetabinom.ab(y, size = size, shape1 = kap1, shape2 = kap2, log = FALSE)

      bounds <- safe_cdf_bounds(pdf, cdf,family, df)
      cbind(bounds$lower, bounds$upper)
    }
  )

  class(obj) <- "marginal.gctsc"
  obj
}



#' Zero-Inflated Beta-Binomial with separate seasonal structure for pi0 and BB mean
#'
#' @rdname marginal.gctsc
#'@export
#' @usage zibb.marg(link = "logit", size,  lambda.lower = NULL, lambda.upper = NULL)
zibb.marg <- function(link = "logit", size, lambda.lower = NULL, lambda.upper = NULL) {

  if (missing(size) || is.null(size)) {
    stop("zibb.marg(): 'size' must be provided.", call. = FALSE)
  }
  
  ok_size <- is.numeric(size) && length(size) == 1L && is.finite(size) && size > 1 &&
    size == as.integer(size)
  
  if (!ok_size) {
    stop("zibb.marg(): 'size' must be a single integer greater than 1.", call. = FALSE)
  }
  
  size <- as.integer(size)
  
  if (link != "logit") {
    stop("zibb.marg(): currently only link = 'logit' is supported.",
         call. = FALSE)
  }
  
  
  invlink <- plogis
  
  obj <- list(
    name = "zibb",
    zero_inflated = TRUE,

    start = function(y, x) {
      X_mu <- x$mu
      X_pi0 <- x$pi0
      
      size_vec <- rep(size, length(y))
     
      
      if (any(y < 0 | y > size_vec)) {
        stop(sprintf("%s(): 'y' must be between 0 and 'size'.", "zibb.marg"),
             call. = FALSE)
      }
      
      if (any(abs(y - round(y)) > sqrt(.Machine$double.eps))) {
        stop(sprintf("%s(): 'y' must contain integer-valued counts.", "zibb.marg"),
             call. = FALSE)
      }
      
      nonzero_idx <- (y > 0)
      
      if (any(nonzero_idx)) {
        fit_mu <- glm.fit( x =X_mu[nonzero_idx, , drop = FALSE], 
                        y = cbind(y[nonzero_idx], size_vec[nonzero_idx]
                         - y[nonzero_idx]), family = binomial(link = link))
        beta <- coef(fit_mu)
        } else {
          beta <- rep(0, ncol(X_mu))
        }
        
      if (anyNA(beta)) {
        stop(sprintf(
          "%s(): NA in mean-component coefficient estimates. Check covariates for collinearity or data issues.","zibb.marg"
        ), call. = FALSE)
      }
      
      prob_hat <- as.vector(invlink(X_mu %*% beta))
      
      # Crude moment estimate of beta-binomial overdispersion.
      # Var(Y_t) = m p_t(1-p_t) {1 + (m - 1) rho}
      var_binom <- size_vec * prob_hat * (1 - prob_hat)
      resid2 <- (y - size_vec * prob_hat)^2
      
      rho_raw <- mean((resid2 / pmax(var_binom, .Machine$double.eps) - 1) /
                        (size_vec - 1), na.rm = TRUE)
      
      if (!is.finite(rho_raw) || rho_raw < 0 || rho_raw > 0.9) {
        rho_hat <- 0.2
      } else {
        rho_hat <- pmin(pmax(rho_raw, 1e-3), 0.9)
      }
      
      # rho_hat <- pmin(pmax(rho_hat, 1e-6), 1 - 1e-6)
      rho_logit <- qlogis(rho_hat)
    
      
      kap1 <- prob_hat * (1 - rho_hat) / rho_hat
      kap2 <- (1 - prob_hat) * (1 - rho_hat) / rho_hat
      
      f0 <- VGAM::dbetabinom.ab( 0,size = size_vec, shape1 = kap1,
        shape2 = kap2 )
      
      # Crude estimate of structural-zero probability.
      pi0_est <- (as.numeric(y == 0) - f0) /
        pmax(1 - f0, .Machine$double.eps)
      
      pi0_est <- pmin(pmax(pi0_est, 1e-6), 1 - 1e-6)
      
      
      # Fit logistic regression for π₀
      fit_pi0 <- glm.fit(X_pi0, y =pi0_est, family = quasibinomial(link = "logit"))
      alpha <- coef(fit_pi0)
      
      
      if (anyNA(alpha)) {
        stop(sprintf(
          "%s(): NA in zero-inflation coefficient estimates. Check covariates for collinearity or data issues.","zibb.marg"
        ), call. = FALSE)
      }
      
      
      lambda <- c(beta, rho_logit, alpha)
      names(lambda) <- c( prefixed_names(X_mu, "mu_"), "logit_rho",
                          prefixed_names(X_pi0, "pi0_") )

      # Attach bounds if provided
      lambda <- attach_bounds(lambda, lambda.lower, lambda.upper, "zibb.marg()")
      lambda
    }

    ,

    npar = function(x) {
      if (!is.list(x) || is.null(x$mu) || is.null(x$pi0)) {
        stop("x must be a list with elements 'mu' and 'pi0'")
      }
      ncol(x$mu)+1L + ncol(x$pi0)
    },

    bounds = function(y, x, lambda,family = "gaussian", df=NULL) {
      X_mu <- x$mu
      X_pi0 <-x$pi0
      
      p_mu <- ncol(X_mu)
      p_pi0 <- ncol(X_pi0)
      
      beta <- lambda[1:p_mu]
      rho_logit <-lambda[p_mu + 1L]
      alpha <- lambda[(p_mu + 2L):(p_mu + 1L + p_pi0)]


      prob <- as.vector(invlink(X_mu %*% beta))
      rho <- plogis(rho_logit)
      pi0 <- as.vector(plogis(X_pi0 %*% alpha))

      size_vec <- rep(size, length(y))
      kap1 <- prob * (1 - rho) / rho
      kap2 <- (1 - prob) * (1 - rho) / rho

      f0 <- VGAM::dbetabinom.ab(0, size = size_vec, shape1 = kap1, shape2 = kap2)
      pmf_bbinom <- VGAM::dbetabinom.ab(y, size = size_vec, shape1 = kap1, shape2 = kap2)
      cdf_bbinom <- VGAM::pbetabinom.ab(y, size = size_vec, shape1 = kap1, shape2 = kap2)

      pdf <- ifelse(y == 0, pi0 + (1 - pi0) * f0, (1 - pi0) * pmf_bbinom)
      cdf <- pi0 + (1 - pi0) * cdf_bbinom

      bounds <- safe_cdf_bounds(pdf, cdf,family, df)
      cbind(bounds$lower, bounds$upper)
    }
  )

  class(obj) <- "marginal.gctsc"
  obj
}




################ utilities##########################
attach_bounds <- function(lambda, lower, upper, where = "marginal") {
  if (!is.null(lower) && length(lower) != length(lambda)) {
    stop(sprintf("%s: lambda.lower must match length of coefficients (got %d vs %d).",
                 where, length(lower), length(lambda)), call. = FALSE)
  }
  if (!is.null(upper) && length(upper) != length(lambda)) {
    stop(sprintf("%s: lambda.upper must match length of coefficients (got %d vs %d).",
                 where, length(upper), length(lambda)), call. = FALSE)
  }
  if (!is.null(lower)) attr(lambda, "lower") <- lower
  if (!is.null(upper)) attr(lambda, "upper") <- upper
  lambda
}


get_qfun <- function(family, df) {
  switch(family,
         "t" = function(p) qt(p, df),
         "gaussian" = qnorm,
         stop("Unknown family in get_qfun.")
  )
}

safe_cdf_bounds <- function(pdf, cdf, family, df) {
  
  EPS <- sqrt(.Machine$double.eps)
  EPS1 <- 1 - EPS
  cdf_l <- pmax(EPS, pmin(EPS1, cdf - pdf))
  cdf_u <- pmax(EPS, pmin(EPS1, cdf))
  cdf_l <- pmin(cdf_l, pmax(0, cdf_u - EPS))
  qfun <- get_qfun(family, df)
  list(lower = qfun(cdf_l), upper = qfun(cdf_u))
}


prefixed_names <- function(X, prefix, base = "X") {
  cn <- colnames(X)
  
  if (is.null(cn)) {
    cn <- paste0(base, seq_len(ncol(X)))
  }
  
  cn <- gsub("^\\(Intercept\\)$", "intercept", cn)
  
  paste0(prefix, cn)
}

