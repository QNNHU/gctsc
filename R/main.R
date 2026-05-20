#' Fit a Copula-Based Count Time Series Model
#'
#' Fits a Gaussian or Student--t copula model to a univariate count time
#' series with flexible discrete marginal distributions and latent ARMA
#' dependence.
#'
#' Supported marginal distributions include Poisson, negative binomial,
#' binomial, beta-binomial, and their zero-inflated variants. The latent
#' dependence is specified through a correlation model such as
#' \code{\link{arma.cormat}}.
#'
#' The copula likelihood involves a high-dimensional rectangle probability.
#' This probability is approximated using one of the following methods:
#' \itemize{
#'   \item \code{"TMET"}: Time Series Minimax Exponential Tilting,
#'   \item \code{"GHK"}: Geweke--Hajivassiliou--Keane simulation,
#'   \item \code{"CE"}: Continuous Extension,
#'   \item \code{"GHK_mvt"}: GHK approximation for the multivariate
#'         Student--t rectangle probability. This option is
#'         experimental and is mainly intended for comparisons.
#' }
#'
#' The model interface follows the usual R formula convention. An intercept
#' is included by default and can be removed using \code{-1} or \code{0 +}.
#' For non-zero-inflated marginals, \code{formula} may be a standard formula,
#' such as \code{y ~ x1 + x2}, or a named list \code{list(mu = y ~ x1 + x2)}.
#' For zero-inflated marginals, \code{formula} must be a named list with
#' components \code{mu} and \code{pi0}, for example
#' \code{list(mu = y ~ x1 + x2, pi0 = ~ z1)}.
#'
#' @param formula A model formula or a named list of formulas. For
#'   non-zero-inflated marginals, this may be a formula such as
#'   \code{y ~ x1 + x2} or \code{list(mu = y ~ x1 + x2)}. For zero-inflated
#'   marginals, this must be a named list with both \code{mu} and \code{pi0}
#'   components, e.g. \code{list(mu = y ~ x1, pi0 = ~ z1)}. The \code{mu}
#'   formula must include the response variable.
#'
#' @param data A data frame containing the response and all covariates
#'   referenced in \code{formula}. This argument is required.
#'
#' @param marginal A marginal model object inheriting class
#'   \code{"marginal.gctsc"}, such as \code{\link{poisson.marg}},
#'   \code{\link{negbin.marg}}, \code{\link{binom.marg}},
#'   \code{\link{bbinom.marg}}, \code{\link{zip.marg}},
#'   \code{\link{zib.marg}}, or \code{\link{zibb.marg}}.
#'
#' @param cormat A correlation model object inheriting class
#'   \code{"cormat.gctsc"}, such as \code{\link{arma.cormat}}.
#'
#' @param method Character string specifying the likelihood approximation
#'   method. One of \code{"TMET"}, \code{"GHK"}, \code{"CE"}, or
#'   \code{"GHK_mvt"}.
#'
#' @param c Numeric smoothing constant used by the CE method. Must be a
#'   single value between 0 and 1. Ignored by TMET and GHK methods.
#'
#' @param QMC Logical; if \code{TRUE}, quasi-Monte Carlo sampling is used
#'   for simulation-based methods.
#'
#' @param pm Positive integer specifying the truncated AR order used by TMET
#'   when approximating ARMA(\eqn{p,q}) dependence. This is mainly relevant
#'   when \eqn{q > 0}. Default is \code{30}.
#'
#' @param start Optional numeric vector of starting values. The vector should
#'   contain the marginal parameters followed by the dependence parameters.
#'   If \code{NULL}, starting values are constructed from the marginal and
#'   correlation objects.
#'
#' @param family Copula family. One of \code{"gaussian"} or \code{"t"}.
#'   Default is \code{"gaussian"}.
#'
#' @param df Degrees of freedom for the Student--t copula. Must be a single
#'   finite numeric value greater than 2 when \code{family = "t"}.
#'
#' @param options A list of computational options, usually created by
#'   \code{\link{gctsc.opts}}. Important components include:
#'   \itemize{
#'   \item \code{M}: positive integer or vector of two positive integers specifying
#'   the number of Monte Carlo or quasi-Monte Carlo samples used by
#'   simulation-based methods. If two values are supplied, staged optimization is
#'   performed, using the first value for an initial fit and the second value for
#'   refinement. This option is ignored for \code{method = "CE"}.
#'     \item \code{seed}: optional integer seed used to make the simulated
#'           likelihood approximation reproducible;
#'     \item \code{opt}: optimization function used for likelihood
#'           maximization.
#'   }
#'   Supplying \code{options$seed} is strongly recommended for TMET and GHK
#'   methods because it uses common random numbers across likelihood
#'   evaluations, making the objective function more stable and improving
#'   numerical Hessian estimation.
#'
#' @details
#' \strong{Formula interface.}
#' For non-zero-inflated marginals, users may write either
#' \code{formula = y ~ x1 + x2} or
#' \code{formula = list(mu = y ~ x1 + x2)}. Internally, both are represented
#' as a list with component \code{mu}.
#'
#' For zero-inflated marginals, users must supply both the mean/count
#' component and the zero-inflation component:
#' \preformatted{
#' formula = list(mu  = y ~ x1 + x2,pi0 = ~ z1 + z2)
#' }
#' The response variable is taken from \code{formula$mu}. The \code{pi0}
#' formula should be one-sided. Intercepts are handled by
#' \code{\link[stats]{model.matrix}} following standard R formula rules.
#'
#' \strong{Missing values.}
#' Missing values are not handled automatically. Users should remove or
#' impute missing values before calling \code{gctsc}. This avoids ambiguity
#' in the time series dependence structure.
#'
#' \strong{Dependence.}
#' The dependence parameters are encoded through \code{cormat}. ARMA(0,0)
#' is not supported. For ARMA dependence, admissible starting values should
#' satisfy the usual causality and invertibility conditions.
#'
#' \strong{Seed and numerical stability.}
#' Simulation-based likelihood approximations are random unless a seed is
#' supplied. If \code{options$seed} is provided, the same random stream is
#' used across likelihood evaluations, which can make optimization and
#' standard error estimation more stable. If no seed is supplied, the model
#' can still be fitted, but the approximate likelihood and numerical Hessian
#' may be less stable.
#'
#' @return An object of class \code{"gctsc"} containing, among others:
#' \itemize{
#'   \item \code{coef}: parameter estimates;
#'   \item \code{maximum}: approximate maximized log-likelihood;
#'   \item \code{se}: standard errors, when available;
#'   \item \code{formula}: normalized model formula list;
#'   \item \code{terms}: model terms for each component;
#'   \item \code{model}: model frames for each component;
#'   \item \code{call}: matched function call.
#' }
#'
#' @references
#' Nguyen, Q. N. and De Oliveira, V. (2026), Approximating Gaussian Copula
#' Models for Count Time Series: Connecting the Distributional Transform and
#' a Continuous Extension, \emph{Journal of Applied Statistics},
#' \strong{53}: 1--22.
#'
#' Nguyen, Q. N. and De Oliveira, V. (2026), Likelihood Inference in Gaussian
#' Copula Models for Count Time Series via Minimax Exponential Tilting,
#' \emph{Computational Statistics & Data Analysis}, \strong{218}: 108344.
#'
#' Nguyen, Q. N. and De Oliveira, V. (2026), Scalable Likelihood Inference
#' for Student--\eqn{t} Copula Count Time Series, \emph{Stats},
#' \strong{9}: 1--49.
#'
#' @examples
#' ## Example 1: Gaussian copula, Poisson marginal, AR(1)
#' set.seed(42)
#' n <- 500
#' sim_dat <- sim_poisson(mu = 10, tau = 0.3, arma_order = c(1, 0),
#'                        nsim = n, family = "gaussian")
#'
#' dat <- data.frame(y = sim_dat$y)
#'
#' fit_gauss <- gctsc(
#'   y ~ 1,
#'   data = dat,
#'   marginal = poisson.marg(lambda.lower = 0),
#'   cormat = arma.cormat(p = 1, q = 0), family = "gaussian",
#'   method = "CE",
#'   options = gctsc.opts(M = 1000, seed = 42)
#' )
#' summary(fit_gauss)
#'
#' ## Example 2: Student--t copula
#' sim_dat_t <- sim_poisson(mu = 10, tau = 0.3, arma_order = c(1, 0),
#'                          nsim = 500, family = "t", df = 10)
#'
#' dat_t <- data.frame(y = sim_dat_t$y)
#'
#' fit_t <- gctsc(
#'   y ~ 1,
#'   data = dat_t,
#'   marginal = poisson.marg(lambda.lower = 0),
#'   cormat = arma.cormat(p = 1, q = 0), family ="t",
#'   df= 10, method = "CE",
#'   options = gctsc.opts(M = 1000, seed = 42)
#' )
#' summary(fit_t)
#'
#'
#' @seealso \code{\link{gctsc.opts}}, \code{\link{arma.cormat}},
#'   \code{\link{poisson.marg}}, \code{\link{zip.marg}},
#'   \code{\link{zib.marg}}, \code{\link{zibb.marg}}
#'
#' @export
#'
gctsc <- function(formula=NULL, data, marginal, cormat,
                  method = c("TMET", "GHK", "CE","GHK_mvt"),
                  c = 0.5, QMC = TRUE, pm = 30, start = NULL,
                  family ="gaussian",df=10,
                  options = gctsc.opts()) {

  
  method <- match.arg(method, c("TMET", "GHK", "CE", "GHK_mvt"))
  
  if (!family %in% c("gaussian", "t")) {
    stop(sprintf("%s(): 'family' must be either 'gaussian' or 't'.", "gctsc"),
         call. = FALSE)
  }
  
  
  ## ---- t copula checks ----
  if (family == "t") {
    if (is.null(df))
      stop("For a Student-t copula, 'df' must be provided.")
    
    if (!is.numeric(df) || length(df) != 1 || !is.finite(df)){
      stop(sprintf("%s(): 'df' must be a single finite numeric value.",  "gctsc"),
           call. = FALSE)}
    
    if (df <= 2){
      stop(sprintf("%s(): 'df' must be greater than 2 for the t copula.",  "gctsc"),
           call. = FALSE)
    }
  }
  
  
  # checking marginal
  
  if (!inherits(marginal, "marginal.gctsc"))
    stop(sprintf("%s(): 'marginal' must be a 'marginal.gctsc' object.", "gctsc"), call. = FALSE)
  
  # checking correlation
  if (!inherits(cormat, "cormat.gctsc"))
    stop(sprintf("%s(): 'cormat' must be a 'cormat.gctsc' object (e.g., arma.cormat()).", "gctsc"), call. = FALSE)
  
  # checking formula
  if (missing(formula) || is.null(formula)) {
    stop(sprintf("%s(): 'formula' must be supplied.", "gctsc"), call. = FALSE)
  }
  
  # Checking for formula type
  if (inherits(formula, "formula")) {
    formula <- list(mu = formula)
  } else if (!is.list(formula)) {
    stop(sprintf("%s(): 'formula' must be a formula or a named list of formulas.", "gctsc"), call. = FALSE)
  }
  
  # checking for formula for zero inflated
  if (isTRUE(marginal$zero_inflated)) {
    
    # ZI needs mu and pi0; 
    if (is.null(formula$mu) || is.null(formula$pi0)) {
      stop(sprintf("%s(): for zero-inflated marginals, supply both formula$mu and formula$pi0, 
                   e.g., formula = list(mu = y ~ 1, pi0 = ~ 1).",  "gctsc"), call. = FALSE)
    } 
    
    if (!inherits(formula$mu,  "formula")) stop(sprintf("%s(): formula$mu must be a formula.",  "gctsc"), call. = FALSE)
    
    if (!inherits(formula$pi0, "formula")) stop(sprintf("%s(): formula$pi0 must be a formula.", "gctsc"), call. = FALSE)
    
  } else {
    # non-ZI: only mu
    if (is.null(formula$mu)) {
      if (length(formula) == 1L && inherits(formula[[1L]], "formula")) {formula$mu <- formula[[1L]]
      }else{
        stop(sprintf("%s(): provide a formula for non-zero-inflated marginals,
                     e.g., y ~ 1 or list(mu = y ~ 1).", "gctsc"), call. = FALSE)}
    }
    
    if (!inherits(formula$mu, "formula")) {stop(sprintf("%s(): must be a formula.", "gctsc"), call. = FALSE)}
  }
    
  
  
  # checking for data
  if (missing(data) || is.null(data)) {
    stop(sprintf("%s():'data' must be supplied.","gctsc"), call. = FALSE)
  }
  
  if (!is.data.frame(data)) {
    stop(sprintf("%s(): 'data' must be a data frame.", "gctsc"), call. = FALSE)
  }
  
  
  # === Build model frames ===
  if (isTRUE(marginal$zero_inflated)) {
    mf_mu  <- model.frame(formula$mu, data = data, na.action = na.pass)
    y      <- model.response(mf_mu)
    
    if (is.null(y)) {
      stop(sprintf("%s(): formula$mu must include a response, e.g., y ~ x1.", "gctsc"),
           call. = FALSE)
    }
    
    if (!is.numeric(y)) {
      stop(sprintf("%s(): response must be numeric counts.", "gctsc"), call. = FALSE)
    }
    
    
    mf_pi0 <- model.frame(formula$pi0, data = data, na.action = na.pass)
   
    X_mu   <- model.matrix(formula$mu,   data = mf_mu)
    X_pi0  <- model.matrix(formula$pi0, data = mf_pi0)
    
    x = list(mu = X_mu, pi0 = X_pi0)
    terms = list(mu = terms(formula$mu), pi0 = terms(formula$pi0))
    model = list(mu = mf_mu,pi0 = mf_pi0)
  } else {
    mf_mu <- model.frame(formula$mu, data = data, na.action = na.pass)
    y  <- model.response(mf_mu)
    if (is.null(y)) {
      stop(sprintf("%s(): formula must include a response, e.g., y ~ x1.", "gctsc"),
           call. = FALSE)
    }
    
    if (!is.numeric(y)) {
      stop(sprintf("%s(): response must be numeric counts.", "gctsc"), call. = FALSE)
    }
    
    X_mu <- model.matrix(formula$mu, data = mf_mu)
    
    x <- list(mu = X_mu)
    terms <- list(mu = terms(formula$mu))
    model = list(mu = mf_mu)
    
  }
  
  
  # Validate option
  
  if (is.null(options)) {
    options <- list()
  }
  
  if (is.null(options) || is.null(options$seed)) {
    message("gctsc(): 'options$seed' is strongly recommended for reproducible and stable likelihood approximation.")
  }
  

  
  # GHK / TMET: only validate M/seed if present
  if (!is.null(options$M)) {
    M <- options$M
    
    ok <- is.numeric(M) && length(M) %in% c(1L, 2L) && all(is.finite(M)) &&
      all(M == as.integer(M)) && all(M > 0)
    
    if (!ok) {
      stop(sprintf( "%s(): options$M must be a positive integer or a vector 
                    of two positive integers.","gctsc"),call. = FALSE)
    }
    
    options$M <- as.integer(M)
  }
  
  # seed checking
  if (!is.null(options$seed)) {
    s <- options$seed
    ok <- is.numeric(s) && length(s) == 1L && is.finite(s) && s == as.integer(s)
    if (!ok) {
      stop(sprintf("%s(): options$seed must be NULL or a single integer.", "gctsc"), call. = FALSE)
    }
  }
  
  # QMC must be logical scalar (for all methods)
  if (!is.logical(QMC) || length(QMC) != 1L || is.na(QMC)) {
    stop(sprintf("%s(): 'QMC' must be TRUE/FALSE.", "gctsc"), call. = FALSE)
  }
  
  # checking c for CE
  if(method=="CE"){
    if (!is.null(c)) {
      ok <- is.numeric(c) && length(c) == 1L && is.finite(c) && c >= 0 && c <= 1
      
      if (!ok) {
        stop(sprintf("%s(): 'c' must be a single numeric value between 0 and 1.", "gctsc"),
             call. = FALSE)
      }
    }
  }
  
  # checking pm for TMET
  
  if(method =="TMET"){
    if (!is.null(pm)) {
      ok <- is.numeric(pm) && length(pm) == 1L && is.finite(pm) && pm == as.integer(pm) && pm > 0
      if (!ok) {
        stop(sprintf("%s(): 'pm' must be a single positive integer.", "gctsc"), call. = FALSE)
      }
    }
  }
  
  
  # checking for length mismatch of y and x
  n <- length(y)
  nm <- vapply(x, nrow, integer(1))
  
  if (any(nm != n)) {
    bad <- names(nm)[nm != n]
    
    stop(sprintf(
      "%s(): nrow mismatch: y has length %d; got %s.", "gctsc", n,
      paste(sprintf("nrow(%s) = %d", bad, nm[bad]), collapse = ", ")
    ), call. = FALSE)
  }
    
  # save family and df to marginal to compute bounds
  marginal$family <- family
  marginal$df <- df
  
  # fitting
  fit <- gctsc.fit(x = x, y = y, marginal = marginal, cormat = cormat,
                   method = method, c = c, QMC = QMC, pm = pm,df=df,
                   start = start, options = options, family = family)
  
  
  
  fit$call <- match.call(expand.dots = FALSE)
  fit$formula <- formula
  fit$terms <- terms
  fit$model <- model
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
gctsc.fit <- function(x, y, marginal, cormat,
                      method, c = 0.5, QMC = TRUE, pm = 30, df=10,
                      start = NULL, options = gctsc.opts(),
                      family = c("gaussian","t")) {
  
  # Missing handling
  missing_y <- anyNA(y)
  
  missing_x <- any(vapply(x, function(X) anyNA(X), logical(1)))
  
  if (missing_y || missing_x) {
    stop(sprintf(
      "%s(): missing values detected in the response or design matrix. Please remove or impute missing values before fitting the model.",
      "gctsc"
    ), call. = FALSE)
  }
  
  
  nbeta <- marginal$npar(x)
  ntau  <- cormat$npar
  
  # Starting values & bounds 
  beta_tmpl <- marginal$start(y, x)
  tau_tmpl  <- cormat$start(y)
  
  
  lb <- c(if(is.null(attr(beta_tmpl, "lower"))){rep(-Inf, length(beta_tmpl))} else attr(beta_tmpl,"lower"),
          if(is.null(attr(tau_tmpl,  "lower"))){  rep(-Inf, length(tau_tmpl))} else attr(tau_tmpl,  "lower") )
  ub <- c(if(is.null(attr(beta_tmpl, "upper"))){rep(Inf, length(beta_tmpl))} else attr(beta_tmpl,"upper"),
          if(is.null(attr(tau_tmpl,  "upper"))){  rep(Inf, length(tau_tmpl))} else attr(tau_tmpl,  "upper") )
  
  if (is.null(start)) {
    init_eta <- c(beta_tmpl, tau_tmpl)
  } else {
    if (!is.numeric(start) || length(start) != (nbeta + ntau) || any(!is.finite(start)))
      stop("gctsc.fit(): 'start' must be numeric, finite, length nbeta + ntau.", call. = FALSE)
    init_eta <- start
  }
  
  f <- structure(list(
    y = y, x = x, c = c, n = length(y), method = method,
    marginal = marginal, cormat = cormat,
    ibeta = 1:nbeta, itau = (nbeta + 1):(nbeta + ntau),
    nbeta = nbeta, ntau = ntau, QMC = QMC, pm = pm,
    call = match.call(), init_eta = init_eta, coef = init_eta,
    lower = lb, upper = ub, options = options,family =family, df=df
  ), class = "gctsc")
  
  gctsc.estimate(f)
}



#' @keywords internal
#' @noRd
gctsc.estimate <- function(cf) {
  start <- cf$init_eta
  low <- cf$lower
  up <- cf$upper
  penalty <- -sqrt(.Machine$double.xmax)
  
  M_vec <- cf$options$M
  
  if (cf$method == "CE") {
    M_vec <- NA_integer_
  }
  
  ans <- NULL
  
  if (cf$method != "CE" && !is.null(cf$options$seed)) {
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    
    if (has_seed) {
      seed.keep <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
  }
  
  for (Mi in M_vec) {
    
    if (cf$method != "CE") {
      cf$options$M_current <- Mi
      log.lik <- llk(cf, Mi, penalty)
    } else {
      log.lik <- llk(cf, cf$options$M[1], penalty)
    }
    
    ans <- suppressWarnings(
      cf$options$opt(start, log.lik, low, up)
    )
    
    start <- ans$estimate
  }
 
  eta <- ans$estimate
  names(eta) <- names(cf$coef)
  
  cf$coef <- eta
  cf$maximum <- ans$maximum
  cf$convergence <- ans$convergence
  cf$M_used <- if (cf$method == "CE") NA_integer_ else M_vec
  system.time({
  xlik <- llk(cf)
  log.lik <- function(th) xlik(pmax(low,pmin(up,th)))
  eps <- .Machine$double.eps^(1/4)
  relStep <- 0.1
  maxtry <- 10
  delta <- ifelse(abs(eta)<1, eps, eps*eta)
  di <- function(i,delta) {
    x1 <- x2 <- eta
    x1[i] <- x1[i] - delta[i]
    x2[i] <- x2[i] + delta[i]
    (log.lik(x2)-log.lik(x1))/(2*delta[i])
  }
  while (1) {
    cf$jac <- sapply(seq_along(eta),di,delta)
    if( all(is.finite(cf$jac)) ) break
    delta <- delta/2
    maxtry <- maxtry - 1
    if (maxtry<0) stop("impossible to compute a finite jacobian")
  }
  
  a <- svd(cf$jac)
  a$d <- pmax(a$d,sqrt(.Machine$double.eps)*a$d[1])
  cf$hessian <- nlme::fdHess(rep(0,length(eta)),
                            function(tx) sum(log.lik(eta+a$v%*%(tx/a$d))),
                            minAbsPar=1,.relStep=relStep)$Hessian
  cf$hessian = (cf$hessian+t(cf$hessian))/2
  cf$hessian <- a$v%*%(outer(a$d,a$d)*(cf$hessian))%*%t(a$v)
  })
  h <- svd(cf$hessian)
  idx <- h$d > sqrt(.Machine$double.eps)*h$d[1]
  vcov <- h$u[,idx,drop=FALSE]%*%( (1/h$d[idx])*t(h$u[,idx,drop=FALSE]))
  
  if (!inherits(vcov, "try-error")) {
    cf$se <- sqrt(diag(vcov))
  }
  
  
  cf
}


#' Set Options for Gaussian and Student t Copula Time Series Model
#'
#' Creates a control list for simulation and likelihood approximation in the
#' Gaussian and Student t copula model, including the random seed, Monte Carlo
#' settings, and optimization controls.
#'
#' @param seed Integer specifying the random seed used for Monte Carlo or
#'   quasi-Monte Carlo simulation during likelihood evaluation. Setting a seed is
#'   recommended for simulation-based methods because it makes the objective
#'   function reproducible across optimization steps.
#' @param M Positive integer or vector of two positive integers. Number of Monte
#'   Carlo or quasi-Monte Carlo samples used in the likelihood approximation. If
#'   a single value is supplied, that value is used throughout the optimization.
#'   If a vector of length two is supplied, staged optimization is used: the
#'   model is first fitted using the first value of \code{M}, and the resulting
#'   estimates are then used as starting values for a second fit using the second
#'   value of \code{M}. This option is used only for simulation-based methods
#'   such as \code{"GHK"} and \code{"TMET"} and is ignored for
#'   \code{method = "CE"}.
#' @param ... Additional control arguments passed to \code{\link[stats]{optim}}.
#'
#' @return
#' A list with components:
#' \item{\code{seed}}{Integer. The random seed used.}
#' \item{\code{M}}{Positive integer or vector of two positive integers specifying
#'   the Monte Carlo or quasi-Monte Carlo sample sizes.}
#' \item{\code{opt}}{A function used internally by \code{gctsc()} to optimize
#'   the approximate log-likelihood.}
#'
#' @export
gctsc.opts <- function(seed = 1, M = c(100, 1000), ...) {
  control <- list(...)
  
  ok_M <- is.numeric(M) &&
    length(M) %in% c(1L, 2L) &&
    all(is.finite(M)) &&
    all(M == as.integer(M)) &&
    all(M > 0)
  
  if (!ok_M) {
    stop(
      "gctsc.opts(): 'M' must be a positive integer or a vector of two positive integers.",
      call. = FALSE
    )
  }
  
  M <- as.integer(M)
  
  opt <- function(start, llk_fn, lower, upper) {
    fn.opt <- function(x) {
      if (any(x <= lower | x >= upper)) return(1e10)
      
      val <- -sum(llk_fn(x))
      
      if (!is.finite(val)) return(1e10)
      val
    }
    
    ans <- stats::optim(
      start,
      fn.opt,
      method = "BFGS",
      hessian = FALSE,
      control = control
    )
    
    if (ans$convergence) {
      warning(paste("optim exits with code", ans$convergence), call. = FALSE)
    }
    
    list(
      estimate = ans$par,
      maximum = ans$value,
      convergence = ans$convergence
    )
  }
  
  list(seed = seed, M = M, opt = opt)
}


#' @keywords internal
#' @noRd
llk <- function(md, M, penalt=NA){
  n <- md$n
  y <- md$y
  x <- md$x
  is.int <- !(md$method %in% c("CE"))
  bounds <- md$marginal$bounds
  ibeta <- md$ibeta
  itau<- md$itau
  od <- md$cormat$od
  p <- md$cormat$od[1]  # AR order
  q <- md$cormat$od[2]  # MA order
  if (missing(M)) {M <- max(md$options$M)}
  seed <- md$options$seed
  c <- md$c
  family <- md$family
  df<- md$df
  method <- md$method
  QMC <- md$QMC
  pm <- md$pm
  cfg <- list(method = method, arg2 = if (is.int) M else c, ret_llk = TRUE,
    pm = pm, od=od, QMC=QMC, df=df )
  cache <- new.env()
  function( eta ) {
    beta <- eta[ibeta]
    if (!identical(cache$beta,beta)) {
      ab <-  bounds(y,x,beta, family, df)
      assign("beta",beta,envir=cache)
      assign("ab",ab,envir=cache)
    } else {
      ab <- get("ab",envir=cache)
    }
    if (is.null(ab) || any(is.nan(ab))) {
      return(penalt)
    }
    tau <- eta[itau]
    if (!identical(cache$tau,tau)) {
      if (p > 0) {
        ar_coefs <- tau[1:p]  # First p elements are AR coefficients
        ar_roots <- polyroot(c(1, -ar_coefs))  # Note the negation for AR polynomial
        if (any(Mod(ar_roots) <= 1.01)) return(NA)  # Penalize invalid AR
      }

      # MA roots check (if q > 0)
      if (q > 0) {
        ma_coefs <- tau[(p + 1):(p + q)]  # Next q elements are MA coefficients
        ma_roots <- polyroot(c(1, ma_coefs))
        if (any(Mod(ma_roots) <= 1.01)) return(NA)  # Penalize invalid MA
      }
      assign("tau",tau,envir=cache)
      
    }
    if (is.int) set.seed(seed)
    lk <-  llk.fn(cfg, ab, tau,family)
    if ( all(is.finite(lk)) ) lk else penalt
  }
}


#' @keywords internal
#' @noRd
llk.fn <- function(cfg, ab, tau, family) {
  
  method  <- cfg$method
  arg2    <- cfg$arg2
  ret_llk <- cfg$ret_llk
  od      <- cfg$od
  QMC     <- cfg$QMC
  pm      <- cfg$pm
  df      <- cfg$df
  
  result <- switch(method,
                   "CE" = loglik_ce( ab = ab, tau = tau, c = arg2, od = od,
                                     ret_llk = ret_llk, df = df, family = family),
                   
                   "GHK" = loglik_ghk( ab = ab, tau = tau, M = arg2, od = od,
                                       QMC = QMC, ret_llk = ret_llk,  df = df,family = family,
                                       engine = "mvmn"),
                   
                   "GHK_mvt" = loglik_ghk(ab = ab, tau = tau,M = arg2, od = od,
                                          QMC = QMC,ret_llk = ret_llk,df = df,family = family,
                                          engine = "mvt"),
                   
                   "TMET" = loglik_tmet(ab = ab, tau = tau,M = arg2, od = od,
                                        QMC = QMC,ret_llk = ret_llk,df = df,family = family),
                   
                   stop("Unknown method")
  )
  
  
  return(result)
}

