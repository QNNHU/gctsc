#' @name predict.gctsc
#' @title One-Step-Ahead Predictive Distribution for Copula Count Time Series Models
#'
#' @description
#' Computes the one-step-ahead predictive distribution for a fitted Gaussian or
#' Student-\eqn{t} copula count time series model.
#'
#' The predictive probability mass function is evaluated on the count grid
#' \code{0:y_max}. For bounded marginal distributions, \code{y_max} is set to
#' the size parameter at the prediction time. For unbounded marginal
#' distributions, \code{y_max} may be supplied by the user; otherwise, it is
#' selected automatically as \code{ceiling(max(y) + k * sd(y))}, where \code{y}
#' is the fitted response series.
#'
#' The function returns summary statistics of the predictive distribution. If
#' the observed response at the prediction time is included in \code{newdata},
#' the Continuous Ranked Probability Score (CRPS) and Logarithmic Score (LOGS)
#' are also computed.
#'
#' @param object A fitted model object of class \code{"gctsc"}, as returned by
#'   \code{\link{gctsc}}.
#'
#' @param newdata Optional one-row \code{data.frame} containing the covariate
#'   values at the prediction time point. The variables in \code{newdata} should
#'   match those used in the fitted model formula. If the fitted model is
#'   intercept-only, \code{newdata} may be omitted.
#'
#'   If the observed response at the prediction time is available, it may also
#'   be included in \code{newdata} using the same response variable name as in
#'   the fitted model. In this case, CRPS and LOGS are computed and returned.
#'
#' @param y_max Optional nonnegative integer specifying the largest count value
#'   included in the predictive grid \code{0:y_max}. For bounded marginal
#'   distributions, this value is determined by the size parameter at the
#'   prediction time. For unbounded marginal distributions, if \code{y_max} is
#'   \code{NULL}, it is selected automatically as
#'   \code{ceiling(max(y) + k * sd(y))}, where \code{y} is the fitted response
#'   series.
#'
#' @param k Nonnegative numeric multiplier used to choose \code{y_max} when
#'   \code{y_max = NULL}. The default is \code{k = 3}.
#'
#' @param ... Ignored. Included for S3 method compatibility.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{mean}: Predictive mean.
#'   \item \code{median}: Predictive median.
#'   \item \code{mode}: Predictive mode.
#'   \item \code{variance}: Predictive variance.
#'   \item \code{p_y}: Predictive probability mass function evaluated on
#'         \code{y_grid}.
#'   \item \code{y_grid}: Count grid \code{0:y_max} over which the predictive
#'         distribution is evaluated.
#'   \item \code{lower}, \code{upper}: Lower and upper endpoints of the 95\%
#'         predictive interval.
#'   \item \code{CRPS}: Continuous Ranked Probability Score, returned only if
#'         the observed response is supplied in \code{newdata}.
#'   \item \code{LOGS}: Logarithmic Score, returned only if the observed response
#'         is supplied in \code{newdata}.
#'   \item \code{y_true}: Observed response value, returned only if supplied in
#'         \code{newdata}.
#' }
#'
#' @details
#' The function constructs prediction design matrices from the stored model
#' terms using \code{\link[stats]{model.matrix}}. Therefore, the same formula
#' convention used in model fitting is used for prediction, including automatic
#' intercept handling and factor-variable expansion.
#'
#' For zero-inflated marginals, \code{newdata} is used to construct both the
#' mean-component design matrix and the zero-inflation design matrix. The column
#' names of the new design matrices must match those from the fitted model.
#'
#' For Gaussian copulas, the predictive distribution is computed using the 
#' approximation method stored in the fitted object.
#' For Student-\eqn{t} copulas, the predictive distribution is currently computed using the
#' GHK approximation.
#'
#' If the observed response is included in \code{newdata}, it must be a single
#' nonnegative integer count and should not exceed \code{y_max}.
#'
#' @references
#' Nguyen, Q. N. and De Oliveira, V. (2026), Likelihood Inference in Gaussian
#' Copula Models for Count Time Series via Minimax Exponential Tilting,
#' \emph{Computational Statistics & Data Analysis}, \strong{218}: 108344.
#'
#' Nguyen, Q. N. and De Oliveira, V. (2026), Scalable Likelihood Inference
#' for Student--\eqn{t} Copula Count Time Series, \emph{Stats},
#' \strong{9}(2): 43.
#'
#' @examples
#' # Simulate Poisson AR(1) data
#' set.seed(1)
#' y_sim <- sim_poisson(
#'   mu = 10,
#'   tau = 0.2,
#'   arma_order = c(1, 0),
#'   nsim = 200,
#'   family = "gaussian"
#' )$y
#'
#' dat <- data.frame(y = y_sim)
#'
#' # Fit Gaussian copula model
#' fit <- gctsc(
#'   formula = y ~ 1,
#'   data = dat,
#'   marginal = poisson.marg(link = "log"),
#'   cormat = arma.cormat(p = 1, q = 0),
#'   method = "GHK",
#'   family = "gaussian",
#'   options = gctsc.opts(M = 1000, seed = 42)
#' )
#'
#' # One-step-ahead prediction for an intercept-only model
#' pred <- predict(fit, y_max = 30)
#'
#' # If the future observed value is available, include it in newdata
#' pred_score <- predict(fit, newdata = data.frame(y = 8), y_max = 30)
#'
#' @seealso \code{\link{gctsc}}, \code{\link{arma.cormat}},
#'   \code{\link{gctsc.opts}}
#'
#' @method predict gctsc
#' @export
#' 
predict.gctsc <- function(object, newdata = NULL, y_max = NULL, k = 3, ...) {
  fn <- "predict.gctsc"
  
  if (!inherits(object, "gctsc")) {
    stop(sprintf("%s(): object must be of class 'gctsc'.", fn),
         call. = FALSE)
  }
  
  if (is.null(newdata)) {
    newdata <- data.frame(.dummy = 1)
  }
  
  if (!is.data.frame(newdata)) {
    stop(sprintf("%s(): 'newdata' must be a data frame.", fn),
         call. = FALSE)
  }
  
  if (nrow(newdata) != 1L) {
    stop(sprintf("%s(): 'newdata' must have exactly one row for one-step prediction.", fn),
         call. = FALSE)
  }
  
  bounds <- object$marginal$bounds
  y <- object$y
  x <- object$x
  seed <- object$options$seed
  family <- object$family
  method <- object$method
  
  # Determine y_max for predictive distribution
  if(!object$marginal$bounded){
    if (is.null(y_max)) {
      y_sd <- stats::sd(y)
      y_max <- ceiling(max(y) + k * y_sd)
    }
    
    if (!is.numeric(y_max) || length(y_max) != 1L || !is.finite(y_max) || y_max < 0 ||
        y_max != as.integer(y_max)) {
      stop(sprintf("%s(): 'y_max' must be a single nonnegative integer.", fn),
           call. = FALSE)
    } 
    
    
    y_max <- as.integer(y_max)
    
    if (!is.finite(y_max) || y_max < max(y)) {
      y_max <- max(y)
    }
  } else {
    y_max <- object$marginal$size
    
  }
  
  y_grid <- 0:y_max
  n_grid <- length(y_grid)
  
  # Check whether observed response is supplied in newdata
  y_name <- as.character(object$formula$mu[[2L]])
  has_y <- y_name %in% names(newdata)
  y_true <- NULL
  
  if (isTRUE(has_y)) {
    y_true <- newdata[[y_name]]
    
    if (length(y_true) != 1L || !is.numeric(y_true) ||
        !is.finite(y_true) || y_true < 0 ) {
      stop(sprintf("%s(): observed response in 'newdata' must be a single nonnegative integer count.", fn),
           call. = FALSE)
    }
  }
  
  # Mean component design matrix
  terms_mu <- delete.response(object$terms$mu)
  X_mu_new <- model.matrix(terms_mu, data = newdata)
  
  if (!identical(colnames(X_mu_new), colnames(object$x$mu))) {
    missing_cols <- setdiff(colnames(object$x$mu), colnames(X_mu_new))
    extra_cols <- setdiff(colnames(X_mu_new), colnames(object$x$mu))
    
    stop(sprintf(
      "%s(): prediction design matrix does not match fitted mean design matrix. Missing: %s. Extra: %s.",
      fn,
      paste(missing_cols, collapse = ", "),
      paste(extra_cols, collapse = ", ")
    ), call. = FALSE)
  }
  
  X_mu_rep <- X_mu_new[rep(1L, n_grid), , drop = FALSE]
  
  if (isTRUE(object$marginal$zero_inflated)) {
    terms_pi0 <- object$terms$pi0
    X_pi0_new <- model.matrix(terms_pi0, data = newdata)
    
    if (!identical(colnames(X_pi0_new), colnames(object$x$pi0))) {
      missing_cols <- setdiff(colnames(object$x$pi0), colnames(X_pi0_new))
      extra_cols <- setdiff(colnames(X_pi0_new), colnames(object$x$pi0))
      
      stop(sprintf(
        "%s(): prediction design matrix does not match fitted zero-inflation design matrix. Missing: %s. Extra: %s.",
        fn,
        paste(missing_cols, collapse = ", "),
        paste(extra_cols, collapse = ", ")
      ), call. = FALSE)
    }
    
    X_pi0_rep <- X_pi0_new[rep(1L, n_grid), , drop = FALSE]
    
    x_rep <- list(mu = X_mu_rep, pi0 = X_pi0_rep)
  } else {
    x_rep <- list(
      mu = X_mu_rep
    )
  }
  
  # Bounds for observed values and prediction grid
  ab <- bounds(y,x,object$coef[object$ibeta],family = family,df = object$df)
  ab_p <- bounds(y_grid,x_rep,object$coef[object$ibeta],family = family,
    df = object$df)
  
  pred_input <- list(a = ab[, 1],b = ab[, 2],ap = ab_p[, 1],bp = ab_p[, 2],
    y_max = y_max,M = max(object$options$M), QMC = object$QMC )
  
  od <- object$cormat$od
  tau <- object$coef[object$itau]
  pm <- object$pm
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (family == "gaussian") {
    if (method == "TMET") {
      p_y <- pred_tmet_mvn(pred_input, tau = tau, od = od, pm =pm )$p_y
    } else{
      p_y <- pred_ghk_mvn(pred_input, tau = tau, od = od)$p_y
      } 
    } else {
      pred_input$df <- object$df
      p_y <- pred_mvt(args = pred_input, tau = tau, od = od)$p_y
  }
  
  # Normalize just in case of small numerical error
  p_y <- as.vector(p_y)
  p_y <- pmax(p_y, 0)
  p_y <- p_y / sum(p_y)
  
  # Summary statistics
  cdf <- cumsum(p_y)
  
  mean_pred <- sum(y_grid * p_y)
  var_pred <- sum((y_grid^2) * p_y) - mean_pred^2
  median_pred <- y_grid[which(cdf >= 0.5)[1]]
  mode_pred <- y_grid[which.max(p_y)]
  
  alpha <- 0.05
  
  lower <- y_grid[which(cdf >= alpha / 2)[1]]
  upper <- y_grid[which(cdf >= 1 - alpha / 2)[1]]
  
  out <- list(mean = mean_pred, median = median_pred, mode = mode_pred,
    variance = var_pred, p_y = p_y, y_grid = y_grid, 
    lower = lower,upper = upper )
  
  if (isTRUE(has_y)) {
    indicator <- as.numeric(y_grid >= y_true)
    
    out$CRPS <- sum((cdf - indicator)^2)
    out$LOGS <- -log(pmax(p_y[y_true + 1L], .Machine$double.eps))
    out$y_true <- y_true
  }
  
  out
}




pred_ghk_mvn <- function(args, tau, od) {
  fn <- "pred_ghk_mvn"
  
  if (anyNA(args$a) || any(is.nan(args$a))) {
    return(list(p_y = NA))
  }
  
  n <- length(args$a)
  
  model <- arma_model(tau = tau,od = od,n = n,fn = fn)
  
  predmvn_ghk(args, model)
}





pred_tmet_mvn <- function(args, tau, od, pm){
  fn <- "pred_tmet_mvn"
  
  if (anyNA(args$a) || any(is.nan(args$a))) {
    return(list(p_y = NA))
  }
  
  n <- length(args$a)
  
  arma <- arma_model(tau = tau,od = od,n = n,fn = fn
  )
  
  p0 <- od[1]
  q0 <- od[2]
  pm <- if (q0 == 0) arma$p else pm
  
  NN <- build_NN(n, pm)
  
  tmet_obj <- cond_mv_tmet(NN, tau, od)
  
  lower <- args$a
  upper <- args$b
  
  z0 <- truncnorm::etruncnorm(lower, upper)
  z0_delta0 <- c(z0, rep(0, n))
  
  solv_delta <- stats::optim(
    z0_delta0,
    fn = function(x, ...) {
      ret <- grad_jacprod(x, ..., retProd = FALSE)
      0.5 * sum(ret$grad^2)
    },
    gr = function(x, ...) {
      ret <- grad_jacprod(x, ..., retProd = TRUE)
      ret$jac_grad
    },
    method = "L-BFGS-B", Condmv_Obj = tmet_obj,
    a = lower,  b = upper,
    lower = c(lower, rep(-Inf, n)),
    upper = c(upper, rep(Inf, n)),
    control = list(maxit = 500)
  )
  
  if (any(solv_delta$par[seq_len(n)] < lower) ||
      any(solv_delta$par[seq_len(n)] > upper)) {
    warning(
      sprintf("%s(): optimal x is outside the integration region during minimax tilting.", fn),
      call. = FALSE
    )
  }
  
  delta <- solv_delta$par[(n + 1L):(2L * n)]
  
  model <- c(arma,
    list(
      delta = delta, Theta = rbind(tmet_obj$Theta),
      condSd = sqrt(tmet_obj$cond_var),
      v = tmet_obj$cond_var
    )
  )
  
  predmvn_tmet(args, model)  
}




pred_mvt <- function(args, tau, od) {
  fn <- "pred_mvt"
  
  if (anyNA(args$a) || any(is.nan(args$a))) {
    return(list(p_y = NA))
  }
  
  n <- length(args$a)
  
  model <- arma_model(tau = tau,od = od,n = n,fn = fn)
  
  predmvt(args, model)
}


arma_model <- function(tau, od, n, fn = "prediction") {
  if (length(tau) != sum(od)) {
    stop(sprintf("%s(): length of 'tau' must match ARMA order.", fn),
         call. = FALSE)
  }
  
  if (all(od == 0)) {
    stop(sprintf("%s(): ARMA(0,0) is not supported.", fn),
         call. = FALSE)
  }
  
  p0 <- od[1]
  q0 <- od[2]
  
  iar <- if (p0 > 0) seq_len(p0) else integer(0)
  ima <- if (q0 > 0) (p0 + 1L):(p0 + q0) else integer(0)
  
  phi <- if (p0 > 0) tau[iar] else 0
  theta <- if (q0 > 0) tau[ima] else 0
  
  p <- if (p0 == 0) 1L else p0
  q <- if (q0 == 0) 1L else q0
  m <- max(p, q)
  
  Tau <- list(phi = phi, theta = theta)
  
  list(phi = phi, theta = theta, theta_r = c(1, theta, numeric(n)),
       n = n, p = p, q = q,  m = m, sigma2 = 1 / sum(ma.inf(Tau)^2),
       gamma = aacvf(Tau, n - 1)
  )
}


