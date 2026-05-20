#' @name residuals.gctsc
#' @title Randomized Quantile Residuals for Copula Count Time Series Models
#'
#' @description
#' Computes randomized quantile residuals for a fitted Gaussian or
#' Student--t copula count time series model.
#'
#' For discrete responses, residuals are constructed using the
#' randomized probability integral transform (PIT) of Dunn and Smyth
#' (1996). The conditional probabilities required for the PIT are
#' approximated according to the fitted copula family and likelihood
#' method. For Gaussian copula models fitted by TMET or GHK, the same
#' simulation-based method is used in the residual computation. For
#' Gaussian copula models fitted by CE, the required conditional
#' probabilities are approximated by GHK. For Student--t copula models,
#' the required conditional probabilities are approximated by GHK.
#'
#' @param object A fitted model object of class \code{"gctsc"},
#'   as returned by \code{\link{gctsc}}.
#' @param ... Ignored. Included for S3 method compatibility.
#'
#' @return A list of class \code{"gctsc.residuals"} containing:
#' \itemize{
#'   \item \code{residuals}: Numeric vector of randomized quantile residuals.
#'   \item \code{pit}: Numeric vector of probability integral transform values.
#' }
#'
#' @details
#' For observation \eqn{y_t}, let \eqn{F_t(y_t^-|y)} and \eqn{F_t(y_t|y)}
#' denote the conditional CDF evaluated at \eqn{y_t - 1} and \eqn{y_t} given past observations,
#' respectively. The PIT value is computed as
#' \deqn{
#' e_t = F_t(y_t^-|y) + u_t \{F_t(y_t|y) - F_t(y_t^-|y)\},
#' }
#' where \eqn{u_t \sim \mathrm{Uniform}(0,1)}.
#'
#'For Gaussian copulas, residuals are obtained as
#'\eqn{r_t = \Phi^{-1}(e_t)}.
#'
#'For Student--t copulas with degrees of freedom \code{df},
#'the residuals are defined as \eqn{r_t = t_{\nu}^{-1}(e_t)},
#'where \eqn{t_{\nu}^{-1}} denotes the quantile function of the
#'Student--t distribution.
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
#' # Simulate Poisson AR(1) data under a Gaussian copula
#' set.seed(1)
#' y <- sim_poisson(mu = 5, tau = 0.7,
#'                  arma_order = c(1, 0),
#'                  nsim = 500,
#'                  family = "gaussian")$y
#'
#' fit <- gctsc(
#'   y ~ 1,
#'   data = data.frame(y),
#'   marginal = poisson.marg(),
#'   cormat = arma.cormat(1, 0),
#'   family = "gaussian",
#'   method = "CE",
#'   options = gctsc.opts(seed = 1, M = 1000)
#' )
#'
#' res <- residuals(fit)
#' hist(res$residuals, main = "Randomized Quantile Residuals")
#' hist(res$pit, main = "PIT Histogram")
#' @seealso \code{\link{gctsc}}, \code{\link{sim_gctsc}},
#'   \code{\link{pmvn}}, \code{\link{pmvt}}, \code{\link{predict.gctsc}}
#' @method residuals gctsc
#' @export
residuals.gctsc <- function(object, ...) {
  fn <- "residuals.gctsc"
  
  if (!inherits(object, "gctsc")) {
    stop(sprintf("%s(): object must be of class 'gctsc'.", fn),
         call. = FALSE)
  }
  
  seed <- object$options$seed
  
  if (!is.null(seed)) {
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    
    if (has_seed) {
      seed.keep <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv),
              add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    
    set.seed(seed)
  }
  
  bounds <- object$marginal$bounds
  family <- object$family
  
  ab <- bounds(
    object$y,
    object$x,
    object$coef[object$ibeta],
    family = family,
    df = object$df
  )
  
  res_method <- object$method
  
  if (family == "gaussian" && res_method == "CE") {
    res_method <- "GHK"
  }
  
  if (family == "t") {
    res_method <- "GHK"
  }
  
  cfg <- list(
    method = res_method,
    arg2 = max(object$options$M),
    ret_llk = FALSE,
    pm = object$pm,
    od = object$cormat$od,
    QMC = object$QMC,
    df = object$df
  )
  
  tau <- object$coef[object$itau]
  
  res <- llk.fn(cfg, ab, tau, family)$summary_stats
  
  u <- stats::runif(object$n)
  
  pit_vals <- res[, 1] + u * (res[, 2] - res[, 1])
  
  pit_vals <- pmin(pmax(pit_vals, .Machine$double.eps),
                   1 - .Machine$double.eps)
  
  if (family == "gaussian") {
    res_vals <- stats::qnorm(pit_vals)
  } else {
    res_vals <- stats::qt(pit_vals, df = object$df)
  }
  
  out <- list(
    residuals = res_vals,
    pit = pit_vals
  )
  
  class(out) <- "gctsc.residuals"
  out
}



#' @title Diagnostic Plots for Fitted Copula Count Time Series Models
#'
#' @description
#' Produces diagnostic plots for a fitted Gaussian or Student--t copula
#' count time series model of class \code{"gctsc"}.
#'
#' The diagnostics are based on randomized quantile residuals and
#' probability integral transform (PIT) values.
#'
#' @param x A fitted model object of class \code{"gctsc"}.
#' @param caption Optional character vector of length 5 providing captions
#'   for the plots.
#' @param main Optional main titles for the plots.
#' @param level Confidence level for the Q--Q envelope (default 0.95).
#' @param col.lines Color used for reference lines.
#' @param ... Additional graphical arguments passed to plotting functions.
#'
#' @details
#' The following diagnostic plots are produced:
#' \enumerate{
#'   \item Time series of randomized quantile residuals.
#'   \item Q--Q plot against the reference distribution.
#'   \item Histogram of PIT values.
#'   \item Autocorrelation function (ACF) of residuals.
#'   \item Partial autocorrelation function (PACF) of residuals.
#' }
#'
#' For Gaussian copulas, residuals are compared against the standard
#' normal distribution. For Student--t copulas, residuals are compared
#' against a Student--t distribution with degrees of freedom obtained from fitted model.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @seealso \code{\link{residuals.gctsc}}
#' @export
#'
#' @return Invisibly returns \code{NULL}. The function is called for its
#'   side effect of producing diagnostic plots.
#' @examples
#' # Simulate data from a Poisson AR(1) model
#' set.seed(123)
#' n <- 2000
#' mu <- 5
#' phi <- 0.5
#' arma_order <- c(1, 0)
#' y <- sim_poisson(mu = mu, tau = phi, arma_order = arma_order, nsim = n)$y
#'
#' # Fit the model using the CE method
#' fit <- gctsc(y~1, data = data.frame(y),
#'   marginal = poisson.marg(link = "identity", lambda.lower = 0),
#'   cormat = arma.cormat(p = 1, q = 0), family ="gaussian",
#'   method = "CE",
#'   options = gctsc.opts(seed = 1, M = 1000),
#'   c = 0.5
#' )
#'
#' # Produce diagnostic plots
#' par(mfrow = c(2, 3))
#' plot(fit)

#' @seealso \code{\link{residuals.gctsc}} for computing the residuals used in the plots.
#'
#' @method plot gctsc
#' @export
plot.gctsc <- function(x,caption = rep("", 5),main = rep("", 5),
                       level = 0.95,col.lines = "gray", ...) {
  
  if (!inherits(x, "gctsc"))
    stop("plot.gctsc() is only for objects of class 'gctsc'.")
  
  object  <- x
  family  <- object$family
  res_out <- residuals(object)
  res     <- res_out$residuals
  pit     <- res_out$pit
  
  has_res_na <- anyNA(res)
  has_pit_na <- anyNA(pit)
  
  op <- par(no.readonly = TRUE)
  par(mfrow = c(2, 3))
  on.exit({
    try(par(op), silent = TRUE)
  }, add = TRUE)
  
  ## 1. Time series of residuals
  if (!has_res_na) {
    plot(res, type = "l",xlab = "Time",
         ylab = "Quantile residual", main = main[1], ...)
    abline(h = 0, col = col.lines)
    mtext(caption[1], 3, 0.25)
  } else {
    frame()
  }
  
  ## 2. Q-Q plot
  if (!has_res_na) {
    if (family == "t") {
      
      df  <- object$df                 # t degrees of freedom
      B   <- 1000                      # number of Monte Carlo replicates
      n   <- length(res)
      
      # 1. Theoretical quantiles
      p  <- ppoints(n)
      tq <- qt(p, df = df)
      
      # 2. Simulate t samples to estimate the 2.5% and 97.5% envelopes
      sim_q <- replicate(B, sort(rt(n, df = df)))
      lo <- apply(sim_q, 1, quantile, probs = 0.025)
      hi <- apply(sim_q, 1, quantile, probs = 0.975)
      
      # 3. Sort residuals
      res_sorted <- sort(res)
      
      # 4. Plot QQ with envelope
      plot(tq, res_sorted, pch = 19, col = "gray90",
           main = "QQ Plot vs t-distribution with 95% envelope",
           xlab = "t quantiles", ylab = "Sorted residuals")
      
      # envelope as gray ribbon
      polygon(c(tq, rev(tq)), c(lo, rev(hi)),
              col = adjustcolor("gray90", 0.6), border = NA)
      
      points(tq, res_sorted, pch = 5, col = "gray20", cex = 0.5, lwd = 0.5)
      abline(0, 1, col = "gray70", lwd = 2)
      
      
    } else {
      car::qqPlot(res, envelope = level, grid = FALSE,xlab = "Normal quantiles", 
                  ylab = "Sorted quantile residuals", col.lines = col.lines)
    }

    mtext(caption[2], 3, 0.25)
  } else {
    frame()
  }
  
  ## 3. PIT histogram
  if (!has_pit_na) {
    hist(pit, breaks = 20, col = "skyblue", border = "white",freq = FALSE,
         xlim = c(0, 1),main = main[3], xlab = "PIT values")
    abline(h = 1, col = col.lines, lty = 2)
    mtext(caption[3], 3, 0.25)
  } else {
    frame()
  }
  
  ## 4. ACF
  if (!has_res_na) {
    acf(res, main = main[4],ci.col = col.lines, na.action = na.pass,plot = TRUE)
    mtext(caption[4], 3, 0.25)
  } else {
    frame()
  }
  
  ## 5. PACF
  if (!has_res_na) {
    pacf(res,main = main[5],ci.col = col.lines,na.action = na.pass,plot = TRUE)
    mtext(caption[5], 3, 0.25)
  } else {
    frame()
  }
  
  ## 6th panel left empty
  frame()
  
  invisible(NULL)
}


