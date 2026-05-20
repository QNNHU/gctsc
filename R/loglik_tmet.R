#' @keywords internal
#' @noRd
loglik_tmet <- function(ab, tau, od,
                        family,
                        pm = 30,
                        M = 1000,
                        QMC = TRUE,
                        ret_llk = TRUE,
                        df = NULL) {
  
  if (any(is.na(ab)))
    return(-1e20)
  
  tryCatch(
    tmet_core(ab[,1], ab[,2], tau, od,
              family = family,
              pm = pm,
              M = M,
              QMC = QMC,
              ret_llk = ret_llk,
              df = df),
    error = function(e) -1e20
  )
}


tmet_core <- function(lower, upper, tau, od, family, pm = 30, M = 1000,
                      QMC = TRUE, ret_llk = TRUE, df = NULL) {
  
  if (any(upper < lower)) {
    stop("Invalid bounds: upper < lower.")
  }
  
  n <- length(lower)
  p <- od[1]
  q <- od[2]
  
  if (!family %in% c("gaussian", "t")) {
    stop("family must be either 'gaussian' or 't'.")
  }
  
  if (family == "t" && (is.null(df) || !is.finite(df) || df <= 2)) {
    stop("For the t copula, 'df' must be a finite value greater than 2.")
  }
  
  pm <- if (q == 0) p else pm
  NN <- build_NN(n, pm)
  
  tmet_obj <- cond_mv_tmet(NN, tau, od)
  
  if (family == "gaussian") {
    
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
      method = "L-BFGS-B",
      Condmv_Obj = tmet_obj,
      a = lower,
      b = upper,
      lower = c(lower, rep(-Inf, n)),
      upper = c(upper, rep(Inf, n)),
      control = list(maxit = 500)
    )
    
    delta <- solv_delta$par[(n + 1):(2 * n)]
    
    sampler_input <- list(
      a = lower,
      b = upper,
      delta = delta,
      condSd = sqrt(tmet_obj$cond_var),
      phi = tmet_obj$phi,
      q = tmet_obj$q,
      m = tmet_obj$m,
      Theta = tmet_obj$Theta,
      M = M,
      QMC = QMC
    )
    
    if (ret_llk) {
      return(ptmvn_tmet(sampler_input))
    } else {
      return(rtmvn_tmet(sampler_input))
    }
  }
  
  if (family == "t") {
    
    trunc_expect <- truncnorm::etruncnorm(lower, upper)
    z0 <- c(trunc_expect, rep(0, n))
    
    z0[n] <- 1
    z0[2 * n] <- 0
    
    sol <- lm_sparse_solver(
      z0,
      Condmv_Obj = tmet_obj,
      a = lower,
      b = upper,
      nu = df,
      maxit = 500,
      tol = 1e-1,
      verbose = FALSE
    )
    
    # if (isFALSE(sol$converged)) {
    #   stop("TMET tilting solver failed: ", sol$message, call. = FALSE)
    # }
    
    z <- sol$x
    delta <- z[(n + 1):(2 * n)]
    kap <- delta[n]
    
    const <- log(2 * pi) / 2 -
      lgamma(df / 2) -
      (df / 2 - 1) * log(2) +
      TruncatedNormal::lnNpr(-kap, Inf) +
      0.5 * kap^2
    
    sampler_input <- list(
      a = lower,
      b = upper,
      condSd = sqrt(tmet_obj$cond_var),
      phi = tmet_obj$phi,
      q = tmet_obj$q,
      m = tmet_obj$m,
      delta = delta,
      df = df,
      Theta = tmet_obj$Theta,
      M = M,
      QMC = QMC
    )
    
    if (ret_llk) {
      return(ptmvt_tmet(sampler_input) + const / n)
    } else {
      return(rtmvt(sampler_input))
    }
  }
}


build_NN <- function(n, pm) {
  NN <- matrix(NA_integer_, nrow = n, ncol = pm + 1)
  NN[1, 1] <- 1
  
  if (n > 1) {
    row_idx <- 2:n
    col_idx <- 0:pm
    
    idx_mat <- outer(row_idx - 1, col_idx, FUN = function(i, j) i - j)
    valid_mask <- outer(row_idx - 1, col_idx, FUN = function(i, j) j <= i)
    idx_mat[!valid_mask] <- NA_integer_
    
    NN[2:n, ] <- idx_mat + 1
  }
  
  NN
}