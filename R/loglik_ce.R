#' @keywords internal
#' @noRd
loglik_ce <- function(ab, tau, od, family, c = 0.5, ret_llk = TRUE, df = NULL) {
  
  if (length(tau) != sum(od))
    stop("Length of 'tau' must equal p+q.")
  
  if (all(od == 0))
    stop("ARMA(0,0) not supported.")
  
  tryCatch(
    ce_core(ab[,1], ab[,2], tau, od, family = family, c = c,
            ret_llk = ret_llk, df = df),
    error = function(e) -1e20
  )
}



ce_core <- function(lower, upper, tau, od, family, c = 0.5, ret_llk = TRUE,
                    df = NULL) {
  
  EPS <- sqrt(.Machine$double.eps)
  EPS1 <- 1 - EPS
  
  if (anyNA(lower) || anyNA(upper) || any(upper < lower - EPS))
    return(-1e20)
  
  n <- length(lower)
  if(ret_llk){
  # Copula CDF and density difference
  if (family == "gaussian") {
    cdf <- pnorm(upper)
    pdf <- cdf - pnorm(lower)
    r <- qnorm(pmin(EPS1, pmax(EPS, cdf - c * pdf)))
  } else {
    cdf <- pt(upper, df)
    pdf <- cdf - pt(lower, df)
    r <- t_qt_safe(cdf - c * pdf, df = df)
  }
  }
  
  # Extract ARMA
  p0 <- od[1]
  q0 <- od[2]
  
  if (length(tau) != sum(od)) {
    stop(sprintf("%s(): length of 'tau' must match sum(od).", "CE"),
         call. = FALSE)
  }
  
  
  iar <- if (p0 > 0) seq_len(p0) else integer(0)
  ima <- if (q0 > 0) (p0 + 1L):(p0 + q0) else integer(0)
  
  phi <- if (p0 > 0) tau[iar] else 0
  theta <- if (q0 > 0) tau[ima] else 0
  
  # Innovations algorithm assumes p >= 1 and q >= 1.
  p <- if (p0 == 0) 1L else p0
  q <- if (q0 == 0) 1L else q0
  
  m <- max(p, q)
  Tau <- list(phi = phi, theta = theta)
  sigma2 <- 1 / sum(ma.inf(Tau)^2)
  gamma  <- aacvf(Tau, n - 1)
  theta_r <- c(1, theta, numeric(n))
  
  if (ret_llk) {
    
    model <- list(
      phi = phi, theta = theta, r = r, theta_r = theta_r, n = n,
      p = p, q = q, m = m, sigma2 = sigma2,  a = pdf
    )
    
    if (family == "gaussian") {
      results <- ptmvn_ce(gamma, model)
    } else {
      model$df <- df
      results <- ptmvt_ce(gamma, model)
    }
    
    return(results$llk)
    
  } else {
    cond_mv <- cond_mv_ghk(n, tau, od)
    sampler_input <- list(a = lower,b = upper,condSd = sqrt(cond_mv$cond_var),
      M = 1000,phi = phi,q = q,m = m,Theta = cond_mv$Theta,QMC = TRUE
    )
    return(rtmvn_ghk(sampler_input)$summary_stats)
  }
}


t_qt_safe <- function(p, df, eps = 1e-8) {
  # Clip probabilities
  p <- pmin(1 - eps, pmax(eps, p))
  
  if (is.infinite(df) || df >= 100) {
    # Gaussian fallback for large df
    return(qnorm(p))
  }
  
  # For extreme tails, normal approx is safer
  if (any(p < 1e-6 | p > 1 - 1e-6)) {
    return(qnorm(p))  # fast, stable fallback
  }
  
  # Otherwise use qt
  return(qt(p, df))
}
