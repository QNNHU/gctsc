#' @keywords internal
#' @noRd
loglik_ghk <- function(ab, tau, od, family, M = 1000,QMC = TRUE,
                       ret_llk = TRUE, df = NULL, engine = "mvmn") {
  
  if (length(tau) != sum(od))
    stop("Length of 'tau' must equal p+q.")
  
  if (all(od == 0))
    stop("ARMA(0,0) not supported.")
  
  if (any(is.na(ab)))
    return(-1e20)
  
  a <- ab[, 1]
  b <- ab[, 2]
  
  out <- tryCatch(
    ghk_core(a, b, tau, od, family = family, M = M, QMC = QMC,
             ret_llk = ret_llk, df = df, engine = engine),
    error = function(e) {
      message("GHK failed: ", e$message)
      -1e20
    }
  )

  return(out)
}


ghk_core <- function(lower, upper, tau, od, family, M = 1000, QMC = TRUE,
                     ret_llk = TRUE, df = NULL, engine = "mvmn") {
  
  if (any(upper < lower)) {
    stop("Invalid bounds: some upper bounds are smaller than lower bounds.")
  }
  
  if (!family %in% c("gaussian", "t")) {
    stop("family must be either 'gaussian' or 't'.")
  }
  
  if (family == "t" && (is.null(df) || !is.finite(df) || df <= 2)) {
    stop("For the t copula, 'df' must be a finite value greater than 2.")
  }
  
  n <- length(lower)
  
  ghk_obj <- cond_mv_ghk(n, tau, od)
  
  sampler_input <- list(a = lower, b = upper, condSd = sqrt(ghk_obj$cond_var),
    M = M, phi = ghk_obj$phi, q = ghk_obj$q, m = ghk_obj$m, Theta = ghk_obj$Theta,
    QMC = QMC
  )
  
  if (family == "gaussian") {
    if (ret_llk) {
      return(ptmvn_ghk(sampler_input))
    } else {
      return(rtmvn_ghk(sampler_input))
    }
  }
  
  if (family == "t") {
    sampler_input$df <- df
    
    if (ret_llk) {
      if (engine == "mvt") {
        return(ptmvt_ghk(sampler_input))
      } else {
        return(ptmvmn_ghk(sampler_input))
      }
    } else {
        return( rtmvt(sampler_input))
    }
  }
}