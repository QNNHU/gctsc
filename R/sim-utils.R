## -----------------------------------------------------------
## simulation Validate helpers 
## -----------------------------------------------------------

.scalar_or_nsim <- function(x, nsim, name, fn) {
  len <- length(x)
  if (len %in% c(1L, nsim)) {
    if (!is.numeric(x) || any(!is.finite(x))) {
      stop(sprintf("%s(): '%s' must be numeric and finite.", fn, name), call. = FALSE)
    }
    return(invisible(TRUE))
  }
  stop(sprintf("%s(): '%s' must be scalar or length nsim = %d. Got length %d.",
               fn, name, nsim, len), call. = FALSE)
}

.check_size <- function(size, fn, scalar_only = TRUE, nsim = NULL) {
  if (scalar_only) {
    if (length(size) != 1L)
      stop(sprintf("%s(): 'size' must be a single positive integer.", fn), call. = FALSE)
  } else {
    .scalar_or_nsim(size, nsim, "size", fn)
  }
  if (!is.numeric(size) || !is.finite(size) || size <= 0 || size != as.integer(size)) {
    stop(sprintf("%s(): 'size' must be a positive integer.", fn), call. = FALSE)
  }
}

.check_common <- function(nsim, tau, arma_order, seed, family, df, fn) {
  # nsim
  if (!is.numeric(nsim) || length(nsim) != 1L || !is.finite(nsim) ||
      nsim < 1 || nsim != as.integer(nsim)) {
    stop(sprintf("%s(): 'nsim' must be a single positive integer. Got %s.",
                 fn, paste(nsim, collapse = ",")), call. = FALSE)
  }
  # arma_order
  if (!is.numeric(arma_order) || length(arma_order) != 2L ||
      any(!is.finite(arma_order)) || any(arma_order < 0) ||
      any(arma_order != as.integer(arma_order))) {
    stop(sprintf("%s(): 'arma_order' must be integer c(p, q) with p,q >= 0. Got %s.",
                 fn, paste(arma_order, collapse = ",")), call. = FALSE)
  }
  if (all(arma_order == 0L)) {
    stop(sprintf("%s(): ARMA(0,0) is not supported.", fn), call. = FALSE)
  }
  # tau
  if (!is.numeric(tau) || any(!is.finite(tau))) {
    stop(sprintf("%s(): 'tau' must be numeric and finite.", fn), call. = FALSE)
  }
  if (length(tau) != sum(arma_order)) {
    stop(sprintf("%s(): length(tau) must equal p+q = %d. Got %d.",
                 fn, sum(arma_order), length(tau)), call. = FALSE)
  }
  # seed
  if (!is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed != as.integer(seed))) {
    stop(sprintf("%s(): 'seed' must be NULL or a single integer.", fn), call. = FALSE)
  }
  
  # Stationarity/invertibility checks
  p <- arma_order[1]; q <- arma_order[2]
  phi   <- if (p > 0) tau[seq_len(p)] else numeric(0)
  theta <- if (q > 0) tau[p + seq_len(q)] else numeric(0)
  
  if (p > 0) {
    ar_poly <- c(1, -phi)                     # 1 - sum phi_i z^i
    r <- polyroot(ar_poly)
    if (any(Mod(r) <= 1)) {
      stop(sprintf("%s(): AR polynomial is not stationary (a root has |z| <= 1). Coefficients: %s",
                   fn, paste(phi, collapse = ",")), call. = FALSE)
    }
  }
  if (q > 0) {
    ma_poly <- c(1, theta)                     # 1 + sum theta_j z^j
    r <- polyroot(ma_poly)
    if (any(Mod(r) <= 1)) {
      stop(sprintf("%s(): MA polynomial is not invertible (a root has |z| <=1). Coefficients: %s",
                   fn, paste(theta, collapse = ",")), call. = FALSE)
    }
  }
  # family
  if (!is.character(family) || length(family) != 1L ||
      !(family %in% c("gaussian", "t"))) {
    stop(sprintf("%s(): 'family' must be 'gaussian' or 't'. Got %s.",
                 fn, paste(family, collapse = ",")), call. = FALSE)
  }
  
  # df (only for t copula)
  if (family == "t") {
    if (is.null(df) || !is.numeric(df) || length(df) != 1L ||
        !is.finite(df) || df <= 2) {
      stop(sprintf("%s(): 'df' must be a single numeric value > 2 when family = 't'. Got %s.",
                   fn, paste(df, collapse = ",")), call. = FALSE)
    }
  }
  invisible(TRUE)
}

.check_mu <- function(mu, nsim, fn) {
  .scalar_or_nsim(mu, nsim, "mu", fn)
  if (any(!(mu > 0))) {
    stop(sprintf("%s(): 'mu' must be strictly > 0.", fn), call. = FALSE)
  }
}

.check_prob <- function(prob, nsim, fn) {
  .scalar_or_nsim(prob, nsim, "prob", fn)
  if (any(!(prob > 0 & prob < 1))) {
    stop(sprintf("%s(): 'prob' must lie in (0, 1).", fn), call. = FALSE)
  }
}

.recyclen <- function(x, n, name = deparse(substitute(x))) {
  if (length(x) == n) return(x)
  if (length(x) == 1L) return(rep(x, n))
  stop(sprintf("Parameter '%s' must be length 1 or %d (got %d).", name, n, length(x)))
}


.check_pi0 <- function(pi0, nsim, fn) {
  .scalar_or_nsim(pi0, nsim, "pi0", fn)
  if (any(!(pi0 >= 0 & pi0 < 1))) {
    stop(sprintf("%s(): 'pi0' must lie in [0, 1).", fn), call. = FALSE)
  }
}

.check_dispersion <- function(dispersion, nsim, fn) {
  .scalar_or_nsim(dispersion, nsim, "dispersion", fn)
  if (any(!(dispersion > 0))) {
    stop(sprintf("%s(): 'dispersion' must be strictly > 0.", fn), call. = FALSE)
  }
}


.check_rho <- function(rho, nsim, fn) {
  .scalar_or_nsim(rho, nsim, "rho", fn)
  if (any(!(rho > 0 & rho < 1))) {
    stop(sprintf("%s(): 'rho' must lie in (0, 1).", fn), call. = FALSE)
  }
}



