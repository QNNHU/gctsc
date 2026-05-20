#' @keywords internal
#' @noRd
cond_mv_base <- function(n, tau, od, fn = "cond_mv_base") {
  p0 <- od[1]
  q0 <- od[2]
  
  if (length(tau) != sum(od)) {
    stop(sprintf("%s(): length of 'tau' must match sum(od).", fn),
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
  gamma <- aacvf(Tau, n - 1)
  
  model <- list(
    phi = phi,
    theta_r = c(1, theta, numeric(n)),
    p = p,
    q = q,
    m = m,
    sigma2 = sigma2,
    gamma = gamma,
    n = n
  )
  
  out_cpp <- compute_cond_var(gamma, model)
  
  list(
    cond_var = out_cpp$v * sigma2,
    Theta = out_cpp$Theta,
    phi = phi,
    p = p0,
    q = q0,
    m = m,
    Tau = Tau,
    gamma = gamma
  )
}

#' @keywords internal
#' @noRd
cond_mv_ghk <- function(n, tau, od) {
  base <- cond_mv_base(n, tau, od, fn = "cond_mv_ghk")
  
  list(
    cond_var = base$cond_var,
    Theta = base$Theta,
    m = base$m,
    phi = base$phi,
    p = base$p,
    q = base$q
  )
}


#' @keywords internal
#' @noRd
cond_mv_tmet <- function(NN, tau, od) {
  n <- nrow(NN)
  
  if (!all(NN[, 1] == seq_len(n))) {
    stop("cond_mv_tmet(): unexpected NN; first column must be 1:n.",
         call. = FALSE)
  }
  
  pm <- ncol(NN) - 1L
  
  base <- cond_mv_base(n, tau, od, fn = "cond_mv_tmet")
  
  arma_blup_coef <- -ar.inf(base$Tau, pm)[-1]
  
  cond_mean_coeff <- matrix(0, nrow = n, ncol = pm)
  
  if (pm > 2L) {
    DL <- durbin_levinson(base$gamma, pm, base$cond_var[-1])
    
    for (i in 2:(pm - 1L)) {
      cond_mean_coeff[i, seq_len(i - 1L)] <- DL[i - 1L, seq_len(i - 1L)]
    }
  }
  
  start_i <- max(2L, pm)
  
  if (start_i <= n) {
    for (i in start_i:n) {
      cond_mean_coeff[i, seq_len(pm)] <- arma_blup_coef[seq_len(pm)]
    }
  }
  
  B <- sparse_B(NN, cond_mean_coeff, n, pm)
  
  list(
    cond_var = base$cond_var,
    cond_mean_coeff = cond_mean_coeff,
    B = B,
    NN = NN,
    Theta = base$Theta,
    phi = base$phi,
    p = base$p,
    q = base$q,
    m = base$m
  )
}

#' @keywords internal
#' @noRd
sparse_B <- function(NNarray, cond_mean_coeff, n, pm) {
  nnz_per_row <- pmin(pm + 1, 1:n)
  total_nnz <- sum(nnz_per_row)
  
  B_row_inds <- integer(total_nnz)
  B_col_inds <- integer(total_nnz)
  B_vals     <- numeric(total_nnz)
  
  pos <- 1
  for (i in 1:n) {
    k <- nnz_per_row[i]
    idx <- pos:(pos + k - 1)
    
    B_row_inds[idx] <- i
    B_col_inds[idx] <- NNarray[i, 1:k]
    
    if (k > 1) {
      B_vals[idx[-1]] <- cond_mean_coeff[i, seq_len(k - 1)]
    }
    
    pos <- pos + k
  }
  
  Matrix::sparseMatrix(i = B_row_inds, j = B_col_inds, x = B_vals, dims = c(n, n))
}



