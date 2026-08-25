#' @useDynLib gctsc, .registration = TRUE
#' @importFrom Rcpp sourceCpp
# core_utils.R
# Internal utilities for TMET and GHK simulation, gradient, and likelihood computation

#' @keywords internal
#' @noRd
sample_mvn <- function(obj, lower, upper, delta = NULL, M = 1000, QMC = TRUE, ret_llk = TRUE,
                          method = c("TMET", "GHK")) {
  method <- match.arg(method)

  sampler_input <- list(
    a = lower,
    b = upper,
    condSd = sqrt(obj$cond_var),
    M = M,
    phi = obj$phi,
    q = obj$q,
    m = obj$m,
    Theta = obj$Theta,
    QMC = QMC
  )

  # Include delta only for TMET


  if (ret_llk) {
    if (method == "TMET") {
      sampler_input$delta <- delta
    }
    switch(method,
           "TMET" = ptmvn_tmet(sampler_input),
           "GHK"  = ptmvn_ghk(sampler_input),
           stop("Unknown method")
    )
  } else {
    if (method == "TMET") {
      sampler_input$delta <- delta
    }
    switch(method,
           "TMET" = rtmvn_tmet(sampler_input),
           "GHK"  = rtmvn_ghk(sampler_input),
           stop("Unknown method")
    )
  }
}



#' @keywords internal
#' @noRd
sample_mvt <- function(obj, lower, upper, delta = NULL, M = 1000, QMC = TRUE, ret_llk = TRUE,
                          method = c("TMET", "GHK","GHK_alg1"), df) {
  method <- match.arg(method)
  
  sampler_input <- list(
    a = lower,
    b = upper,
    condSd = sqrt(obj$cond_var),
    M = M,
    phi = obj$phi,
    q = obj$q,
    m = obj$m,
    Theta = obj$Theta,
    QMC = QMC,
    df= df
    
  )
  
  # Include delta only for TMET
  
  
  if (ret_llk) {
    if (method == "TMET") {
      sampler_input$delta <- delta
    }
    switch(method,
           "TMET" = ptmvt_tmet(sampler_input),
           "GHK"  = ptmvt_ghk(sampler_input),
           "GHK_alg1" = ptmvt_ghk_alg1(sampler_input),
           stop("Unknown method")
    )
  } else {
    if (method == "TMET") {
      sampler_input$delta <- delta
    }
    switch(method,
           "TMET" = rtmvt_tmet(sampler_input),
           "GHK"  = rtmvt_ghk(sampler_input),
           stop("Unknown method")
    )
  }
  
  
}


# Compute gradient --------------------------------------------

#' @keywords internal
#' @noRd
grad_jacprod<- function(z_delta, Condmv_Obj, a, b, retProd = TRUE) {
  n <- length(a)
  z <- z_delta[1:n]
  delta <- z_delta[(n + 1):(2 * n)]
  D <- sqrt(Condmv_Obj$cond_var)
  B <- Condmv_Obj$B
  mu_c <- as.vector(B %*% z)
  a_t_shift <- (a - mu_c) / D - delta
  b_t_shift <- (b - mu_c) / D - delta
  log_diff_cdf <- TruncatedNormal::lnNpr(a_t_shift, b_t_shift)
  pl <- exp(-0.5 * a_t_shift^2 - log_diff_cdf)  / sqrt(2 * pi)
  pu <- exp(-0.5 * b_t_shift^2 - log_diff_cdf) / sqrt(2 * pi)
  ld <- pl - pu

  # compute grad ------------------------------------------------
  dpsi_dz <- as.vector(Matrix::t(B) %*%
                         (delta / D + ld / D)) - delta / D
  dpsi_ddelta <- delta - (z - mu_c) / D + ld

  # build return list ----------------------------------------------
  rslt <- list(
    grad = c(dpsi_dz, dpsi_ddelta)
  )
  if (retProd) {

    a_t_shift[is.infinite(a_t_shift)] <- 0
    b_t_shift[is.infinite(b_t_shift)] <- 0
    LD <- (-ld^2) + a_t_shift * pl - b_t_shift * pu

    H11_dpsi_dz <- as.vector(Matrix::t(LD / D / D * as.vector(B %*% dpsi_dz)) %*%B)
    H12_dpsi_ddelta <- as.vector(Matrix::t((dpsi_ddelta + LD * dpsi_ddelta) / D) %*% B) - dpsi_ddelta / D
    H21_dpsi_dz <-  as.vector(B %*% dpsi_dz) / D * (1 + LD) - dpsi_dz / D
    H22_dpsi_ddelta  <- ((1 + LD) * dpsi_ddelta)
    rslt$jac_grad <- c(
      H11_dpsi_dz + H12_dpsi_ddelta,
      H21_dpsi_dz + H22_dpsi_ddelta
    )
  }
  return(rslt)
}



grad_jac_psiT  <- function(z_delta, Condmv_Obj, a, b, nu, deriv = c("grad","jac","both")) {
  deriv <- match.arg(deriv)

  n <- length(a)
  z <- delta <- numeric(n)
  z[-n] <- z_delta[1:(n - 1)]
  w <- z_delta[n]
  delta[-n] <- z_delta[(n + 1):(2 * n - 1)]
  kap <- z_delta[2 * n]

  a <- a / sqrt(nu)
  b <- b / sqrt(nu)
  B <- Condmv_Obj$B
  D <- sqrt(Condmv_Obj$cond_var)

  mu_c <- as.vector(B %*% z)
  a_tilde_shift <- (a * w - mu_c) / D - delta
  b_tilde_shift <- (b * w - mu_c) / D - delta
  
  # if (any(!(a_tilde_shift < b_tilde_shift)) || !is.finite(w) || w <= 0) {
  #   bad_grad <- rep(NaN, 2 * n)
  #   bad_jac  <- Matrix::Diagonal(x = rep(NaN, 2 * n))
  #   if (deriv == "grad") return(bad_grad)
  #   if (deriv == "jac")  return(bad_jac)
  #   return(list(grad = bad_grad, jac = bad_jac))
  # }
  

  log_diff_cdf <- TruncatedNormal::lnNpr(a_tilde_shift, b_tilde_shift)
  
  pl <- exp(-0.5 * a_tilde_shift^2 - log_diff_cdf) / sqrt(2 * pi)
  pu <- exp(-0.5 * b_tilde_shift^2 - log_diff_cdf) / sqrt(2 * pi)
  ld <- pl - pu
  
  
  ## gradient ------------------------------------------------------
  if (deriv %in% c("grad","both")) {
    dpsi_dz <- as.vector(Matrix::t(B) %*% (delta / D + ld / D)) - delta / D
    dpsi_ddelta <- delta - (z - mu_c) / D + ld
    dpsi_dw <- (nu - 1) / w - kap + sum((b/D) * pu - (a/D) * pl)
    dpsi_dk <- kap - w + exp(stats::dnorm(kap, log = TRUE) -
                               TruncatedNormal::lnNpr(-kap, Inf))
    grad <- c(dpsi_dz[-n], dpsi_dw, dpsi_ddelta[-n], dpsi_dk)
  } else {
    grad <- NULL
  }

  ## Jacobian ------------------------------------------------------
  if (deriv %in% c("jac","both")) {
    dld_dw <- (-(a/D) * a_tilde_shift * pl) + (b/D) * b_tilde_shift * pu +
      (a/D) * pl^2 + (b/D) * pu^2 -
      ((a + b)/D) * (pl * pu)

    LD <- (-ld^2) + a_tilde_shift * pl - b_tilde_shift * pu

    # sparse helpers
    LD_diag   <- Matrix::Diagonal(x = LD)
    Dinv_diag <- Matrix::Diagonal(x = 1/D)
    LDI_diag  <- Matrix::Diagonal(x = 1 + LD)

    # blocks
    H11 <- Matrix::t(B) %*% Dinv_diag %*% LD_diag %*% Dinv_diag %*% B
    H12 <- Matrix::t(B) %*% Dinv_diag %*% (dld_dw)
    H13 <- Matrix::t(B) %*% (LDI_diag %*% Dinv_diag) - Dinv_diag
    H14 <- Matrix(0, nrow = n, ncol = 1, sparse = TRUE)
    H11 <- H11[-n, -n]
    H12 <- H12[-n, , drop = FALSE]
    H13 <- H13[-n, -n]
    H14 <- H14[-n, , drop = FALSE]
    first_row <- cbind(H11, H12, H13, H14)

    H21 <- Matrix(Matrix::t(H12), sparse = TRUE)
    H22 <- Matrix( -(nu - 1) / (w^2) +
                     sum((a/D)^2 * a_tilde_shift * pl -
                           (b/D)^2 * b_tilde_shift * pu -
                           (b/D * pu - a/D * pl)^2),
                   nrow = 1, ncol = 1, sparse = TRUE)
    H23 <- Matrix(dld_dw, nrow = 1, sparse = TRUE)
    H23 <- H23[, -n, drop = FALSE]
    H24 <- Matrix(-1, nrow = 1, ncol = 1, sparse = TRUE)
    second_row <- cbind(H21, H22, H23, H24)

    H31 <- Matrix(Matrix::t(H13), sparse = TRUE)
    H32 <- Matrix(dld_dw, nrow = n, ncol = 1, sparse = TRUE)
    H33 <- LDI_diag
    H34 <- Matrix(0, nrow = n, ncol = 1, sparse = TRUE)
    H32 <- H32[-n, , drop = FALSE]
    H33 <- H33[-n, -n]
    H34 <- H34[-n, , drop = FALSE]
    third_row <- cbind(H31, H32, H33, H34)

    H41 <- Matrix(0, nrow = 1, ncol = n-1, sparse = TRUE)
    H42 <- Matrix(-1, nrow = 1, ncol = 1, sparse = TRUE)
    H43 <- Matrix(0, nrow = 1, ncol = n-1, sparse = TRUE)
    H44 <- Matrix(-kap * exp(stats::dnorm(kap, log = TRUE) -
                               TruncatedNormal::lnNpr(-kap, Inf)) -
                    exp(stats::dnorm(kap, log = TRUE) -
                          TruncatedNormal::lnNpr(-kap, Inf))^2 + 1,
                  nrow = 1, ncol = 1, sparse = TRUE)
    fourth_row <- cbind(H41, H42, H43, H44)

    J <- rbind(first_row, second_row, third_row, fourth_row)
  } else {
    J <- NULL
  }

  ## return
  if (deriv == "grad") return(grad)
  if (deriv == "jac")  return(J)
  list(grad = grad, jac = J)
}


cg_solve <- function(Hfun, b, Mdiag = NULL, tol = 1e-6, maxiter = 1000, verbose = FALSE) {
  n <- length(b)
  x <- rep(0, n)        # initial guess
  r <- b - Hfun(x)      # residual

  # Jacobi preconditioner: diagonal of A
  if (is.null(Mdiag)) {
    z <- r              # no preconditioning
  } else {
    z <- r / Mdiag
  }

  p <- z
  rzold <- sum(r * z)

  for (i in seq_len(maxiter)) {
    Ap <- Hfun(p)
    alpha <- rzold / sum(p * Ap)

    x <- x + alpha * p
    r <- r - alpha * Ap

    if (sqrt(sum(r * r)) < tol) {
      if (verbose) cat(sprintf("   CG converged at iter %d: |r|=%.3e\n", i, sqrt(sum(r*r))))
      break
    }

    if (is.null(Mdiag)) {
      z <- r
    } else {
      z <- r / Mdiag
    }

    rznew <- sum(r * z)
    beta <- rznew / rzold
    p <- z + beta * p
    rzold <- rznew

    if (verbose) {
      cat(sprintf("   CG iter %d: |r|=%.3e\n", i, sqrt(sum(r*r))))
    }
  }
  x
}


###better solver
# lm_sparse_solver <- function(x0, Condmv_Obj, a, b, nu,
#                              maxit = 50, tol = 1e-8,
#                              lambda0 = 1e-3, verbose = FALSE) {
#   x <- x0
#   lambda <- lambda0
# 
#   # --- Scaling block (clipped) ------------------------------------
# 
#   # scale_vec <- c(1 / D_clip[-length(D_clip)], 1,
#   #                1 / D_clip[-length(D_clip)], 1)
#   # scale_vec[!is.finite(scale_vec)] <- 1
#   # x <- x * scale_vec
# 
#   # --- Feasibility check ------------------------------------------
#   check_feasible <- function(z_delta) {
#     n <- length(a)
#     z <- delta <- rep(0, n)
#     z[-n] <- z_delta[1:(n-1)]
#     w <- z_delta[n]
#     delta[-n] <- z_delta[(n+1):(2*n-1)]
#     kap <- z_delta[2*n]
# 
#     a_scaled <- a / sqrt(nu)
#     b_scaled <- b / sqrt(nu)
#     B <- Condmv_Obj$B
#     D <- sqrt(Condmv_Obj$cond_var)
#     mu_c <- as.vector(B %*% z)
# 
#     a_tilde_shift <- (a_scaled * w - mu_c) / D - delta
#     b_tilde_shift <- (b_scaled * w - mu_c) / D - delta
# 
#     bad <- which(a_tilde_shift >= b_tilde_shift - 1e-6 * abs(D) | w < 0)
#     list(ok = length(bad) == 0, bad = bad)
#   }
# 
#   # --- Main loop ---------------------------------------------------
#   for (k in seq_len(maxit)) {
#     gj <- grad_jac_psiT(x, Condmv_Obj, a, b, nu, deriv = "both")
#     f  <- gj$grad
#     J  <- gj$jac
#     g  <- as.vector(Matrix::t(J) %*% f)
# 
#     grad_norm <- sqrt(sum(g^2))
#     if (grad_norm < tol) {
#       if (verbose) cat("Converged (small grad) at iter", k, "\n")
#       # x <- x / scale_vec
#       return(list(x = x, fval = f, iter = k, lambda = lambda, converged = TRUE))
#     }
# 
#     # (JᵀJ + λI) step = -Jᵀf  solved by CG
#     Hfun <- function(v) as.numeric(Matrix::t(J) %*% (J %*% v) + lambda * v)
#     Mdiag <- Matrix::colSums(J^2) + lambda
#     rhs <- -g
#     system.time({
#     cg_tol_iter <- max(3e-4 * (0.9)^(k - 1), 1e-7)
#     step <- cg_solve(Hfun, rhs, Mdiag = Mdiag, tol = cg_tol_iter, maxiter = 500)})
# 
# 
# 
# 
#     # --- Feasibility-only line search -------------------------------
#     t <- 1
#     repeat {
#       chk <- check_feasible(x + t * step)
#       if (chk$ok || t < 1e-8) break
#       t <- t / 2
#     }
# 
#     ts <- t * step
# 
#     # LM-consistent predicted & actual decrease
#     pred <- 0.5 * sum(ts * (lambda * ts - g))        # always ≥0 for LM
#     f_new <- grad_jac_psiT(x + ts, Condmv_Obj, a, b, nu, deriv = "grad")
#     phi_old <- 0.5 * sum(f * f)
#     phi_new <- 0.5 * sum(f_new * f_new)
#     act <- phi_old - phi_new
#     rho <- if (pred > 0) act / pred else -Inf
# 
#     # --- Accept/reject & λ update ----------------------------------
#     if (rho > 0) {
#       x <- x + ts
#       f <- f_new
#       if (rho > 0.75)      lambda <- lambda * 0.5
#       else if (rho < 0.25) lambda <- lambda * 2.0
#     } else {
#       lambda <- lambda * 2.0
#     }
# 
#     if (verbose) {
#       cat(sprintf("Iter %d: |g|=%.3e, step=%.2e, t=%.2f, λ=%.3e, rho=%.3f\n",
#                   k, grad_norm, norm(step, "2"), t, lambda, rho))
#     }
#   }
# 
#   # --- Return if not converged ------------------------------------
#   # x <- x / scale_vec
#   list(x = x, fval = f, iter = maxit, lambda = lambda, converged = FALSE)
# }
# 



lm_sparse_solver <- function(x0, jac, Condmv_Obj, a, b, nu,
                             maxit=50, tol=1e-8,
                             lambda0=1e-3, verbose=FALSE) {
  x <- x0
  lambda <- lambda0


  check_feasible <- function(z_delta) {
    n <- length(a)
    z <- delta <- rep(0, n)
    z[-n] <- z_delta[1:(n-1)]
    w <- z_delta[n]
    delta[-n] <- z_delta[(n+1):(2*n-1)]
    kap <- z_delta[2*n]

    a_scaled <- a / sqrt(nu)
    b_scaled <- b / sqrt(nu)
    B <- Condmv_Obj$B
    D <- sqrt(Condmv_Obj$cond_var)
    mu_c <- as.vector(B %*% z)

    a_tilde_shift <- (a_scaled * w - mu_c) / D - delta
    b_tilde_shift <- (b_scaled * w - mu_c) / D - delta

    bad <- which(a_tilde_shift >= b_tilde_shift - 1e-6 * abs(D) | w < 0)
    list(ok = length(bad) == 0, bad = bad)
  }




  for (k in seq_len(maxit)) {
    gj <- grad_jac_psiT(x, Condmv_Obj, a, b, nu, deriv="both")
    f  <- gj$grad      # gradient vector (∇Φ = Jᵀψ)
    J  <- gj$jac       # Jacobian (sparse)



    # g <- as.vector(crossprod(J, f))
    g <- as.vector(Matrix::t(J) %*% f)

    if (sqrt(sum(g^2)) < tol) {
      if (verbose) cat("Converged at iter", k, "\n")
      return(list(x=x, fval=f, iter=k, lambda=lambda, converged=TRUE))
    }

    # define Hfun for CG: (J^T J + lambda I) v
      Hfun <- function(v) { as.numeric(Matrix::t(J) %*% (J %*% v) + lambda * v) }
#
    Mdiag <- Matrix::colSums(J^2) + lambda

    rhs <- -g
    step <- cg_solve(Hfun, rhs, Mdiag = Mdiag, tol = 1e-6, maxiter = 500, verbose = FALSE)

    # backtracking line search with feasibility
    t <- 1
    repeat {
      chk <- check_feasible(x + t*step)
      if (chk$ok) break
      if (verbose) cat("Infeasible indices at iter", k, ":", chk$bad, "\n")
      t <- t/2
      if (t < 1e-8) break
    }

    x_new <- x + t*step
    f_new <- grad_jac_psiT(x_new, Condmv_Obj, a, b, nu, deriv="grad")
    rho <- sum(f^2) - sum(f_new^2)

    if (rho > 0) {
      x <- x_new
      lambda <- lambda / 2
    } else {
      lambda <- lambda * 2
    }

    if (verbose) {
      cat(sprintf("Iter %d: |grad|=%.3e, step=%.2e, lambda=%.3e\n",
                  k, sqrt(sum(g^2)), norm(step, "2"), lambda))
    }
  }

  list(x=x, fval=fn(x), iter=maxit, lambda=lambda, converged=FALSE)
}
#
#



