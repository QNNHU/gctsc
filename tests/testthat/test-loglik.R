test_that(" Poisson AR(1)", {
  mu <- 100
  phi <- 0.2
  theta <- 0
  tau <- c(phi)
  arma_order <- c(1, 0)
  n <- 500
  seed <- 1

  sim_data <- sim_poisson(mu, tau, arma_order, nsim = n, seed = 1)
  y <- sim_data$y
  X <- as.matrix(rep(1, n))
  marginal <- poisson.marg(link = "identity")
  ab <- marginal$bounds(y, list(mu = as.matrix(X)),lambda = mu)
  llk_tmet_qmc <- pmvn(lower = ab[,1], upper= ab[,2], tau=tau, 
                        od = arma_order, method = "TMET", QMC = TRUE)
  llk_tmet_mc  <- pmvn(lower = ab[,1], upper= ab[,2], tau=tau,
                       od = arma_order,method = "TMET", QMC = FALSE)
  llk_ghk_mc   <- pmvn(lower = ab[,1], upper= ab[,2], tau=tau, 
                       od = arma_order, method = "GHK", QMC = TRUE)
  llk_ce       <- pmvn(lower = ab[,1], upper= ab[,2], tau = tau,
                       c = 0.5, od = arma_order, method = "CE")

})

test_that(" negative binomial ARMA(1,1)", {
  mu <- 10
  dispersion <- 2
  phi <- 0.8
  theta <- 0.2
  tau <- c(phi,theta)
  arma_order <- c(1, 1)
  n <- 500
  seed <- 42

  sim_data <-sim_negbin(mu, dispersion, tau, arma_order, nsim = n, seed = 1)

  y <- sim_data$y
  X <- as.matrix(rep(1, n))
  marginal <- negbin.marg(link = "identity")
  ab <- marginal$bounds(y, list(mu = as.matrix(X)),lambda = mu)
  llk_tmet_qmc <- pmvn(lower = ab[,1], upper= ab[,2], tau=tau, 
                       od = arma_order, method = "TMET", QMC = TRUE)
  llk_tmet_mc  <- pmvn(lower = ab[,1], upper= ab[,2], tau=tau,
                       od = arma_order,method = "TMET", QMC = FALSE)
  llk_ghk_mc   <- pmvn(lower = ab[,1], upper= ab[,2], tau=tau, 
                       od = arma_order, method = "GHK", QMC = TRUE)
  llk_ce       <- pmvn(lower = ab[,1], upper= ab[,2], tau = tau,
                       c = 0.5, od = arma_order, method = "CE")
  


})


test_that(" Poisson AR(1) with covariate", {
  n <- 1000
  beta <- c(1,3,0.5,1)
  arma_order <- c(1,0)
  xi <- numeric(n)
  zeta <- rnorm(n)
  seed <-1
  phi <-0.3
  tau <- c(phi)
  # Generate the covariate xi
  for (j in 3:n) {
    xi[j] <- 0.6 * xi[j - 1] - 0.4 * xi[j - 2] + zeta[j]
  }

  X <- as.matrix(data.frame(x1 = rep(1,n) , x2 =  sin(2 * pi * (1:n) / 12),
                            x3 = cos(2 * pi * (1:n) / 12), x4 = xi))

  # Generate the response variable Y based on the model
  mu <- exp(X%*%beta)

  sim_data <- sim_poisson(mu, tau, arma_order, nsim = n, seed = 1)

  Y<- sim_data$y
  marginal <- poisson.marg(link="log")
  ab <- marginal$bounds(Y, list(mu = as.matrix(X)),lambda = beta)
  
  data <- data.frame(cbind(y = Y, X))
  llk_tmet_qmc <- pmvn(lower = ab[,1], upper= ab[,2], tau=tau, 
                       od = arma_order, method = "TMET", QMC = TRUE)
  llk_tmet_mc  <- pmvn(lower = ab[,1], upper= ab[,2], tau=tau,
                       od = arma_order,method = "TMET", QMC = FALSE)
  llk_ghk_mc   <- pmvn(lower = ab[,1], upper= ab[,2], tau=tau, 
                       od = arma_order, method = "GHK", QMC = TRUE)
  llk_ce       <- pmvn(lower = ab[,1], upper= ab[,2], tau = tau,
                       c = 0.5, od = arma_order, method = "CE")
})



