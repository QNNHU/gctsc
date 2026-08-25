
test_that("(Poisson AR(1)) ", {
  mu <- 10
  phi <- 0.2
  theta <- 0
  tau <- c(phi)
  arma_order <- c(1, 0)
  n <- 500
  seed <- 42

  expect_no_error({
  sim_data <- sim_poisson(mu, tau, arma_order, nsim = n, seed = 1)
  y <- sim_data$y
  fit <- gctsc(formula = y~1, data = data.frame(y),
                   marginal = gctsc::poisson.marg(link ="log"),
                   cormat = gctsc::arma.cormat(p = 1, q = 0),
                   method = "CE",  
                   options = gctsc.opts(seed = 42, M = 1000),
                   QMC = TRUE)
  summary(fit)
  par(mfrow= c(2,3))
  plot(fit)
  predict(fit)
})
})

#
test_that("(NB ARMA(1,1)) ", {
  mu <- 10
  dispersion <- 1
  phi <- 0.8
  theta <- 0.2
  tau <- c(phi,theta)
  arma_order <- c(1, 1)
  n <- 500
  seed <- 42

  expect_no_error({
  sim_data <-sim_negbin(mu, dispersion, tau, arma_order, nsim = n, seed = 1)

  y <- sim_data$y
  x <- as.matrix(rep(1, n))
  fit <- gctsc(formula = y~1, data = data.frame(y),
                   marginal = negbin.marg(lambda.lower = c(0,0), link="identity"),
                   cormat = arma.cormat(p = 1, q = 1),
                   method = "CE",  
                   options = gctsc.opts(seed = 1, M = 1000),
                   QMC = TRUE)
  summary(fit)
  par(mfrow= c(2,3))
  plot(fit)
  predict(fit)
  })
})
#
test_that("(Poisson AR(1)) with covariate", {
  n <- 500
  beta <- c(0.1,0.3,1,1)
  arma_order <- c(1,0)
  xi <- numeric(n)
  zeta <- rnorm(n)
  seed <-1
  phi <-0.8
  tau <- c(phi)
  # Generate the covariate xi
  for (j in 3:n) {
    xi[j] <- 0.6 * xi[j - 1] - 0.4 * xi[j - 2] + zeta[j]
  }

  X <- as.matrix(data.frame(x1 = rep(1,n) , x2 =  sin(2 * pi * (1:n) / 12),
                            x3 = cos(2 * pi * (1:n) / 12), x4 = xi))

  # Generate the response variable Y based on the model
  mu <- exp(X%*%beta)
  
  expect_no_error({
  sim_data <- sim_poisson(mu, tau, arma_order, nsim = n, seed = 1)

  Y<- sim_data$y
  data <- data.frame(cbind(Y,X))
  dta_train <- data[1:499,]
  fit <- gctsc(formula = Y~x2+x3+x4, data=data,
                   marginal = poisson.marg(link  = "log"),
                   cormat = arma.cormat(p = 1, q = 0),
                   method = "CE",  
                   options = gctsc.opts(seed = 1))
  summary(fit)
  par(mfrow= c(2,3))
  plot(fit)
  predict(fit, newdata = data[500,], y_max = max(Y)+1)

})
  
})
