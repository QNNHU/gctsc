## -------------------------------
## Example: Poisson AR(1) model
## -------------------------------

## --- Parameter setup ---
n <- 2000
mu <- 2
phi <- 0.8
arma_order <- c(1, 0)
tau <- c(phi)
family <- "gaussian"
## --- Simulate data ---
set.seed(7)
sim_data <- sim_poisson(mu = mu,
                        tau = tau,
                        arma_order = arma_order,
                        family = "gaussian",
                        nsim = n)
y <- sim_data$y
X <- matrix(1, nrow = n)
x <- list(mu = X)

## --- Compute truncation bounds ---
marg <- poisson.marg(link = "identity")

ab <- marg$bounds(y, x, mu, family = family)

## --- Likelihood approximation ---
llk_tmet <- pmvn(lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order, method= "TMET",
                      pm = 30, QMC = TRUE)

llk_ghk  <- pmvn( lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order, method ="GHK",
                      QMC = TRUE)

llk_ce  <- pmvn( lower = ab[,1], upper = ab[,2],
                  tau = tau, od = arma_order, method ="CE",
                  QMC = TRUE)

c(TMET = llk_tmet, GHK = llk_ghk, CE = llk_ce)

## --- Fit Gaussian copula model using CE ---
system.time(
fit <- gctsc(
  formula = y ~ 1, data= data.frame(y),
  marginal = poisson.marg(link ="log"),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "CE",
  family   = "gaussian",
  QMC      = TRUE,
  options = gctsc.opts(seed=1)
))


plot(fit)       # residual diagnostics
predict(fit)    # one-step forecasting

## --- Fit Gaussian copula model using GHK ---
system.time(
fit_GHK <- gctsc(
  formula = y ~ 1, data= data.frame(y),
  marginal = poisson.marg(),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "GHK",
  family   = "gaussian",
  QMC      = TRUE,
  options = gctsc.opts(seed=1)
))

plot(fit_GHK)

## -------------------------------
## Example: Poisson AR(1) model with covariates
## -------------------------------

n <- 1000
phi <- 0.8
tau <- c(phi)
arma_order <- c(1, 0)

## --- Generate covariates (seasonal + autoregressive) ---
set.seed(1)
zeta <- rnorm(n)
xi <- numeric(n)
for (j in 3:n) {
  xi[j] <- 0.6 * xi[j - 1] - 0.4 * xi[j - 2] + zeta[j]
}

X <- as.matrix(data.frame(
  x1 = 1,
  x2 = sin(2 * pi * (1:n) / 12),
  x3 = cos(2 * pi * (1:n) / 12),
  x4 = xi
))


beta <- c(0.1, 0.3, 1, 3)
mu <- exp(X %*% beta)

## --- Simulate Poisson response ---
sim_data <- sim_poisson(mu = mu, tau = tau,
                        arma_order = arma_order, family   = "gaussian",
                        nsim = n, seed = 1)
y <- sim_data$y

## --- Compute bounds and log-likelihood approximations ---
marginal <- gctsc::poisson.marg(link = "log")
ab <- marginal$bounds(y, list(mu = X), beta, family   = "gaussian")
system.time(
llk_tmet_qmc <- pmvn(
  lower = ab[, 1],
  upper = ab[, 2],
  tau   = tau,
  od    = arma_order,
  method = "TMET"
))

## --- Fit Gaussian copula model ---
data_df <- data.frame(Y = y, X)
data_train <- data_df[1:1000,]
system.time(
fit1 <- gctsc(
  formula  = Y ~ x2 + x3 + x4,
  data     = data_train,
  marginal = gctsc::poisson.marg(link = "log"),
  cormat   = gctsc::arma.cormat(p = 1, q = 0),
  family   = "gaussian",
  method   = "CE",
  options  = gctsc.opts(seed = 1, M = c(100,1000))
))



summary(fit1)
plot(fit1)
residuals(fit1)
predict(fit1, newdata= data_df[1000,], y_max = max(y))
