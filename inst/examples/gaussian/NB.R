## -------------------------------
## Example: negative binomial ARMA(1,1) model
## -------------------------------

## --- Parameter setup ---
n <- 1000
mu <- 30
phi <- 0.5
theta <- 0
arma_order <- c(1, 0)
tau <- c(phi)
dispersion <- 0.15

## --- Simulate data ---
set.seed(7)
sim_data <- sim_negbin(mu = mu, dispersion,
                        tau = tau,
                        arma_order = arma_order, family ="gaussian",
                        nsim = n)
y <- sim_data$y
X <- matrix(1, nrow = n)

## --- Compute truncation bounds ---
marg <- negbin.marg(link = "identity")
ab <- marg$bounds(y, x = list(mu=X), c(mu,dispersion),family ="gaussian")

## --- Likelihood approximation ---
llk_tmet <- pmvn(lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order, 
                      pm = 30, QMC = TRUE, method = "TMET")

llk_ghk  <- pmvn( lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order,
                      QMC = TRUE, method ="GHK")

llk_ce  <- pmvn( lower = ab[,1], upper = ab[,2],
                  tau = tau, od = arma_order,
                  QMC = TRUE, method ="CE")


c(TMET = llk_tmet, GHK = llk_ghk, CE = llk_ce)


## --- Fit Gaussian copula model using GHK ---
system.time({
fit_ghk <- gctsc(
  formula = y ~ 1,data= as.data.frame(y),
  marginal = negbin.marg(),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "GHK",
  family = "gaussian",
  QMC      = TRUE
)})


## --- Fit Gaussian copula model using CE ---
system.time({
  fit_CE <- gctsc(
    formula = y ~ 1,data= as.data.frame(y),
    marginal = negbin.marg(),
    cormat   = arma.cormat(p = 1, q = 0),
    method   = "CE",
    family = "gaussian",
    QMC      = TRUE
  )})


plot(fit_CE)



## -------------------------------
## Example: Negative Binomial AR(1) model with covariates
## -------------------------------

library(gctsc)

n <- 500
phi <- 0.5
tau <- c(phi)
arma_order <- c(1, 0)

## Covariate regression coefficients
beta <- c(1, 0.3, 1, 0.5)

## Negative binomial dispersion (gctsc parameterization: variance = mu + mu^2/dispersion)
dispersion <- 2

## Generate covariates (seasonal + autoregressive component)
set.seed(1)
zeta <- rnorm(n)
xi <- numeric(n)
for (j in 3:n) {
  xi[j] <- 0.6 * xi[j - 1] - 0.4 * xi[j - 2] + zeta[j]
}

X <- as.matrix(data.frame(
  x1 = rep(1, n),
  x2 = sin(2 * pi * (1:n) / 12),
  x3 = cos(2 * pi * (1:n) / 12),
  x4 = xi
))

## Compute mean function
mu <- as.vector(exp(X %*% beta))

## Simulate Negative Binomial counts with AR(1) latent process
sim_data <- sim_negbin(
  mu        = mu,
  dispersion = dispersion,
  tau       = tau,
  arma_order = arma_order,
  nsim      = n,
  family = "gaussian",
  seed      = 10
)
y <- sim_data$y



## Assemble data for model fitting
data_df <- data.frame(Y = y, X)

## Fit Gaussian Copula model using CE
fit_nb <- gctsc(
  formula  = Y ~ x2 + x3 + x4,
  data     = data_df[1:499,],
  marginal = negbin.marg(link = "log"),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "CE",
  family = "gaussian",
  options  = gctsc.opts(seed = 1, M = 1000)
)

summary(fit_nb)
plot(fit_nb)

## One-step-ahead prediction
predict(fit_nb,newdata = data_df[500,])
