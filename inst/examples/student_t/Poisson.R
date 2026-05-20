## -------------------------------
## Example: Poisson AR(1) model
## -------------------------------

## --- Parameter setup ---
n <- 500
mu <- 10
phi <- 0.2
arma_order <- c(1, 0)
tau <- c(phi)
family <- "t"
df= 15
## --- Simulate data ---
set.seed(7)
sim_data <- sim_poisson(mu = mu,
                        tau = tau,
                        arma_order = arma_order,
                        family = "t",
                        nsim = n, df= df)
y <- sim_data$y
X <- matrix(1, nrow = n)
x <- list(mu = X)

## --- Compute truncation bounds ---
marg <- poisson.marg(link ="log")
ab <- marg$bounds(y, x, log(mu), family = family, df= df)

## --- Likelihood approximation ---
llk_tmet <- pmvt(lower = ab[,1], upper = ab[,2],
                 tau = tau, od = arma_order, method= "TMET",
                 pm = 30, QMC = TRUE, df= df)

llk_ghk  <- pmvt( lower = ab[,1], upper = ab[,2],
                  tau = tau, od = arma_order, method ="GHK",
                  QMC = TRUE,df= df)

llk_ce  <- pmvt( lower = ab[,1], upper = ab[,2],
                 tau = tau, od = arma_order, method ="CE",
                 QMC = TRUE,df= df)

c(TMET = llk_tmet, GHK = llk_ghk, CE =llk_ce )


## --- Fit t copula model using TMET ---
system.time(
  fit_TMET <- gctsc(
    formula = y ~ 1,data= data.frame(y),
    marginal = poisson.marg(link ="log"),
    cormat   = arma.cormat(p = 1, q = 0),
    method   = "TMET",
    family   = "t",
    df= 15,
    QMC      = TRUE
  ))



plot(fit_TMET)       # residual diagnostics
predict(fit_TMET)    # one-step forecasting

## --- Fit t copula model using GHK ---
system.time(
fit_ghk <- gctsc(
  formula = y ~ 1,data= data.frame(y),
  marginal = poisson.marg(lambda.lower = 0),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "GHK",
  family   = "t",
  df= 15,
  QMC      = TRUE
))

plot(fit_ghk)

## -------------------------------
## Example: Poisson AR(1) model with covariates
## -------------------------------

n <- 500
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
                        arma_order = arma_order, family   = "t", df= df,
                        nsim = n, seed = 1)
y <- sim_data$y

## --- Compute bounds and log-likelihood approximations ---
marginal <- poisson.marg(link = "log")
ab <- marginal$bounds(y, x = list(mu= X), beta, family   = "t", df= df)

llk_tmet_qmc <- pmvt(
  lower = ab[, 1],
  upper = ab[, 2],
  tau   = tau,
  od    = arma_order,
  method ="TMET",
  df = 10
)

## --- Fit t copula model ---
data_df <- data.frame(Y = y, X)
data_train <- data_df[1:400,]
fit <- gctsc(
  formula  = Y ~ x2 + x3 + x4,
  data     = data_train,
  marginal = poisson.marg(link = "log"),
  cormat   = arma.cormat(p = 1, q = 0),
  family   = "t",
  method   = "TMET",
  df= 15,
  options  = gctsc.opts(seed = 1, M = 1000)
)

summary(fit)
plot(fit)
predict(fit, newdata = data_df[401,])
