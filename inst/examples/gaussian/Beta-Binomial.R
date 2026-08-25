## -------------------------------
## Example: Beta–Binomial AR(1) model
## -------------------------------

library(gctsc)

## --- Parameter setup ---
n    <- 500
size <- 24
beta0 <- 0.2                  # logit-scale intercept
prob <- plogis(beta0)         # simulation-scale probability
rho  <- 0.18                  # intra-class correlation
phi  <- 0.8                   # AR(1) parameter
tau  <- c(phi)
arma_order <- c(1, 0)

## --- Simulate Beta–Binomial time series ---
set.seed(1)
sim_data <- sim_bbinom(
  prob       = rep(prob, n),
  rho        = rho,
  size       = size,
  tau        = tau,
  arma_order = arma_order,
  family = "gaussian",
  nsim       = n
)
y <- sim_data$y

X <- matrix(1, nrow = n)

## --- Compute truncation bounds ---
marg <- bbinom.marg(link = "logit", size= size)
ab <- marg$bounds(y, list(mu=X), c(beta0,qlogis(rho)),family ="gaussian")

## --- Likelihood approximation ---
llk_tmet <- pmvn(lower = ab[,1], upper = ab[,2],
                 tau = tau, od = arma_order, 
                 pm = 30, QMC = TRUE,  method ="TMET")

llk_ghk  <- pmvn( lower = ab[,1], upper = ab[,2],
                  tau = tau, od = arma_order,
                  QMC = TRUE,  method = "GHK")

c(TMET = llk_tmet, GHK = llk_ghk)


## --- Fit Gaussian copula Beta–Binomial model using TMET ---
fit_bbinom <- gctsc(
  formula  = y ~ 1, data= data.frame(y),
  marginal = bbinom.marg(link = "logit", size = size),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "CE", family = "gaussian",
  options  = gctsc.opts(seed = 1)
)

summary(fit_bbinom)
plot(fit_bbinom)
predict(fit_bbinom)

## -------------------------------
## Example: Beta–Binomial AR(1) with covariates
## -------------------------------

library(gctsc)
## --- Parameter setup ---
n    <- 500
size <- 24
phi  <- 0.5
tau  <- c(phi)
arma_order <- c(1, 0)

## Overdispersion / intra-class correlation parameterization
rho <- 1 / (1 + 5)

## --- Construct covariates (seasonal + AR structure) ---
zeta <- rnorm(n)
xi   <- numeric(n)
for (j in 3:n) {
  xi[j] <- 0.6 * xi[j-1] - 0.4 * xi[j-2] + zeta[j]
}

X <- data.frame(
  x1 = 1,
  x2 = sin(2*pi*(1:n)/12),
  x3 = cos(2*pi*(1:n)/12),
  x4 = xi
)

## True logit(prob)
beta_true <- c(0.2, 0.3, 0.5, 0.3)
logit_prob <- as.matrix(X) %*% beta_true
prob <- plogis(logit_prob)

## --- Simulate Beta–Binomial time series ---
set.seed(1)
sim_data <- sim_bbinom(
  prob       = prob,
  rho        = rho,
  size       = size,
  tau        = tau,
  arma_order = arma_order,
  family = "gaussian",
  nsim       = n
)
y <- sim_data$y


## --- Compute truncation bounds ---
marg <- bbinom.marg(link = "logit", size= size)
ab <- marg$bounds(y, list(mu = as.matrix(X)), c(beta_true, qlogis(rho)),family ="gaussian")

## --- Likelihood approximation ---
llk_tmet <- pmvn(lower = ab[,1], upper = ab[,2],
                 tau = tau, od = arma_order, 
                 pm = 30, QMC = TRUE,  method = "TMET")

llk_ghk  <- pmvn( lower = ab[,1], upper = ab[,2],
                  tau = tau, od = arma_order,
                  QMC = TRUE, method ="GHK")

c(TMET = llk_tmet, GHK = llk_ghk)

## --- Fit Gaussian copula Beta–Binomial model ---
data_df <- data.frame(y = y, X)
fit_bbinom_cov <- gctsc(
  formula  = y ~ x2 + x3 + x4,
  data     = data_df,
  marginal = bbinom.marg(link = "logit", size = size),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "TMET",family = "gaussian",
  options  = gctsc.opts(seed = 1)
)

summary(fit_bbinom_cov)
plot(fit_bbinom_cov)
predict(fit_bbinom_cov, newdata =  data_df[500, ])

