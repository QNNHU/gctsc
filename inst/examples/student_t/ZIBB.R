## -------------------------------
## Example: Zero-Inflated Beta-Binomial AR(1) model
## -------------------------------

library(gctsc)
## --- Parameter setup ---
n    <- 500
size <- 24
phi  <- 0.5
tau  <- c(phi)
arma_order <- c(1, 0)

## True parameters
beta_mu  <- 1.2                 
prob     <- plogis(beta_mu)   
rho      <- 0.1                
pi0      <- 0.2                
df= 10

## --- Simulate ZIBB time series ---
set.seed(7)
sim_data <- sim_zibb(
  prob       = rep(prob, n),
  rho        = rho,
  pi0        = rep(pi0, n),
  size       = size,
  tau        = tau,
  arma_order = arma_order,
  family = "t",
  df= 10,
  nsim       = n
)
y <- sim_data$y

X <- list(mu = matrix(1, nrow = n), pi0 = matrix(1, nrow = n))

## --- Compute truncation bounds ---
lambda <- c(qlogis(prob), qlogis(rho), qlogis(pi0))

marg <- zibb.marg(link = "logit", size= size)
ab <- marg$bounds(y, X, lambda, family ="t", df= df)

## --- Likelihood approximation ---
llk_tmet <- pmvt(lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order, 
                      pm = 30, QMC = TRUE, df= df, method = "TMET")

llk_ghk  <- pmvt( lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order,
                      QMC = TRUE, df= df, method = "GHK")

c(TMET = llk_tmet, GHK = llk_ghk)
## --- Fit ZIBB copula model using TMET ---
system.time(
fit_zibb <- gctsc(
  formula  = list(mu = y ~ 1, pi0 = ~ 1), data= data.frame(y),
  marginal = zibb.marg(link = "logit", size = size),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "TMET", family = "t", df= 10,
  options  = gctsc.opts(seed = 7)
))

summary(fit_zibb)
plot(fit_zibb)
predict(fit_zibb)

## -------------------------------
## Example: ZIBB AR(1) with seasonal zero-inflation π₀(t)
## -------------------------------
library(gctsc)

## --- Parameter setup ---
n    <- 1000
size <- 24
phi  <- 0.85
tau  <- c(phi)
arma_order <- c(1, 0)
df= 10

prob <- plogis(0.2)        # constant prob(t)
rho  <- 0.16               # ICC for BB component

## Seasonal covariates for π₀(t)
X_pi <- cbind(
  1,
  sin(2*pi*(1:n)/12),
  cos(2*pi*(1:n)/12)
)
colnames(X_pi) <- c("int", "sin", "cos")

beta_pi <- c(1, 2, 5)
pi0 <- plogis(X_pi %*% beta_pi)

## --- Simulate ZIBB data ---
set.seed(1)
sim_data <- sim_zibb(
  prob       = rep(prob, n),
  rho        = rho,
  pi0        = pi0,
  size       = size,
  tau        = tau,
  arma_order = arma_order,
  family = "t",
  df= 10,
  nsim       = n
)
y <- sim_data$y

X <- list(mu = matrix(1, nrow = n), pi0 = as.matrix(X_pi))

## --- Compute truncation bounds ---
lambda <- c(qlogis(prob), qlogis(rho), beta_pi)

marg <- zibb.marg(link = "logit", size= size)
ab <- marg$bounds(y, X, lambda, family ="t", df= df)

## --- Likelihood approximation ---
llk_tmet <- pmvt(lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order, 
                      pm = 30, QMC = TRUE, df= df, method  = "TMET")

llk_ghk  <- pmvt( lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order,
                      QMC = TRUE, df= df, method = "GHK")

c(TMET = llk_tmet, GHK = llk_ghk)
## --- Fit using formula interface ---
data <- data.frame(y = y, X_pi)

system.time(
fit_zibb_seasonal <- gctsc(
  formula = list(mu = y ~ 1,
                 pi0 = ~ sin + cos),
  data     = data,
  marginal = zibb.marg(link = "logit", size = size),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "TMET",family = "t", df= 10,
  options  = gctsc.opts(seed = 1)
))

summary(fit_zibb_seasonal)
plot(fit_zibb_seasonal)
predict(fit_zibb_seasonal, newdata  = data[200, ])


## -------------------------------
## Example: ZIBB AR(1) with covariates in μ(t) and π₀(t)
## -------------------------------


library(gctsc)

## --- Parameter setup ---
n    <- 1500
size <- 10
phi  <- 0.5
tau  <- c(phi)
arma_order <- c(1, 0)
rho  <- 0.3
df= 10

## --- Covariates for μ(t) ---
X_mu <- cbind(
  1,
  sin(2*pi*(1:n)/365),
  cos(2*pi*(1:n)/365)
)
colnames(X_mu) <- c("int", "sin", "cos")

beta_mu <- c(1.2, 2, 3)
prob <- plogis(X_mu %*% beta_mu)

## --- Covariates for π₀(t) ---
day_of_year <- rep(1:365, length.out = n)
season <- factor(ifelse(day_of_year < 100, "Winter",
                        ifelse(day_of_year < 180, "Spring",
                               ifelse(day_of_year < 270, "Summer", "Fall"))))
X_pi <- model.matrix(~ season)

alpha_pi <- c(0.2, 0.3, 0.5, 0.7)
pi0 <- plogis(X_pi %*% alpha_pi)

## --- Simulate ZIBB data ---
set.seed(1)
sim_data <- sim_zibb(
  prob       = prob,
  rho        = rho,
  pi0        = pi0,
  size       = size,
  tau        = tau,
  arma_order = arma_order,
  nsim       = n,
  family     = "t",
  df         = df
  
)
y <- sim_data$y

X <- list(mu = as.matrix(X_mu), pi0 = as.matrix(X_pi))

## --- Compute truncation bounds ---
lambda <- c(beta_mu, qlogis(rho), alpha_pi)

marg <- zibb.marg(link = "logit", size= size)
ab <- marg$bounds(y, X, lambda, family ="t", df= df)

## --- Likelihood approximation ---
llk_tmet <- pmvt(lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order, 
                      pm = 30, QMC = TRUE, df= df, method = "TMET")

llk_ghk  <- pmvt( lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order,
                      QMC = TRUE, df= df, method = "GHK")

c(TMET = llk_tmet, GHK = llk_ghk)


## --- Fit ZIBB copula model using formula interface ---
data_df <- data.frame(y = y, X_mu, X_pi)
system.time(
fit_zibb_cov <- gctsc(
  formula = list(mu  = y ~ sin + cos,
                 pi0 = ~ seasonSpring + seasonSummer + seasonWinter),
  data     = data_df,
  marginal = zibb.marg(link = "logit", size = size),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "TMET",family = "t", df= 10,
  options  = gctsc.opts(seed = 1, M = 1000),
  start = c(beta_mu, qlogis(rho), alpha_pi, phi)
)
)
summary(fit_zibb_cov)
plot(fit_zibb_cov)

predict(fit_zibb_cov, X_test = data_df[500, ])

