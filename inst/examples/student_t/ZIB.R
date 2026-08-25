## -------------------------------
## Example: Zero-Inflated Binomial AR(1) model
## -------------------------------

library(gctsc)


## --- Parameter setup ---
n    <- 500
size <- 24
prob <- 0.2           # Binomial success probability
pi0  <- 0.2           # Zero-inflation probability
phi  <- 0.8           # AR(1) parameter
tau  <- c(phi)
arma_order <- c(1, 0)
df <- 10

## --- Simulate ZIB count time series ---
set.seed(1)
sim_data <- sim_zib(
  prob       = rep(prob, n),
  pi0        = rep(pi0, n),
  size       = size,
  tau        = tau,
  arma_order = arma_order,
  family     = "t", 
  df         = 10,
  nsim       = n
)
y <- sim_data$y



X <- list(mu = matrix(1, nrow = n), pi0 = matrix(1, nrow = n))

## --- Compute truncation bounds ---
lambda <- c(qlogis(prob), qlogis(pi0))

marg <- zib.marg(link = "logit", size= size)
ab <- marg$bounds(y, X, lambda, family ="t", df= df)

## --- Likelihood approximation ---
llk_tmet <- pmvt(lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order, 
                      pm = 30, QMC = TRUE, df= df, method = "TMET")

llk_ghk  <- pmvt( lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order,
                      QMC = TRUE, df= df, method ="GHK")

c(TMET = llk_tmet, GHK = llk_ghk)

## --- Fit t copula ZIB model using TMET ---
system.time(
fit_zib <- gctsc(
  formula  = list(
    mu  = y ~ 1,  
    pi0 = ~ 1      
  ), data= data.frame(y),
  marginal = zib.marg(link = "logit", size = size),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "TMET", family = "t", df= 10,
  options  = gctsc.opts(seed = 1, M = c(100,1000))
))

summary(fit_zib)
plot(fit_zib)
predict(fit_zib)


## -------------------------------
## Example: Zero-Inflated Binomial AR(1) with covariates
## -------------------------------

library(gctsc)


n    <- 2000
size <- 24
prob <- 0.2
phi  <- 0.8
tau  <- c(phi)
arma_order <- c(1, 0)
df <- 10
## --- Construct seasonal covariates for π0(t) ---
day_of_year <- rep(1:365, length.out = n)
season <- factor(ifelse(day_of_year < 100, "Winter",
                        ifelse(day_of_year < 180, "Spring",
                               ifelse(day_of_year < 270, "Summer", "Fall"))))

X_pi <- model.matrix(~ season)
colnames(X_pi) <- make.names(colnames(X_pi))

## True zero-inflation function
beta_pi  <- c(0.2, 1, 0.3, 0.5)
logit_pi <- X_pi %*% beta_pi
pi0      <- plogis(logit_pi)

## --- Simulate ZIB time series ---
set.seed(1)
sim_data <- sim_zib(
  prob       = rep(prob, n),
  pi0        = pi0,
  size       = size,
  tau        = tau,
  arma_order = arma_order, family = "t", df= 10,
  nsim       = n
)
y <- sim_data$y

dta <- data.frame(y = y, X_pi)


X <- list(mu =  matrix(1, nrow = n), pi0 = X_pi)

## --- Compute truncation bounds ---
lambda <- c(qlogis(prob), beta_pi)

marg <- zib.marg(link = "logit", size= size)
ab <- marg$bounds(y, X, lambda, family ="t", df= df)

## --- Likelihood approximation ---
llk_tmet <- pmvt(lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order, 
                      pm = 30, QMC = TRUE, df= df, method ="TMET")

llk_ghk  <- pmvt( lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order,
                      QMC = TRUE, df= df, method ="GHK")

c(TMET = llk_tmet, GHK = llk_ghk)


## --- Fit t copula ZIB model (TMET) ---
system.time(
fit_zib_cov <- gctsc(
  formula  = list(
    mu  = y ~ 1,
    pi0 = ~ seasonSpring + seasonSummer + seasonWinter
  ),
  data     = dta,
  marginal = zib.marg(link = "logit", size = size),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "TMET", family = "t", df= 10,
  options  = gctsc.opts(seed = 1, M = c(1000))
))

summary(fit_zib_cov)
plot(fit_zib_cov)
predict(fit_zib_cov, newdata  = dta[2000, ])
