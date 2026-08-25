## -------------------------------
## Example: Zero-Inflated Poisson AR(1) model
## -------------------------------

library(gctsc)

## --- Parameter setup ---
n   <- 500
mu  <- 10               # Poisson mean
pi0 <- 0.5              # zero-inflation probability (probability scale)
phi <- 0.5
tau <- c(phi)
arma_order <- c(1, 0)
df <- 10

## --- Simulate ZIP data (pi0 on probability scale) ---
set.seed(1)
sim_data <- sim_zip(
  mu         = rep(mu, n),
  pi0        = rep(pi0, n),
  tau        = tau,
  arma_order = arma_order,
  family     = "t",
  df         = df,
  nsim       = n
)
y <- sim_data$y

X <- list(mu = as.matrix(rep(1,n)), pi0 = as.matrix(rep(1,n)))

## --- Compute truncation bounds ---
marg <- zip.marg(link = "identity")
ab <- marg$bounds(y, X, c(mu,pi0),family ="t", df= df)

## --- Likelihood approximation ---
llk_tmet <- pmvt(lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order, 
                      pm = 30, QMC = TRUE, df= df, method ="TMET")

llk_ghk  <- pmvt( lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order,
                      QMC = TRUE, df= df, method = "GHK")

c(TMET = llk_tmet, GHK = llk_ghk)

## --- Fit ZIP t copula model using TMET ---

fit_zip <- gctsc(
  formula  = list(mu = y ~ 1, pi0 = ~ 1), data= data.frame(y),
  marginal = zip.marg(link = "identity"),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "TMET",
  family = "t",
  df= df,
  options  = gctsc.opts(seed = 1, M = 1000)
)

summary(fit_zip)
plot(fit_zip)
predict(fit_zip)
