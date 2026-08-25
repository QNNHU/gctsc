## ---------- Poisson AR(1) ----------
library(VeccTMVN)
n <- 2
mu <- 1


phi <- 0.9
theta <- 0
arma_order <- c(1,0)
tau <- c(phi)
seed=1
df = 3

phi <- 0
theta <- -0.2
arma_order <- c(0,1)
tau <- c(theta)
seed=6
df = 3

phi <- 0.8
theta <- 0.8
arma_order <- c(1,1)
tau <- c(phi, theta)
seed=81
df = 10



sim_data <- sim_poisson(mu, tau, arma_order, nsim = n, seed = seed, family ="t", df= df)
y <- sim_data$y
# mu=0.85
# tau=0.1615 
y <- c(0,1)
x <- matrix(1, nrow = n)
marginal <- poisson.marg()
ab <- marginal$bounds(y, x, mu, family ="t", df=df)
covM <- toeplitz(stats::ARMAacf(ar = phi, ma = theta, lag.max = n - 1, pacf = FALSE))
pm =30

lower = ab[, 1]
upper = ab[, 2]
od = arma_order


# Likelihood approximation

system.time({
llk_tmet_qmc <- pmvt_tmet_optim(lower = lower, upper = upper , tau, od = od, df=df, pm =pm )})

system.time({
  llk_tmet_qmc <- pmvt_tmet(lower = lower, upper = upper , tau, od = od, pm =pm, df=df )})


 

system.time({
llk_ghk_qmc  <- pmvt_ghk(lower = ab[, 1], upper = ab[, 2], tau, od = arma_order, df= df)})
system.time({
llk_ghk_mc_alg1   <- pmvt_ghk_mvt(ab[, 1], ab[, 2], tau, od = arma_order, df= df)})
system.time({
llk_vmet     <- pmvt_vmet(lower = ab[, 1], upper = ab[, 2], tau = tau, od = arma_order, df= df)})
llk_ce       <- pmvt_ce(lower = ab[, 1], upper = ab[, 2], tau = tau, od = arma_order, df= df)
system.time({
log_likelihood <-   TruncatedNormal::pmvt( sigma = covM, df=df, lb = ab[, 1], ub =   ab[, 2], type = c("qmc"), 
                                               B = 1e4)})



## ---------- Negative Binomial ARMA(1,1) ----------
mu <- c(10)   # Mean of NegBin
phi <- c(0.8)   # AR parameter
theta <- c(0.8) # MA parameter
df <- c(10)       # df for t copula
dispersion<- c(2) # dispersion parameter for NegBin
tau <- c(phi,theta)
arma_order <- c(1, 1)
seed=82
n= 1000
sim_data <- sim_negbin(phi, dispersion, tau, arma_order, nsim = n, seed = seed, family ="t", df=df)
y <- sim_data$y
x <- matrix(1, nrow = n)

marginal <- negbin.marg()
ab <- marginal$bounds(y, x, c(mu,dispersion), family="t", df=df)
covM <- toeplitz(stats::ARMAacf(ar = phi, ma = theta, lag.max = n - 1, pacf = FALSE))


lower = ab[, 1]
upper = ab[, 2]
od = arma_order

set.seed(seed)
VeccTMVN::pmvt(lower = lower, upper = upper, sigma = covM,
               reorder = 0, NLevel1 = 1, delta = rep(0, n),
               NLevel2 = 1000, verbose = FALSE, retlog = TRUE, m = 30, df= df)



system.time({
  llk_ghk_qmc  <- pmvt_ghk(lower=lower, upper=upper, tau, od = arma_order, df= df, M=1000)})



# Likelihood approximation
system.time({
llk_tmet_optim <- pmvt_tmet_optim(lower= lower, upper = upper, tau, od = arma_order,df=df,pm =30)})
system.time({
  llk_tmet_qmc <- pmvt_tmet(lower= lower, upper = upper, tau, od = arma_order,df=df,pm =30, M=1000)})


system.time({
llk_vmet     <- pmvt_vmet(lower = lower, upper = upper, tau = tau, od = arma_order,  df= df)})
llk_ce       <- pmvt_ce(ab[, 1], ab[, 2], tau, c = 0.5, od = arma_order,  df= df)
llk_ghk_alg1 <-  pmvt_ghk_mvt(lower, upper, tau, od = arma_order, df= df, M=1000)

system.time({
  log_likelihood <-   (TruncatedNormal::pmvt( sigma = covM, df=df, lb = ab[, 1],
                                              ub =   ab[, 2], type = c("qmc"), 
                                                 B = 1e4))})



## ---------- Poisson AR(1) with Covariates ----------
n <- 1000
phi <- 0.8
tau <- c(phi)
arma_order <- c(1, 0)
df=7
# Covariate generation (seasonal + autoregressive)
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

# Simulate Poisson response
sim_data <- sim_poisson(mu, tau, arma_order, nsim = n, seed = 1, family ="t", df= df)
y <- sim_data$y

# Compute bounds and log-likelihood approximations
marginal <- poisson.marg(link = "log")
ab <- marginal$bounds(y, X, beta)

llk_tmet_qmc <- pmvt_tmet(ab[, 1], ab[, 2], tau, od = arma_order, df= df)
llk_ghk_qmc  <- pmvt_ghk(ab[, 1], ab[, 2], tau, od = arma_order, df= df)
llk_ghk_alg1   <- pmvt_ghk_mvt(ab[, 1], ab[, 2], tau, od = arma_order, QMC = FALSE, df= df)
llk_ce       <- pmvt_ce(ab[, 1], ab[, 2], tau, c = 0.5, od = arma_order, df= df)


## ---------- ZIP AR(1) ----------
n <- 1000
mu <- 10
pi0 <- 0.1
phi <- 0.5
tau <- c(phi)
arma_order <- c(1, 0)
df= 5
sim_data <- sim_zip(mu, pi0 = pi0, tau, arma_order, nsim = n, seed = 1, df= df, family ="t")
y <- sim_data$y

x_mu <- matrix(1, nrow = n)
x_pi <- matrix(1, nrow = n)
lambda <- c(mu, pi0)

marginal <- zip.marg(link = "identity")
ab <- marginal$bounds(y, x = list(mu = x_mu, pi0 = x_pi), lambda)

llk_tmet_qmc <- pmvt_tmet(ab[, 1], ab[, 2], tau, od = arma_order,df= df)
llk_ghk_qmc  <- pmvt_ghk(ab[, 1], ab[, 2], tau, od = arma_order,df= df)
llk_vmet     <- pmvt_vmet(ab, tau = tau, od = arma_order,df= df)
llk_ce       <- pmvt_ce(ab[, 1], ab[, 2], tau, c = 0.5, od = arma_order,df= df)
llk_ghk_alg1   <- pmvt_ghk_mvt(ab[, 1], ab[, 2], tau, od = arma_order, QMC = FALSE, df= df)

## ---------- ZIP AR(1) with Covariates in π₀ ----------
n <- 5000
mu <- 15
phi <- 0.5
tau <- c(phi)
arma_order <- c(1, 0)

# Seasonal factor for zero-inflation
day_of_year <- rep(1:365, length.out = n)
season <- factor(ifelse(day_of_year < 100, "Winter",
                        ifelse(day_of_year < 180, "Spring",
                               ifelse(day_of_year < 270, "Summer", "Fall"))))
x_pi <- model.matrix(~ season)

beta_pi <- c(0.2, 1, 4, 6)
logit_pi <- as.vector(x_pi %*% beta_pi)
pi0 <- plogis(logit_pi)
lambda <- c(mu, beta_pi)

# Simulate ZIP response
sim_data <- sim_zip(mu, pi0, tau, arma_order, nsim = n, seed = 1)
y <- sim_data$y

marginal <- zip.marg(link = "identity")
ab <- marginal$bounds(y, x = list(mu = matrix(1, nrow = n), pi0 = x_pi), lambda)

llk_tmet_qmc <- pmvn_tmet(ab[, 1], ab[, 2], tau, od = arma_order)
llk_tmet_mc  <- pmvn_tmet(ab[, 1], ab[, 2], tau, od = arma_order, QMC = FALSE)
llk_ghk_qmc  <- pmvn_ghk(ab[, 1], ab[, 2], tau, od = arma_order)
llk_ghk_mc   <- pmvn_ghk(ab[, 1], ab[, 2], tau, od = arma_order, QMC = FALSE)
llk_vmet     <- loglik_vmet(ab, tau = tau, od = arma_order)
llk_ce       <- pmvn_ce(ab[, 1], ab[, 2], tau, c = 0.5, od = arma_order)


## ---------- ZIB AR(1) with Covariates ----------
n <- 5000
phi <- 0.8
tau <- c(phi)
arma_order <- c(1, 0)
size <- 24
pi0 <- 0.3

# Covariates for logit(mu)
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
prob <- plogis(X %*% beta)


# Simulate ZIB response
sim_data <- sim_zib(prob, pi0*rep(1,n), size, tau, arma_order, nsim = n, seed = 1)
y <- sim_data$y

lambda <- c(beta, pi0)
marginal <- zib.marg(link = "logit", size = size)
ab <- marginal$bounds(y, list(mu = X, pi0 = matrix(1, nrow = n)), lambda)

llk_tmet_qmc <- pmvn_tmet(ab[, 1], ab[, 2], tau, od = arma_order)
llk_tmet_mc  <- pmvn_tmet(ab[, 1], ab[, 2], tau, od = arma_order, QMC = FALSE)
llk_ghk_qmc  <- pmvn_ghk(ab[, 1], ab[, 2], tau, od = arma_order)
llk_ghk_mc   <- pmvn_ghk(ab[, 1], ab[, 2], tau, od = arma_order, QMC = FALSE)
llk_vmet     <- loglik_vmet(ab, tau = tau, od = arma_order)
llk_ce       <- pmvn_ce(ab[, 1], ab[, 2], tau, c = 0.5, od = arma_order)


## ---------- ZIBB AR(1) ----------
n <- 1000
phi <- 0.8
tau <- c(phi)
arma_order <- c(1, 0)
size <- 24
mu <- 0.5
rho <- 0.15
pi0 <- 0.7

sim_data <- sim_zibb(0.5, rho, pi0, size, tau, arma_order, nsim = n, seed = 1)
y <- sim_data$y

lambda <- c(mu, qlogis(rho), qlogis(pi0))
marginal <- zibb.marg(link = "identity", size = size)
ab <- marginal$bounds(y, x = list(mu = matrix(1, nrow = n), pi0 = matrix(1, nrow = n)), lambda)

llk_tmet_qmc <- pmvn_tmet(ab[, 1], ab[, 2], tau, od = arma_order)
llk_tmet_mc  <- pmvn_tmet(ab[, 1], ab[, 2], tau, od = arma_order, QMC = FALSE)
llk_ghk_qmc  <- pmvn_ghk(ab[, 1], ab[, 2], tau, od = arma_order)
llk_ghk_mc   <- pmvn_ghk(ab[, 1], ab[, 2], tau, od = arma_order, QMC = FALSE)
llk_vmet     <- loglik_vmet(ab, tau = tau, od = arma_order)
llk_ce       <- pmvn_ce(ab[, 1], ab[, 2], tau, c = 0.5, od = arma_order)


## ---------- ZIBB AR(1) with Covariates in μ(t) ----------
n <- 1000
phi <- 0.5
tau <- c(phi)
arma_order <- c(1, 0)
size <- 24
rho <- 0.15
pi0 <- 0.3
beta <- c(0.2, 0.3, 0.5)

X <- as.matrix(data.frame(
  x1 = 1,
  x2 = sin(2 * pi * (1:n) / 12),
  x3 = cos(2 * pi * (1:n) / 12)
))
prob <- plogis(X %*% beta)

sim_data <- sim_zibb(prob, rho, pi0, size, tau, arma_order, nsim = n, seed = 1)
y <- sim_data$y

lambda <- c(beta, qlogis(rho), qlogis(pi0))
marginal <- zibb.marg(link = "logit", size = size)
ab <- marginal$bounds(y, x = list(prob = X, pi0 = rep(1, nrow = n)), lambda)

llk_tmet_qmc <- pmvn_tmet(ab[, 1], ab[, 2], tau, od = arma_order)
llk_tmet_mc  <- pmvn_tmet(ab[, 1], ab[, 2], tau, od = arma_order, QMC = FALSE)
llk_ghk_qmc  <- pmvn_ghk(ab[, 1], ab[, 2], tau, od = arma_order)
llk_ghk_mc   <- pmvn_ghk(ab[, 1], ab[, 2], tau, od = arma_order, QMC = FALSE)
llk_vmet     <- loglik_vmet(ab, tau = tau, od = arma_order)
llk_ce       <- pmvn_c


## ---------- ZIBB AR(1) with Seasonal π₀(t) ----------
n <- 1000
phi <- 0.85
tau <- c(phi)
arma_order <- c(1, 0)
size <- 24
prob <- 0.5
rho <- 0.16

# Seasonal design matrix for π₀
X_pi <- cbind(
  1,
  sin(2 * pi * (1:n) / 12),
  cos(2 * pi * (1:n) / 12)
)
colnames(X_pi) <- c("x1", "x2", "x3")

beta_pi <- c(1, 2, 5)
pi0 <- plogis(X_pi %*% beta_pi)

# Simulate ZIBB response with constant mu and varying pi0
sim_data <- sim_zibb(prob, rho, pi0, size, tau, arma_order, nsim = n, seed = 1)
y <- sim_data$y

lambda <- c(mu, qlogis(rho), beta_pi)
marginal <- zibb.marg(link = "logit", size = size)
ab <- marginal$bounds(y, x = list(mu = matrix(1, nrow = n), pi0 = X_pi), lambda)

llk_tmet_qmc <- pmvn_tmet(ab[, 1], ab[, 2], tau, od = arma_order)
llk_tmet_mc  <- pmvn_tmet(ab[, 1], ab[, 2], tau, od = arma_order, QMC = FALSE)
llk_ghk_qmc  <- pmvn_ghk(ab[, 1], ab[, 2], tau, od = arma_order)
llk_ghk_mc   <- pmvn_ghk(ab[, 1], ab[, 2], tau, od = arma_order, QMC = FALSE)
llk_vmet     <- loglik_vmet(ab, tau = tau, od = arma_order)
llk_ce       <- pmvn_ce(ab[, 1], ab[, 2], tau, c = 0.5, od = arma_order)


## ---------- ZIBB AR(1) with Covariates in μ(t) and π₀(t) ----------
n <- 5000
phi <- 0.9
tau <- c(phi)
arma_order <- c(1, 0)
size <- 24
rho <- 0.333

# Covariates for μ(t)
X_beta <- cbind(
  1,
  sin(2 * pi * (1:n) / 365),
  cos(2 * pi * (1:n) / 365)
)
colnames(X_beta) <- c("x1", "x2", "x3")
beta_mu <- c(1.2, 2, 3)
prob <- plogis(X_beta %*% beta_mu)

# Seasonal factor covariates for π₀(t)
day_of_year <- rep(1:365, length.out = n)
season <- factor(ifelse(day_of_year < 100, "Winter",
                        ifelse(day_of_year < 180, "Spring",
                               ifelse(day_of_year < 270, "Summer", "Fall"))))
X_pi <- model.matrix(~ season)
alpha_pi <- c(0.2, 1, 4, 6)
pi0 <- plogis(X_pi %*% alpha_pi)

# Simulate data
sim_data <- sim_zibb(prob, rho, pi0, size, tau, arma_order, nsim = n, seed = 1)
y <- sim_data$y

lambda <- c(beta_mu, qlogis(rho), alpha_pi)
marginal <- zibb.marg(link = "logit", size = size)
ab <- marginal$bounds(y, x = list(mu = X_beta, pi0 = X_pi), lambda)

llk_tmet_qmc <- pmvn_tmet(ab[, 1], ab[, 2], tau, od = arma_order)
llk_tmet_mc  <- pmvn_tmet(ab[, 1], ab[, 2], tau, od = arma_order, QMC = FALSE)
llk_ghk_qmc  <- pmvn_ghk(ab[, 1], ab[, 2], tau, od = arma_order)
llk_ghk_mc   <- pmvn_ghk(ab[, 1], ab[, 2], tau, od = arma_order, QMC = FALSE)
llk_vmet     <- loglik_vmet(ab, tau = tau, od = arma_order)
llk_ce       <- pmvn_ce(ab[, 1], ab[, 2], tau, c = 0.5, od = arma_order)

