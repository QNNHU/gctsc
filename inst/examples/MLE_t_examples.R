
# Example: Gaussian Copula Time Series Estimation for Various Marginals

## ---------- Poisson AR(1) ----------
n <- 5000
mu <- 10
phi <- -0.8
theta <- 0
arma_order <- c(1, 0)
tau <- c(phi)
seed= 5
df<- 3
sim_data <- sim_poisson(mu, tau, arma_order, nsim = n, seed =seed, family ="t", df= df)
y <- sim_data$y
x <- matrix(1, nrow = n)
marginal <- poisson.marg()
ab <- marginal$bounds(y, x, mu, family ="t", df=df)
system.time({
fit_tmet <- gctsc(y~1,
                  marginal = poisson.marg(lambda.lower = c(0)),
                  cormat = arma.cormat(p = 1, q = 0),
                  method = "TMET", family="t", df=df, 
                  options = gctsc.opts(seed= 1,M=1000),
                  QMC = TRUE)})

plot(fit_tmet)

## ---------- Negative Binomial ARMA(1,1) ----------
mu_values <- c(10)   # Mean of NegBin
phi_values <-0.2  # AR parameter
theta_values <- c(0.2) # MA parameter
df <- 10      # df for t copula
dispersion<- c(2) # dispersion parameter for NegBin
tau <- c(phi_values, theta_values)
arma_order <- c(1,1)
seed=10
n=1000
sim_data <- sim_negbin(mu_values, dispersion, tau, arma_order, nsim = n, seed = seed,family ="t", df= df)
y <- sim_data$y
x <- matrix(1, nrow = n)
set.seed(seed)

fit_ce2 <- gctsc(y ~ 1,
            marginal = negbin.marg(lambda.lower = c(0,0)),
            cormat = arma.cormat(p = 1, q = 1),
            method = "TMET", family ="t", df = df,
            QMC = TRUE,
            options = gctsc.opts(seed = 1, M = 1000))
plot(fit_ce2)


profile_llk <- function(log_nu) {
  nu <- 2 + exp(log_nu)   # enforce nu > 2
  
  fit <- gctsc(y ~ 1,
               marginal = negbin.marg(lambda.lower = c(0,0)),
               cormat = arma.cormat(p = 1, q = 1),
               method = "CE", family ="t", df = df,
               QMC = TRUE,
               options = gctsc.opts(seed = 1, M = 1000))
  
  return(fit$maximum)  # for minimizer
}

opt <- optim(par = log(60 - 2),
             fn = profile_llk,
             method = "Brent",
             lower = log(1),
             upper = log(300))

nu_hat <- 2 + exp(opt$par)


## ---------- Poisson AR(1) with Covariates ----------
n <- 1000
phi <- 0.5
tau <- c(phi)
arma_order <- c(1, 0)
beta <- c(1, 0.3, 1, 0.5)

xi <- numeric(n)
zeta <- rnorm(n)
for (j in 3:n) xi[j] <- 0.6 * xi[j - 1] - 0.4 * xi[j - 2] + zeta[j]

X <- as.matrix(data.frame(
  x1 = rep(1, n),
  x2 = sin(2 * pi * (1:n) / 12),
  x3 = cos(2 * pi * (1:n) / 12),
  x4 = xi
))

df= 3
mu <- exp(X %*% beta)
sim_data <- sim_poisson(mu, tau, arma_order, nsim = n, seed = 1, family = "t",
                        df=df )
y <- sim_data$y

# fit_fit <- gctsc.fit(X, y,
#                      marginal = poisson.marg(link = "log"),
#                      cormat = arma.cormat(p = 1, q = 0),
#                      method = "CE",c=0.6,
#                      options = gctsc.opts(seed = 1))

data_df <- data.frame(Y = y, X)

fit_fml <- gctsc(Y ~ x2 + x3 + x4, data = data_df,
                 marginal = poisson.marg(link = "log"),
                 cormat = arma.cormat(p = 1, q = 0), family ="t", df=df,
                 method = "TMET",
                 options = gctsc.opts(seed = 1, M=1000))
predict(fit_fml, X_test = ( data_df[401,3:5]))




## ---------- NB AR(1) with Covariates ----------
n <- 500
phi <- 0.5
theta <- 0.3
tau <- c(phi)
arma_order <- c(1, 0)
beta <- c(1, 0.3, 1, 0.5)
dispersion <- 2
xi <- numeric(n)
zeta <- rnorm(n)
dispersion =2
for (j in 3:n) xi[j] <- 0.6 * xi[j - 1] - 0.4 * xi[j - 2] + zeta[j]

X <- as.matrix(data.frame(
  x1 = rep(1, n),
  x2 = sin(2 * pi * (1:n) / 12),
  x3 = cos(2 * pi * (1:n) / 12),
  x4 = xi
))

df= 3
mu <- exp(X %*% beta)
sim_data <- sim_negbin(mu, dispersion = dispersion, tau, arma_order, nsim = n, seed = 1, family = "t",
                        df=df )
y <- sim_data$y

# fit_fit <- gctsc.fit(X, y,
#                      marginal = poisson.marg(link = "log"),
#                      cormat = arma.cormat(p = 1, q = 0),
#                      method = "CE",c=0.6,
#                      options = gctsc.opts(seed = 1))

data_df <- data.frame(Y = y, X)

fit_fml <- gctsc(Y ~ x2 + x3 + x4, data = data_df,
                 marginal = negbin.marg(link = "log"),
                 cormat = arma.cormat(p = 1, q = 0), family ="gaussian",
                 method = "GHK",
                 options = gctsc.opts(seed = 1, M=1000))
predict(fit_fml, X_test = ( data_df[401,3:5]))


## ---------- ZIP AR(1) ----------
n <- 500
mu <- 10
pi0 <- 0.4
logit_pi0 <- qlogis(0.4)
phi <- 0.5
tau <- c(phi)
arma_order <- c(1, 0)

# Simulate data from ZIP model with AR(1) latent correlation
sim_data <- sim_zip(mu = mu, pi0 = pi0, tau = tau, arma_order = arma_order, nsim = n, seed = 1)
y <- sim_data$y

# Design matrices for both mean and zero-inflation components
x_mu <- matrix(1, nrow = n)
x_pi <- matrix(1, nrow = n)

# # Fit CE method using design matrix input
# fit_ghk <- gctsc.fit(
#   x = list(mu = x_mu, pi0 = x_pi), y = y,
#   marginal = zip.marg(link = "identity", lambda.lower = c(0, -Inf)),
#   cormat = arma.cormat(p = 1, q = 0),
#   method = "GHK", c = 0.6,
#   options = gctsc.opts(seed = 1, M = 1000)
# )

# Fit GHK method using formula input
fit_tmet <- gctsc(
  formula = list(mu = y ~ 1, pi0 = ~1),
  marginal = zip.marg(link = "identity", lambda.lower = c(0,-Inf)),
  cormat = arma.cormat(p = 1, q = 0),
  method = "TMET",
  options = gctsc.opts(seed = 1, M = 1000)
)

predict(fit_tmet)

## ---------- ZIP AR(1) with Covariates ----------
n <- 1000
mu <- 5
phi <- 0.5
tau <- c(phi)
arma_order <- c(1, 0)

# Create seasonal covariate for pi0 using factor
day_of_year <- rep(1:365, length.out = n)
season <- factor(ifelse(day_of_year < 100, "Winter",
                        ifelse(day_of_year < 180, "Spring",
                               ifelse(day_of_year < 270, "Summer", "Fall"))))
x_pi <- model.matrix(~ season)

# True pi0 value
beta_pi <- c(0.2, 1, 0.3, 0.5)
logit_pi <- as.vector(x_pi %*% beta_pi)
pi0 <- plogis(logit_pi)

# Simulate data
sim_data <- sim_zip(mu = mu, pi0 = pi0, tau = tau, arma_order = arma_order, nsim = n, seed = 1)
y <- sim_data$y

# Fit CE method using formula interface
data_df <- data.frame(y = y, x_pi)
colnames(data_df) <- make.names(colnames(data_df))  # ensure syntactically valid names

fit_tmet <- gctsc(
  formula = list(mu = y ~ 1, pi0 = ~ seasonSpring + seasonSummer + seasonWinter),
  data = data_df,
  marginal = zip.marg(link = "identity"),
  cormat = arma.cormat(p = 1, q = 0),
  method = "TMET", 
  options = gctsc.opts(seed = 1, M = 1000)
)

# # Fit GHK method using design matrices
# x_mu <- matrix(1, nrow = n)
# fit_ghk <- gctsc.fit(
#   x = list(mu = x_mu, pi0 = x_pi), y = y,
#   marginal = zip.marg(link = "identity"),
#   cormat = arma.cormat(p = 1, q = 0),
#   method = "GHK", c = 0.5,
#   options = gctsc.opts(seed = 1, M = 1000)
# )

predict(fit_tmet, X_test = data_df[499,])
## ---------- Binomial AR(1) ----------
n <- 500
size <- 24
prob <- rep(0.3, n)
phi <- 0.8
tau <- c(phi)
arma_order <- c(1, 0)

# Simulate Binomial count time series
sim_data <- sim_binom(prob = 0.3, size = size, tau = tau, arma_order = arma_order, nsim = n, seed = 1)
y <- sim_data$y
x <- matrix(1, nrow = n)

# # Fit CE method using x-input
fit_ce <- gctsc.fit(
  x = x, y = y,
  marginal = binom.marg(link = "identity", size = size),
  cormat = arma.cormat(p = 1, q = 0),
  method = "CE", c = 0.5,
  options = gctsc.opts(seed = 1, M = 1000)
)

# Fit CE method using formula interface
fit_fml <- gctsc(
  formula = y ~ 1,
  marginal = binom.marg(link = "logit", size = size),
  cormat = arma.cormat(p = 1, q = 0),
  method = "GHK",
  options = gctsc.opts(seed = 1, M = 1000)
)

predict(fit_fml)

##---------zero inflated binomial---------
n <- 500
size <- 24
prob <- 0.2
rho <-0.18    
phi <- 0.8
tau <- c(phi)
arma_order <- c(1, 0)
pi0 <- 0.2
# Simulate data from Beta-Binomial with AR(1) structure
sim_data <- sim_zib(prob = prob, pi0 = pi0, size = size, tau = tau,
                    arma_order = arma_order, nsim = n, seed = 1)
y <- sim_data$y
x <- matrix(1, nrow = n)

# Fit TMET method using formula interface
fit_tmet <- gctsc(
  formula = list(mu = y ~ 1, pi0 = ~ 1),
  marginal =  zib.marg(link = "logit", size = size),
  cormat = arma.cormat(p = 1, q = 0),
  method = "TMET",
  options = gctsc.opts(seed = 1, M = 1000)
)

predict(fit_tmet)

##---------zero inflated binomial with covariate---------
n <- 500
size <- 24
prob <- 0.2
rho <-0.18    
phi <- 0.8
tau <- c(phi)
arma_order <- c(1, 0)
# Create seasonal covariate for pi0 using factor
day_of_year <- rep(1:365, length.out = n)
season <- factor(ifelse(day_of_year < 100, "Winter",
                        ifelse(day_of_year < 180, "Spring",
                               ifelse(day_of_year < 270, "Summer", "Fall"))))
x_pi <- model.matrix(~ season)

# True pi0 value
beta_pi <- c(0.2, 1, 0.3, 0.5)
logit_pi <- as.vector(x_pi %*% beta_pi)
pi0 <- as.vector(plogis(logit_pi))

# Simulate data from Beta-Binomial with AR(1) structure
sim_data <- sim_zib(prob = prob, pi0 = pi0, size = size, tau = tau,
                    arma_order = arma_order, nsim = n, seed = 1)
y <- sim_data$y

df <- data.frame(y,x_pi)
# Fit TMET method using formula interface
fit <- gctsc(
  formula = list(mu = y ~ 1, pi0 = ~ seasonSpring + seasonSummer + seasonWinter), data=df,
  marginal =  zib.marg(link = "logit", size = size),
  cormat = arma.cormat(p = 1, q = 0),
  method = "TMET",
  options = gctsc.opts(seed = 1, M = 1000)
)
summary(fit)

fit <- gctsc.fit(
  x = list(mu = as.matrix(rep(1,n)), pi0 = x_pi), y=y,
  marginal =  zib.marg(link = "logit", size = size),
  cormat = arma.cormat(p = 1, q = 0),
  method = "GHK",
  options = gctsc.opts(seed = 1, M = 1000)
)

predict(fit, X_test = df)


## ---------- Beta-Binomial AR(1) ----------
n <- 500
size <- 24
beta_0 <- 0.2
prob <- plogis(beta_0)
rho <-0.18    
phi <- 0.8
tau <- c(phi)
arma_order <- c(1, 0)

# Simulate data from Beta-Binomial with AR(1) structure
sim_data <- sim_bbinom(prob = prob, rho = rho, size = size, tau = tau,
                       arma_order = arma_order, nsim = n, seed = 1)
y <- sim_data$y
x <- matrix(1, nrow = n)

# Fit GHK method with x-input
fit_ghk <- gctsc.fit(
  x = x, y = y,
  marginal = bbinom.marg( size = size),
  cormat = arma.cormat(p = 1, q = 0),
  method = "GHK",
  options = gctsc.opts(seed = 1, M = 1000)
)

# Fit TMET method using formula interface
fit_tmet <- gctsc(
  formula = y ~ 1,
  marginal = bbinom.marg(link = "logit", size = size),
  cormat = arma.cormat(p = 1, q = 0),
  method = "TMET",
  options = gctsc.opts(seed = 1, M = 1000)
)
predict(fit_tmet)
## ---------- Beta-Binomial AR(1) with Covariates ----------
n <- 500
size <- 24
phi <- 0.5
tau <- c(phi)
arma_order <- c(1, 0)
rho <- 5
rho_vgam <- 1 / (1 + rho)

# Simulate covariate: seasonal + AR(2)-type xi
zeta <- rnorm(n)
xi <- numeric(n)
for (j in 3:n) {
  xi[j] <- 0.6 * xi[j - 1] - 0.4 * xi[j - 2] + zeta[j]
}

# Design matrix for logit(prob)
X <- data.frame(
  x1 = 1,
  x2 = sin(2 * pi * (1:n) / 12),
  x3 = cos(2 * pi * (1:n) / 12),
  x4 = xi
)

# Simulate true logit probabilities
beta <- c(0.2, 0.3, 0.5, 0.3)
logit_prob <- as.vector(as.matrix(X) %*% beta)
prob <- plogis(logit_prob)

# Simulate Beta-Binomial count data
sim_data <- sim_bbinom(prob = prob, rho = rho_vgam, size = size,
                       tau = tau, arma_order = arma_order, nsim = n, seed = 1)
y <- sim_data$y

# Fit using x-based interface
fit_ghk <- gctsc.fit(
  x = as.matrix(X), y = y,
  marginal = bbinom.marg(link = "logit", size = size),
  cormat = arma.cormat(p = 1, q = 0),
  method = "GHK",
  options = gctsc.opts(seed = 1, M = 1000)
)

# Fit using formula interface
data_df <- data.frame(y = y, X)
fit_ce <- gctsc(
  formula = y ~ x2 + x3 + x4, data = data_df,
  marginal = bbinom.marg(link = "logit", size = size),
  cormat = arma.cormat(p = 1, q = 0),
  method = "CE",
  options = gctsc.opts(seed = 1, M = 1000)
)

predict(fit_ce, X_test = data_df[500,])
## ---------- ZIBB AR(1) ----------
n <- 500
size <- 24
phi <- 0.5
tau <- c(phi)
arma_order <- c(1, 0)
rho <- 0.18
logit_rho <-qlogis(rho)
pi0 <- 0.7
beta_0 <- 0.2
prob <- plogis(0.2)

# Simulate ZIBB data with constant mean and zero-inflation
sim_data <- sim_zibb(prob, rho = rho_vgam, pi0 = pi0,
                     size = size, tau = tau, arma_order = arma_order, nsim = n, seed = 1)
y <- sim_data$y

# Design matrices
x_mu <- matrix(1, nrow = n)
x_pi <- matrix(1, nrow = n)

#Fit GHK method (x-input)
fit_ghk <- gctsc.fit(
  x = list(mu = x_mu, pi0 = x_pi), y = y,
  marginal = zibb.marg(link = "logit", size = size),
  cormat = arma.cormat(p = 1, q = 0),
  method = "GHK",
  options = gctsc.opts(seed = 1, M = 1000)
)

# Fit TMET method (formula interface)
fit_tmet <- gctsc(
  formula = list(mu = y ~ 1, pi0 = ~ 1),
  data = data.frame(y),
  marginal = zibb.marg(link = "logit", size = size),
  cormat = arma.cormat(p = 1, q = 0),
  method = "CE",
  options = gctsc.opts(seed = 1, M = 1000)
)

predict(fit_tmet)

## ---------- ZIBB AR(1) with Seasonal π₀(t) ----------
n <- 1500
size <- 24
phi <- 0.85
tau <- c(phi)
arma_order <- c(1, 0)
rho <- 0.16
logit_rho <-qlogis(rho)
beta_0 <- 0.2
prob <- plogis(0.2)
# Covariate matrix for π₀(t): seasonal (sin, cos)
X_pi <- cbind(
  1,
  sin(2 * pi * (1:n) / 12),
  cos(2 * pi * (1:n) / 12)
)
colnames(X_pi) <- c("x1", "x2", "x3")

# Generate time-varying π₀
beta_pi <- c(1, 2, 5)
pi0 <- as.vector(plogis(X_pi %*% beta_pi))

# Simulate ZIBB data with constant mu and seasonal π₀
sim_data <- sim_zibb(prob, rho = rho, pi0 = pi0,family = "t", df=10,
                     size = size, tau = tau, arma_order = arma_order, nsim = n, seed = 1)
y <- sim_data$y

# Fit TMET method using x-input
X_mu <- matrix(1, nrow = n)
fit_tmet <- gctsc.fit(
  x = list(mu = X_mu, pi0 = X_pi), y = y,
  marginal = zibb.marg(link = "logit", size = size),
  cormat = arma.cormat(p = 1, q = 0), family = "t",df= 10,
  method = "GHK"
)

# # Fit GHK method using formula interface
data_df <- data.frame(y = y, X_pi)
fit_ghk <- gctsc(
  formula = list(mu = y ~ 1, pi0 = ~ x2 + x3),
  data = data_df,
  marginal = zibb.marg(link = "logit", size = size),
  cormat = arma.cormat(p = 1, q = 0),
  method = "GHK"
)

predict(fit_ghk,X_test = data_df[500,])

## ---------- ZIBB AR(1) with Covariates in μ(t) and π₀(t) ----------
n <- 500
size <- 10
phi <- 0.5
tau <- c(phi)
arma_order <- c(1, 0)
rho <- 0.3
logit_rho <-qlogis(rho)
# Covariates for μ(t): seasonal sine/cosine
X_mu <- cbind(
  1,
  sin(2 * pi * (1:n) / 365),
  cos(2 * pi * (1:n) / 365)
)
colnames(X_mu) <- c("x1", "x2", "x3")
beta_mu <- c(1.2, 2, 3)
prob <- as.vector(plogis(X_mu %*% beta_mu))

# Covariates for π₀(t): seasonal categories
day_of_year <- rep(1:365, length.out = n)
season <- factor(ifelse(day_of_year < 100, "Winter",
                        ifelse(day_of_year < 180, "Spring",
                               ifelse(day_of_year < 270, "Summer", "Fall"))))
X_pi <- model.matrix(~ season)
alpha_pi <- c(0.2, 0.3, 0.5, 0.7)
pi0 <- as.vector(plogis(X_pi %*% alpha_pi))

# Simulate ZIBB data
sim_data <- sim_zibb(prob, rho = rho, pi0 = pi0,
                     size = size, tau = tau, arma_order = arma_order, nsim = n, seed = 1)
y <- sim_data$y

# Fit CE method using x-input
fit_ce <- gctsc.fit(
  x = list(mu = X_mu, pi0 = X_pi), y = y,
  marginal = zibb.marg(link = "logit", size = size),
  cormat = arma.cormat(p = 1, q = 0),
  method = "CE", c = 0.6,
  options = gctsc.opts(seed = 1, M = 1000)
)

## Fit TMET using formula interface
data_df <- data.frame(y = y, X_mu, X_pi)
fit_tmet <- gctsc(
  formula = list(mu = y ~ x2 + x3, pi0 = ~ seasonSpring + seasonSummer  + seasonWinter),
  data = data_df,
  marginal = zibb.marg(link = "logit", size = size),
  cormat = arma.cormat(p = 1, q = 0),
  method = "CE",
  options = gctsc.opts(seed = 1, M = 1000)
)

predict(fit_tmet,X_test = data_df[500,])
