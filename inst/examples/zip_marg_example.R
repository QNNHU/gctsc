# Simulate zero-inflated Poisson data with seasonal covariates
library(gctsc)
set.seed(1)
n <- 500
mu <- 10
phi <- 0.5
tau <- c(phi)
arma_order <- c(1, 0)

# Seasonal factor for zero-inflation component
day_of_year <- rep(1:365, length.out = n)
season <- factor(ifelse(day_of_year < 100, "Winter",
                        ifelse(day_of_year < 180, "Spring",
                               ifelse(day_of_year < 270, "Summer", "Fall"))))

# Design matrix for pi0
X_pi <- model.matrix(~ season)
colnames(X_pi) <- c("Intercept", "Spring", "Summer", "Winter")
beta_pi <- c(0.2, 1, 0.3, 0.5)
logit_pi <- as.vector(X_pi %*% beta_pi)
pi0 <- plogis(logit_pi)

# Simulate data
sim_data <- sim_zip(mu = mu, pi0 = pi0, tau = tau, arma_order = arma_order, nsim = n, seed = 1)
y <- sim_data$y
X_mu <- matrix(1, n, 1)

# Prepare training data
train_idx <- 1:499
data_train <- data.frame(cbind(y, X_mu,X_pi ))

# Fit the Gaussian copula time series model
fit <- gctsc(
  formula = list(mu = y ~ 1, pi0 = ~ Spring + Summer + Winter),
  data = data_train,
  marginal = zip.marg(link = "identity"),
  cormat = arma.cormat(p = 1, q = 0),
  method = "GHK",
  options = gctsc.opts(seed = 1, M = 1000)
)

# Predict for a new observation (e.g., time point 991)
X_test <- list(mu = matrix(1, 1, 1), pi0 = X_pi[500, , drop = FALSE])
predict(fit, X_test = X_test)
