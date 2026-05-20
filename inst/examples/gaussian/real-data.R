## -------------------------------------------
## Example: Real Data – ZIBB Gaussian Copula Model
## Kickapoo Downtown Airport “Hot Hours”
## -------------------------------------------

library(gctsc)

## --- Load data ---
data("KCWC", package = "copTSC")
KCWC$date <- as.Date(KCWC$date)
y <- KCWC$hot

n <- length(y)
time <- 1:n

## Seasonal covariates for π0(t)
X_pi <- cbind(
  Intercept = 1,
  x_sin = sin(2 * pi * time / 365),
  x_cos = cos(2 * pi * time / 365)
)


data <- data.frame(cbind(y, X_pi))


## --- Train/Test split ---
n_train <- 1000
train_data <- data[1:n_train,]

## ===========================================
##   Fit ZIBB Marginal + AR(1) Copula (TMET)
## ===========================================
fit_tmet <- gctsc(
  formula  = list(mu = y ~ 1, pi0 = ~ x_sin + x_cos),
  data     = train_data,
  marginal = zibb.marg(link = "logit", size = 24),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "TMET", family ="gaussian",
  options  = gctsc.opts(seed = 1, M = c(100,1000))
)

summary(fit_tmet)
plot(fit_tmet)         # PIT, residuals, ACF plot

## --- One-step-ahead prediction ---
t_pred <- n_train + 1
pred_tmet <- predict(
  fit_tmet, 
  newdata = data[1001,]
)

pred_tmet

## ===========================================
##   Fit ZIBB Marginal (GHK for comparison)
## ===========================================
fit_ghk <- gctsc(
  formula  = list(mu = y~ 1, pi0 = ~ x_sin + x_cos),
  data     = train_data,
  marginal = zibb.marg(link = "logit", size = 24),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "GHK", family ="gaussian",
  options  = gctsc.opts(seed = 1, M = 1000)
)

pred_ghk <- predict(
  fit_ghk, 
  newdata = data[1001,]
)

pred_ghk

## ===========================================
##   Optional: Fit Zero-Inflated Binomial (ZIB)
## ===========================================
fit_zib <- gctsc(
  formula  = list(mu = y ~ 1, pi0 = ~ x_sin + x_cos),
  data     = train_data,
  marginal = zib.marg(link = "logit", size = 24),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "TMET",
  options  = gctsc.opts(seed = 1, M = 1000)
)

summary(fit_zib)
