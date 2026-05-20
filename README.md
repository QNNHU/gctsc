---
output: github_document
---

# gctsc

gctsc provides fast and scalable likelihood inference for Gaussian and Student–t copula models for count time series.

The package supports a wide range of discrete marginals, including:

- Poisson
- Negative Binomial
- Binomial
- Beta-Binomial
- Zero-Inflated Poisson (ZIP)
- Zero-Inflated Binomial (ZIB)
- Zero-Inflated Beta-Binomial (ZIBB)

The latent dependence structure is modeled through ARMA(p, q) processes.

Likelihood evaluation is available through the following approximation methods:

- TMET — Time Series Minimax Exponential Tilting
- GHK — Geweke–Hajivassiliou–Keane simulator
- CE — Continuous Extension

The implementation exploits ARMA structure for efficient high-dimensional computation.

Additional features include:

- Simulation utilities
- Randomized quantile residual diagnostics
- Predictive distributions and scoring rules (CRPS, LOGS)

## Installation

From CRAN:

```r
install.packages("gctsc")
```

From Github:
remotes::install_github("QNNHU/gctsc")

## Quick Example:
```r
library(gctsc)

# Simulate Poisson AR(1) data under a Gaussian copula
set.seed(1)
y <- sim_poisson(
  mu = 5,
  tau = 0.5,
  arma_order = c(1, 0),
  nsim = 300,
  family = "gaussian"
)$y

# Fit model
fit <- gctsc(
  y ~ 1,
  data = data.frame(y = y),
  marginal = poisson.marg(link = "log"),
  cormat = arma.cormat(p = 1, q = 0),
  method = "TMET",
  family = "gaussian",
  options = gctsc.opts(seed =1, M = c(100,1000))
)

summary(fit)

# Diagnostic plots
plot(fit)

# One-step prediction
predict(fit)
```
# What Makes gctsc Different?

Compared to existing implementations, `gctsc` added:

- Exploits ARMA structure for scalable likelihood evaluation in time series settings

- Supports zero-inflated marginals with flexible covariate specification, including seasonal components

- Implements scalable minimax exponential tilting (TMET) for efficient likelihood approximation

- Provides a linear-cost GHK importance sampling implementation

- Implements fast continuous extension method

- Supports Student–t copulas for modeling heavy-tailed dependence

- Computes full predictive distributions for discrete time series

# References

If you use this package in published work, please cite:


Nguyen, Q. N. and De Oliveira, V. (2026). Approximating Gaussian Copula Models for Count Time Series: Connecting the Distributional Transform and a Continuous Extension. Journal of Applied Statistics, 53, 1--22.

Nguyen, Q. N. and De Oliveira, V. (2026). Likelihood Inference in Gaussian Copula Models for Count Time Series via Minimax Exponential Tilting. Computational Statistics and Data Analysis, 218, 108344.

Nguyen, Q. N. and De Oliveira, V. (2026). Scalable Likelihood Inference for Student--t Copula Count Time Series. Stats, 9, 1--49.




