# Load required libraries
library(bayesmeta)
library(ggplot2)
library(VGAM)

# Define PDF and CDF of the mixture distribution
pdf1 <- function(x) {
  0.2 * dweibull(x, shape = 3, scale = 1) + 0.8 * dbeta(x, shape1 = 0.5, shape2 = 0.5)
}

cdf1 <- function(x) {
  0.2 * pweibull(x, shape = 3, scale = 1) + 0.8 * pbeta(x, shape1 = 0.5, shape2 = 0.5)
}

# Define double gamma PDF and CDF with shape parameter = 0.2
# Double gamma distribution: symmetric gamma distribution with scale parameter b
# PDF: f(x) = (1/(2*b*gamma(shape))) * |x/b|^(shape-1) * exp(-|x/b|)
# For shape = 0.2, this gives heavy-tailed behavior

ddoublegamma <- function(x, shape = 0.2, scale = 1) {
  (1/(2 * scale * gamma(shape))) * (abs(x/scale))^(shape-1) * exp(-abs(x/scale))
}

pdoublegamma <- function(x, shape = 0.2, scale = 1) {
  ifelse(x < 0, 
         (1/2) * pgamma(-x, shape = shape, scale = scale, lower.tail = FALSE),
         1 - (1/2) * pgamma(x, shape = shape, scale = scale, lower.tail = FALSE))
}

# Define noise-to-signal ratios (NSRs)
nsr_values <- c(0.1, 0.2, 0.5)

# the SD of the mixture
sigma_x <- sqrt(0.1458)

# For double gamma distribution with shape = 0.2
# Variance = 2 * shape * scale^2 (for shape > 0)
# For shape = 0.2: variance = 2 * 0.2 * scale^2 = 0.4 * scale^2
# Therefore: scale = sqrt(var_noise / 0.4)

shape_parameter <- 2

b_values_doubleneg <- sapply(nsr_values, function(nsr) {
  var_noise <- (nsr * sigma_x)^2  # Noise variance based on NSR
  sqrt(var_noise / (2 * shape_parameter))
})

# Define double gamma PDF and CDF with shape = 0.2
pdf2 <- function(x, b) ddoublegamma(x, shape = shape_parameter, scale = b)
cdf2 <- function(x, b) pdoublegamma(x, shape = shape_parameter, scale = b)

# Convolution function
convolve_cdf <- Vectorize(function(t, b) {
  integrate(function(x) pdf2(x, b) * cdf1(t - x), -Inf, Inf, subdivisions = 1000)$value
})

# Generate evaluation points
y <- seq(-1, 2, length.out = 1000)
d <- cdf1(y)

# Compute convoluted CDFs for different NSRs using double gamma error (shape=0.2)
d_list <- lapply(b_values_doubleneg, function(b) convolve_cdf(y, b))
names(d_list) <- c("10%", "20%", "50%")

# Compute hazard functions
H <- -log(1 - d)
CH_values <- lapply(d_list, function(d) -log(1 - d))

# Remove NA values for consistency
valid_indices <- complete.cases(H, CH_values[[1]], CH_values[[2]], CH_values[[3]])
y <- y[valid_indices]
H <- H[valid_indices]
CH_values <- lapply(CH_values, function(CH) CH[valid_indices])

# Create dataframe for plotting
Data2 <- data.frame(
  x = rep(y, times = 4),
  y = c(H, CH_values[[1]], CH_values[[2]], CH_values[[3]]),
  group = rep(c("0%", "10%", "20%", "50%"), each = length(y))
)

# Plot hazard functions
ggplot(Data2, aes(x = x, y = y, linetype = group)) +
  geom_line() +
  labs(y = "Hazard function", x = "x", linetype = "NSR",
       title = "Double Gamma Error (Shape = 0.2)") +
  scale_linetype_manual(values = c("0%" = "solid", "10%" = "dashed", "20%" = "twodash", "50%" = "dotted")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))




# ------------------------------------------
# 1. Parameters
# ------------------------------------------
set.seed(123)  # For reproducibility; optional

n   <- 100
r   <- 4
k   <- 2
p   <- 1
sim <- 5000
a   <- 1.5
alpha <- 0.1

# Theoretical reference values
TR <- sapply(seq_len(r), function(i) {
  1 - exp(-i / k)
})

# ------------------------------------------
# 2. Define Krec
# ------------------------------------------
# Mathematica:
# 1 - GammaRegularized[i, -k*Log[1-u]]
#
# R's pgamma() gives the regularized lower
# incomplete gamma function.

Krec <- function(i, u) {
  pgamma(-k * log1p(-u), shape = i, rate = 1)
}

# ------------------------------------------
# 3. Calculate the coefficients coe[i, j]
# ------------------------------------------
u_grid <- (0:n) / n

# Rows correspond to i = 1,...,r
# Columns correspond to j = 1,...,n

coe <- t(vapply(
  seq_len(r),
  function(i) diff(Krec(i, u_grid)),
  numeric(n)
))

# ------------------------------------------
# 4. Define the test statistic
# ------------------------------------------
test_statistic <- function(x) {
  
  # Sort the sample
  x <- sort(x)
  
  # Estimated expected generalized order statistics
  EXu <- drop(coe %*% x)
  
  # Empirical CDF evaluated at each estimate
  EMP <- vapply(
    EXu,
    function(value) mean(x <= value),
    numeric(1)
  )
  
  # Positive-part differences
  z <- pmax(0, TR - EMP)
  
  # L_p norm; here p = 1
  if (p == 1) {
    return(sum(abs(z)))
  } else if (is.infinite(p)) {
    return(max(abs(z)))
  } else {
    return(sum(abs(z)^p)^(1 / p))
  }
}

# ------------------------------------------
# 5. Simulate under the null hypothesis
#    X ~ Exponential(rate = 1)
# ------------------------------------------
SI <- replicate(sim, {
  x <- rexp(n, rate = 1)
  test_statistic(x)
})

# ------------------------------------------
# 6. Simulate under the Weibull alternative
#    X ~ Weibull(shape = a, scale = 1)
# ------------------------------------------
s <- replicate(sim, {
  x <- rweibull(n, shape = a, scale = 1)
  test_statistic(x)
})

# ------------------------------------------
# 7. Calculate empirical upper-tail p-values
# ------------------------------------------
null_cdf <- ecdf(SI)

DD <- 1 - null_cdf(s)

# ------------------------------------------
# 8. Estimated rejection rate
# ------------------------------------------
rejection_rate <- mean(DD <= alpha)

# Display results
cat("Theoretical reference values:\n")
print(TR)

cat("\nEstimated rejection rate:\n")
print(rejection_rate)





set.seed(1)

n   <- 100
r   <- 4
k   <- 2
p   <- 1
sim <- 5000
a   <- 1.5
TR  <- 1 - exp(-(1:r) / k)

# Krec[i, u] = 1 - GammaRegularized[i, -k*log(1-u)]  (= pgamma lower tail; 1 at u = 1)
Krec <- function(i, u) pgamma(-k * log1p(-u), shape = i)   # u = 1 gives Inf -> 1

# n x r coefficient matrix: coe[jj, i] = Krec(i, jj/n) - Krec(i, (jj-1)/n)
u   <- (0:n) / n
coe <- sapply(1:r, function(i) diff(Krec(i, u)))

lp_norm <- function(x, p) sum(abs(x)^p)^(1 / p)

# Statistic for one sample: || max(0, TR - ECDF(EXu)) ||_p
stat_one <- function(x) {
  x   <- sort(x)
  EXu <- drop(x %*% coe)            # length r
  EMP <- ecdf(x)(EXu)               # step ECDF, like CDF[EmpiricalDistribution[...]]
  lp_norm(pmax(0, TR - EMP), p)
}

# Null: Exp(1)
SI <- replicate(sim, stat_one(rexp(n, 1)))

# Alternative: Weibull(shape = a, scale = 1)
s <- replicate(sim, stat_one(rweibull(n, shape = a, scale = 1)))

# p-values: DD = 1 - ECDF_SI(s) = P(SI > s)
DD <- vapply(s, function(si) mean(SI > si), numeric(1))

# Estimated power at level 0.1
mean(DD <= 0.1)
