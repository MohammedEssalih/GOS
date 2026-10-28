library(pracma)
library(flexsurv)
library(purrr)
library(tidyverse)



lp_norm <- function(x, p = 2) {
  if (length(p) != 1L || is.na(p) || p < 1) {
    stop("p must be a single number greater than or equal to 1.")
  }
  
  if (p == Inf) {
    return(max(abs(x), na.rm = FALSE))
  }
  
  sum(abs(x)^p)^(1 / p)
}


# ---------------------------------------------------------
# 2. Interpolate the empirical CDF
# ---------------------------------------------------------
interp <- function(data, x_out) {
  if (!length(data) || anyNA(data) || any(!is.finite(data))) {
    stop("data must be non-empty and contain finite values.")
  }
  
  sorted_data <- sort(data)
  unique_data <- unique(sorted_data)
  unique_cdf <- ecdf(sorted_data)(unique_data)
  
  approx(
    x = unique_data,
    y = unique_cdf,
    xout = x_out,
    rule = 2,
    ties = "ordered"
  )$y
}


# ---------------------------------------------------------
# 3. CDF of generalized order statistics
# ---------------------------------------------------------
cdf_gos <- function(u, gamma) {
  if (!length(gamma) ||
      anyNA(gamma) ||
      any(!is.finite(gamma)) ||
      any(gamma <= 0) ||
      anyDuplicated(gamma)) {
    stop("gamma must contain distinct, finite, positive values.")
  }
  
  out <- numeric(length(u))
  inside <- which(u > 0 & u < 1)
  
  out[u >= 1] <- 1
  
  if (!length(inside)) {
    return(out)
  }
  
  # Compute coefficients only once per gamma vector
  r <- length(gamma)
  
  if (r == 1L) {
    coeff <- 1
  } else {
    coeff <- vapply(seq_len(r), function(i) {
      prod(gamma[-i] / (gamma[-i] - gamma[i]))
    }, numeric(1))
  }
  
  t <- -log1p(-u[inside])
  
  survival <- vapply(t, function(tt) {
    sum(coeff * exp(-gamma * tt))
  }, numeric(1))
  
  out[inside] <- 1 - survival
  
  # Protect against small floating-point errors
  pmin(1, pmax(0, out))
}


# ---------------------------------------------------------
# 4. Compute coefficients
# ---------------------------------------------------------
compute_coefficients_gos <- function(n, m) {
  if (length(n) != 1L || is.na(n) || n < 1 ||
      n != as.integer(n)) {
    stop("n must be a positive integer.")
  }
  
  if (length(m) != 1L || is.na(m) ||
      m < 1 || m > n || m != as.integer(m)) {
    stop("m must be an integer between 1 and n.")
  }
  
  coe <- matrix(0, nrow = n, ncol = m)
  u <- (0:n) / n
  
  for (i in seq_len(m)) {
    gamma <- n - seq_len(i) + 1
    cdf <- cdf_gos(u, gamma)
    coe[, i] <- diff(cdf)
  }
  
  coe
}

Test_IHR <- function(a,n,m,p) {
  gamma_list <- lapply(
    seq_len(m),
    function(i) n - seq_len(i) + 1
  )
  
  S <- function(i, n) {
    sum(1 / seq(n - i + 1, n))  # Harmonic series for the order statistics
  }
  
  # True "pi" based on exponential distribution
  TR <- sapply(1:m, function(i) 1 - exp(-S(i, n)))
  
  coe <- compute_coefficients_gos(n, m)
  
  # Storage for results
  EMP <- matrix(0, R, m)
  EXu <- matrix(0, R, m)
  
  # 5000 Simulate samples from Exp and compute statistics
  for (ii in 1:R) {
    aau <- sort(rexp(n,1))  # Sample from Exp(1) to have the critical values
    
    # Estimate the expected order statistics
    EXu[ii, ] <- sapply(1:m, function(i) sum(aau * coe[, i]))
    
    # Empirical CDF of expected order statistics
    #EMP[ii, ] <- sapply(1:m, function(i) ecdf(aau)(EXu[ii, i]))
    EMP[ii, ] <- sapply(1:m, function(i) interp(aau,EXu[ii, i]))
    
  }
  
  # To obtain the critical value for testing convexity and concavity
  if (p == Inf) {
    cx <- sapply(1:R, function(ii) max(pmax(0, (TR - EMP[ii, ])*(1))))  # Convexity test
    cv <- sapply(1:R, function(ii) max(pmax(0, (-TR + EMP[ii, ])*(1)))) # Concavity test
  } else {
    cx <- sapply(1:R, function(ii) lp_norm(pmax(0, (TR - EMP[ii, ])*(1)),p))  # Convexity test
    cv <- sapply(1:R, function(ii) lp_norm(pmax(0,(-TR + EMP[ii, ])*(1)),p)) # Concavity test
  }
  # Critical values (90th percentile)
  qcv <- quantile(cv, 0.9, na.rm = TRUE)
  qcx <- quantile(cx, 0.9, na.rm = TRUE)
  for (ii in 1:R) {
    aau <- sort(rweibull(n,a,1))  # Sample from Weibull distribution
    
    # Estimate the expected order statistics
    EXu[ii, ] <- sapply(1:m, function(i) sum(aau * coe[, i]))
    
    # Empirical CDF of expected order statistics
    EMP[ii, ] <- sapply(1:m, function(i) interp(aau,EXu[ii, i]))
  }
  
  # Power calculation for convexity and concavity
  # the positve patr of TR - EMP[ii, ]
  cx=0
  cv=0
  if (p == Inf) {
    cx <- sapply(1:R, function(ii) max(pmax(0, (TR - EMP[ii, ])*(1))))  # Convexity test
    cv <- sapply(1:R, function(ii) max(pmax(0, (-TR + EMP[ii, ])*(1)))) # Concavity test
  } else {
    cx <- sapply(1:R, function(ii) lp_norm(pmax(0, (TR - EMP[ii, ])*(1)),p))  # Convexity test
    cv <- sapply(1:R, function(ii) lp_norm(pmax(0,(-TR + EMP[ii, ])*(1)),p)) # Concavity test
  }
  v1 <- mean(cx >= qcx, na.rm = TRUE)  # Power for convexity
  v2 <- mean(cv >= qcv, na.rm = TRUE)  # Power for concavity
  
  return(data.frame(a=a,IHR_power = v1, DHR_power = v2))
}


R <- 500
a <- seq(0.5, 2, by = 0.1)
p <- 1
n <- 25


m=2

Test_IHR(1.5,n=50,m=2,1)



# Define the different 'm' values
m_values <- c(2,5,10,15,20,25,30,35,40,45,50,100)
#m_values <- c(2,5,10,15,20,25)
# Initialize an empty data frame for results
results <- data.frame()

# Loop through each value of m, generating data and appending to results
for (m in m_values) {
  temp_data <- do.call(rbind, lapply(a, function(i) Test_IHR(i, n, m, p)))
  temp_data <- as.data.frame(temp_data)  # Ensure it's a data frame
  temp_data <- mutate(temp_data, source = paste("m =", m))  # Add source column for each m value
  colnames(temp_data) <- colnames(results)  # Ensure column names are consistent
  results <- rbind(results, temp_data)  # Append to results
}


# Convert columns
colnames(results) <- c("a", "IHR_power", "DHR_power", "source")

# Plot the power of convexity and concavity
ggplot(results, aes(x = a, y = DHR_power, shape = factor(source))) +
  geom_line() +
  geom_point() +
  labs(title = "Power of DHR Test for n=100",
       x = "Shape parameter 'a'",
       y = "Power") +
  theme_minimal() +
  ylim(0, 1) +
  xlim(0.5, 2) +
  theme(legend.position = "right") +
  scale_shape_manual(name = "Source", values = c(1, 3, 4, 5, 16,25,18,24,78,10,11,14))
