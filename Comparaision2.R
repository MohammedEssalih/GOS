library(pracma)
library(flexsurv)
library(purrr)
library(tidyverse)

interp <- function(data, x_out) {
  # Sort data and calculate the empirical CDF
  sorted_data <- sort(data)
  cdf_values <- ecdf(sorted_data)(sorted_data) # Empirical CDF values
  
  # Use approx to perform linear interpolation
  interp_result <- approx(x = sorted_data, y = cdf_values, xout = x_out)
  
  return(interp_result$y) # Return interpolated CDF values
}



# Function to compute the coefficients for order statistics
compute_coefficients <- function(n, m) {
  coe <- matrix(0, n, m)
  for (j in 1:n) {
    for (i in 1:m) {
      coe[j, i] <- pbeta(j/n, i, 1 - i + m) - pbeta((j - 1)/n, i, 1 - i + m)
    }
  }
  return(coe)
}

# Function to compute the power for a given shape parameter 'a' and value of m
Lando_IHR <- function(a,n,m,p) {
  S <- function(i, m) {
    sum(1 / seq(m - i + 1, m))  # Harmonic series for the order statistics
  }
  
  # True "pi" based on exponential distribution
  TR <- sapply(1:m, function(i) 1 - exp(-S(i, m)))
  
  coe <- compute_coefficients(n, m)
  
  # Storage for results
  EMP <- matrix(0, 5000, m)
  EXu <- matrix(0, 5000, m)
  
  # 5000ulate samples from Exp and compute statistics
  for (ii in 1:5000) {
    aau <- sort(rexp(n,1))  # Sample from Gamma(a, 1)
    
    # Estimate the expected order statistics
    EXu[ii, ] <- sapply(1:m, function(i) sum(aau * coe[, i]))
    
    # Empirical CDF of expected order statistics
    # EMP[ii, ] <- sapply(1:m, function(i) ecdf(aau)(EXu[ii, i]))
    EMP[ii, ] <- sapply(1:m, function(i) interp(aau,EXu[ii, i]))
    
  }
  
  # To obtain the critical value for testing convexity and concavity
  if (p == Inf) {
    cx <- sapply(1:5000, function(ii) max(pmax(0, TR - EMP[ii, ])))  # Convexity test
    cv <- sapply(1:5000, function(ii) max(pmax(0, -TR + EMP[ii, ]))) # Concavity test
  } else {
    cx <- sapply(1:5000, function(ii) lp_norm(pmax(0, TR - EMP[ii, ]),p))  # Convexity test
    cv <- sapply(1:5000, function(ii) lp_norm(pmax(0,-TR + EMP[ii, ]),p)) # Concavity test
  }
  # Critical values (90th percentile)
  qcv <- quantile(cv, 0.9)
  qcx <- quantile(cx, 0.9)
  for (ii in 1:5000) {
    aau <- sort(rgamma(n,a,1))  # Sample from Weibull distribution
    
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
    cx <- sapply(1:5000, function(ii) max(pmax(0, TR - EMP[ii, ])))
    cv <- sapply(1:5000, function(ii) max(pmax(0, -TR + EMP[ii, ])))
  } else {
    cx <- sapply(1:5000, function(ii) lp_norm(pmax(0, TR - EMP[ii, ]),p))
    cv <- sapply(1:5000, function(ii) lp_norm(pmax(0, -TR + EMP[ii, ]),p))
  }
  v1 <- mean(cx >= qcx)  # Power for convexity
  v2 <- mean(cv >= qcv)  # Power for concavity
  
  return(data.frame(a=a,IHR_power = v1, DHR_power = v2))
}



Lando_IHR_W <- function(a,n,m,p) {
  S <- function(i, m) {
    sum(1 / seq(m - i + 1, m))  # Harmonic series for the order statistics
  }
  
  # True "pi" based on exponential distribution
  TR <- sapply(1:m, function(i) 1 - exp(-S(i, m)))
  
  coe <- compute_coefficients(n, m)
  
  # Storage for results
  EMP <- matrix(0, 5000, m)
  EXu <- matrix(0, 5000, m)
  
  # 5000ulate samples from Exp and compute statistics
  for (ii in 1:5000) {
    aau <- sort(rexp(n,1))  # Sample from Gamma(a, 1)
    
    # Estimate the expected order statistics
    EXu[ii, ] <- sapply(1:m, function(i) sum(aau * coe[, i]))
    
    # Empirical CDF of expected order statistics
    # EMP[ii, ] <- sapply(1:m, function(i) ecdf(aau)(EXu[ii, i]))
    EMP[ii, ] <- sapply(1:m, function(i) interp(aau,EXu[ii, i]))
    
  }
  
  # To obtain the critical value for testing convexity and concavity
  if (p == Inf) {
    cx <- sapply(1:5000, function(ii) max(pmax(0, TR - EMP[ii, ])*(1-TR)))  # Convexity test
    cv <- sapply(1:5000, function(ii) max(pmax(0, -TR + EMP[ii, ])*(1-TR))) # Concavity test
  } else {
    cx <- sapply(1:5000, function(ii) lp_norm(pmax(0, TR - EMP[ii, ])*(1-TR),p))  # Convexity test
    cv <- sapply(1:5000, function(ii) lp_norm(pmax(0,-TR + EMP[ii, ])*(1-TR),p)) # Concavity test
  }
  # Critical values (90th percentile)
  qcv <- quantile(cv, 0.9)
  qcx <- quantile(cx, 0.9)
  for (ii in 1:5000) {
    aau <- sort(rgamma(n,a,1))  # Sample from Weibull distribution
    
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
    cx <- sapply(1:5000, function(ii) max(pmax(0, TR - EMP[ii, ])*(1-TR)))
    cv <- sapply(1:5000, function(ii) max(pmax(0, -TR + EMP[ii, ])*(1-TR)))
  } else {
    cx <- sapply(1:5000, function(ii) lp_norm(pmax(0, TR - EMP[ii, ]),p)*(1-TR))
    cv <- sapply(1:5000, function(ii) lp_norm(pmax(0, -TR + EMP[ii, ]),p)*(1-TR))
  }
  v1 <- mean(cx >= qcx)  # Power for convexity
  v2 <- mean(cv >= qcv)  # Power for concavity
  
  return(data.frame(a=a,IHR_power = v1, DHR_power = v2))
}
Lando_IHR(1.2,100,8,1)
Lando_IHR_W(1.2,100,8,1)

