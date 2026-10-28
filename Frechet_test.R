
library(VGAM)
# ---------------------------------------------------------
# 1. Lp norm
# ---------------------------------------------------------
lp_norm <- function(x, p = 2) {
  if (length(p) != 1L || is.na(p) || p < 1) {
    stop("p must be a single number greater than or equal to 1.")
  }
  if (is.infinite(p)) return(max(abs(x)))
  sum(abs(x)^p)^(1 / p)
}

# ---------------------------------------------------------
# 2. Interpolated empirical CDF
# ---------------------------------------------------------
interp <- function(data, x_out) {
  if (!length(data) || anyNA(data) || any(!is.finite(data))) {
    stop("data must be non-empty and contain finite values.")
  }
  sorted_data <- sort(data)
  unique_data <- unique(sorted_data)
  unique_cdf  <- ecdf(sorted_data)(unique_data)
  approx(x = unique_data, y = unique_cdf, xout = x_out,
         rule = 2, ties = "ordered")$y
}

# ---------------------------------------------------------
# 3. CDF of the r-th GOS (equal gamma = k)
# ---------------------------------------------------------
cdf_gos_equal <- function(u, r, k) {
  if (length(r) != 1L || is.na(r) || !is.finite(r) || r < 1 || r != floor(r)) {
    stop("r must be a positive integer.")
  }
  if (length(k) != 1L || is.na(k) || !is.finite(k) || k <= 0) {
    stop("k must be a finite positive number.")
  }
  out <- numeric(length(u))
  out[u >= 1] <- 1
  inside <- which(u > 0 & u < 1)
  if (length(inside)) {
    z <- -k * log1p(-u[inside])
    out[inside] <- pgamma(z, shape = r)
  }
  out
}

# ---------------------------------------------------------
# 4. Coefficients: n x m matrix
# ---------------------------------------------------------
compute_coefficients_gos <- function(n, m, k) {
  if (length(n) != 1L || is.na(n) || n < 1 || n != as.integer(n)) {
    stop("n must be a positive integer.")
  }
  if (length(m) != 1L || is.na(m) || m < 1 || m > n || m != as.integer(m)) {
    stop("m must be an integer between 1 and n.")
  }
  u   <- (0:n) / n
  coe <- matrix(0, nrow = n, ncol = m)
  for (i in seq_len(m)) coe[, i] <- diff(cdf_gos_equal(u, i, k))
  coe
}

# ---------------------------------------------------------
# 5. Log-logistic LL(a, 1) generator: cdf x^a / (1 + x^a)
# ---------------------------------------------------------
rll <- function(n, a = 1) {
  u <- runif(n)
  (u / (1 - u))^(1 / a)
}

# ---------------------------------------------------------
# 6. Simulated statistics (convexity / concavity of the odds-type transform)
# ---------------------------------------------------------
simulate_stats <- function(rgen, R, n, m, coe, TR, p) {
  out <- matrix(0, nrow = R, ncol = 2, dimnames = list(NULL, c("cx", "cv")))
  for (ii in seq_len(R)) {
    aau <- sort(rgen(n))
    EXu <- drop(aau %*% coe)
    EMP <- interp(aau, EXu)
    d   <- (TR - EMP) * (1 - TR)
    out[ii, "cx"] <- lp_norm(pmax(0,  d), p)
    out[ii, "cv"] <- lp_norm(pmax(0, -d), p)
  }
  out
}

# ---------------------------------------------------------
# 7. Test of IOR / DOR
#    G = LL(1,1);  critical values under LL(1,1);  power under LL(a,1)
# ---------------------------------------------------------
Test_IOR <- function(a, n, m, p, k, R = 5000, alpha = 0.1, seed = NULL) {
  if (k <= 1) stop("k must be > 1 for the LL(1,1) baseline (finite mean of G^{-1}(U_r)).")
  if (!is.null(seed)) set.seed(seed)
  
  # TR_r = G( E[G^{-1}(U_r)] ) = 1 - ((k-1)/k)^r  for G = LL(1,1)
  TR  <- 1 - ((k - 1) / k)^seq_len(m)
  coe <- compute_coefficients_gos(n, m, k)
  
  # Critical values (computed once): null LL(1,1)
  null <- simulate_stats(function(n) rll(n, 1), R, n, m, coe, TR, p)
  qcx  <- quantile(null[, "cx"], 1 - alpha, na.rm = TRUE)
  qcv  <- quantile(null[, "cv"], 1 - alpha, na.rm = TRUE)
  
  # Power: alternative LL(a, 1)
  res <- lapply(a, function(ai) {
    alt <- simulate_stats(function(n) rfrechet(n, 0,1,1), R, n, m, coe, TR, p)
    data.frame(a         = ai,
               IOR_power = mean(alt[, "cx"] > qcx, na.rm = TRUE),
               DOR_power = mean(alt[, "cv"] > qcv, na.rm = TRUE))
  })
  do.call(rbind, res)
}

# ---------------------------------------------------------
# Example
# ---------------------------------------------------------
Test_IOR(a = 1, n = 100, m = 10, p = 1, k = 15, R = 5000, seed = 1)

# Power curve (a < 1 DOR, a = 1 null, a > 1 IOR)
# Test_IOR(a = seq(0.5, 2, by = 0.1), n = 100, m = 1, p = 1, k = 10, R = 5000, seed = 1)

library(flexsurv)

library(purrr)
# Function to compute the coefficients for order statistics

interp <- function(data, x_out) {
  # Sort data and calculate the empirical CDF
  sorted_data <- sort(data)
  cdf_values <- ecdf(sorted_data)(sorted_data) # Empirical CDF values
  
  # Use approx to perform linear interpolation
  interp_result <- approx(x = sorted_data, y = cdf_values, xout = x_out)
  
  return(interp_result$y) # Return interpolated CDF values
}

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
Lando_IOR <- function(a,n,m,p) {
  
  # True "pi" based on exponential distribution
  TR <- sapply(1:m, function(i) i/ m)
  
  coe <- compute_coefficients(n, m)
  
  # Storage for results
  EMP <- matrix(0, 5000, m)
  EXu <- matrix(0, 5000, m)
  # the log logistic distribution 
  
  # 5000ulate samples from Exp and compute statistics
  for (ii in 1:5000) {
    aau <- sort(rllogis(n,1,1))
    
    # Estimate the expected order statistics
    EXu[ii, ] <- sapply(1:m, function(i) sum(aau * coe[, i]))
    
    # Empirical CDF of expected order statistics
    EMP[ii, ] <- sapply(1:m, function(i) interp(aau,EXu[ii, i]))
  }
  
  # To obtain the critical value for testing convexity and concavity
  if (p == Inf) {
    CX <- sapply(1:5000, function(ii) max(pmax(0, TR - EMP[ii, ])))  # Convexity test
    CV <- sapply(1:5000, function(ii) max(pmax(0, -TR + EMP[ii, ]))) # Concavity test
  } else {
    CX <- sapply(1:5000, function(ii) lp_norm(pmax(0, TR - EMP[ii, ]),p))  # Convexity test
    CV <- sapply(1:5000, function(ii) lp_norm(pmax(0,-TR + EMP[ii, ]),p)) # Concavity test
  }
  # Critical values (90th percentile)
  qcv <- quantile(CV, 0.9)
  qcx <- quantile(CX, 0.9)
  for (ii in 1:5000) {
    aau <- sort(rfrechet(n,0,a,1))  # Sample from Weibull distribution
    
    # Estimate the expected order statistics
    EXu[ii, ] <- sapply(1:m, function(i) sum(aau * coe[, i]))
    
    # Empirical CDF of expected order statistics
    EMP[ii, ] <- sapply(1:m, function(i) interp(aau,EXu[ii, i]))
  }
  
  # Power calculation for convexity and concavity
  # the positve patr of TR - EMP[ii, ]
  
  if (p == Inf) {
    cx <- sapply(1:5000, function(ii) max(pmax(0, TR - EMP[ii, ])))
    cv <- sapply(1:5000, function(ii) max(pmax(0, -TR + EMP[ii, ])))
  } else {
    cx <- sapply(1:5000, function(ii) lp_norm(pmax(0, TR - EMP[ii, ]),p))
    cv <- sapply(1:5000, function(ii) lp_norm(pmax(0, -TR + EMP[ii, ]),p))
  }
  v1 <- mean(cx >=qcx)#+ sum(rbernoulli(sum(cx == qcx),0.1))/5000  # Power for convexity
  v2 <- mean(cv >= qcv)  # Power for concavity
  
  return(data.frame(a=a,IOR_power = v1, DOR_power = v2))
}



Lando_IOR(1,100,10,1)
