# Piecewise Weibull-Exponential distribution

pWE <- function(x, b = 2, t = 0.8) {
  x <- pmax(x, 0)
  H <- ifelse(x < t,
              x^(b+1)/((b+1)*t^b),
              x - b*t/(b+1))
  1 - exp(-H)
}

dWE <- function(x, b = 2, t = 0.8) {
  h <- ifelse(x < t, (pmax(x, 0)/t)^b, 1)
  ifelse(x < 0, 0, h * (1 - pWE(x, b, t)))
}

qWE <- function(u, b = 2, t = 0.8) {
  z <- -log1p(-u)
  ifelse(z <= t/(b+1),
         t*((b+1)*z/t)^(1/(b+1)),
         z + b*t/(b+1))
}

rWE <- function(n, b = 2, t = 0.8) {
  qWE(runif(n), b, t)
}

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
# 3. CDF of the r-th generalized order statistic (equal gamma = k)
#    -log(1-U_(r)) ~ Gamma(r, rate = k)  =>  P(U_(r) <= u) = pgamma(-k*log(1-u)... )
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
    out[inside] <- pgamma(z, shape = r)      # = 1 - exp(-z) * sum_{j<r} z^j/j!
  }
  out
}

# ---------------------------------------------------------
# 4. Coefficients: n x m matrix, column i = weights for the i-th GOS
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
  for (i in seq_len(m)) {
    coe[, i] <- diff(cdf_gos_equal(u, i, k))   # FIX: i, not m
  }
  coe
}

# ---------------------------------------------------------
# 5. Simulate the test statistics (convexity / concavity) from a sampler
# ---------------------------------------------------------
simulate_stats <- function(rgen, R, n, m, coe, TR, p) {
  out <- matrix(0, nrow = R, ncol = 2,
                dimnames = list(NULL, c("cx", "cv")))
  for (ii in seq_len(R)) {
    aau <- sort(rgen(n))
    EXu <- drop(aau %*% coe)                   # expected GOS, length m
    EMP <- interp(aau, EXu)                    # ecdf at expected GOS
    d   <- (TR - EMP) * (1)
    out[ii, "cx"] <- lp_norm(pmax(0,  d), p)   # convexity test
    out[ii, "cv"] <- lp_norm(pmax(0, -d), p)   # concavity test
  }
  out
}

# ---------------------------------------------------------
# 6. Test of IHR / DHR: critical values under Exp(1), power under Weibull(a)
# ---------------------------------------------------------
Test_IHR <- function(a, n, m, p, k, R = 5000, alpha = 0.1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # True "pi" based on the exponential distribution
  TR  <- 1 - exp(-seq_len(m) / k)
  coe <- compute_coefficients_gos(n, m, k)
  
  # Critical values (computed once): sample from Exp(1)
  null <- simulate_stats(function(n) rexp(n, 1), R, n, m, coe, TR, p)
  qcx  <- quantile(null[, "cx"], 1 - alpha, na.rm = TRUE)
  qcv  <- quantile(null[, "cv"], 1 - alpha, na.rm = TRUE)
  
  # Power: sample from Weibull(a, 1) for each a
  res <- lapply(a, function(ai) {
    alt <- simulate_stats(function(n) rWE(n, 2, 0.8), R, n, m, coe, TR, p)
    data.frame(a         = ai,
               IHR_power = mean(alt[, "cx"] > qcx, na.rm = TRUE),
               DHR_power = mean(alt[, "cv"] > qcv, na.rm = TRUE))
  })
  do.call(rbind, res)
}

# ---------------------------------------------------------
# Example
# ---------------------------------------------------------
Test_IHR(a = 1, n = 100, m = 8, p = 1, k = 24, R = 5000, seed = 1)

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
    aau <- sort(rWE(n,2,0.8))  # Sample from Weibull distribution
    
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
Lando_IHR(1,100,24,1)
