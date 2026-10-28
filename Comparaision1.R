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
    d   <- (TR - EMP) 
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
    alt <- simulate_stats(function(n) rgamma(n, ai, 1), R, n, m, coe, TR, p)
    data.frame(a         = ai,
               IHR_power = mean(alt[, "cx"] > qcx, na.rm = TRUE),
               DHR_power = mean(alt[, "cv"] > qcv, na.rm = TRUE))
  })
  do.call(rbind, res)
}




simulate_stats_W <- function(rgen, R, n, m, coe, TR, p) {
  out <- matrix(0, nrow = R, ncol = 2,
                dimnames = list(NULL, c("cx", "cv")))
  for (ii in seq_len(R)) {
    aau <- sort(rgen(n))
    EXu <- drop(aau %*% coe)                   # expected GOS, length m
    EMP <- interp(aau, EXu)                    # ecdf at expected GOS
    d   <- (TR - EMP) * (1-TR)
    out[ii, "cx"] <- lp_norm(pmax(0,  d), p)   # convexity test
    out[ii, "cv"] <- lp_norm(pmax(0, -d), p)   # concavity test
  }
  out
}

Test_IHR_W <- function(a, n, m, p, k, R = 5000, alpha = 0.1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # True "pi" based on the exponential distribution
  TR  <- 1 - exp(-seq_len(m) / k)
  coe <- compute_coefficients_gos(n, m, k)
  
  # Critical values (computed once): sample from Exp(1)
  null <- simulate_stats_W(function(n) rexp(n, 1), R, n, m, coe, TR, p)
  qcx  <- quantile(null[, "cx"], 1 - alpha, na.rm = TRUE)
  qcv  <- quantile(null[, "cv"], 1 - alpha, na.rm = TRUE)
  
  # Power: sample from Weibull(a, 1) for each a
  res <- lapply(a, function(ai) {
    alt <- simulate_stats_W(function(n) rgamma(n, ai, 1), R, n, m, coe, TR, p)
    data.frame(a         = ai,
               IHR_power = mean(alt[, "cx"] > qcx, na.rm = TRUE),
               DHR_power = mean(alt[, "cv"] > qcv, na.rm = TRUE))
  })
  do.call(rbind, res)
}

# ---------------------------------------------------------
# Example
# ---------------------------------------------------------
Test_IHR(a = 1.2, n = 100, m = 8, p = 1, k = 1, R = 5000, seed = 1)
Test_IHR_W(a = 1.2, n = 100, m = 8, p = 1, k = 1, R = 5000, seed = 1)
a=seq(0.5,2,0.1)
n=100
m=10
k=5
p=1
R=5000
y=lapply(a,function (i) Test_IHR_W(i,n,m,p,k,R)$IHR_power)
plot(a,y,type="l",col="red")
t=lapply(a, function(i) Test_IHR(i,n,m,k,p,R)$IHR_power)
lines(a,t)
