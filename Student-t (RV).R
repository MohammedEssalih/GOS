set.seed(12345)

# ============================================================
# 1. Parameters
# ============================================================
n       <- 100
kList   <- c(1, 2, 3, 5)
rList   <- 2:floor(sqrt(n))
rmax    <- max(rList)
p       <- 2
alpha   <- 0.05
B0      <- 3000
B1      <- 1000
df_vals <- seq(0.1, 1, by = 0.1)

# Lp norm
lp_norm <- function(x, p) {
  sum(abs(x)^p)^(1 / p)
}

# ============================================================
# 2. Construct and run the test for one k
# ============================================================
make_test <- function(k) {
  
  # Grid: the final point is Inf
  z <- -log1p(-(0:n) / n)
  
  # Weight matrix: rmax x n
  W <- t(vapply(
    seq_len(rmax),
    function(j) {
      diff(pgamma(k * z, shape = j))
    },
    numeric(n)
  ))
  
  # Theoretical record probabilities
  pi_rec <- 1 - exp(-(seq_len(rmax)) / k)
  
  # Empirical CDF grid
  Fgrid <- (0:(n - 1)) / (n - 1)
  
  # Test statistic
  stat <- function(x) {
    
    sx <- sort(x)
    
    # Weighted estimates
    mu <- drop(W %*% sx)
    
    # Estimate the CDF at the weighted estimates
    cdf_est <- approx(
      x = sx,
      y = Fgrid,
      xout = mu,
      rule = 2,
      ties = "ordered"
    )$y
    
    # Positive deviations, as in the original test
    d <- pmax(0, cdf_est - pi_rec)
    
    # Maximize over the selected number of records
    max(vapply(
      rList,
      function(r) lp_norm(d[seq_len(r)], p),
      numeric(1)
    ))
  }
  
  # ==========================================================
  # 3. Null simulations: exponential distribution
  # ==========================================================
  null_stats <- replicate(
    B0,
    stat(rexp(n, rate = 1))
  )
  
  critical <- unname(
    quantile(null_stats, probs = 1 - alpha, type = 1)
  )
  
  # ==========================================================
  # 4. Power simulations: Student-t alternatives
  # ==========================================================
  power <- vapply(
    df_vals,
    function(df) {
      
      alt_stats <- replicate(
        B1,
        stat(rt(n, df = df))
      )
      
      mean(alt_stats > critical)
    },
    numeric(1)
  )
  
  list(
    critical = critical,
    power = power
  )
}

# ============================================================
# 5. Run simulations for all k
# ============================================================
results <- lapply(kList, make_test)
names(results) <- paste0("k=", kList)

critical_values <- vapply(
  results,
  function(x) x$critical,
  numeric(1)
)

power_matrix <- sapply(
  results,
  function(x) x$power
)

power_table <- data.frame(
  df = df_vals,
  power_matrix,
  check.names = FALSE
)

names(power_table) <- c("df", paste0("k=", kList))

# ============================================================
# 6. Display results
# ============================================================
cat("\nCritical values:\n")
print(round(critical_values, 4))

cat("\nEstimated rejection probabilities:\n")
print(power_table, row.names = FALSE)

# ============================================================
# 7. Plot empirical power
# ============================================================
cols <- seq_along(kList)

matplot(
  x = power_table$df,
  y = as.matrix(power_table[, -1]),
  type = "b",
  pch = 19,
  lty = 1,
  lwd = 2,
  col = cols,
  ylim = c(0, 1),
  xlab = "Student-t degrees of freedom",
  ylab = "Estimated rejection probability"
)

abline(h = alpha, lty = 2, col = "grey50")

legend(
  "bottomright",
  legend = paste("k =", kList),
  col = cols,
  pch = 19,
  lty = 1,
  lwd = 2,
  bty = "n"
)
