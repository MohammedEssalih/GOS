set.seed(12345)

# ============================================================
# 1. Parameters
# ============================================================
n       <- 100L
mList   <- seq_len(floor(sqrt(n)))   # m = 1, ..., floor(sqrt(n))
mmax    <- max(mList)
p       <- 2
alpha   <- 0.05
B0      <- 3000L
B1      <- 1000L
df_vals <- seq(0.1, 1, by = 0.1)

lp_norm <- function(x, p) {
  sum(abs(x)^p)^(1 / p)
}

# ============================================================
# 2. Weights and theoretical bounds
# ============================================================
tgrid <- (0:n) / n

build_weights <- function(m, tgrid) {
  vapply(
    seq_len(m),
    function(j) diff(pbeta(tgrid, shape1 = j, shape2 = m - j + 1)),
    numeric(length(tgrid) - 1L)
  ) |> t()
}

W <- lapply(mList, build_weights, tgrid = tgrid)

H <- c(0, cumsum(1 / seq_len(mmax)))   # H[m + 1] = H_m
pi_os <- lapply(mList, function(m) {
  j <- seq_len(m)
  1 - exp(-(H[m + 1L] - H[m - j + 1L]))
})

Fgrid <- seq(0, 1, length.out = n)

# ============================================================
# 3. Test statistic: max over m of the L^p norm
# ============================================================
stat <- function(x) {
  sx <- sort(x)
  max(vapply(seq_along(mList), function(i) {
    mu  <- as.numeric(W[[i]] %*% sx)
    cdf <- approx(
      x = sx,
      y = Fgrid,
      xout = mu,
      method = "linear",
      rule = 2,
      ties = "ordered"
    )$y

    lp_norm(pmax(0, cdf - pi_os[[i]]), p)
  }, numeric(1)))
}

# ============================================================
# 4. Null simulations: exponential distribution
# ============================================================
null_stats <- replicate(B0, stat(rexp(n, rate = 1)))
critical   <- unname(stats::quantile(null_stats, probs = 1 - alpha, type = 1))

# ============================================================
# 5. Power simulations: Student-t alternatives
# ============================================================
power <- vapply(df_vals, function(df) {
  mean(replicate(B1, stat(rt(n, df = df))) > critical)
}, numeric(1))

power_table <- data.frame(
  df = df_vals,
  rejection_rate = power
)

# ============================================================
# 6. Display results
# ============================================================
cat("\nCritical value:", round(critical, 4), "\n")
cat("\nEstimated rejection probabilities:\n")
print(power_table, row.names = FALSE)

# ============================================================
# 7. Plot
# ============================================================
plot(
  power_table$df,
  power_table$rejection_rate,
  type = "b",
  pch = 19,
  lwd = 2,
  ylim = c(0, 1),
  xlab = "Student-t degrees of freedom",
  ylab = "Estimated rejection probability"
)
abline(h = alpha, lty = 2, col = "grey50")
