set.seed(12345)

# ============================================================
# 1. Parameters
# ============================================================
n      <- 100
size   <- 2000                       # simulations per shape
B0     <- 5000                       # null simulations (more = more stable critical values)
level  <- 0.1
shapes <- seq(0.3, 1, by = 0.05)     # Weibull shapes in [0.3, 1]

# ============================================================
# 2. Statistic: max distance from the convex minorant
# ============================================================
d_stat <- function(a) {
  n  <- length(a)
  x  <- sort(a)
  Tt <- sum(x)
  k  <- seq_len(n)
  S  <- (cumsum(x) + (n - k) * x) / Tt
  S[n] <- 1
  h  <- (k - 1) / n
  dmax <- 0
  for (j in 2:n) {
    cq <- h[j]
    if (j < n) {
      i  <- seq_len(j - 1)
      kk <- (j + 1):n
      L <- outer(i, kk, function(ii, kk_)
        h[ii] + (S[j] - S[ii]) * (h[kk_] - h[ii]) / (S[kk_] - S[ii]))
      cq <- min(cq, L)
    }
    dmax <- max(dmax, abs(h[j] - cq))
  }
  dmax
}

# ============================================================
# 3. Null distribution under Exp(1) (computed once)
# ============================================================
sim_n <- sort(replicate(B0, d_stat(rexp(n, rate = 1))))

# p-value = P(null > z) = 1 - ECDF_null(z)
pval <- function(z) 1 - findInterval(z, sim_n) / length(sim_n)

# ============================================================
# 4. Power for each shape
# ============================================================
power <- vapply(shapes, function(a) {
  zz <- replicate(size, d_stat(rweibull(n, shape = a, scale = 1)))
  mean(pval(zz) <= level)
}, numeric(1))

se <- sqrt(power * (1 - power) / size)     # Monte Carlo standard error

power_table <- data.frame(shape = shapes, power = power, se = se)
print(round(power_table, 4), row.names = FALSE)

# ============================================================
# 5. Plot the power curve
# ============================================================
plot(shapes, power, type = "b", pch = 19, lwd = 2, ylim = c(0, 1),
     xlab = "Weibull shape", ylab = "Rejection rate",
     main = paste0("Power curve, n = ", n, ", level = ", level))
arrows(shapes, pmax(0, power - 1.96 * se), shapes, pmin(1, power + 1.96 * se),
       angle = 90, code = 3, length = 0.03, col = "grey50")
abline(h = level, lty = 2, col = "grey50")   # nominal level
abline(v = 1,     lty = 3, col = "grey50")   # null (exponential)