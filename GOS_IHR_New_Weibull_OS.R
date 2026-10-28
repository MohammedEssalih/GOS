set.seed(12345)

# ---- parameters ------------------------------------------------------------
n      <- 100
mList  <- 2:floor(sqrt(n))
p      <- 1
alpha  <- 0.1
B0     <- 3000     # null simulations
B1     <- 1000     # simulations per alternative
shapes <- c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8)
shapes <- seq(0.3, 1, by = 0.05) 
alt    <- "weibull"   # "weibull" (as labelled) or "t" (as actually coded)

# ---- weights of the L-statistics -------------------------------------------
# W[[k]] is an m x n matrix: row j = diff of Beta cdf of the j-th order stat of m
tgrid <- (0:n) / n
weights <- lapply(mList, function(m) {
  t(sapply(1:m, function(j) diff(pbeta(tgrid, j, m - j + 1))))
})

# ---- theoretical bounds (G = Exp(1)) ---------------------------------------
H  <- function(m) sum(1 / seq_len(m))          # harmonic number, H(0) = 0
pi_list <- lapply(mList, function(m) {
  1 - exp(-sapply(1:m, function(j) H(m) - H(m - j)))
})

lp_norm <- function(x, p) sum(abs(x)^p)^(1 / p)

# ---- test statistic --------------------------------------------------------
stat <- function(x) {
  sx   <- sort(x)
  Fgrid <- (0:(n - 1)) / (n - 1)               # same ecdf values as Mathematica
  vals <- vapply(seq_along(mList), function(k) {
    mu  <- drop(weights[[k]] %*% sx)           # length m
    cdf <- approx(sx, Fgrid, xout = mu, rule = 2, ties = "ordered")$y
    lp_norm(pmax(0, cdf - pi_list[[k]]), p)
  }, numeric(1))
  max(vals)
}

# ---- calibration under G = Exp(1) ------------------------------------------
nullStats <- replicate(B0, stat(rexp(n, 1)))
critical  <- quantile(nullStats, 1 - alpha, type = 1, names = FALSE)  # type 1 ~ Mathematica
cat("Valore critico =", critical, "\n")

# ---- power -----------------------------------------------------------------
rgen <- function(a) {
  if (alt == "weibull") rweibull(n, shape = a, scale = 1) else rt(n, df = a)
}

power1 <- sapply(shapes, function(a) {
  mean(replicate(B1, stat(rgen(a))) > critical)
})

results <- data.frame(shape = shapes, rejection_rate = power1)
print(results)

# ---- plot ------------------------------------------------------------------
plot(results$shape, results$rejection_rate, type = "b", pch = 19,
     ylim = c(0, 1),
     xlab = if (alt == "weibull") "Weibull shape" else "t degrees of freedom",
     ylab = "Rejection rate")
abline(v = 1, h = alpha, lty = 2, col = "grey50")
lines(results$shape,power)
