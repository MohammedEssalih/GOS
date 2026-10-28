set.seed(12345)

# ---- parameters ------------------------------------------------------------
n      <- 100
kList  <- c(1, 2, 3, 5)           # values of k (k-th upper records)
rList  <- 2:floor(sqrt(n))        # number of records used
rmax   <- max(rList)
p      <- 2
alpha  <- 0.1
B0     <- 3000                    # null simulations per k
B1     <- 1000                    # power simulations per (k, shape)
shapes <- c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8)
shapes <- seq(0.3, 1, by = 0.05) 
lp_norm <- function(x, p) sum(abs(x)^p)^(1 / p)

# ---- everything for one value of k -----------------------------------------
make_test <- function(k) {
  z  <- -log1p(-(0:n) / n)                              # Inf at the last point
  W  <- t(sapply(1:rmax, function(j) diff(pgamma(k * z, shape = j))))  # rmax x n
  pi_rec <- 1 - exp(-(1:rmax) / k)                      # G(E[Y_j]), E[Y_j] = j/k
  Fgrid  <- (0:(n - 1)) / (n - 1)
  
  stat <- function(x) {
    sx  <- sort(x)
    mu  <- drop(W %*% sx)
    cdf <- approx(sx, Fgrid, xout = mu, rule = 2, ties = "ordered")$y
    d   <- pmax(0, cdf - pi_rec)
    max(vapply(rList, function(r) lp_norm(d[1:r], p), numeric(1)))
  }
  
  nullStats <- replicate(B0, stat(rexp(n, 1)))
  critical  <- quantile(nullStats, 1 - alpha, type = 1, names = FALSE)
  
  power <- sapply(shapes, function(a)
    mean(replicate(B1, stat(rweibull(n, shape = a, scale = 1))) > critical))
  
  list(critical = critical, power = power)
}

# ---- run for all k ---------------------------------------------------------
res <- lapply(kList, make_test)
names(res) <- paste0("k=", kList)

critical_values <- sapply(res, `[[`, "critical")
power_table <- data.frame(shape = shapes,
                          sapply(res, `[[`, "power"), check.names = FALSE)

cat("Valori critici:\n"); print(round(critical_values, 4))
cat("\nRejection rates:\n"); print(power_table)

# ---- plot ------------------------------------------------------------------
cols <- seq_along(kList)
matplot(power_table$shape, power_table[, -1], type = "b", pch = 19, lty = 1,
        col = cols, ylim = c(0, 1),
        xlab = "Weibull shape", ylab = "Rejection rate")
lines(power_table$shape,power)
abline(v = 1, h = alpha, lty = 2, col = "grey50")
legend("topright", legend = paste("k =", kList), col = cols, pch = 19, lty = 1)
