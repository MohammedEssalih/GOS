
# ============================================================
# 1. CDF of the generalized uniform order statistic
#    U_{r, gamma} = 1 - prod(U_i^(1/gamma_i))
# ============================================================




# ============================================================
# 2. Set parameters
# ============================================================

n <- 10
r <- 4

u <- seq(0, 1, length.out = 501)

# Classical order statistic:
# gamma_i = n - i + 1
gamma_classical <- n - seq_len(r) + 1



# ============================================================
# 3. Calculate the CDFs
# ============================================================

# Generalized order statistic, classical parameters
cdf_classical_gos <- cdf_gos(u, gamma_classical)

# Classical Beta CDF
cdf_beta <- pbeta(u, shape1 = r, shape2 = n - r + 1)



# Panel 1: classical GOS versus Beta CDF
plot(u, cdf_beta, type = "l", lwd = 3, col = "blue",
     xlab = "u", ylab = "CDF",
     main = "Classical order statistic",
     ylim = c(0, 1))
lines(u, cdf_classical_gos, col = "red", lwd = 2, lty = 2)



legend("topleft",
       legend = c("Beta CDF", "GOS CDF"),
       col = c("blue", "red"),
       lty = c(1, 2), lwd = c(3, 2), bty = "n")


# ============================================================
# 6. Simulation check
# ============================================================

set.seed(123)
N <- 200000

simulate_gos <- function(N, gamma) {
  r <- length(gamma)
  U <- matrix(runif(N * r), nrow = N, ncol = r)
  1 - exp(-rowSums(-log(U) / rep(gamma, each = N)))
}

# Same underlying construction, two parameter choices
sim_classical <- simulate_gos(N, gamma_classical)
sim_general <- simulate_gos(N, gamma_general)

# Empirical CDF versus theoretical CDF
ecdf_classical <- ecdf(sim_classical)
ecdf_general <- ecdf(sim_general)

plot(u, cdf_beta, type = "l", lwd = 3, col = "blue",
     xlab = "u", ylab = "CDF",
     main = "Simulation validation", ylim = c(0, 1))
lines(u, cdf_gos(u, gamma_classical),
      col = "red", lwd = 2, lty = 2)
lines(u, ecdf_classical(u),
      col = "black", lwd = 1, lty = 3)

legend("topleft",
       legend = c("Beta theoretical", "GOS theoretical",
                  "GOS empirical"),
       col = c("blue", "red", "black"),
       lty = c(1, 2, 3), lwd = c(3, 2, 1), bty = "n")

cat("Simulation check: classical GOS vs Beta\n")
cat("Empirical CDF at u = 0.5:",
    ecdf_classical(0.5), "\n")
cat("Beta CDF at u = 0.5:",
    pbeta(0.5, r, n - r + 1), "\n")








# Generalized coefficient matrix
compute_coefficients_gos <- function(n, m, gamma_list) {
  # gamma_list[[i]] contains (gamma_1, ..., gamma_i)
  # for the i-th generalized order statistic.
  stopifnot(length(gamma_list) >= m)
  
  coe <- matrix(0, nrow = n, ncol = m)
  
  for (i in seq_len(m)) {
    gamma <- gamma_list[[i]]
    
    if (length(gamma) != i)
      stop("gamma_list[[i]] must have length i.")
    
    for (j in seq_len(n)) {
      coe[j, i] <-
        cdf_gos(j / n, gamma) -
        cdf_gos((j - 1) / n, gamma)
    }
  }
  
  coe
}





cdf_gos <- function(u, gamma) {
  r <- length(gamma)
  
  # Formula below assumes distinct positive gamma values
  if (any(gamma <= 0)) stop("All gamma values must be positive.")
  if (anyDuplicated(gamma)) {
    stop("This formula requires distinct gamma values.")
  }
  
  out <- numeric(length(u))
  
  for (k in seq_along(u)) {
    x <- u[k]
    
    if (x <= 0) {
      out[k] <- 0
    } else if (x >= 1) {
      out[k] <- 1
    } else {
      t <- -log1p(-x)
      
      # Survival function of a sum of exponentials
      coeff <- sapply(seq_len(r), function(i) {
        if (r == 1) return(1)
        prod(gamma[-i] / (gamma[-i] - gamma[i]))
      })
      
      out[k] <- 1 - sum(coeff * exp(-gamma * t))
    }
  }
  
  out
}





#########################################

Value=Test_IHR(2,25,5,1,R=500)
n=50
m_values <- c(1, 5, 10, 20)

# Initialize an empty data frame for results
results <- data.frame()

# Loop through each value of m, generating data and appending to results
for (m in m_values) {
  temp_data <- do.call(rbind, lapply(a, function(i) Test_IHR(i, n, m, p,R=500)))
  temp_data <- as.data.frame(temp_data)  # Ensure it's a data frame
  temp_data <- mutate(temp_data, source = paste("m =", m))  # Add source column for each m value
  colnames(temp_data) <- colnames(results)  # Ensure column names are consistent
  results <- rbind(results, temp_data)  # Append to results
}


# Convert columns
colnames(results) <- c("a", "IHR_power", "DHR_power", "source")

# Plot the power of convexity and concavity
ggplot(results, aes(x = a, y = IHR_power, shape = factor(source))) +
  geom_line() +
  geom_point() +
  labs(title = "Power of IHR Test for n=25",
       x = "Shape parameter 'a'",
       y = "Power") +
  theme_minimal() +
  ylim(0, 1) +
  xlim(0, 2) +
  theme(legend.position = "right") +
  scale_shape_manual(name = "Source", values = c(1, 3, 4, 5, 16))



compute_coefficients_gos <- function(n, m, gamma_list) {
  stopifnot(length(gamma_list) >= m)
  
  coe <- matrix(0, nrow = n, ncol = m)
  
  for (i in seq_len(m)) {
    gamma <- gamma_list[[i]]
    
    if (length(gamma) != i) {
      stop("gamma_list[[i]] must have length i.")
    }
    
    for (j in seq_len(n)) {
      coe[j, i] <-
        cdf_gos(j / n, gamma) -
        cdf_gos((j - 1) / n, gamma)
    }
  }
  
  coe
}


# compute_coefficients_gos <- function(n, m, gamma_list) {
#   stopifnot(length(gamma_list) >= m)
#   
#   coe <- matrix(0, nrow = n, ncol = m)
#   
#   for (i in seq_len(m)) {
#     gamma <- gamma_list[[i]]
#     
#     stopifnot(length(gamma) == i)
#     
#     for (j in seq_len(n)) {
#       coe[j, i] <-
#         cdf_gos(j / n, gamma) -
#         cdf_gos((j - 1) / n, gamma)
#     }
#   }
#   
#   coe
# }

n <- 20
m <- 5

gamma_list <- lapply(
  seq_len(m),
  function(i) m - seq_len(i) + 1
)

gamma_list <- lapply(
  seq_len(m),
  function(i) n - seq_len(i) + 1
)

A <- compute_coefficients(n, m)
B <- compute_coefficients_gos(n, m, gamma_list)

all.equal(A, B, tolerance = 1e-10)
