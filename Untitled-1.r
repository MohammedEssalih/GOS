
gos_generate <- function(gamma, Finv) {
  r <- length(gamma)
  U <- runif(r)
  
  x <- numeric(r)
  prod_term <- 1
  
  for(i in 1:r){
    prod_term <- prod_term * U[i]^(1/gamma[i])
    x[i] <- Finv(1 - prod_term)
  }
  
  return(x)
}

Finv <- function(u) -log(1-u)

gos_generate(c(5,4,3,2,1), Finv)
