A <- matrix(c(
  -3, -3, -6,  6,  1,
  -1, -1, -1,  1, -2,
  0,  0, -1,  1, -1,
  0,  0,  1,  1, -1
), nrow = 4, byrow = TRUE)

B <- matrix(c(
  4,  0, -1, -2,  0,
  -5,  0,  3,  5,  0,
  2,  0, -1, -2,  0,
  6,  0, -3, -6,  0
), nrow = 4, byrow = TRUE)

b <- matrix(c(6, -1, -4, 6), nrow = 4)

pinv_cours <- function(M, tol = 1e-8) {
  MtM <- t(M) %*% M                
  eg  <- eigen(MtM, symmetric = TRUE)
  
  lambda <- eg$values
  V      <- eg$vectors
  
  sigma <- sqrt(pmax(lambda, 0))    
  r <- sum(sigma > tol)             
  
  Vr <- V[, 1:r, drop = FALSE]
  Dinv <- diag(1 / sigma[1:r], nrow = r, ncol = r)
  
  Ur <- M %*% Vr %*% Dinv          
  
  Mplus <- Vr %*% Dinv %*% t(Ur)   
  
  list(pinv = Mplus, rank = r, singular_values = sigma)
  
}

resA <- pinv_cours(A, tol = 1e-8)
resB <- pinv_cours(B, tol = 1e-8)

Aplus <- resA$pinv
Bplus <- resB$pinv

cat("rang(A) =", resA$rank, "\n")
cat("sigma(A) =", signif(resA$singular_values, 10), "\n\n")

cat("rang(B) =", resB$rank, "\n")
cat("sigma(B) =", signif(resB$singular_values, 10), "\n\n")

xA <- Aplus %*% b 
xB <- Bplus %*% b

cat("Solution x pour Ax=b :\n")
print(signif(xA, 10))

cat("\nSolution x pour Bx=b :\n")
print(signif(xB, 10))