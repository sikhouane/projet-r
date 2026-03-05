A <- matrix(c(
  -18,  13,  -4,   4,
  2,  19,  -4,  12,
  -14,  11, -12,   8,
  -2,  21,   4,   8
), nrow = 4, byrow = TRUE)

B <- matrix(c(
  6, -8, -4,  5, -4,
  2,  7, -5, -6,  4,
  0, -1, -8,  2,  2,
  -1, -2,  4,  4, -8
), nrow = 4, byrow = TRUE)

svd_cours <- function(M, tol = 1e-10){
  
  AtA <- t(M) %*% M
  
  eg <- eigen(AtA)
  lambda <- eg$values
  V <- eg$vectors
  
  sigma <- sqrt(pmax(lambda,0))
  
  r <- sum(sigma > tol)
  
  Vr <- V[,1:r]
  
  D <- diag(sigma[1:r])
  
  Dinv <- diag(1/sigma[1:r])
  
  Ur <- M %*% Vr %*% Dinv
  
  M_rec <- Ur %*% D %*% t(Vr)
  
  list(
    U = Ur,
    D = D,
    V = Vr,
    singular_values = sigma,
    rank = r,
    reconstruction_error = norm(M - M_rec,"F")
  )
}

resA <- svd_cours(A)
resB <- svd_cours(B)

cat("----- A -----\n")
print(resA$singular_values)
cat("rang =",resA$rank,"\n")
cat("erreur =",resA$reconstruction_error,"\n\n")

cat("----- B -----\n")
print(resB$singular_values)
cat("rang =",resB$rank,"\n")
cat("erreur =",resB$reconstruction_error,"\n")