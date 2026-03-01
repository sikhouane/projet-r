############################################################
# Méthode:
# 1) AtA symétrique
# 2) Jacobi => valeurs/vecteurs propres de AtA
# 3) sigma_i = sqrt(lambda_i)
# 4) u_i = A v_i / ||A v_i||  (exactement comme dans le cours)
# 5) A = U_r D V_r^T   (SVD réduite)
############################################################

# (p,q) max hors diagonale
max_offdiag <- function(S) {
  n <- nrow(S)
  max_val <- 0
  p <- 1; q <- 2
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      v <- abs(S[i, j])
      if (v > max_val) {
        max_val <- v
        p <- i; q <- j
      }
    }
  }
  list(max_val = max_val, p = p, q = q)
}

# Jacobi pour matrice symétrique: S = V diag(vals) V^T
jacobi_eigen <- function(S, tol = 1e-12, max_iter = 20000) {
  if (nrow(S) != ncol(S)) stop("S doit être carrée.")
  if (max(abs(S - t(S))) > 1e-8) stop("S doit être symétrique.")
  
  n <- nrow(S)
  A <- S
  V <- diag(1, n)
  
  for (iter in 1:max_iter) {
    info <- max_offdiag(A)
    if (info$max_val < tol) break
    
    p <- info$p; q <- info$q
    app <- A[p,p]; aqq <- A[q,q]; apq <- A[p,q]
    if (abs(apq) < tol) next
    
    # angle phi = 1/2 atan2(2 apq, aqq - app)
    phi <- 0.5 * atan2(2 * apq, (aqq - app))
    c <- cos(phi); s <- sin(phi)
    
    # mise à jour A (symétrique)
    for (k in 1:n) {
      if (k != p && k != q) {
        aik <- A[k,p]
        akq <- A[k,q]
        A[k,p] <- c * aik - s * akq
        A[p,k] <- A[k,p]
        A[k,q] <- s * aik + c * akq
        A[q,k] <- A[k,q]
      }
    }
    
    A[p,p] <- c*c*app - 2*s*c*apq + s*s*aqq
    A[q,q] <- s*s*app + 2*s*c*apq + c*c*aqq
    A[p,q] <- 0
    A[q,p] <- 0
    
    # mise à jour V
    for (k in 1:n) {
      vkp <- V[k,p]
      vkq <- V[k,q]
      V[k,p] <- c * vkp - s * vkq
      V[k,q] <- s * vkp + c * vkq
    }
  }
  
  list(vals = diag(A), vecs = V)
}

# Norme euclidienne d'un vecteur
norm2 <- function(x) sqrt(sum(x * x))

############################################################
# SVD réduite selon le cours (sans Schmidt)
############################################################
svd_cours_sans_schmidt <- function(A, tol_rank = 1e-10,
                                   jacobi_tol = 1e-14, jacobi_max_iter = 50000) {
  m <- nrow(A); n <- ncol(A)
  AtA <- t(A) %*% A
  
  eig <- jacobi_eigen(AtA, tol = jacobi_tol, max_iter = jacobi_max_iter)
  V <- eig$vecs
  
  # Normaliser les colonnes de V (au cas où)
  for (j in 1:n) {
    nj <- sqrt(sum(V[, j] * V[, j]))
    if (nj > 0) V[, j] <- V[, j] / nj
  }
  
  # IMPORTANT (conforme au cours): sigma_i = ||A v_i||
  sigma <- rep(0, n)
  Av_list <- vector("list", n)
  for (i in 1:n) {
    avi <- A %*% V[, i, drop = FALSE]
    Av_list[[i]] <- avi
    sigma[i] <- sqrt(sum(avi * avi))
  }
  
  # Trier décroissant
  ord <- order(sigma, decreasing = TRUE)
  sigma <- sigma[ord]
  V <- V[, ord, drop = FALSE]
  Av_list <- Av_list[ord]
  
  # rang numérique
  r <- sum(sigma > tol_rank)
  if (r == 0) stop("Rang nul (toutes valeurs singulières ~ 0).")
  
  Vr <- V[, 1:r, drop = FALSE]
  Dr <- diag(sigma[1:r], r, r)
  
  # Construire U_r : u_i = (A v_i) / ||A v_i||  (exact cours)
  Ur <- matrix(0, m, r)
  for (i in 1:r) {
    Ur[, i] <- (Av_list[[i]] / sigma[i])[, 1]
  }
  
  # Reconstruction + pseudo-inverse
  A_recons <- Ur %*% Dr %*% t(Vr)
  A_pinv <- Vr %*% diag(1 / sigma[1:r], r, r) %*% t(Ur)
  
  list(
    U = Ur, D = Dr, V = Vr,
    singular_values = sigma[1:r],
    rank = r,
    A_recons = A_recons,
    recons_error_fro = sqrt(sum((A - A_recons)^2)),
    A_pinv = A_pinv
  )
}

############################################################
# Données de l'exercice 4.2 (page 5 du PDF)
############################################################

A <- matrix(c(
  -18, 13, -4,  4,
  2, 19, -4, 12,
  -14, 11, -12,  8,
  -2, 21,  4,  8
), nrow = 4, byrow = TRUE)

B <- matrix(c(
  6, -8, -4,  5, -4,
  2,  7, -5, -6,  4,
  0, -1, -8,  2,  2,
  -1, -2,  4,  4, -8
), nrow = 4, byrow = TRUE)

resA <- svd_cours_sans_schmidt(A)
resB <- svd_cours_sans_schmidt(B)

cat("===== SVD réduite de A =====\n")
cat("rang(A) =", resA$rank, "\n")
cat("valeurs singulières(A) =\n"); print(resA$singular_values)
cat("Erreur Frobenius ||A - UDV^T|| = ", resA$recons_error_fro, "\n\n")

cat("U_A =\n"); print(resA$U)
cat("D_A =\n"); print(resA$D)
cat("V_A =\n"); print(resA$V)

cat("\n===== SVD réduite de B =====\n")
cat("rang(B) =", resB$rank, "\n")
cat("valeurs singulières(B) =\n"); print(resB$singular_values)
cat("Erreur Frobenius ||B - UDV^T|| = ", resB$recons_error_fro, "\n\n")

cat("U_B =\n"); print(resB$U)
cat("D_B =\n"); print(resB$D)
cat("V_B =\n"); print(resB$V)
