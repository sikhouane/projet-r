# fonctions communes

#newton multi-dim
newton <- function(x0, grad_fun, hess_fun, maxit=50, tol=1e-10) {
  x <- x0
  for (it in 1:maxit) {
    g <- grad_fun(x)
    ng <- sqrt(dot(g,g))
    if (ng < tol) return(list(x=x, it=it, grad_norm=ng, converged=TRUE))
    H <- hess_fun(x)
    step <- gauss_solve(H, g) #H * step=g
    x <- x - step
  }
  list(x=x, it=maxit, grad_norm=sqrt(dot(grad_fun(x),grad_fun(x))), converged=FALSE)
}

#produit scalaire
#elimination gauss

# EXERCICE 2.2.1

f_val <- function(x, y) exp(x*y) - exp(x) + 2

f_grad <- function(u) {
  x <- u[1]; y <- u[2]
  c(y*exp(x*y) - exp(x),
    x*exp(x*y))
}

f_hess <- function(u) {
  x <- u[1]; y <- u[2]
  fxx <- y*y*exp(x*y) - exp(x)
  fyy <- x*x*exp(x*y)
  fxy <- (1 + x*y)*exp(x*y)
  matrix(c(fxx, fxy,
           fxy, fyy), nrow=2, byrow=TRUE)
}

#newton depuis un point de départ (proche de (0,1))
res_f <- newton(c(0.2, 0.8), f_grad, f_hess)

#classification via hessienne 2x2 (sans eigen)
classify_hess_2d <- function(H) {
  a <- H[1,1]; b <- H[1,2]; c <- H[2,1]; d <- H[2,2]
  detH <- a*d - b*c
  if (detH > 0 && a > 0) return("Minimum local")
  if (detH > 0 && a < 0) return("Maximum local")
  if (detH < 0) return("Point selle")
  "Test inconclusif"
}
class_f <- classify_hess_2d(f_hess(res_f$x))


# EXERCICE 2.2.2

# EXERCICE 2.2.3



