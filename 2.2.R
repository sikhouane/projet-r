############################# EXERCICE 2.2.1

f <- function(x, y) {
  exp(x*y) - exp(x) + 2
}

grad_f <- function(x, y) {
  gx <- y*exp(x*y) - exp(x)
  gy <- x*exp(x*y)
  c(gx, gy)
}

hess_f <- function(x, y) {
  fxx <- y*y*exp(x*y) - exp(x)
  fyy <- x*x*exp(x*y)
  fxy <- (1 + x*y)*exp(x*y)
  
  matrix(c(fxx, fxy,
           fxy, fyy), nrow=2, byrow=TRUE)
}

x0 <- 0
y0 <- 1

grad_f(x0, y0)
H <- hess_f(x0, y0)

# test hessienne (det)
detH <- H[1,1]*H[2,2] - H[1,2]^2


############################# EXERCICE 2.2.2

g <- function(x, y, z) {
  (x + z^2) * exp(x*(y^2 + z^2 + 1))
}

grad_g <- function(x, y, z) {
  A <- y^2 + z^2 + 1
  ex <- exp(A*x)
  gx <- ex * (1 + A*(x + z^2))
  gy <- 2*x*y*(x + z^2)*ex
  gz <- 2*z*ex * (1 + x*(x + z^2))
  c(gx, gy, gz)
}

x0 <- -1
y0 <- 0
z0 <- 0
grad_g(x0, y0, z0)
g(x0, y0, z0)


############################EXERCICE 2.2.3

h <- function(x,y,z) {
  log(x*y*z) - log(x)*log(y)*log(z)
}

grad_h <- function(x,y,z) {
  lx <- log(x)
  ly <- log(y)
  lz <- log(z)
  gx <- (1 - ly*lz)/x
  gy <- (1 - lx*lz)/y
  gz <- (1 - lx*ly)/z

  c(gx, gy, gz)
}

p1 <- c(exp(1), exp(1), exp(1))
p2 <- c(exp(-1), exp(-1), exp(-1))
grad_h(p1[1], p1[2], p1[3])
grad_h(p2[1], p2[2], p2[3])
h(p1[1], p1[2], p1[3])
h(p2[1], p2[2], p2[3])


hess_h <- function(x,y,z) {
  
  lx <- log(x)
  ly <- log(y)
  lz <- log(z)
  
  hxx <- (ly*lz - 1)/(x^2)
  hyy <- (lx*lz - 1)/(y^2)
  hzz <- (lx*ly - 1)/(z^2)
  
  hxy <- -(lz)/(x*y)
  hxz <- -(ly)/(x*z)
  hyz <- -(lx)/(y*z)
  
  matrix(c(hxx,hxy,hxz,
           hxy,hyy,hyz,
           hxz,hyz,hzz),
         3,3,byrow=TRUE)
}

H1 <- hess_h(p1[1], p1[2], p1[3])
H2 <- hess_h(p2[1], p2[2], p2[3])