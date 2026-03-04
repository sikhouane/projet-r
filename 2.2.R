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


x <- seq(0, 2, length.out = 80)
y <- seq(0, 2, length.out = 80)

Z <- matrix(0, nrow = length(x), ncol = length(y))
for (i in 1:length(x)) {
  for (j in 1:length(y)) {
    Z[i, j] <- f(x[i], y[j])
  }
}

#3D
persp(x, y, Z,
      theta = 35, phi = 25,
      xlab = "x", ylab = "y", zlab = "f(x,y)",
      main = "Surface de f(x,y) sur [0,2]^2")

#courbes de niveau + point critique
contour(x, y, Z,
        xlab = "x", ylab = "y",
        main = "Courbes de niveau de f(x,y) sur [0,2]^2")
points(0, 1, pch = 19)
text(0, 1, "(0,1)", pos = 4)

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

####coupes 1D de g
xg <- seq(-4, 2, length.out = 500)

g_x_00 <- rep(0, length(xg))
g_x_01 <- rep(0, length(xg))

for (i in 1:length(xg)) {
  g_x_00[i] <- g(xg[i], 0, 0)
  g_x_01[i] <- g(xg[i], 0, 1)
}

plot(xg, g_x_00, type = "l",
     xlab = "x", ylab = "g(x,0,z)",
     main = "Coupes : g(x,0,0) (trait plein) et g(x,0,1) (pointillé)")
lines(xg, g_x_01, lty = 2)

# Marquer le minimum global (-1,0,0)
points(-1, g(-1,0,0), pch = 19)
text(-1, g(-1,0,0), "(-1,0,0)", pos = 4)
abline(v = -1, lty = 3)

####contour de g(x,y,0)
x <- seq(-3, 1, length.out = 90)
y <- seq(-2, 2, length.out = 90)
Zg <- matrix(0, nrow = length(x), ncol = length(y))

for (i in 1:length(x)) {
  for (j in 1:length(y)) {
    Zg[i, j] <- g(x[i], y[j], 0)
  }
}

contour(x, y, Zg,
        xlab = "x", ylab = "y",
        main = "Courbes de niveau de g(x,y,0)")
points(-1, 0, pch = 19)
text(-1, 0, "(-1,0,0)", pos = 4)

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

####coupe h(t,t,t)
t <- seq(0.05, 6, length.out = 600)
ht <- rep(0, length(t))

for (i in 1:length(t)) {
  ht[i] <- h(t[i], t[i], t[i])
}

plot(t, ht, type = "l",
     xlab = "t", ylab = "h(t,t,t)",
     main = "Coupe : h(t,t,t)")

abline(v = exp(1), lty = 2)
abline(v = exp(-1), lty = 2)

points(exp(1), h(exp(1),exp(1),exp(1)), pch = 19)
text(exp(1), h(exp(1),exp(1),exp(1)), "(e,e,e)", pos = 4)

points(exp(-1), h(exp(-1),exp(-1),exp(-1)), pch = 19)
text(exp(-1), h(exp(-1),exp(-1),exp(-1)), "(1/e,1/e,1/e)", pos = 4)

####contour h(x,y,z0) avec z0 = e
x <- seq(0.2, 5, length.out = 120)
y <- seq(0.2, 5, length.out = 120)
z0 <- exp(1)

Zh <- matrix(0, nrow = length(x), ncol = length(y))

for (i in 1:length(x)) {
  for (j in 1:length(y)) {
    Zh[i, j] <- h(x[i], y[j], z0)
  }
}

contour(x, y, Zh,
        xlab = "x", ylab = "y",
        main = "Courbes de niveau de h(x,y,e)")

# projection du point (e,e,e) sur ce plan z=e
points(exp(1), exp(1), pch = 19)
text(exp(1), exp(1), "(e,e,e)", pos = 4)
