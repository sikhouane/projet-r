# ============================================================
#  EXERCICE 3.2 — Optimisation sous contraintes (Lagrange)
# ============================================================
# Methode (points 1, 2, 3) :
#
#  POINT 1 — Lagrangien et condition necessaire
#    L(X,λ) = f(X) - λ1*g1(X) - ... - λp*gp(X)
#    => ∇L(X0, λ0) = 0
#
#  POINT 2 — Critere via D²L seule
#    Si D²L(X0,λ0) definie POSITIVE  => minimum local
#    Si D²L(X0,λ0) definie NEGATIVE  => maximum local
#    (Sinon : indetermine, on passe au point 3)
#
#  POINT 3 — Critere via hessienne bordee H(X0)
#    H(X0) = [ 0         Dg(X0)        ]  ordre n+p
#            [ Dg(X0)^T  D²L(X0,λ0)   ]
#    - Les (n-p) derniers mineurs principaux alternent en signe
#      ET det(H) meme signe que (-1)^n  =>  MAXIMUM local
#    - Les (n-p) derniers mineurs principaux alternent en signe
#      ET det(H) meme signe que (-1)^p  =>  MINIMUM local
# ============================================================

cat("========================================================\n")
cat("  OPTIMISATION SOUS CONTRAINTES - Methode de Lagrange\n")
cat("  (Points 1, 2, 3)\n")
cat("========================================================\n\n")

# ============================================================
# FONCTION CENTRALE — applique EXACTEMENT les points 2 et 3
# pour conclure sur la nature du point critique
# ============================================================
analyser_point <- function(label, X0, lambda0, Dg, D2L) {

  n <- ncol(Dg)  # nombre de variables
  p <- nrow(Dg)  # nombre de contraintes

  cat(sprintf("\n=== %s ===\n", label))
  cat(sprintf("X0 = (%s)\n", paste(round(X0,4), collapse=", ")))
  cat(sprintf("lambda0 = (%s)\n", paste(round(lambda0,4), collapse=", ")))

  # --- Verification qualification des contraintes
  rang <- qr(Dg)$rank
  cat(sprintf("Dg(X0) : rang = %d  (p = %d)  => contraintes %s\n",
              rang, p, ifelse(rang==p,"QUALIFIEES","NON qualifiees")))

  # ----------------------------------------------------------
  # POINT 2 du cours : critere via D²L seule
  # ----------------------------------------------------------
  cat("\n-- POINT 2 : D²L(X0, lambda0) --\n")
  print(round(D2L, 6))

  vp <- eigen(D2L, symmetric=TRUE)$values
  cat(sprintf("Valeurs propres de D²L : %s\n",
              paste(round(vp,6), collapse=", ")))

  conclusion_D2L <- "indeterminee"
  if (all(vp > 1e-10)) {
    conclusion_D2L <- "MINIMUM local (D²L definie positive)"
    cat("=> D²L definie POSITIVE => MINIMUM local\n")
  } else if (all(vp < -1e-10)) {
    conclusion_D2L <- "MAXIMUM local (D²L definie negative)"
    cat("=> D²L definie NEGATIVE => MAXIMUM local\n")
  } else {
    cat("=> D²L indefinie : le critere du point 2 ne conclut pas.\n")
    cat("   => On passe au POINT 3 (hessienne bordee).\n")
  }

  # ----------------------------------------------------------
  # POINT 3 du cours : critere via hessienne bordee H(X0)
  # ----------------------------------------------------------
  cat("\n-- POINT 3 : Hessienne bordee H(X0) --\n")

  taille <- n + p
  H <- matrix(0, taille, taille)
  H[1:p, (p+1):taille]          <- Dg
  H[(p+1):taille, 1:p]          <- t(Dg)
  H[(p+1):taille, (p+1):taille] <- D2L

  cat(sprintf("H(X0) de taille (%d+%d) = %d x %d :\n", n, p, taille, taille))
  print(round(H, 6))

  # (n-p) mineurs principaux bordes
  # Mineur k = det de la sous-matrice de taille (2p+k) x (2p+k)
  nb_mineurs <- n - p
  cat(sprintf("Nombre de mineurs a examiner : n-p = %d-%d = %d\n", n, p, nb_mineurs))

  if (nb_mineurs == 0) {
    cat("n = p : point isole, pas de mineur additionnel.\n")
    cat(sprintf("\n>>> CONCLUSION FINALE : %s\n", conclusion_D2L))
    return(invisible(NULL))
  }

  mineurs <- numeric(nb_mineurs)
  for (k in 1:nb_mineurs) {
    sz         <- 2*p + k
    mineurs[k] <- det(H[1:sz, 1:sz])
    cat(sprintf("  Mineur %d (sous-matrice %dx%d) : det = %+.6f\n",
                k, sz, sz, mineurs[k]))
  }

  # Signes calcules automatiquement depuis les resultats
  signes_obs <- sign(mineurs)

  # Signes attendus selon le cours :
  # MAX : det(H) meme signe que (-1)^n, mineurs alternent
  # MIN : det(H) meme signe que (-1)^p, mineurs alternent
  # L'alternance pour MAX : signe du k-eme mineur = (-1)^(n) * (-1)^(k-1)
  # L'alternance pour MIN : signe du k-eme mineur = (-1)^(p) * (-1)^(k-1)
  signe_max_attendu <- (-1)^n * (-1)^(0:(nb_mineurs-1))
  signe_min_attendu <- (-1)^p * (-1)^(0:(nb_mineurs-1))

  est_max <- all(signes_obs == sign(signe_max_attendu))
  est_min <- all(signes_obs == sign(signe_min_attendu))

  cat(sprintf("\n  Signes observes     : %s\n",
              paste(ifelse(signes_obs>0,"+","-"), collapse=" ")))
  cat(sprintf("  Signes attendus MAX : %s  [(-1)^n=%+.0f]\n",
              paste(ifelse(sign(signe_max_attendu)>0,"+","-"), collapse=" "),
              (-1)^n))
  cat(sprintf("  Signes attendus MIN : %s  [(-1)^p=%+.0f]\n",
              paste(ifelse(sign(signe_min_attendu)>0,"+","-"), collapse=" "),
              (-1)^p))

  # Conclusion automatique basee sur les vrais resultats calcules
  cat("\n>>> CONCLUSION (Point 3 du cours) : ")
  if (est_max && !est_min) {
    cat("MAXIMUM LOCAL sous les contraintes\n")
    cat(sprintf("    Mineurs alternent en signe ET det(H)=%+.4f meme signe que (-1)^n=%+.0f\n",
                mineurs[nb_mineurs], (-1)^n))
  } else if (est_min && !est_max) {
    cat("MINIMUM LOCAL sous les contraintes\n")
    cat(sprintf("    Mineurs alternent en signe ET det(H)=%+.4f meme signe que (-1)^p=%+.0f\n",
                mineurs[nb_mineurs], (-1)^p))
  } else {
    cat("Critere des mineurs non concluant (point selle ou cas limite)\n")
  }

  # Coherence entre point 2 et point 3
  if (conclusion_D2L != "indeterminee") {
    nat2 <- if(grepl("MINIMUM", conclusion_D2L)) "MINIMUM" else "MAXIMUM"
    nat3 <- if(est_min) "MINIMUM" else if(est_max) "MAXIMUM" else "indetermine"
    if (nat2 == nat3) {
      cat("    [Coherence Point 2 / Point 3 : OK]\n")
    } else {
      cat(sprintf("    [ATTENTION : Point 2 dit %s, Point 3 dit %s]\n", nat2, nat3))
    }
  }
}


# ============================================================
# EXERCICE 1
# f(x,y) = xy,  g(x,y) = x - y + 6 = 0
# n = 2,  p = 1
# ============================================================
cat("\n##########################################################\n")
cat("## EXERCICE 1 : f(x,y) = xy,  g = x - y + 6 = 0\n")
cat("##########################################################\n")

# --- Lagrangien : L(x,y,λ) = xy - λ(x - y + 6)
# --- ∇L = 0 :
#   dL/dx = y - λ = 0  =>  λ = y
#   dL/dy = x + λ = 0  =>  λ = -x
#   => y = -x
#   contrainte x - (-x) + 6 = 0 => 2x = -6 => x = -3

x0_1 <- -3;  y0_1 <- 3;  l0_1 <- 3

cat(sprintf("Resolution ∇L=0 :\n"))
cat(sprintf("  λ=y et λ=-x => y=-x\n"))
cat(sprintf("  contrainte : x-(-x)+6=0 => x=-3, y=3, λ=3\n"))
cat(sprintf("Point critique X0 = (%.4f, %.4f)\n", x0_1, y0_1))
cat(sprintf("Valeur f(X0) = %.4f\n", x0_1 * y0_1))
cat(sprintf("Verification g(X0) = %.6f (doit etre 0)\n", x0_1 - y0_1 + 6))

# Jacobienne Dg = [dg/dx, dg/dy] = [1, -1]
Dg_1 <- matrix(c(1, -1), nrow=1)

# D²L :
# L = xy - λ(x-y+6)
# d²L/dx²   = 0
# d²L/dy²   = 0
# d²L/dx dy = 1  (car d/dy(y) = 1)
D2L_1 <- matrix(c(0, 1,
                   1, 0), nrow=2, byrow=TRUE)

analyser_point("Exo 1 -- X0=(-3, 3)", c(x0_1, y0_1), l0_1, Dg_1, D2L_1)

# Trace sur [-8,4] x [-4,12] pour voir le point et la contrainte
x_g <- seq(-8, 4,  length.out=300)
y_g <- seq(-4, 12, length.out=300)
Z1  <- outer(x_g, y_g, function(x, y) x * y)
xc1 <- seq(-8, 4, length.out=200)
yc1 <- xc1 + 6  # contrainte : y = x + 6

png("~/Documents/develop/esgi/5a/t2/maths/outputs/ex1_trace.png", width=900, height=700, res=110)
par(mar=c(4,4,3.5,1))
image(x_g, y_g, Z1, col=hcl.colors(60,"YlOrRd"),
      xlab="x", ylab="y",
      main="Exo 1 : f(x,y)=xy  |  contrainte x-y+6=0")
contour(x_g, y_g, Z1, add=TRUE, col="grey40", lwd=0.8, nlevels=20)
lines(xc1, yc1, col="blue", lwd=2.5)
points(x0_1, y0_1, pch=21, bg="green4", col="black", cex=2.5, lwd=2)
text(x0_1+0.3, y0_1+0.6,
     sprintf("MINIMUM\n(%.0f, %.0f)\nf=%.0f", x0_1, y0_1, x0_1*y0_1),
     col="green4", cex=0.85, font=2)
legend("topright",
       legend=c("f(x,y)=xy (niveaux)","contrainte y=x+6","Minimum local"),
       col=c("grey40","blue","green4"),
       lty=c(1,1,NA), pch=c(NA,NA,21), pt.bg=c(NA,NA,"green4"),
       lwd=c(1,2.5,2), bg="white", cex=0.85)
dev.off()
cat("-> ex1_trace.png genere\n")


# ============================================================
# EXERCICE 2
# f(x,y,z) = x²+y²+z²
# g1 = x+y+z-1 = 0,  g2 = z = 0
# n = 3,  p = 2
# ============================================================
cat("\n##########################################################\n")
cat("## EXERCICE 2 : f=x2+y2+z2\n")
cat("##   g1=x+y+z-1=0,  g2=z=0\n")
cat("##########################################################\n")

# --- Lagrangien : L = x²+y²+z² - λ1(x+y+z-1) - λ2*z
# --- ∇L = 0 :
#   dL/dx = 2x - λ1 = 0          => λ1 = 2x
#   dL/dy = 2y - λ1 = 0          => λ1 = 2y  =>  x = y
#   dL/dz = 2z - λ1 - λ2 = 0     => λ2 = 2z - λ1
#   g2 : z = 0
#   g1 : x + y + 0 = 1 et x=y  =>  x = y = 1/2

x0_2 <- 1/2;  y0_2 <- 1/2;  z0_2 <- 0
l1_2 <- 1;    l2_2 <- -1

cat("Resolution ∇L=0 :\n")
cat("  g2 => z=0 ; x=y ; g1 => 2x=1 => x=y=1/2\n")
cat(sprintf("  λ1 = 2x = 1,  λ2 = 2z-λ1 = -1\n"))
cat(sprintf("Point critique X0 = (%.4f, %.4f, %.4f)\n", x0_2, y0_2, z0_2))
cat(sprintf("Valeur f(X0) = %.4f\n", x0_2^2+y0_2^2+z0_2^2))

Dg_2 <- matrix(c(1,1,1,
                  0,0,1), nrow=2, byrow=TRUE)
# L = x²+y²+z² - λ1(x+y+z-1) - λ2*z
# d²L/dx² = 2,  d²L/dy² = 2,  d²L/dz² = 2,  croises = 0
D2L_2 <- 2 * diag(3)

analyser_point("Exo 2 -- X0=(1/2, 1/2, 0)",
               c(x0_2, y0_2, z0_2), c(l1_2, l2_2), Dg_2, D2L_2)


# ============================================================
# EXERCICE 3
# f(x,y,z) = x²+y²+z²
# g1 = x²+y²+2z-6 = 0,  g2 = x-y-z = 0
# n = 3,  p = 2
# ============================================================
cat("\n##########################################################\n")
cat("## EXERCICE 3 : f=x2+y2+z2\n")
cat("##   g1=x2+y2+2z-6=0,  g2=x-y-z=0\n")
cat("##########################################################\n")

# --- Lagrangien : L = x²+y²+z² - λ1(x²+y²+2z-6) - λ2(x-y-z)
# --- ∇L = 0 :
#   dL/dx = 2x - 2λ1*x - λ2 = 0   => x(2-2λ1) = λ2       (i)
#   dL/dy = 2y - 2λ1*y + λ2 = 0   => y(2-2λ1) = -λ2      (ii)
#   dL/dz = 2z - 2λ1 + λ2 = 0     => 2z = 2λ1 - λ2       (iii)
#
# De (i) et (ii) : x(2-2λ1) = -y(2-2λ1)
# CAS A : 2-2λ1 ≠ 0  =>  x = -y
#   (v) z = x-y = 2x
#   g1 : x²+x²+4x-6=0  =>  2x²+4x-6=0  =>  (x+3)(x-1)=0
#   => (x,y,z) = (-3,3,-6) ou (1,-1,2)
# CAS B : λ1 = 1, λ2 = 0
#   (iii) z=1 ; (v) y=x-1 ; g1 : 2x²-2x-3=0 => x=(1±√7)/2

sqrt7 <- sqrt(7)
pts3 <- list(
  list(x=-3,          y= 3,          z=-6, l1=4.5, l2= 21),
  list(x= 1,          y=-1,          z= 2, l1=1.5, l2= -1),
  list(x=(1+sqrt7)/2, y=(1+sqrt7)/2-1, z=1, l1=1,  l2=  0),
  list(x=(1-sqrt7)/2, y=(1-sqrt7)/2-1, z=1, l1=1,  l2=  0)
)

cat("Points critiques trouves :\n")
for (i in seq_along(pts3)) {
  pt <- pts3[[i]]
  fv <- pt$x^2 + pt$y^2 + pt$z^2
  cat(sprintf("  P%d : (%.4f, %.4f, %.4f)  f=%.4f  λ=(%.4f,%.4f)\n",
              i, pt$x, pt$y, pt$z, fv, pt$l1, pt$l2))
}

for (pt in pts3) {
  fv <- pt$x^2 + pt$y^2 + pt$z^2
  # Jacobienne : g1 non-lineaire => depend du point
  Dg_3 <- matrix(c(2*pt$x, 2*pt$y, 2,
                    1,      -1,    -1), nrow=2, byrow=TRUE)
  # D²L : d²L/dx² = 2-2λ1,  d²L/dy² = 2-2λ1,  d²L/dz² = 2
  # termes croises = 0 car L ne contient pas xy, xz, yz
  D2L_3 <- diag(c(2-2*pt$l1, 2-2*pt$l1, 2))

  analyser_point(
    sprintf("Exo 3 -- (%.3f, %.3f, %.3f) f=%.3f",
            pt$x, pt$y, pt$z, fv),
    c(pt$x, pt$y, pt$z), c(pt$l1, pt$l2), Dg_3, D2L_3)
}


# ============================================================
# EXERCICE 4
# f(x,y) = x²*y,  g = 2x²+y²-3 = 0
# n = 2,  p = 1
# ============================================================
cat("\n##########################################################\n")
cat("## EXERCICE 4 : f(x,y)=x2*y,  g=2x2+y2-3=0\n")
cat("##########################################################\n")

# --- Lagrangien : L = x²y - λ(2x²+y²-3)
# --- ∇L = 0 :
#   dL/dx = 2xy - 4λx = 0  =>  2x(y-2λ) = 0
#   dL/dy = x² - 2λy = 0
#   g     : 2x² + y² = 3
#
# CAS A : x = 0  =>  g: y²=3 => y=±√3,  λ=0  (de dL/dy: 0=0)
# CAS B : x≠0  =>  y = 2λ
#   dL/dy : x² = 2λy = 4λ² => λ=±x/2
#   g     : 2x² + 4λ² = 3 => 2x²+x² = 3 => x²=1 => x=±1
#   y = 2λ = ±1

pts4 <- list(
  list(x= 0,  y= sqrt(3), l= 0),
  list(x= 0,  y=-sqrt(3), l= 0),
  list(x= 1,  y= 1,       l= 0.5),
  list(x=-1,  y= 1,       l= 0.5),
  list(x= 1,  y=-1,       l=-0.5),
  list(x=-1,  y=-1,       l=-0.5)
)

cat("Points critiques trouves :\n")
for (i in seq_along(pts4)) {
  pt <- pts4[[i]]
  fv <- pt$x^2 * pt$y
  cat(sprintf("  P%d : (%.4f, %.4f)  f=%.4f  λ=%.4f\n",
              i, pt$x, pt$y, fv, pt$l))
}

for (pt in pts4) {
  fv <- pt$x^2 * pt$y
  Dg_4  <- matrix(c(4*pt$x, 2*pt$y), nrow=1)
  # L = x²y - λ(2x²+y²-3)
  # d²L/dx²   = 2y - 4λ
  # d²L/dy²   = -2λ
  # d²L/dx dy = 2x
  D2L_4 <- matrix(c(2*pt$y - 4*pt$l,  2*pt$x,
                     2*pt$x,           -2*pt$l), nrow=2, byrow=TRUE)

  analyser_point(
    sprintf("Exo 4 -- (%.3f, %.3f) f=%.3f", pt$x, pt$y, fv),
    c(pt$x, pt$y), pt$l, Dg_4, D2L_4)
}

# --- Trace sur [-2,2]²
x4 <- seq(-2, 2, length.out=300)
y4 <- seq(-2, 2, length.out=300)
Z4 <- outer(x4, y4, function(x,y) x^2*y)
theta4 <- seq(0, 2*pi, length.out=300)
xc4 <- sqrt(3/2)*cos(theta4); yc4 <- sqrt(3)*sin(theta4)

pts4_plot <- data.frame(
  x  = c(0,      0,       1, -1,  1, -1),
  y  = c(sqrt(3),-sqrt(3),1,  1, -1, -1),
  f  = c(0,      0,       1,  1, -1, -1),
  nat= c("f=0","f=0","MAX","MAX","MIN","MIN")
)
png("~/Documents/develop/esgi/5a/t2/maths/outputs/ex4_trace.png", width=900, height=700, res=110)
par(mar=c(4,4,3.5,1))
image(x4, y4, Z4, col=hcl.colors(60,"RdYlBu"),
      xlab="x", ylab="y",
      main="Exo 4 : f(x,y)=x2y  |  contrainte 2x2+y2=3")
contour(x4, y4, Z4, add=TRUE, col="grey40", lwd=0.7, nlevels=15)
lines(xc4, yc4, col="navy", lwd=2.5)
for (i in 1:nrow(pts4_plot)) {
  ci <- if(pts4_plot$f[i]>0.5)"red" else if(pts4_plot$f[i]< -0.5)"green4" else "orange"
  points(pts4_plot$x[i], pts4_plot$y[i], pch=21,
         bg=ci, col="black", cex=2.2, lwd=1.5)
  text(pts4_plot$x[i]+0.1, pts4_plot$y[i]+0.13,
       sprintf("%s\nf=%.0f", pts4_plot$nat[i], pts4_plot$f[i]),
       cex=0.75, font=2, col=ci)
}
legend("topright",
       legend=c("f=x2y","contrainte 2x2+y2=3",
                "Maximum local (f=1)","Minimum local (f=-1)","f=0 (neutre)"),
       col=c("grey40","navy","red","green4","orange"),
       lty=c(1,1,NA,NA,NA), pch=c(NA,NA,21,21,21),
       pt.bg=c(NA,NA,"red","green4","orange"),
       lwd=c(1,2.5,2,2,2), bg="white", cex=0.8)
dev.off()
cat("-> ex4_trace.png genere\n")


# ============================================================
# EXERCICE 5
# f = (x²y²)/(x²+y²) + 2x²+2y²+4xy+5x+5y+20
# g = 2x²+2y²+4xy - 8 = 3  <=>  2(x+y)² = 11
# n = 2,  p = 1
# ============================================================
cat("\n##########################################################\n")
cat("## EXERCICE 5 : f complexe,  g=2(x+y)^2-11=0\n")
cat("##########################################################\n")

f5      <- function(x,y){ r2<-x^2+y^2; if(r2<1e-14)return(NA);
                          (x^2*y^2)/r2+2*x^2+2*y^2+4*x*y+5*x+5*y+20 }
g5      <- function(x,y) 2*x^2+2*y^2+4*x*y - 11
df5_dx  <- function(x,y){ r2<-x^2+y^2; 2*x*y^4/r2^2+4*x+4*y+5 }
df5_dy  <- function(x,y){ r2<-x^2+y^2; 2*y*x^4/r2^2+4*y+4*x+5 }
dg5_dx  <- function(x,y) 4*x+4*y
dg5_dy  <- function(x,y) 4*y+4*x

# Jacobien numerique (diff. finies centrales)
jac_num <- function(F, v, eps=1e-7){
  n<-length(v); J<-matrix(0,n,n)
  for(j in 1:n){ vp<-vm<-v; vp[j]<-vp[j]+eps; vm[j]<-vm[j]-eps
    J[,j]<-(F(vp)-F(vm))/(2*eps) }; J
}
# Systeme KKT : [df/dx - λdg/dx,  df/dy - λdg/dy,  g] = 0
F5_kkt <- function(v){
  x<-v[1]; y<-v[2]; lam<-v[3]
  c(df5_dx(x,y)-lam*dg5_dx(x,y),
    df5_dy(x,y)-lam*dg5_dy(x,y), g5(x,y))
}
# Newton from scratch
newton5 <- function(v0, tol=1e-12, maxit=200){
  v<-v0
  for(i in 1:maxit){
    Fv<-tryCatch(F5_kkt(v),error=function(e)rep(NA,3))
    if(any(is.na(Fv))||max(abs(Fv))<tol) break
    J<-jac_num(F5_kkt,v)
    if(abs(det(J))<1e-15) break
    v<-v-solve(J)%*%Fv }
  list(sol=v, res=tryCatch(F5_kkt(v),error=function(e)rep(NA,3)))
}

# Balayage sur les deux droites contrainte x+y = ±√5.5
s5 <- sqrt(5.5); sols5 <- list()
for(sg in c(1,-1)){
  for(t in seq(-4,4,by=0.15)){
    x_i<-(sg*s5+t)/2; y_i<-(sg*s5-t)/2
    if(is.na(f5(x_i,y_i))) next
    gy<-dg5_dy(x_i,y_i)
    l_i<-if(abs(gy)>1e-8) df5_dy(x_i,y_i)/gy else 0
    res<-tryCatch(newton5(c(x_i,y_i,l_i)),error=function(e)NULL)
    if(is.null(res)||any(is.na(res$res))) next
    if(max(abs(res$res))>1e-8) next
    sol<-round(res$sol,6)
    dup<-any(sapply(sols5,function(s) max(abs(s-sol))<1e-4))
    if(!dup) sols5[[length(sols5)+1]]<-sol
  }
}
cat(sprintf("Points critiques trouves : %d\n", length(sols5)))

# Derivees secondes analytiques de f (pour D²L)
d2f5_dxx <- function(x,y){ r2<-x^2+y^2; 4+(2*y^4*r2^2-8*x^2*y^4*r2)/r2^4 }
d2f5_dyy <- function(x,y){ r2<-x^2+y^2; 4+(2*x^4*r2^2-8*y^2*x^4*r2)/r2^4 }
d2f5_dxy <- function(x,y){ r2<-x^2+y^2; 4+8*x*y^3/r2^2-8*x*y^5/r2^3 }
# d²g/dx² = 4,  d²g/dy² = 4,  d²g/dxdy = 4

for(sol in sols5){
  xs<-sol[1]; ys<-sol[2]; ls<-sol[3]; fs<-f5(xs,ys)
  cat(sprintf("\nPoint (%.4f, %.4f)  lambda=%.4f  f=%.6f\n",xs,ys,ls,fs))
  Dg_5  <- matrix(c(dg5_dx(xs,ys), dg5_dy(xs,ys)), nrow=1)
  # D²L = D²f - λ D²g
  H11 <- d2f5_dxx(xs,ys) - ls*4
  H22 <- d2f5_dyy(xs,ys) - ls*4
  H12 <- d2f5_dxy(xs,ys) - ls*4
  D2L_5 <- matrix(c(H11,H12,H12,H22), 2, 2)
  analyser_point(sprintf("Exo 5 -- (%.3f,%.3f) f=%.4f",xs,ys,fs),
                 c(xs,ys), ls, Dg_5, D2L_5)
}
if(length(sols5)>0){
  fv5<-sapply(sols5,function(s) f5(s[1],s[2]))
  cat(sprintf("\nMinimum Exo 5 : f = %.6f\n", min(fv5)))
  cat(sprintf("Maximum Exo 5 : f = %.6f\n",  max(fv5)))
}

# --- Trace sur [-3,3]²
x5<-seq(-3,3,length.out=300); y5<-seq(-3,3,length.out=300)
Z5<-outer(x5,y5,function(x,y){
  r2<-x^2+y^2
  ifelse(r2<1e-10,NA,(x^2*y^2)/r2+2*x^2+2*y^2+4*x*y+5*x+5*y+20)})
xc5<-seq(-3,3,length.out=400)
yc5p<-s5-xc5; yc5n<- -s5-xc5
png("~/Documents/develop/esgi/5a/t2/maths/outputs/ex5_trace.png", width=900, height=700, res=110)
par(mar=c(4,4,3.5,1))
image(x5,y5,Z5,col=hcl.colors(60,"Spectral",rev=TRUE),
      xlab="x",ylab="y",
      main="Exo 5 : f(x,y)  |  contrainte 2(x+y)2=11")
contour(x5,y5,Z5,add=TRUE,col="grey30",lwd=0.6,nlevels=25)
ok_p<-yc5p>=-3&yc5p<=3; ok_n<-yc5n>=-3&yc5n<=3
if(any(ok_p)) lines(xc5[ok_p],yc5p[ok_p],col="navy",lwd=2.5)
if(any(ok_n)) lines(xc5[ok_n],yc5n[ok_n],col="navy",lwd=2.5,lty=2)
if(length(sols5)>0){
  fv5<-sapply(sols5,function(s) f5(s[1],s[2]))
  for(i in seq_along(sols5)){
    xs<-sols5[[i]][1]; ys<-sols5[[i]][2]
    if(xs< -3||xs>3||ys< -3||ys>3) next
    ci<-if(fv5[i]==min(fv5))"green4" else "red"
    ni<-if(fv5[i]==min(fv5))"MIN" else "MAX"
    points(xs,ys,pch=21,bg=ci,col="black",cex=2.2,lwd=1.5)
    text(xs+0.15,ys+0.15,sprintf("%s\nf=%.2f",ni,fv5[i]),
         cex=0.75,font=2,col=ci)
  }
}
legend("topright",
       legend=c("f(x,y)","x+y=+sqrt(5.5)","x+y=-sqrt(5.5)","Minimum","Maximum"),
       col=c("grey30","navy","navy","green4","red"),
       lty=c(1,1,2,NA,NA),pch=c(NA,NA,NA,21,21),
       pt.bg=c(NA,NA,NA,"green4","red"),
       lwd=c(1,2.5,2.5,2,2),bg="white",cex=0.8)
dev.off()
cat("-> ex5_trace.png genere\n")

# ============================================================
# TABLEAU RECAPITULATIF
# ============================================================
cat("\n========================================================\n")
cat("  TABLEAU RECAPITULATIF — Exercice 3.2\n")
cat("========================================================\n")
cat(sprintf("%-4s %-22s %7s %-25s\n","Exo","X0","f(X0)","Conclusion"))
cat(strrep("-",62),"\n")
rows <- list(
  list("1", "(-3, 3)",           -9,   "MINIMUM local"),
  list("2", "(1/2, 1/2, 0)",      0.5, "MINIMUM local"),
  list("3", "(-3, 3, -6)",       54,   "MAXIMUM local"),
  list("3", "(1, -1, 2)",         6,   "MINIMUM local"),
  list("3", "((1+v7)/2,...)",    NA,   "voir console"),
  list("3", "((1-v7)/2,...)",    NA,   "voir console"),
  list("4", "(1,1) et (-1,1)",   1,   "MAXIMUM local"),
  list("4", "(1,-1) et (-1,-1)",-1,   "MINIMUM local"),
  list("4", "(0, +/-sqrt(3))",   0,   "f=0 neutre"),
  list("5", "Newton numerique", NA,   "voir console")
)
for(r in rows) cat(sprintf("%-4s %-22s %7s %-25s\n",
  r[[1]], r[[2]],
  ifelse(is.na(r[[3]]),"—",sprintf("%.4g",r[[3]])),
  r[[4]]))

cat("\nTracés generes : ex1_trace.png | ex4_trace.png | ex5_trace.png\n")
cat("Script termine.\n")
