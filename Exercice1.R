# =========================
# EXO 1 - Moindres carrés 
# =========================

# Données (X en centaines d'euros)
t <- c(0, 1, 2, 3, 4, 5, 6)
X <- c(90, 73.8, 60, 49.5, 40.5, 33, 27)

# 1) Transformation
Z <- log(X)

# 2) Calcul des coefficients a et b (moindres carrés)
n <- length(t)
t_bar <- sum(t) / n
Z_bar <- sum(Z) / n

S_tt <- sum( (t - t_bar)^2 )
S_tZ <- sum( (t - t_bar) * (Z - Z_bar) )

a <- S_tZ / S_tt
b <- Z_bar - a * t_bar

cat("Droite d'ajustement sur Z = ln(X):\n")
cat("Z(t) = a*t + b avec\n")
cat("a =", a, "\n")
cat("b =", b, "\n\n")

# 3) Modèle exponentiel X(t) = A * exp(a*t), avec A = exp(b)

A <- exp(b)
cat("Modele exponentiel:\n")
cat("X(t) = A * exp(a*t)\n")
cat("A =", A, "\n\n")

# 4) Estimation à 10 ans
t_pred <- 10
X10 <- A * exp(a * t_pred)

cat("Estimation X(10) (en centaines d'euros) =", X10, "\n")
cat("Donc en euros ~", 100 * X10, "€\n\n")

# 5) Graphique : points (X) + courbe ajustée
t_grid <- seq(0, 10, by = 0.1)
X_fit <- A * exp(a * t_grid)

plot(t, X,
     pch = 19,
     xlab = "t (annees)",
     ylab = "X(t) (centaines d'euros)",
     main = "Ajustement exponentiel par moindres carres (via ln(X))")
lines(t_grid, X_fit, lwd = 2)
abline(v = 10, lty = 2)

