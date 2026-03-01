# nbre de clients observes par jour
obs <- c(Lundi = 30, Mardi = 14, Mercredi = 34, Jeudi = 45, Vendredi = 57, Samedi = 20)

# distribution en pourcentage
p <- c(Lundi = 0.10, Mardi = 0.10, Mercredi = 0.15, Jeudi = 0.20, Vendredi = 0.30, Samedi = 0.15)

# taille totale par semaine
n <- sum(obs)

exp <- n * p

cat("Taille n =", n, "\n\n")
cat("Effectifs observés :\n"); print(obs); cat("\n")
cat("Probabilités théoriques :\n"); print(p); cat("\n")
cat("Effectifs attendus :\n"); print(exp); cat("\n")

# test du chi-deux dadéquation
test <- chisq.test(x = obs, p = p, rescale.p = TRUE)

cat("\n--- Résultat du test Chi-deux ---\n")
print(test)

# decision au seuil alpha = 5%
alpha <- 0.05
cat("\nSeuil alpha =", alpha, "\n")
if (test$p.value < alpha) {
  cat("Décision : on REJETTE l'hypothèse d'adéquation (la distribution théorique ne colle pas aux données).\n")
  cat("Conclusion : l'étude de marché n'est PAS 'bonne' par rapport aux observations.\n")
} else {
  cat("Décision : on NE REJETTE PAS l'hypothèse d'adéquation (cohérent avec la distribution théorique).\n")
  cat("Conclusion : l'étude de marché est compatible avec les observations.\n")
}