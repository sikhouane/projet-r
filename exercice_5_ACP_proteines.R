library(FactoMineR)
library(factoextra)
library(corrplot)

donnees1 <- read.table("~/Documents/develop/esgi/5a/t2/maths/donnees1.txt", 
                       header = TRUE, row.names = 1, sep = "\t")

actif <- donnees1[, 1:9]
EST   <- factor(donnees1$EST, levels = c(0, 1), labels = c("Ouest", "Est"))


# ==============================================================================
#  PARTIE 1 — premiere approche
# ==============================================================================


# Q1 — Matrice de corrélation ---------------------------------------------------

R_mat <- cor(actif)
round(R_mat, 2)

corrplot(R_mat,
         method      = "color",
         type        = "upper",
         addCoef.col = "black",
         number.cex  = 0.65,
         tl.col      = "black",
         tl.srt      = 45,
         title       = "Matrice de corrélation — variables alimentaires",
         mar         = c(0, 0, 2, 0))

# Q2 — Identifier les corrélations fortes et les quasi-nulles -------------------

# On extrait toutes les paires du triangle supérieur pour avoir un tableau propre
paires <- data.frame(V1 = character(), V2 = character(),
                     r  = numeric(), stringsAsFactors = FALSE)
p <- ncol(actif)
for (i in 1:(p - 1))
  for (j in (i + 1):p)
    paires <- rbind(paires, data.frame(
      V1 = colnames(actif)[i],
      V2 = colnames(actif)[j],
      r  = round(R_mat[i, j], 3),
      stringsAsFactors = FALSE))

# Seuils choisis : fort si |r| >= 0.5, quasi nul si |r| <= 0.15
SEUIL_FORT   <-  0.50
SEUIL_FAIBLE <-  0.15

pos <- paires[paires$r >=  SEUIL_FORT, ]
neg <- paires[paires$r <= -SEUIL_FORT, ]
nul <- paires[abs(paires$r) <= SEUIL_FAIBLE, ]

cat("Corrélations positives fortes :\n")
pos <- pos[order(-pos$r), ]
for (i in seq_len(nrow(pos)))
  cat(sprintf("  %-6s – %-6s : r = %+.3f\n", pos$V1[i], pos$V2[i], pos$r[i]))

cat("\nCorrélations négatives fortes :\n")
if (nrow(neg) == 0) cat("  (aucune)\n") else {
  neg <- neg[order(neg$r), ]
  for (i in seq_len(nrow(neg)))
    cat(sprintf("  %-6s – %-6s : r = %+.3f\n", neg$V1[i], neg$V2[i], neg$r[i]))
}

cat("\nVariables quasi décorrélées :\n")
if (nrow(nul) == 0) cat("  (aucune)\n") else
  for (i in seq_len(nrow(nul)))
    cat(sprintf("  %-6s – %-6s : r = %+.3f\n", nul$V1[i], nul$V2[i], nul$r[i]))

# Remarque : 
# - on voit déjà que CERE et OEUF sont très négativement corrélées,
#   ce qui laisse pressentir une opposition entre régimes céréaliers et animaux.
# - On voit également que VBLA et OEUF sont fortement corrélés positivement,
#   ce qui veut dire que la viande blanche et l'oeuf on tendance à être consommé
#   conjointement


# Q3 — Nuages de points sur 4 couples -------------------------------------------

couples <- list(c("VBLA","OEUF"), c("POISS","FRLEG"),
                c("OEUF","CERE"), c("OEUF","FRLEG"))

op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

for (cv in couples) {
  xv <- cv[1]; yv <- cv[2]
  x  <- actif[[xv]]; y <- actif[[yv]]
  r  <- round(cor(x, y), 3)
  col_pts <- ifelse(EST == "Est", "firebrick", "steelblue3")
  
  plot(x, y,
       xlab = xv, ylab = yv,
       main = sprintf("%s vs %s  (r = %+.3f)", xv, yv, r),
       pch = 19, col = col_pts, cex = 1.1)
  text(x, y, labels = rownames(actif), pos = 3, cex = 0.58, col = col_pts)
  abline(lm(y ~ x), col = "grey40", lty = 2, lwd = 1.5)
  legend("topleft", legend = c("Ouest", "Est"),
         col = c("steelblue3", "firebrick"), pch = 19, cex = 0.7, bty = "n")
  
  mod       <- lm(y ~ x)
  res_abs   <- abs(residuals(mod))
  atypiques <- rownames(actif)[res_abs > mean(res_abs) + 2 * sd(res_abs)]
  cat(sprintf("(%s, %s)  r = %+.3f", xv, yv, r))
  if (length(atypiques) > 0)
    cat(sprintf("  →  atypiques : %s", paste(atypiques, collapse = ", ")))
  cat("\n")
}

par(op)

# Commentaires sur les 4 nuages :
#
# (VBLA, OEUF) r ≈ +0.62 : tendance positive modérée bien visible.
#   Les pays qui consomment beaucoup de viande blanche tendent aussi à
#   consommer beaucoup d'oeufs. Pas d'individu vraiment atypique.
#
# (POISS, FRLEG) r ≈ +0.27 : très faible corrélation. il n'y a pas vraiment de 
# tendance. Le Portugal (pt) est
# clairement atypique, très éloigné du nuage principal en raison de
# sa consommation de poisson exceptionnellement élevée.
#
# (OEUF, CERE) r ≈ -0.71 : corrélation négative forte et bien visible.
#   Les pays de l'Est (rouge) se regroupent à droite (beaucoup de céréales,
#   peu d'oeufs) et les pays de l'Ouest à gauche. Cette opposition sera
#   l'axe principal de l'ACP.
#
# (OEUF, FRLEG) r ≈ -0.05 : quasi nulle, le nuage est totalement dispersé.
#   Ces deux variables sont indépendantes l'une de l'autre.


# ==============================================================================
#  PARTIE 2 — ACP
# ==============================================================================

# Q4 — Reproduction des résultats

acp <- PCA(donnees1, scale.unit = TRUE, quali.sup = 10, ncp = 9, graph = FALSE)

# Graphiques pour visualiser les résultats de Q4
print(fviz_screeplot(acp, addlabels = TRUE, ylim = c(0, 55),
                     main    = "Éboulis des valeurs propres",
                     barfill = "steelblue3", barcolor = "steelblue4",
                     linecolor = "firebrick"))

print(fviz_pca_var(acp, axes = c(1, 2), col.var = "cos2",
                   gradient.cols = c("grey70", "steelblue2", "firebrick"),
                   repel = TRUE, title = "Cercle des corrélations — Axes 1 x 2"))

print(fviz_pca_var(acp, axes = c(1, 3), col.var = "cos2",
                   gradient.cols = c("grey70", "steelblue2", "firebrick"),
                   repel = TRUE, title = "Cercle des corrélations — Axes 1 x 3"))

print(fviz_pca_ind(acp, axes = c(1, 2), habillage = EST,
                   addEllipses = TRUE, ellipse.type = "convex",
                   palette = c("steelblue3", "firebrick"), repel = TRUE,
                   title = "Individus — Axes 1 x 2  |  couleur = EST"))

print(fviz_pca_biplot(acp, axes = c(1, 2), habillage = EST,
                      palette = c("steelblue3", "firebrick"),
                      addEllipses = FALSE, col.var = "black", repel = TRUE,
                      title = "Biplot — Axes 1 x 2  |  couleur = EST"))


# Q5 — Répartition de l'inertie — combien d'axes retenir ? ---------------------
#
# Lecture des valeurs propres du prof :
#   Axe 1 : 4.01  →  4.01/9 = 44.6% d'inertie
#   Axe 2 : 1.63  →  1.63/9 = 18.1%
#   Axe 3 : 1.13  →  1.13/9 = 12.6%
#   Axe 4 : 0.95  →  inférieur à 1, on s'arrête ici
#
# Critère de Kaiser : on retient les axes dont la valeur propre est > 1.
# → 3 axes retenus (4.01 > 1, 1.63 > 1, 1.13 > 1, 0.95 < 1)
#
# Inertie cumulée sur 3 axes : 44.6 + 18.1 + 12.6 = 75.2%
# C'est la qualité globale de représentation de l'espace retenu.
# On considère ce niveau satisfaisant pour une analyse exploratoire.


# Q6 — Quels aliments déterminent les axes ? Effet de taille ? -----------------
#
# Critère : une variable est notable sur un axe si |r| > 1/sqrt(9) ≈ 0.33
# On observant les résultats fournis:
#
# Axe 1 (44.6% d'inertie) :
#   CERE  : r = +0.88  ->  0.88 > 0.33  corrélé positivement
#   NOIX  : r = +0.84  ->  0.84 > 0.33  corrélé positivement
#   OEUF  : r = -0.85  ->  0.85 > 0.33  corrélé negativement
#   LAIT  : r = -0.76  ->  0.76 > 0.33  corrélé negativement
#   VBLA  : r = -0.62  ->  0.62 > 0.33  corrélé negativement
#   VROU  : r = -0.61  ->  0.61 > 0.33  corrélé negativement
#   AMID  : r = -0.59  ->  0.59 > 0.33  corrélé negativement

#   POISS : r = -0.27  ->  0.27 < 0.33  ne contribue pas
#   FRLEG : r = +0.22  ->  0.22 < 0.33  ne contribue pas
#
# Axe 2 (18.1% d'inertie) :
#   POISS : r = -0.83  ->  0.83 > 0.33  corrélé negativement
#   FRLEG : r = -0.69  ->  0.69 > 0.33  corrélé negativement
#   AMID  : r = -0.45  ->  0.45 > 0.33  corrélé negativement
#   Les autres ont |r| < 0.33 -> ne contribuent pas
#
# Axe 3 (12.6% d'inertie) :
#   VBLA  : r = -0.66  ->  0.66 > 0.33  corrélé negativement
#   LAIT  : r = +0.41  ->  0.41 > 0.33  corrélé positivement
#   FRLEG : r = -0.43  ->  0.43 > 0.33  corrélé negativement
#   POISS : r = +0.34  ->  0.34 > 0.33  corrélé positivement
#   Les autres ont |r| < 0.33 -> ne contribuent pas
#
# Effet de taille :
#   Sur l'axe 1, CERE et NOIX sont positifs tandis que OEUF, LAIT, VBLA,
#   VROU, AMID sont négatifs -> les signes sont MIXTES.
#   Il n'y a donc pas d'effet de taille au sens strict.
#   Le centrage-réduction appliqué dans l'ACP a neutralisé l'effet de niveau
#   global. Pour l'obtenir malgré tout, il faudrait travailler sur des données
#   non centrées-réduites ou en proportions.


# Q7 — Quels pays déterminent les axes ? ----------------------------------------
#
# Critère : on s'appuie sur le cos² pour identifier les pays les plus
# représentatifs de chaque axe (cos² élevé = pays bien porté par cet axe).
# Lecture directe dans le tableau du prof.
#
# Axe 1 — pays avec cos² élevé sur Axis1 :
#   yu : cos²(Ax1) = 86.1%  coordonnée = +3.70  →  côté + (céréales/noix)
#   ro : cos²(Ax1) = 80.1%  coordonnée = +2.81  →  côté +
#   de : cos²(Ax1) = 79.6%  coordonnée = -2.14  →  côté − (protéines animales)
#   ie : cos²(Ax1) = 77.6%  coordonnée = -2.72  →  côté −
#   bg : cos²(Ax1) = 74.1%  coordonnée = +3.20  →  côté +
#   be : cos²(Ax1) = 70.9%  coordonnée = -1.66  →  côté −
#   dk : cos²(Ax1) = 66.3%  coordonnée = -2.41  →  côté −
#   al : cos²(Ax1) = 61.2%  coordonnée = +3.56  →  côté +
#
# Axe 2 — pays avec cos² élevé sur Axis2 :
#   pt : cos²(Ax2) = 80.5%  coordonnée = -4.38  →  côté − (poisson)
#   es : cos²(Ax2) = 65.3%  coordonnée = -2.61  →  côté −
#
# Axe 3 — pays avec cos² élevé sur Axis3 :
#   no : cos²(Ax3) = 48.1%  coordonnée = +1.74  →  côté + (lait)
#   hu : cos²(Ax3) = 47.6%  coordonnée = -1.95  →  côté − (viande blanche)
#   pl : cos²(Ax3) = 46.2%  coordonnée = -1.51  →  côté −
#   fi : cos²(Ax3) = 40.8%  coordonnée = +2.09  →  côté +


# Q8 — Interprétation des axes --------------------------------------------------
#
# En croisant Q6 (variables) et Q7 (pays) :
#
# Axe 1 (44.6% d'inertie) :
#   Pôle + : CERE (+0.88), NOIX (+0.84)
#             -> yu (+3.70), al (+3.56), bg (+3.20), ro (+2.81)  [tous bloc Est]
#   Pôle − : OEUF (-0.85), LAIT (-0.76), VBLA (-0.62), VROU (-0.61), AMID (-0.59)
#             -> ie (-2.72), dk (-2.41), de (-2.14)  [bloc Ouest]
#   C'est l'axe du TYPE DE RÉGIME ALIMENTAIRE dominant :
#   régime céréalier (Est) vs régime riche en protéines animales (Ouest).
#
# Axe 2 (18.1% d'inertie) :
#   Pôle − : POISS (-0.83), FRLEG (-0.69)
#             → pt (-4.38), es (-2.61)
#   C'est l'axe de la TRADITION MARITIME et MÉDITERRANÉENNE.
#   Attention : cet axe est fortement influencé par le Portugal (voir Q10).
#
# Axe 3 (12.6% d'inertie) :
#   Pôle + : LAIT (+0.41), POISS (+0.34)  →  fi (+2.09), no (+1.74)
#   Pôle − : VBLA (-0.66), FRLEG (-0.43)  →  hu (-1.95), pl (-1.51)
#   Axe de nuance au sein des pays européens, opposant les pays nordiques
#   à tradition laitière aux pays d'Europe centrale à forte consommation
#   de viande blanche.


# Q9 — Cinq pays les moins bien représentés -------------------------------------
#
# Critère : on additionne les cos² sur les 3 axes retenus (Axis1 + Axis2 + Axis3).
# Plus ce total est faible, plus le pays est mal représenté.
# Lecture et calcul directement depuis le tableau du prof :
#
#        Axis1   Axis2   Axis3   Total sur 3 axes
# ru  :  12.8  +  0.3  +  2.9  =  16.0%  ← le moins bien représenté
# fr  :  25.4  +  7.1  +  0.0  =  32.5%
# ch  :  20.6  + 14.0  +  0.6  =  35.2%
# it  :  44.6  +  3.0  +  0.3  =  47.9%
# pl  :   0.3  +  6.0  + 46.2  =  52.5%
#
# Les 5 pays les moins bien représentés sur les 3 axes retenus sont :
# ru (16.0%), fr (32.5%), ch (35.2%), it (47.9%), pl (52.5%)
#
# Leur profil alimentaire particulier se manifeste principalement sur les
# axes 4 et 5 non retenus. Par exemple ru a 18.0% sur Axis4 et 58.3%
# sur Axis5 — c'est dans ces dimensions que son originalité se révèle.


# Q10 — Le cas du Portugal ------------------------------------------------------
#
# Lecture directe dans les tableaux du prof :
#
# Coordonnées de pt :
#   Axis1 = +1.74  Axis2 = -4.38  Axis3 = -0.04  Axis4 = -0.91  Axis5 = -0.39
#
# Cos² de pt :
#   Axis1 = 12.7%  Axis2 = 80.5%  Axis3 = 0.0%
#   Total sur 3 axes retenus : 12.7 + 80.5 + 0.0 = 93.2%
#
# Sa coordonnée de -4.38 sur l'axe 2 est la plus extrême de tout le dataset
# — le pays suivant est es à -2.61, soit presque deux fois moins loin.
# Son cos² sur l'axe 2 atteint 80.5% : presque toute son originalité
# tient à cet axe unique, lié à sa consommation de poisson exceptionnelle.
#
# L'axe 2 est donc largement construit autour du seul cas Portugal.
# Un pays unique influence fortement la structure d'un axe entier,
# ce qui pose un problème d'interprétation.
#
# Pour le traiter correctement :
#   1. Le passer en individu supplémentaire : les axes seraient construits
#      sur la tendance générale des autres pays, et le Portugal observerait
#      a posteriori sans influencer les calculs
#   2. Appliquer une transformation log sur POISS pour réduire son effet
#      levier et rééquilibrer le poids de cette variable dans l'analyse


# ==============================================================================
#  PARTIE 3 — Validation par la variable supplémentaire EST
# ==============================================================================


# Q11 — Les données sont-elles suffisantes pour les valeurs test ? --------------
#
# On dénombre dans les données : 25 pays au total, 7 du bloc Est, 18 du bloc Ouest.
# La valeur test repose sur une approximation par la loi N(0,1) via le TCL.
# Le TCL est généralement admis à partir de n = 30.
# Ici n = 25 : on est en dessous du seuil.
# Les données ne sont donc pas tout à fait suffisantes au sens strict —
# les valeurs test restent indicatives mais doivent être interprétées
# avec prudence.


# Q12 — Conditions réunies ? Que nous apportent les valeurs test ? --------------
#
# Conditions d'utilisation :
#   1. EST doit être une variable SUPPLÉMENTAIRE dans l'ACP, pour que les
#      valeurs test mesurent une vraie association a posteriori et non un
#      artefact qu'on aurait introduit soi-même.
#      → Condition respectée : EST est en quali.sup, hors des calculs
#   2. n suffisant pour l'approximation normale.
#      → Condition limite : n = 25 < 30 (voir Q11)
#
# Lecture des valeurs test fournies par le prof, seuil bilatéral 5% = 1.96 :
#
#   Axis1 :  2.77  ->  |2.77| > 1.96  →  SIGNIFICATIF
#   Axis2 :  1.84  ->  |1.84| < 1.96  →  non significatif
#   Axis3 : -1.64  ->  |-1.64| < 1.96 →  non significatif
#   Axis4 : -1.79  ->  |-1.79| < 1.96 →  non significatif
#   Axis5 :  1.54  ->  |1.54| < 1.96  →  non significatif
#
# Seul l'axe 1 est significativement associé à EST.
# La valeur test positive (2.77) indique que les pays du bloc Est ont
# des coordonnées significativement plus élevées sur l'axe 1, c'est-à-dire
# qu'ils se situent du côté céréales/noix — contrairement aux pays de l'Ouest
# qui se situent du côté protéines animales.
#
# En excluant la limite sur n, ce résultat est interprétable et important :
# les habitudes alimentaires révèlent SPONTANÉMENT la partition politique
# Est/Ouest, sans que cette information ait été injectée dans l'analyse.
# La Guerre froide se lisait bien dans les assiettes — les contraintes
# économiques et politiques du bloc Est se traduisaient par des régimes
# alimentaires structurellement différents de ceux de l'Ouest.