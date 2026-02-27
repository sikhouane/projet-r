# =========================
# EXO 7 - AFC (cas 1)
# =========================

donnees2 <- matrix(c(
  20,  9,  9, 27, 10, 16, 20,  4,  8,   # PAYS
  42, 10, 22, 51, 18, 28, 38, 12, 22,   # OUVR
  11,  2,  5, 14,  8,  7,  5,  8,  6,   # VEND
  8,  9, 12, 23, 14, 16, 14, 12, 12,   # COMM
  19, 10, 16, 52, 32, 25, 22, 25, 30,   # EMPL
  10,  5, 12, 23, 20, 13, 11, 13, 10,   # TECH
  2,  8,  7,  6, 15,  6,  6,  9,  4,   # UNIV
  8, 42, 23, 24, 46, 22, 22, 34, 16    # LIBE
), byrow = TRUE, nrow = 8)


rownames(donnees2) <- c("PAYS","OUVR","VEND","COMM","EMPL","TECH","UNIV","LIBE")
colnames(donnees2) <- c("seri","gene","gai","honn","intl","serv","cour","comp","disc")

donnees2
addmargins(donnees2) 

# INTRODUCTION : 
#1-	Pourcentage de chaque qualité + pourcentage des catégories pro 

N <- sum(donnees2)     
tot_lignes <- rowSums(donnees2)
tot_col <- colSums(donnees2)

pct_col <- 100 * tot_col / N       # % par qualité
pct_lignes <- 100 * tot_lignes / N # % par profession

pct_col
pct_lignes


# 2-a : nombre de personnes sondés: 
nb_personnes <- N / 3
nb_personnes

# 2-b: •	Proportion des employés pour qui être honnête rend sympathique : 
p_honn_sachant_empl <- donnees2["EMPL","honn"] / tot_lignes["EMPL"]
p_honn_sachant_empl


# 2.c) “Proportion d’employés parmi ceux qui pensent qu’être honnête rend sympathique” 
p_empl_sachant_honn <- donnees2["EMPL","honn"] / tot_col["honn"]
p_empl_sachant_honn


# ============================================================
# PARTIE 2 : AFC (Analyse des correspondances) via SVD
# ============================================================

# 1) Matrice des fréquences + masses
P <- donnees2 / N
r <- rowSums(P)
c <- colSums(P)

Dr_inv_sqrt <- diag(1 / sqrt(r))
Dc_inv_sqrt <- diag(1 / sqrt(c))

# 2) Résidus standardisés (indépendance)
S <- Dr_inv_sqrt %*% (P - r %*% t(c)) %*% Dc_inv_sqrt

# 3) SVD
sv <- svd(S)
eig <- sv$d^2
pct_inertie <- 100 * eig / sum(eig)

# 4) Coordonnées principales
U <- sv$u
V <- sv$v
Delta <- diag(sv$d)

F <- Dr_inv_sqrt %*% U %*% Delta   # lignes (professions)
G <- Dc_inv_sqrt %*% V %*% Delta   # colonnes (qualités)

rownames(F) <- rownames(donnees2)
rownames(G) <- colnames(donnees2)

# 5) Contributions (%)
ctr_rows <- sweep(F^2, 1, r, "*")
ctr_rows <- sweep(ctr_rows, 2, eig, "/") * 100

ctr_cols <- sweep(G^2, 1, c, "*")
ctr_cols <- sweep(ctr_cols, 2, eig, "/") * 100

# 6) Qualité de représentation cos² (%)
cos2_rows <- (F^2 / rowSums(F^2)) * 100
cos2_cols <- (G^2 / rowSums(G^2)) * 100

# 7) Plan (1,2)
F12 <- F[, 1:2]
G12 <- G[, 1:2]

ctr_rows12 <- ctr_rows[, 1:2]
ctr_cols12 <- ctr_cols[, 1:2]

cos2_rows12 <- cos2_rows[, 1:2]
cos2_cols12 <- cos2_cols[, 1:2]

qual_rows_plan12 <- rowSums(cos2_rows12)
qual_cols_plan12 <- rowSums(cos2_cols12)

# ============================================================
# TABLES
# ============================================================

table_inertie <- data.frame(
  Axe = 1:length(eig),
  Valeur_propre = eig,
  Pourcentage_inertie = pct_inertie
)

table_coord_rows <- data.frame(Modalite = rownames(F12), Axe1 = F12[,1], Axe2 = F12[,2])
table_coord_cols <- data.frame(Modalite = rownames(G12), Axe1 = G12[,1], Axe2 = G12[,2])

table_ctr_rows <- data.frame(Modalite = rownames(ctr_rows12), Axe1_pct = ctr_rows12[,1], Axe2_pct = ctr_rows12[,2])
table_ctr_cols <- data.frame(Modalite = rownames(ctr_cols12), Axe1_pct = ctr_cols12[,1], Axe2_pct = ctr_cols12[,2])

table_cos2_rows <- data.frame(Modalite = rownames(cos2_rows12), Axe1 = cos2_rows12[,1], Axe2 = cos2_rows12[,2],
                              Qualite_plan12 = qual_rows_plan12)

table_cos2_cols <- data.frame(Modalite = rownames(cos2_cols12), Axe1 = cos2_cols12[,1], Axe2 = cos2_cols12[,2],
                              Qualite_plan12 = qual_cols_plan12)

# Affichages arrondis 
print(transform(table_inertie,
                Valeur_propre = round(Valeur_propre, 6),
                Pourcentage_inertie = round(Pourcentage_inertie, 1)))

print(transform(table_coord_rows, Axe1 = round(Axe1,4), Axe2 = round(Axe2,4)))
print(transform(table_coord_cols, Axe1 = round(Axe1,4), Axe2 = round(Axe2,4)))

print(transform(table_ctr_rows, Axe1_pct = round(Axe1_pct,1), Axe2_pct = round(Axe2_pct,1)))
print(transform(table_ctr_cols, Axe1_pct = round(Axe1_pct,1), Axe2_pct = round(Axe2_pct,1)))

print(transform(table_cos2_rows,
                Axe1 = round(Axe1,1), Axe2 = round(Axe2,1), Qualite_plan12 = round(Qualite_plan12,1)))
print(transform(table_cos2_cols,
                Axe1 = round(Axe1,1), Axe2 = round(Axe2,1), Qualite_plan12 = round(Qualite_plan12,1)))

# Modalités les plus mal représentées (plan 1-2)
tmp_rows <- round(table_cos2_rows$Qualite_plan12, 1); names(tmp_rows) <- table_cos2_rows$Modalite
tmp_cols <- round(table_cos2_cols$Qualite_plan12, 1); names(tmp_cols) <- table_cos2_cols$Modalite

sort(tmp_rows)  # lignes (petit = mal représenté)
sort(tmp_cols)  # colonnes (petit = mal représenté)

# ============================================================
# GRAPHIQUE : plan factoriel (axes 1-2)
# ============================================================
par(mar = c(4, 4, 2, 1))

plot(F12, type = "n",
     xlab = paste0("Axe 1 (", round(pct_inertie[1], 1), "%)"),
     ylab = paste0("Axe 2 (", round(pct_inertie[2], 1), "%)"),
     main = "AFC - Plan factoriel (axes 1-2)")
abline(h = 0, v = 0, lty = 2)
text(F12, labels = rownames(F12))  # professions
text(G12, labels = rownames(G12))  # qualités


