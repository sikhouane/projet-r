donnees2 <- matrix(c(
  20,  9,  9, 27, 10, 16, 20,  4,  8,
  42, 10, 22, 51, 18, 28, 38, 12, 22,
  11,  2,  5, 14,  8,  7,  5,  8,  6,
  8,  9, 12, 23, 14, 16, 14, 12, 12,
  19, 10, 16, 52, 32, 25, 22, 25, 30,
  10,  5, 12, 23, 20, 13, 11, 13, 10,
  2,  8,  7,  6, 15,  6,  6,  9,  4,
  8, 42, 23, 24, 46, 22, 22, 34, 16
), byrow=TRUE, nrow=8)

rownames(donnees2) <- c("PAYS","OUVR","VEND","COMM","EMPL","TECH","UNIV","LIBE")
colnames(donnees2) <- c("seri","gene","gai","honn","intl","serv","cour","comp","disc")

Nvotes <- sum(donnees2)           # 1200
Npers  <- Nvotes/3                # 400

pct_cols <- 100*colSums(donnees2)/Nvotes
pct_rows <- 100*rowSums(donnees2)/Nvotes

p_honn_given_empl <- donnees2["EMPL","honn"] / rowSums(donnees2)["EMPL"]
p_empl_given_honn <- donnees2["EMPL","honn"] / colSums(donnees2)["honn"]

Nvotes; Npers
pct_cols; pct_rows
p_honn_given_empl; p_empl_given_honn



# Packages
install.packages(c("FactoMineR","factoextra"))
library(FactoMineR)
library(factoextra)


# AFC
res.ca <- CA(donnees2, graph = FALSE)

# 1) Valeurs propres + % inertie
res.ca$eig
# colonnes typiques: eigenvalue / % of variance / cumulative %

# 2) Coordonnées (lignes + colonnes)
res.ca$row$coord
res.ca$col$coord

# 3) Contributions aux axes (%)
res.ca$row$contrib
res.ca$col$contrib

# 4) Qualité de représentation cos² (%)
res.ca$row$cos2
res.ca$col$cos2

# 5) Graphiques
fviz_screeplot(res.ca, addlabels = TRUE)
fviz_ca_biplot(res.ca, repel = TRUE)             
fviz_ca_row(res.ca, repel = TRUE)                
fviz_ca_col(res.ca, repel = TRUE)                

# (Option) Voir qui définit le plus Dim1 / Dim2
fviz_contrib(res.ca, choice = "row", axes = 1, top = 10)
fviz_contrib(res.ca, choice = "col", axes = 1, top = 10)
fviz_contrib(res.ca, choice = "row", axes = 2, top = 10)
fviz_contrib(res.ca, choice = "col", axes = 2, top = 10)

