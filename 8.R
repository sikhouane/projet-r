########################### AFC - Cas 2 : Election 2007

library(FactoMineR)
library(factoextra)


read_election <- function(path) {
  elec <- read.csv(path,
                   sep = ";",
                   fileEncoding = "latin1",
                   check.names = FALSE)
  
  # 1ère colonne = Région
  names(elec)[1] <- "Region"
  rownames(elec) <- elec$Region
  elec$Region <- NULL
  
  if ("France" %in% rownames(elec)) {
    elec <- elec[rownames(elec) != "France", , drop = FALSE]
  }
  
  for (j in 1:ncol(elec)) elec[, j] <- as.numeric(elec[, j])
  
  elec
}

elec <- read_election("data/election 2007.csv")

##Outils : Chi² + AFC + Plot

chi2_report <- function(tab) {
  test <- chisq.test(tab)
  list(
    stat = unname(test$statistic),
    df = unname(test$parameter),
    pvalue = unname(test$p.value)
  )
}

run_ca <- function(tab) {
  CA(tab, graph = FALSE)
}

plot_ca <- function(res, title, subtitle = "", row_col = "steelblue", col_col = "firebrick") {
  fviz_ca_biplot(
    res,
    repel = TRUE,
    col.row = row_col,
    col.col = col_col,
    title = title,
    subtitle = subtitle
  )
}

######### AFC 1 : 1er tour complet

tour1_cols <- c("Besancenot", "Buffet", "Schivardi", "Bayrou", "Bové",
                "Voynet", "de Villiers", "Royal", "Nihous",
                "Le Pen", "Laguiller", "Sarkozy")

tab_t1 <- elec[, tour1_cols, drop = FALSE]

chi_t1 <- chi2_report(tab_t1)
chi_t1

res_t1 <- run_ca(tab_t1)

plot_ca(res_t1,
        title = "AFC 1 — Présidentielle 2007 — 1er tour (tous candidats)",
        subtitle = paste0("Test χ² : p-value = ", signif(chi_t1$pvalue, 3)),
        row_col = "navy",
        col_col = "firebrick")

######### AFC 2 : 1er tour — 4 grands candidats

tab_t1_4 <- elec[, c("Royal", "Sarkozy", "Bayrou", "Le Pen"), drop = FALSE]

chi_t1_4 <- chi2_report(tab_t1_4)
chi_t1_4

res_t1_4 <- run_ca(tab_t1_4)

plot_ca(res_t1_4,
        title = "AFC 2 — 1er tour — 4 grands candidats",
        subtitle = paste0("Test χ² : p-value = ", signif(chi_t1_4$pvalue, 3)),
        row_col = "darkgreen",
        col_col = "firebrick")

####### AFC 3 : 1er tour — blocs politiques

tab_blocs <- data.frame(
  Gauche = rowSums(elec[, c("Besancenot", "Buffet", "Schivardi", "Bové",
                            "Voynet", "Laguiller", "Royal"), drop = FALSE]),
  Centre = elec[, "Bayrou"],
  Droite = rowSums(elec[, c("de Villiers", "Nihous", "Sarkozy"), drop = FALSE]),
  ExtDroite = elec[, "Le Pen"]
)
rownames(tab_blocs) <- rownames(elec)

chi_blocs <- chi2_report(tab_blocs)
chi_blocs

res_blocs <- run_ca(tab_blocs)

plot_ca(res_blocs,
        title = "AFC 3 — 1er tour — Blocs politiques",
        subtitle = paste0("Test χ² : p-value = ", signif(chi_blocs$pvalue, 3)),
        row_col = "gray20",
        col_col = c("red3", "purple3", "royalblue3", "black"))


