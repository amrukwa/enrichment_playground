library(tmod)
load("data/data_lung_cancer.RData")

# Differentiating genes search
source("source/expression_testing.R")
df <- means_tests(data, metaInfo)
de_genes <- rownames(df[df$corrected_pval < 0.05,])
ordered_genes <- rownames(df[with(df, order(corrected_pval)), ])
rm(df)

# ORA
source("source/ora.R")
ora_result <- ora(data, de_genes, KEGGhsa)
significant_ora <- ora_result[ora_result$corrected_pvals < 0.05,]
rora <- significant_ora[order(significant_ora$corrected_pvals), c("Title", "corrected_pvals")]
head(rora, n=15L)
nrow(rora)
rm(rora)

# CERNO
source("source/cerno.R")
cerno_result <- cerno(ordered_genes, KEGGhsa)
significant_cerno <- cerno_result[cerno_result$corrected_pvals < 0.05,]
rcerno <- significant_cerno[order(significant_cerno$corrected_pvals), c("Title", "corrected_pvals")]
head(rcerno, n=15L)
nrow(rcerno)
rm(rcerno)

# Z-transform
source("source/cerno.R")
z_result <- cerno(ordered_genes, KEGGhsa, combining_method = "z-transform")
significant_z <- z_result[z_result$corrected_pvals < 0.05,]
rz <- significant_z[order(significant_z$corrected_pvals), c("Title", "corrected_pvals")]
head(rz, n=15L)
nrow(rz)
rm(rz)

# GSEA
load("data/results.RData")
significant_s2n <- s2n[s2n$NES_pval < 0.05,]
significant_s2n_abs <- s2n_abs[s2n_abs$NES_pval < 0.05,]
significant_lfc <- lfc[lfc$NES_pval < 0.05,]
significant_lfc_abs <- lfc_abs[lfc_abs$NES_pval < 0.05,]
nrow(significant_s2n_abs)
nrow(significant_s2n)
nrow(significant_lfc_abs)
nrow(significant_lfc)

# GSVA
source("source/gsva.R")
gsva_results <- gsva(data, KEGGhsa, metaInfo)
GSVA_ES <- gsva_results$ES
gsva_result <- gsva_results$pvals
gsva_result$Title <- rownames(gsva_result)
rownames(gsva_result) <- NULL
significant_gsva <- gsva_result[gsva_result$corrected_pvals < 0.05,]
rgsva <- significant_gsva[order(significant_gsva$corrected_pvals), c("Title","corrected_pvals")]
head(rgsva, n=15L)
nrow(rgsva)
rm(rgsva)

# PLAGE
source("source/plage.R")
plage_result <- plage(data, KEGGhsa, metaInfo)
significant_plage <- plage_result[plage_result$corrected_pvals < 0.05,]
rplage <- significant_plage[order(significant_plage$corrected_pvals), c("Title", "corrected_pvals")]
head(rplage, n=15L)
nrow(rplage)
rm(rplage)

# results correlation
library(corrplot)
pval_comp=cbind(ORA=ora_result[,c("ID", "corrected_pvals")], CERNO=cerno_result[,"corrected_pvals"], 
                Z=z_result[,"corrected_pvals"],
                S2N=s2n[,"NES_pval"],  S2N_abs=s2n_abs[,"NES_pval"],  
                LFC=lfc[,"NES_pval"], LFC_abs=lfc_abs[,"NES_pval"],
                PLAGE=plage_result[, "corrected_pvals"], GSVA=gsva_result[, "corrected_pvals"])

colnames(pval_comp)[c(1,2)] <- c("ID", "ORA")
rownames(pval_comp) <- pval_comp[,1]
pval_comp[,1] <- NULL
log_pval_comp <- (log(pval_comp))
M <- cor(log_pval_comp, method="spearman")
n <- ncol(log_pval_comp )
pcorr <- matrix(0, n, n)
for (i in 1:(n-1)){
  for (j in (i+1):n){
    val <- cor.test(log_pval_comp[, i], log_pval_comp[, j], method="spearman")$p.value
    pcorr[i,j] <- val
    pcorr[j,i] <- val
  }
}
colnames(pcorr) <- colnames(log_pval_comp)
rownames(pcorr) <- colnames(pcorr)
png("presentation/plots/correlation.png")
correlo <- corrplot(M, method="circle", type="upper", 
                    col=colorRampPalette(c("blue", "white", "red"))(200), 
                    tl.col="black", p.mat = pcorr, sig.level = 0.05, insig = "blank")
dev.off()
rm(correlo, pcorr, M, log_pval_comp, i, j, n)

# joint GS
joint_names <- Reduce(intersect, list(significant_z$Title, significant_cerno$Title, 
                                      significant_lfc_abs$Title, significant_lfc$Title,
                                      significant_s2n_abs$Title, significant_s2n$Title,
                                      significant_ora$Title, significant_plage$Title,
                                      significant_gsva$Title))
length(joint_names)

# integration and so on
Z_vals <- apply(pval_comp, 2, function(x) qnorm(x))
k <- ncol(pval_comp)
zcombined <- apply(Z_vals, 1, function(x) sum(x)/sqrt(k))
pval_combined <- pnorm(zcombined)
combined_pvals <- data.frame(pval_combined)
combined_pvals$Title <- cerno_result$Title
significant_combined <- combined_pvals[combined_pvals$pval_combined < 0.05,]
rownames(significant_combined) <- NULL
significant_combined <- significant_combined[order(significant_combined$pval_combined), 
                                             c("Title", "pval_combined")]
head(significant_combined, n=15L)
nrow(significant_combined)

# Clearing after first set of methods
rm(significant_cerno, significant_combined, significant_gsva, significant_lfc, significant_lfc_abs,
   significant_ora, significant_plage, significant_s2n, significant_s2n_abs, significant_z)
rm(cerno_result, ora_result, z_result, plage_result, gsva_results, gsva_result,
   lfc, lfc_abs, s2n, s2n_abs, Z_vals, k, zcombined)


# Best pathways visualization
source("source/visualization.R")

rank_pathways <- order(pval_combined)
pathways <- KEGGhsa[rank_pathways[1:10]]

subplots <- visualize_pathways(data, metaInfo, pathways, GSVA_ES)
save(subplots, file = "data/plots.RData")

for (i in 1:length(subplots)){
  fig <- subplots[[i]]
  print(pathways[i]$MODULES$Title)
  print(fig)
}

# single cerno
source("source/cerno.R")
heatmap_all_cerno <- cerno_heatmaps(data, KEGGhsa, color_labels=metaInfo$Group, sort_type="abs", 
                                    with_dendro='none')

heatmaps_best <- cerno_heatmaps(data, pathways, color_labels=metaInfo$Group, 
                        sort_type="abs", with_dendro=TRUE)
save(heatmaps_best, subplots, heatmap_all_cerno, file = "data/plots.RData")
#rm(heatmaps_best, heatmap_all_cerno, subplots)
