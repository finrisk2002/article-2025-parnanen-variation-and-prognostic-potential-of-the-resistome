library(mia)
library(miaViz)
library(ggpubr)
library(patchwork)
library(grid)
library(microViz)
library(viridis)
library(multcomp)
library(phyloseq)

theme_set(theme_bw(20))

# Load functions
source("R_functions.R")

# Read the data (ordination pre-calculated)
TSE <- readRDS("../data/TSE.rds") 

# Figure 3B: associations between ARG and driver taxa

# df_SP_axes <- cbind(reducedDim(TSE), colData(TSE)) %>%
df_SP_axes <- cbind(reducedDim(TSE), DF) %>%
                as.data.frame %>% 
                dplyr::rename(ARG_load = SUM_norm)

lm(log(ARG_load)~PC1+PC2+PC3, data = df_SP_axes) %>% summary
lm(log(ARG_load)~Bacteroidaceae+Enterobacteriaceae+Bifidobacteriaceae+Prevotellaceae+Lachnospiraceae, data = df_SP_axes) %>% summary

tops <- c("Bacteroidaceae","Enterobacteriaceae","Prevotellaceae","Bifidobacteriaceae","Lachnospiraceae")
# tops <- c("Bacteroidaceae","Enterobacteriaceae",""Prevotellaceae")
plots <- lapply(tops, function (i) plotfun(i, df_SP_axes)); names(plots) <- tops
# This was previously added on Prevotellaceae; not sure if really needed
# scale_color_manual(values = pal, "") + scale_fill_manual(values = pal, "") 

fig3B <- cowplot::plot_grid(
                   align = "h", 
                   plots[["Bacteroidaceae"]],
                   plots[["Enterobacteriaceae"]],
                   plots[["Prevotellaceae"]],
                   nrow=1, 
		   labels = c("c", "d", "e"))

# Family-level associations
pdf("../RESULTS/R_results/Fig3B.pdf", width=20, height=10)
print(fig3B)
dev.off()

# ===============================================================================================

# Genus-level associations
# between ARG burden and taxa abundance

# TODO: replace with our final selections
taxa <- c("Prevotella", "Bacteroides", "Escherichia", "Faecalibacterium", "Alistipes", "Ruminococcus", "Roseburia", "Dialister", "Bifidobacterium")

# Line plot of ARG burden classes as a function of taxa abundance quantiles
plots <- lapply(taxa, function (tax) plotfun1(tax))
library(patchwork)
p <- (plots[[1]] + plots[[2]] + plots[[3]]) / 
     (plots[[4]] + plots[[5]] + plots[[6]]) /
     (plots[[7]] + plots[[8]] + plots[[9]])


pdf("../RESULTS/R_results/Fig3Bb.pdf", width=20, height=10)
print(p)
dev.off()

# ===============================================================================================

# check that Top_fam are highlighted as well

# Family-level associations
# between ARG burden and taxa abundance

# TODO: replace with our final selections
taxa <- setdiff(as.character(unique(TSE$Top_fam)), c("Other"))

# Line plot of ARG burden classes as a function of taxa abundance quantiles
plots <- lapply(taxa, function (tax) {plotfun1f(tax)})
p <- cowplot::plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], 
                        plots[[6]], plots[[7]], plots[[8]], 
                        ncol=4)

pdf("../RESULTS/R_results/Fig3Bb.pdf", width=20, height=10)
print(p)
dev.off()

# ===============================================================================================

# NMF-level associations
# between ARG burden and abundance

# Line plot of ARG burden classes as a function of taxa abundance quantiles
esnames <- c(
        "Bifi",
	"Esch",	
	"Prev",
        "Bact",
        "Firm"
	)
plots <- lapply(1:5, function (es) plotfunes(es, nbreaks=10, esname=esnames[[es]]))
library(patchwork)
p <- cowplot::plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], ncol=5)

#pdf("../RESULTS/R_results/Fig3Bb.pdf", width=20, height=10)
print(p)
#dev.off()

# ===============================================================================================

# tax <- "Prevotella"
#A <- assay(altExp(TSE, "Genus"), "relabundance")[tax,]
#prop.table(table(cut(A, breaks=c(0, 0.1, 1, 10, 100)/100, include.lowest=TRUE), colData(TSE)$ARG_burd), 1)
#                     low   med_low  med_high      high
#  [0,0.001]    0.2076811 0.2552703 0.2638280 0.2732206
#  (0.001,0.01] 0.2950192 0.2145594 0.2107280 0.2796935
#  (0.01,0.1]   0.3265027 0.2609290 0.2185792 0.1939891
#  (0.1,1]      0.3508772 0.2326468 0.2257818 0.1906941
