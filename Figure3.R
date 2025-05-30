library(cowplot)
library(mia)
library(miaViz)
library(ggpubr)
library(patchwork)
library(grid)
library(microViz)
library(viridis)
library(multcomp)
library(phyloseq)
library(data.table)

theme_set(theme_bw(20))
library(circlize)
library(mia)

# Load functions
source("R_functions.R")

# Read the data (ordination pre-calculated)
TSE <- readRDS(file="../data/TSE.rds")

# First run to get family-ARG associations
Family_ARG_associations <- readRDS("Family_ARG_associations.rds")

source("FigureES.R") 

# -------------------------------------------------------------------
# Commented out from here and moved to enrichments.R
# # Calculate Top ARG enrichments by family abundance bin
# nbreaks <- 5
# theme_set(theme_bw(20))
# sigfam <- setdiff(levels(TSE$Signif_fam), "Other")
# dfe <- calculateEnrichment(altExp(TSE, "Family", withColData=TRUE), varname="ARG_burd", assay.type="relabundance", features=sigfam, nbreaks=nbreaks) 
# 
# # ARG classes vs. enterosignatures
# # corsf <- cor(t(assay(altExp(TSE, "Family_aggregated"), "counts")), TSE$SUM_norm, method="kendall")
# # Calculate all NMF associations at one go to deal with multiple testing properly
# F <- t(assay(altExp(TSE, "FamilyPrevalent"), "counts"))
# colData(TSE) <- cbind(colData(TSE), F)
# selected.vars <- "SUM_norm"
# 
# # Combine family abundances and colData
# DF <- cbind(as.data.frame(colData(TSE)), as.data.frame(t(assay(altExp(TSE, "Family"), "relabundance"))))

##########################################################

# Add tree
source("doTree.R")

##########################################################

source("pcoas.R")

##########################################################

source("enrichments.R")

##########################################################

# Figure 3a
source("resistome_heatmap.R")

##########################################################

# ARG diversity vs. Species diversity
theme_set(theme_bw(20))
fig.div <- ggplot(as.data.frame(colData(TSE)), aes(x=Species_diversity, y=ARG_div)) +
       geom_point(size=1, alpha=0.25) +
       labs(x="Species diversity", y="ARG diversity")
r2 <- cor(TSE$Species_diversity, TSE$ARG_div)
print(paste("Pearson:", r2))
cor.test(TSE$Species_diversity, TSE$ARG_div)

##########################################################

fig3a <- fig.div
fig3b <- figa / pic.es.associations
fig3b <- plot_grid(fig3a, fig3b, ncol=2, rel_widths=c(0.3, 0.7), labels=c("b", "c"), label_size=30)
fig3tot <- plot_grid(p.resistome.scores,
	             fig3b,
                     # p.es.scores + labs(title="Samples (n=7,095)"),
		     fig.tree, # fig3c
		     nrow=3,
		     rel_heights=c(0.15, 0.25, 0.6),
		     labels=c("a", "", "d"), label_size=30)


li2 <- lapply(dic.nmf, function (es) {myplot(assay(altExp(TSE, "ES"), "signal"), es, altExp(TSE, "arg"))})
names(li2) <- dic.nmf[names(li2)]
li2 <- li2[levels(df2$Level)] # sort
figS6C <- cowplot::plot_grid(
  li2[[1]], li2[[2]], li2[[3]], li2[[4]], li2[[5]], nrow=1) +
  theme(plot.margin  = unit(c(1, 1, 1, 1), "pt"),
	plot.title = element_text(size = base.size),  
        legend.text = element_text(size=0.7*base.size),
	legend.title = element_text(size=base.size),
                                  axis.text.x = element_text(size = base.size),
                                  axis.text.y = element_text(size = base.size),
                                  axis.title.x = element_text(size = base.size),
                                  axis.title.y = element_text(size = base.size)	  
				  ) 

# ----------------------------------------------------------------------------

# Final tables
write_xlsx(ES_ARG_associations, "../RESULTS/R_results/Table_S5_ES_ARG_association.xlsx") # Related to Figure 2

# Final figures

# Main Figure 3; tree panel figure
library(Cairo)
w <- 1500
CairoJPEG("../RESULTS/R_results/Fig3.jpg", width=w, height=1.25*w)
print(fig3tot)
dev.off()

pdf("../RESULTS/R_results/FigS5A_Enterosignature_sampleprofile.pdf", width=25, height=4)
print(p.es.scores)
dev.off()

pdf("../RESULTS/R_results/FigS5B_Enterosignature_genusprofile.pdf", width=30, height=4)
draw(p1horizontal, auto_adjust = FALSE)
dev.off()

pdf("../RESULTS/R_results/FigureS6C.pdf", width=25, height=7)
print(figS6C)
dev.off()

# PCoA + ARG enrichments by family
library(cowplot)
figS7 <- plot_grid(metaphlan.plot, fig3enr, ncol=2, rel_widths=c(13,3), labels="auto", label_size=30)
s <- 800
CairoJPEG("../RESULTS/R_results/FigureS7.jpg", width=2*s, height=1.4*s)
print(figS7)
dev.off()







