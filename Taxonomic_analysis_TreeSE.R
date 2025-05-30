# Ordinate

# ARG load & ARG diversity vs. covariates and regions
# ALUE vs. ARG profiles table

# Load data
source("loadData.R")

# Load libraries
library(vegan)
library(mia)
library("readr")

# Test difference in TAXA and ARG sample dissimilarities
TAX_dist <-
  vegdist(t(assay(TSE, "relabundance")), method = "bray")
ARG_dist <-
  vegdist(t(assay(altExp(TSE, "arg"), "relabundance")), method = "bray")
# Slow and falls. Better as part of batch job than interactive.
#print("Mantel test p-value for Taxonomic vs. ARG dissimilarities:")
#print(mantel(ARG_dist, TAX_dist, method = "kendall"))

# Slow
# Assess the association between ALUE and ARG compositions.
adonis_df <-
  pairwise.adonis(
    x = sqrt(t(assay(altExp(TSE, "arg")))),
    factors = TSE$ALUE,
    sim.function = 'vegdist',
    sim.method = 'bray',
    p.adjust.m = 'fdr'
  )

adonis_EAST_WEST <-
  pairwise.adonis(
    x = sqrt(t(assay(altExp(TSE, "arg")))),
    factors = TSE$EAST,
    sim.function = 'vegdist',
    sim.method = 'bray',
    p.adjust.m = 'fdr'
  )

# Let us use writexl library as it is free of Java dependencies and works better across multiple systems
# TODO: manuscript Tables should be numbered instead of naming sheets
library(writexl)
xs <- list()
x <- adonis_df; rownames(x) <- NULL
xs[["PERMANOVA Geographical Regions and ARGs"]] <- x
x <- adonis_EAST_WEST; rownames(x) <- NULL
xs[["PERMANOVA Eastern and Western Finland and ARGs"]] <- adonis_EAST_WEST
write_xlsx(xs, path = "../RESULTS/R_results/Table_S4_ARG_ALUE_pairwise.xlsx")

# -------------------------------------------------------

# Better to have a single output file with the results, so commenting this out (xlsx done above)
#write_delim(
#  adonis_df,
#  file  = ("../RESULTS/R_results/ARG_ALUE_pairwise_ADONIS_TreeSE.tsv"),
#  delim  = "\t"
#)

# Check correlations between distances roughly
#library(microbiome)
#x <- as.matrix(TAX_dist)
#y <- as.matrix(ARG_dist)
#xvec <- x[lower.tri(x)]
#yvec <- y[lower.tri(y)]
#inds <- sample(length(xvec), 1e3);
#cor(xvec[inds], yvec[inds])
#df <- data.frame(x=xvec[inds], y=yvec[inds]);
#ggplot(df, aes(x,y)) + geom_point(alpha = 0.1)


