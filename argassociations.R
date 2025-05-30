

# Test associations between ARG burden and taxa abundance

source("R_functions.R")
taxa <- c("Prevotella", "Bacteroides", "Escherichia", "Faecalibacterium", "Alistipes", "Ruminococcus", "Roseburia", "Dialister", "Bifidobacterium")

# Line plot of ARG burden classes as a function of taxa abundance quantiles
plots <- lapply(taxa, function (tax) plotfun1(tax))
library(patchwork)
(plots[[1]] + plots[[2]] + plots[[3]]) / 
(plots[[4]] + plots[[5]] + plots[[6]]) /
(plots[[7]] + plots[[8]] + plots[[9]])

# Same but now as stacked barplots
plots <- lapply(taxa, function (tax) plotfun2(tax))
library(patchwork)
(plots[[1]] + plots[[2]] + plots[[3]]) / 
(plots[[4]] + plots[[5]] + plots[[6]]) /
(plots[[7]] + plots[[8]] + plots[[9]])

# Showing that these trends are not so clear from mere jitteplots
plots <- lapply(taxa, function (tax) plotfun3(tax))
library(patchwork)
(plots[[1]] + plots[[2]] + plots[[3]]) / 
(plots[[4]] + plots[[5]] + plots[[6]]) /
(plots[[7]] + plots[[8]] + plots[[9]])

# Same but with scatterplot
plots <- lapply(taxa, function (tax) plotfun4(tax))
library(patchwork)
(plots[[1]] + plots[[2]] + plots[[3]]) / 
(plots[[4]] + plots[[5]] + plots[[6]]) /
(plots[[7]] + plots[[8]] + plots[[9]])

# Instead taxa abundance as a function of ARG burden
# Issues: how to treat samples with zero abundance taxa
# (these distort the picture, and removing them is also questionable plus does not make the
# figure any more clear)
plots <- lapply(taxa, function (tax) plotfun5(tax))
library(patchwork)
(plots[[1]] + plots[[2]] + plots[[3]]) / 
(plots[[4]] + plots[[5]] + plots[[6]]) /
(plots[[7]] + plots[[8]] + plots[[9]])



# boxplot(assay(altExp(TSE, "Genus"), "clr")["Prevotella",] ~ colData(TSE)$ARG_burd)
# boxplot(df$Abundance ~ df$ARG_burd)

# tax <- "Prevotella"
#A <- assay(altExp(TSE, "Genus"), "relabundance")[tax,]
#prop.table(table(cut(A, breaks=c(0, 0.1, 1, 10, 100)/100, include.lowest=TRUE), colData(TSE)$ARG_burd), 1)
#                     low   med_low  med_high      high
#  [0,0.001]    0.2076811 0.2552703 0.2638280 0.2732206
#  (0.001,0.01] 0.2950192 0.2145594 0.2107280 0.2796935
#  (0.01,0.1]   0.3265027 0.2609290 0.2185792 0.1939891
#  (0.1,1]      0.3508772 0.2326468 0.2257818 0.1906941
