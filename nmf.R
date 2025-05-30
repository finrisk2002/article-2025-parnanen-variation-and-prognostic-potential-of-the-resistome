library(NMF)
library(mia)
library(ComplexHeatmap)
library(reshape2)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)
library(mia)

# Merge by prevalence 
#altExp(TSE, "GenusPrevalent") <- mergeFeaturesByPrevalence(TSE[rowData(TSE)$Domain=="Bacteria"], rank="Genus", assay.type="relabundance", detection=0.1/100, prevalence=1/100) #moved to Carpentry

set.seed(33)
vec <- 2:7
x <- t(assay(altExp(TSE, "GenusPrevalent"), "relabundance"))

# Fit NMF with the optimal component number
e <- nmf(x, rank=vec, nrun = 10, seed = 12345) # 10-20 mins on atlas
Ncomp <- vec[[which.min(e$measures$silhouette.consensus)]]

set.seed(3221)
nmf.optimal <- nmf(x, Ncomp)

# Pick NMF components
H <- nmf.optimal@fit@H
W <- nmf.optimal@fit@W

# Sort as in original ES paper
nmf.order <- c(4,5,3,1,2)
H <- H[nmf.order, ]
W <- W[,nmf.order]

colnames(W) <- rownames(H) <- paste("nmf", 1:ncol(W), sep = "")

# Add enterosignature scores to colData
colData(TSE) <- DataFrame(as.data.frame(colData(TSE)) %>% select(!starts_with("nmf"))) # first remove conflicting entries
colData(TSE) <- DataFrame(cbind(colData(TSE), W))

# Add NMF loadings in the metadata slot
metadata(TSE)$NMF_loadings <- H

# Obtain the relative abundance of ES in each sample.
Wrel <- t(apply(W, 1, function (x) {x/sum(x)}))

# Define as primary ES the ES of a sample with the highest relative abundance
# and add this to colData
colData(TSE)$ES_primary <- apply(Wrel, 1, which.max)

##################

# source("FigureES.R")

# ES profile across bacteria (normalized to 1 per ES)
# ESprofiles <- apply(H, 1, function (x) {x/sum(x)})

# Original ESs vs. ours (manual mapping)
#- Bacteroides     # ES4
#- Firmicutes      # ES5 (first 10 genera are from Firmicutes and they account for 89% of the ES)
#- Prevotella      # ES3
#- Bifidobacterium # ES1
#- Escherichia     # ES2
# table(mapTaxonomy(TSE, head(names(rev(sort(H[5,]))), 10), to="Phylum"))

########################


