library(NMF)
library(mia)
library(reshape2)
library(dplyr)

# Read in MetaPhlAn4 results
TAX_met4 <- importMetaPhlAn("../RESULTS/final_output_files/metaphlan4_rel_ab_merged.txt")

# Read in the TSE object created with CapentryTreeSE.R script
TSE <- readRDS("../data/TSE.rds") # This is the old TSE

# Set colnames
colnames(TSE) <- colData(TSE)$Row.names

# Agglomerate to species
TAX_met4<-agglomerateByRank(TAX_met4, rank="species")

# Identify common cols
common_cols <- intersect(colnames(TAX_met4), colnames(TSE))

# Match colnames
TAX_met4 <- TAX_met4[, match(colnames(TSE), colnames(TAX_met4)), drop = FALSE]

# Transform to relative abundance
TAX_met4 <- transformAssay(TAX_met4, assay.type="counts", method="relabundance", MARGIN="samples")

########## PCoA ################

# Perform PCoA with Bray-Curtis on relative abundances
TAX_met4 <- scater::runMDS(TAX_met4,
                      FUN = vegan::vegdist,
                      method = "bray",
                      name = "PCoA-met4",
                      ncomponents=50,
                      exprs_values = "relabundance")

colnames(reducedDim(TAX_met4)) <- paste0("PC-meta4", seq_len(ncol(reducedDim(TAX_met4))))

# Switch (the arbitrary) sign so that the figure is more directly
# comparable with earlier publications
reducedDim(TAX_met4)[,1] <- -reducedDim(TAX_met4)[,1]
reducedDim(TAX_met4)[,2] <- -reducedDim(TAX_met4)[,2]

######## Save as altExp ##############
# Save as altExp to the TSE
altExp(TSE, "metaphlan4") <- TAX_met4


######## metaphlan4 family alt exp #############
# Save family level altAxp also
altExp(TSE, "Family-meta4")  <- mergeFeaturesByRank(altExp(TSE, "metaphlan4"), rank="family")

# Get dominant features for families
doms <- addPerSampleDominantFeatures(altExp(TSE, "Family-meta4"), name="top")$top

tops <- unname(sapply(doms, function (x) {sample(unlist(x),1)}))

# Identify the top taxa
top <- microbiome::top(tops, n=11) # Use 11 here in stead of 6 since Bifido is 11th

# Group the rest into the "Other" category
tops[!tops %in% names(top)] <- "Other"

# Store in the main colData 
TSE$Top_fam_meta4 <- factor(tops)

########### ES #################
# Use the same prevalence cut off as for metaphlan3 to calculate nmf
altExp(TSE, "metaphlan4GenusPrevalent") <-
  mergeFeaturesByPrevalence(
    altExp(TSE, "metaphlan4")[rowData(altExp(TSE, "metaphlan4"))$kingdom == "k__Bacteria"],
    rank =
      "genus",
    assay.type = "relabundance",
    detection = 0.1 / 100,
    prevalence = 1 / 100
  )

set.seed(33)
vec <- 2:7
x <- t(assay(altExp(TSE, "metaphlan4GenusPrevalent")))

# Normalize the relab 100
x <- sweep(x, 1, rowSums(x), FUN = "/") * 100


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

colnames(W) <- rownames(H) <- paste("nmf-met4-", 1:ncol(W), sep = "")

# Add enterosignature scores to colData
colData(TSE) <- DataFrame(as.data.frame(colData(TSE)) %>% dplyr::select(!starts_with("nmf-met4-"))) # first remove conflicting entries
colData(TSE) <- DataFrame(cbind(colData(TSE), W))

# Add NMF loadings in the metadata slot
metadata(TSE)$NMF_met4_loadings <- H

# Obtain the relative abundance of ES in each sample.
Wrel <- t(apply(W, 1, function (x) {x/sum(x)}))

# Define as primary ES the ES of a sample with the highest relative abundance
# and add this to colData
colData(TSE)$ES_met4_primary <- apply(Wrel, 1, which.max)

########### ES with high prevalence ########
# Use the same prevalence cut off as for metaphlan3 to calculate nmf
altExp(TSE, "metaphlan4GenusVeryPrevalent") <-
  mergeFeaturesByPrevalence(
    altExp(TSE, "metaphlan4")[rowData(altExp(TSE, "metaphlan4"))$kingdom == "k__Bacteria"],
    rank =
      "genus",
    assay.type = "relabundance",
    detection = 10 / 100,
    prevalence = 1 / 100
  )

set.seed(33)
vec <- 2:7
x <- t(assay(altExp(TSE, "metaphlan4GenusVeryPrevalent")))

# Normalize the relab 100
x <- sweep(x, 1, rowSums(x), FUN = "/") * 100


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

colnames(W) <- rownames(H) <- paste("nmf-met4-high", 1:ncol(W), sep = "")

# Add enterosignature scores to colData
colData(TSE) <- DataFrame(as.data.frame(colData(TSE)) %>% dplyr::select(!starts_with("nmf-met4-high"))) # first remove conflicting entries
colData(TSE) <- DataFrame(cbind(colData(TSE), W))

# Add NMF loadings in the metadata slot
metadata(TSE)$NMF_met4_high_loadings <- H

# Obtain the relative abundance of ES in each sample.
Wrel <- t(apply(W, 1, function (x) {x/sum(x)}))

# Define as primary ES the ES of a sample with the highest relative abundance
# and add this to colData
colData(TSE)$ES_met4_high_primary <- apply(Wrel, 1, which.max)

# ------------- Save ------------------------------

# Save the joint data object
saveRDS(TSE, file="../data/TSE.rds")

