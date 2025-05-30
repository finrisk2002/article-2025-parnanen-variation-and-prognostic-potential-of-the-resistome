rm(list  = ls())

# To ensure that ordination signs stay the same in different runs etc.
set.seed(252452)

# library(tidyverse)
library(caret)
library(phyloseq)
library(polycor)
library(dplyr)
library(microbiome)
library(readr)
library(stringr)
library(mia)
source("R_functions.R")

# Working dir expected to be FR_metagenomes/R
load("../data/FR_metagenomes.RData") # FR02 - where does this come from?
FR02_orig <- FR02

# -------- TAXONOMIC ABUNDANCES -----------------------------------------

# Metaphlan data preparation
metaphlan_table <- read.table("../RESULTS/final_output_files/metaphlan3_joined_bugs_list_species.tsv",
  header = TRUE, row.names  = 1, check.names  = FALSE, sep = "\t")
metaphlan_table_mod <- as.matrix(metaphlan_table[,2:8089])
row.names(metaphlan_table_mod) <- paste("OTU",1:570, sep = "")
metaphlan_tax <- read.table("../RESULTS/final_output_files/metaphlan_species_tax.tsv", sep  = "\t", header = FALSE)
row.names(metaphlan_tax) <- paste("OTU",1:570, sep = "")
metaphlan_tax <- apply(metaphlan_tax, 2, function(y) (gsub(".__", "", y)))
taxtab <- metaphlan_tax
colnames(taxtab) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# ---------- ARG DATA ----------------------------------------------------------

# Prepare data
ARG_bt_res_mod <- ARG_bt_res[,2:8090]
tmp <- read.csv("../RESULTS/final_output_files/resfinder_phenotypes.txt", header = TRUE, sep = "\t", fill  = TRUE)
tmp$V1 <- tmp$Gene_accession.no
ARG_tax$Seq_id <- row.names(ARG_tax)
tmp2 <- merge(tmp, ARG_tax, by  = "V1")
tmp2<-tmp2[!tmp2$Seq_id %>% duplicated(),]
# Add Seq_id as the row name
rownames(tmp2)<-tmp2$Seq_id

# Create gene names from gene accessions by removing everything after the first "_"
tmp2$Gene<-gsub("(.+?)(\\_.*)", "\\1", tmp2$Gene_accession.no.)
tmp2 <- dplyr::select(tmp2, Class, Gene, Gene_accession.no.)


# -------- PHENODATA -----------------------------------------

# Phenodata preparation
FR02_2 <- read.table("../data/Metagenomics_FR0207_endpoints_FU19_2022-05-10.txt", sep  = "\t", header = TRUE)
FR02_2 <- dplyr::filter(FR02_2, PROJECT == "FR02")

# Remove the samples which don't have metagenome data
FR02_2 <- FR02_2[!(FR02_2$Barcode == ""),]
rownames(FR02_2) <- FR02_2$Barcode

# Add population density
df <- data.frame(readr::read_csv("../RESULTS/final_output_files/FR02_pop_dens_df.csv", show_col_types=FALSE))
FR02 <- dplyr::full_join(df, FR02, by = "FID")

# Remove the samples which don't have metagenome data
FR02_2 <- FR02_2[!(FR02_2$Barcode == ""),]
row.names(FR02_2) <- FR02_2$Barcode

# Add hospital exposure K27 
FR02_K27 <- read.table("../data/Metagenomics_FR02_phenotypes_2025-04-24.txt.gz", sep  = "\t", header = TRUE)
# 
FR02_K27 <- FR02_K27 %>% dplyr::select(matches(c("K27", "FID")))
#  

FR02_K27 <- FR02_K27[!(FR02_K27$FID == ""),] 
FR02 <- dplyr::full_join(FR02_K27 , FR02, by = "FID")
FR02$K27_bin <- FR02$K27>1
# # 

# Add healthy food score score

df <-  data.frame(readr::read_csv("../RESULTS/final_output_files/hfc_score_barcode.csv", show_col_types=FALSE))

FR02 <- dplyr::inner_join(df, FR02, by = "Barcode") 

inds <- match(rownames(FR02_2), FR02$Barcode)
tmp <- FR02[inds, !colnames(FR02) %in% c(colnames(FR02_2))]
rownames(tmp) <- FR02$Barcode[inds]

# Join
FR02 <- cbind(tmp, FR02_2)

# Check the number of samples that match
length(intersect(FR02$Barcode, colnames(ARG_bt_res)))

# Check the samples missing from FR02 metadata
# setdiff(colnames(ARG_bt_res), FR02$Barcode)

# Filter variables with near zero variance, high nr of NAs and extra columns
nzv <- nearZeroVar(FR02, saveMetrics = TRUE )

FR02_select <- dplyr::select(FR02, all_of(colnames(FR02[, !nzv$nzv])))
FR02_select <- dplyr::select(FR02_select, where(is.numeric))
FR02_select <- dplyr::select(FR02_select, !starts_with("RX"))
FR02_select <- dplyr::select(FR02_select, !contains(c("AGEDIFF", "AGE", "ika", "naiset", "miehet",
 "INCIDENT", "DETECTABLE", "FR02", "Matrix.2D.tube.barcode")))

# Remove high NAs
na_count <- sapply(FR02_select, function(y) sum(length(which(is.na(y)))))
high_NAs <- colnames(FR02_select[,na_count>500])

FR02_select <- dplyr::select(FR02_select, !ends_with(high_NAs))
FR02_select$BL_AGE <- FR02$BL_AGE
FR02_select$K27_bin <- FR02$K27_bin
FR02_select$HFC_score <- FR02$HFC_score
FR02_select$DEATH_AGEDIFF <- FR02$DEATH_AGEDIFF
FR02_select$AB1_SEPSIS_BACT_AGEDIFF <- FR02$AB1_SEPSIS_BACT_AGEDIFF
FR02_select$AB1_SEPSIS_BACT_AGE <- FR02$AB1_SEPSIS_BACT_AGE
FR02_select$DEATH_AGE <- FR02$DEATH_AGE
FR02_select$BL_USE_RX_L <- FR02$BL_USE_RX_L
FR02_select$RX_J01_AGEDIFF <- FR02$RX_J01_AGEDIFF
FR02_select$RX_J01_AGE <- FR02$RX_J01_AGE
FR02_select$CURR_SMOKE <- FR02$CURR_SMOKE
FR02_select$PREVAL_RX_J01_NEVT <- FR02$PREVAL_RX_J01_NEVT
FR02_select$INCIDENT_RX_J01_AGEDIFF <- FR02$INCIDENT_RX_J01_AGEDIFF
FR02_select$INCIDENT_RX_J01 <- FR02$INCIDENT_RX_J01
FR02_select$K_TPKS_level1 <- substr(FR02$K_TPKS, 1,1)
FR02_select$K_TPKS_level1[is.na(FR02_select$K_TPKS_level1)] <- 0
FR02_select$K_VKS_level1 <- substr(FR02$K_VKS, 1,1)
FR02_select$K_VKS <- FR02$K_VKS
FR02_select$K_TPKS <- FR02$K_TPKS
FR02_select$K_VKS_level1[is.na(FR02_select$K_KS_level1)] <- 0
FR02_select$KY100_14 <- FR02$KY100_14
FR02_select$KY100_17 <- FR02$KY100_17
FR02_select$KY100_18 <- FR02$KY100_18
FR02_select$KY100_21 <- FR02$KY100_21
FR02_select$FOOD_RECOMMENDED_CHOICES <- FR02$FOOD_RECOMMENDED_CHOICES
FR02_select$PREVAL_RX_J01X <- FR02$PREVAL_RX_J01X
FR02_select$PREVAL_RX_J01X_NEVT <- FR02$PREVAL_RX_J01X_NEVT


#--- Phenodata aggregation ----------------------------------------------------------------------------------------------------------

# Combine small categories manually

## Here combine poultry 5/6 -> 5
FR02_select <- FR02_select %>% mutate(KY100_21=replace(KY100_21, KY100_21==6, 5))
#> table(FR02$KY100_21)
#   1    2    3    4    5    6 
# 671 1823 2343 1814  243   24 

## Here combine vegetables 1/2 -> 1
## and shift other categories
#> table(FR02$KY100_14)
#   1    2    3    4    5    6 
# 208  607  779 1388 2477 1599 
FR02_select <- FR02_select %>% mutate(KY100_14=replace(KY100_14, KY100_14==1, 2)) %>% mutate(KY100_14=KY100_14-1)



# ------------- ARG stuff --------------------------------------------------------------------------------

# Normalize for gene lengths
genelengths <- read.csv("../R/genelengths.csv")
ARG_bt_res_length_norm <- ARG_bt_res_mod/genelengths[,3]
read_counts <- data.frame(read.table("../RESULTS/final_output_files/read_counts_table", row.names  = 1))
read_counts_subs <- dplyr::filter(read_counts, row.names(read_counts) %in% colnames(ARG_bt_res_length_norm))
df <- merge(colSums(ARG_bt_res_mod), read_counts, by = 0) 
rownames(df)<- df$Row.names

# Add RA
phenodata.arg <- FR02_select
df$RA<- df[,2]/df[,3]
df2 <- merge(df[,3], data.frame(phenodata.arg), by=0)
phenodata.arg$RA <- df2$RA

# Save non length normalized ARG relative abundances
# Not used anywhere
# df <- colSums(ARG_bt_res_mod)/read_counts

# Transform to RPKM
new_sample_data <- merge(as.data.frame(phenodata.arg), read_counts_subs, by = 0)
row.names(new_sample_data) <- new_sample_data$Row.names

# Add age categories and other phenodata
new_sample_data <- new_sample_data %>%
  mutate(
    age_class = cut(
      BL_AGE,
      breaks = seq(20, 80, by = 10),
      labels = c(
        "24-\n30",
        "30-\n40",
        "40-\n50",
        "50-\n60",
        "60-\n70",
        "70-\n75"
      ),
      include.lowest = TRUE
    )) %>%
  mutate(across('ALUE', str_replace, '2', 'Oulu')) %>%
  mutate(across('ALUE', str_replace, '3', 'Kuopio')) %>%    
  mutate(across('ALUE', str_replace, '4', 'Turku')) %>%
  mutate(across('ALUE', str_replace, '5', 'Helsinki')) %>%
  mutate(across('ALUE', str_replace, '6', 'Karelia')) %>%  
  mutate(across('ALUE', str_replace, '7', 'Lapland')) %>%
  mutate(Region = factor(ALUE))

# Scale pop density
new_sample_data$vaesto_scaled <- as.vector(scale(new_sample_data$vaesto))

# ------------------------------------------------------------------------------

# Make species TreeSE
library(TreeSummarizedExperiment)
# Create Taxa TreeSE; ensure that the features and samples are shared across data elements
common.samples <- intersect(colnames(metaphlan_table_mod), rownames(new_sample_data))
common.feat <- intersect(rownames(metaphlan_table_mod), rownames(taxtab))
TAX <- TreeSummarizedExperiment(assay=SimpleList(counts=metaphlan_table_mod[common.feat, common.samples]),
                                #colData=DataFrame(FR02[common.samples, ]),
                                colData=DataFrame(new_sample_data[common.samples, ]),				
				rowData=DataFrame(taxtab[common.feat, ]))

# ---------------------------------------------------------------------------------------------------------------

# Make ARG TreeSE
common.samples <- intersect(colnames(ARG_bt_res_length_norm), rownames(new_sample_data))
common.feat <- intersect(rownames(ARG_bt_res_length_norm), rownames(tmp2))
ARG <- TreeSummarizedExperiment(
                                assay=SimpleList(counts=as.matrix(ARG_bt_res_length_norm)[common.feat, common.samples]),
                                colData=DataFrame(new_sample_data[common.samples, ]),
				rowData=DataFrame(as.matrix(tmp2)[common.feat, ]))

# Add sample-wise sums of the ARG signal
ARG$SUM <- colSums(assay(ARG, "counts"))
ARG$SUM_norm <- ARG$SUM/ARG$V2*1e6*1e3

# ARG burden
bins <-
  quantile(ARG$SUM_norm, probs = (c(0, 0.9, 1)))  
ARG_burd <-
  cut(
    ARG$SUM_norm,
    breaks  = bins,
    include.lowest  = TRUE,
    label = c("low/moderate", "high")
  )
ARG$ARG_burd <- ARG_burd


## -------- Filtering ----------------

# Remove empty features
TAX <- TAX[mia::getPrevalence(TAX) > 0,]
ARG <- ARG[mia::getPrevalence(ARG) > 0,]

# Subset samples
ARG <- ARG[, ARG$SUM_norm>0]
TAX <- TAX[, colSums(assay(TAX))>0]

# Ensure samples are in the same order
common.samples <- intersect(colnames(TAX), colnames(ARG))
TAX <- TAX[, common.samples]
ARG <- ARG[, common.samples]

## ----  Transformations ----------------

# Add transformations for the samples based on "counts" assay in the taxa experiment
TAX <- transformAssay(TAX, assay.type="counts", method="relabundance", MARGIN="samples")

# Add relative abundances to ARG data (for barplots; not for statistical analyses)
ARG <- transformAssay(ARG, assay.type="counts", method="relabundance", MARGIN="samples")

# ----------------------Agglomerate taxa--------------------------------------------------

# Add Family level abundance table
altExp(TAX, "Family")  <- mergeFeaturesByRank(TAX, rank="Family")
altExp(TAX, "Genus")   <- mergeFeaturesByRank(TAX, rank="Genus")
altExp(TAX, "Species") <- mergeFeaturesByRank(TAX, rank="Species")

# Merge rare genera instead
altExp(TAX, "GenusPrevalent") <- mergeFeaturesByPrevalence(TAX[rowData(TAX)$Domain=="Bacteria"], rank="Genus", assay.type="relabundance", detection=0.1/100, prevalence=1/100)

# Old R version
#altExp(ARG, "ClassPrevalent") <- mergeFeaturesByPrevalence(ARG, f="Class", assay.type="relabundance", detection=0.1/100, prevalence=1/100)
# New R version
altExp(ARG, "ClassPrevalent") <- agglomerateByPrevalence(ARG, rank="Class", assay.type="relabundance", detection=0.1/100, prevalence=1/100)

# Add information on dominant features to colData
# - this differs for very few samples (3-6 depending on random seed) from that obtained earlier with phyloseq
# - and those are all cases where there are several top features with an identical count i.e. the top feature identification picks one of them
## Merge features by Family in each sample

## Add dominant feature per sample to colData (do this for the merged TAX data to
## keep it clear where it belongs)
## Note that a few samples have multiple equally dominant features; let us pick a random selection in those cases
doms <- addPerSampleDominantFeatures(altExp(TAX, "Family"), name="top")$top
tops <- unname(sapply(doms, function (x) {sample(unlist(x),1)}))
# Identify the top taxa
top <- microbiome::top(tops, n=8)
# Group the rest into the "Other" category
tops[!tops %in% names(top)] <- "Other"
# Store in the main colData 
TAX$Top_fam <- factor(tops)
# Aggregate data by combining other than top families
rowData(TAX)$Family_aggregated <- rowData(TAX)$Family
rowData(TAX)$Family_aggregated[!rowData(TAX)$Family %in% names(top)] <- "Other"
rowData(TAX)$Family_aggregated <- factor(rowData(TAX)$Family_aggregated)
altExp(TAX, "Family_aggregated") <- mergeFeatures(TAX, "Family_aggregated")
altExp(TAX, "FamilyPrevalent") <- mergeFeaturesByPrevalence(TAX[rowData(TAX)$Domain=="Bacteria"], rank="Family", assay.type="relabundance", detection=0.1/100, prevalence=1/100)

# -----------------------------------------------------------------------------------------------

## Add dominant Genus per sample to colData (do this for the merged TAX data to
## keep it clear where it belongs)
## Note that a few samples have multiple equally dominant features; let us pick a random selection in those cases
doms <- addPerSampleDominantFeatures(altExp(TAX, "Genus"), name="top")$top
tops <- unname(sapply(doms, function (x) {sample(unlist(x),1)}))
# Identify the top taxa
top <- microbiome::top(tops, n=20) # Top-10 explains 95%; Top-20 explains 99.8%
# Group the rest into the "Other" category
tops[!tops %in% names(top)] <- "Other"
# Store in the main colData 
TAX$Top_genus <- factor(tops)
#cumsum(rev(sort(table(TAX$Top_genus)))/ncol(TAX))[[20]]

# Aggregate data by combining other than top taxa
rowData(TAX)$Genus_aggregated   <- rowData(TAX)$Genus
rowData(TAX)$Genus_aggregated[!rowData(TAX)$Genus %in% names(top)] <- "Other"
rowData(TAX)$Genus_aggregated   <- factor(rowData(TAX)$Genus_aggregated)
altExp(TAX, "Genus_aggregated") <- mergeFeatures(TAX, "Genus_aggregated")

## -------- Agglomerate ARGs ----------------

# Add information on dominant features to colData
# - this differs for very few samples (3-6 depending on random seed) from that obtained earlier with phyloseq
# - and those are all cases where there are several top features with an identical count
# - i.e. the top feature identification picks one of them

## Merge features by Gene in each sample
altExp(ARG, "arg_merged") <- mergeFeatures(ARG, "Gene") 

## Get dominant feature per sample 
## Note that a few samples have multiple equally dominant features; let us pick a random selection in those cases
doms <- addPerSampleDominantFeatures(altExp(ARG, "arg_merged"), name="top")$top

# Some samples may have multiple features that are equally dominant; in those cases select just one of them at random
tops <- unname(sapply(doms, function (x) {sample(unlist(x),1)}))

# Identify the top ARG features among all
top <- microbiome::top(tops, n=5)

# Group the rest into the "Other" category
tops[!tops %in% names(top)] <- "Other"

# Also store in the main colData 
ARG$Top_ARGgene <- factor(tops)

# ===================================================================0

# Agglomerate ARG data to the Class level (this handles the original and relabundances on one go)
altExp(ARG, "arg_class") <- mergeFeatures(ARG, f="Class")

## Add dominant feature per sample to colData (do this for the merged data to keep it clear where it belongs)
## Note that a few samples have multiple equally dominant features; let us pick a random selection in those cases
doms <- addPerSampleDominantFeatures(altExp(ARG, "arg_class"), name="top")$top
tops <- doms
# Identify the top features
top <- microbiome::top(doms, n=5)
# Group the rest into the "Other" category
tops[!tops %in% names(top)] <- "Other"
# Store in the main colData 
ARG$Top_class <- factor(tops)
# Aggregate data by combining other than top classes
rowData(ARG)$Class_top <- rowData(ARG)$Class
rowData(ARG)$Class_top[!rowData(ARG)$Class %in% names(top)] <- "Other"
rowData(ARG)$Class_top <- factor(rowData(ARG)$Class_top)
altExp(ARG, "Class_top") <- mergeFeatures(ARG, f="Class_top")

## ----------CLR and log10p -------------------------------------

# Better added after agglomerations since this has to be done per rank-specific altExps

# Add transformations for the samples based on "counts" assay in the taxa experiment
TAX <- transformAssay(TAX, assay.type="counts", method="clr", pseudocount=1, MARGIN="samples")
altExp(TAX, "Family") <- transformAssay(altExp(TAX, "Family"), assay.type="counts", method="clr", pseudocount=1, MARGIN="samples")
altExp(TAX, "Genus") <- transformAssay(altExp(TAX, "Genus"), assay.type="counts", method="clr", pseudocount=1, MARGIN="samples")
altExp(TAX, "Species") <- transformAssay(altExp(TAX, "Species"), assay.type="counts", method="clr", pseudocount=1, MARGIN="samples")


# Add log10p assay
altExp(TAX, "Genus_aggregated") <- transformAssay(altExp(TAX, "Genus_aggregated"), assay_name="counts", method="log10", pseudocount=1)
# Add clr assay
altExp(TAX, "Genus_aggregated") <- transformAssay(altExp(TAX, "Genus_aggregated"), assay_name="counts", method="clr", pseudocount=1)
# Z transform of CLR counts i.e. add CLR-Z assay
altExp(TAX, "Genus_aggregated") <- transformAssay(altExp(TAX, "Genus_aggregated"), MARGIN = "features", assay_name = "clr", method = "z", name="clrz")

## --------- Alpha diversities -------------------

# Add taxonomic shannon diversity
# FIXME: this seems feature diversity rather than Species diversity (which is one rank higher than the data itself)
TAX <- mia::estimateDiversity(TAX, assay.type="counts", index = "shannon", name="Species_diversity")

# Add observed ARG diversity
ARG <- mia::estimateDiversity(ARG, assay.type="counts", index = "shannon", name="ARG_div")

# Add observed ARG richness
ARG <- mia::estimateRichness(ARG, assay.type="counts", index = "observed", name="ARG_obs")


## ---------- Ordinations ------------

# Perform PCoA with Bray-Curtis on relative abundances
TAX <- scater::runMDS(TAX,
                      FUN = vegan::vegdist,
                      method = "bray",
                      name = "PCoA",
		      ncomponents=50,
                      exprs_values = "relabundance")
colnames(reducedDim(TAX)) <- paste0("PC", seq_len(ncol(reducedDim(TAX))))

# Switch (the arbitrary) sign so that the figure is more directly
# comparable with earlier publications
reducedDim(TAX)[,1] <- -reducedDim(TAX)[,1]
reducedDim(TAX)[,2] <- -reducedDim(TAX)[,2]

# Perform PCoA for the ARG data as well
ARG <- scater::runMDS(ARG,
                                     FUN = vegan::vegdist,
                      		     method = "bray",
                      		     name = "PCoA_ARG",
		                     ncomponents=3,				     
                      		     exprs_values = "relabundance")
colnames(reducedDim(ARG)) <- paste0("PC", seq_len(ncol(reducedDim(ARG))))				     
reducedDim(ARG)[,2] <- -reducedDim(ARG)[,2]


## --------- Combine and Store -------------------------------------

## Combine TreeSummarizedExperiment objects:
## Use Taxonomic data as the main object
TSE <- TAX
## -> add ARG data sets an alternative experiments in the main data object
altExp(TSE, "arg") <- removeAltExps(ARG)
## -> also add the main ARG data set as an alternative experiment
altExps(TSE) <- c(altExps(TSE), altExps(ARG))

# Combine TAX and ARG sample data 
colData(TSE) <- dplyr::full_join(as.data.frame(colData(TAX)),
                                 as.data.frame(colData(ARG))) %>%
		       DataFrame

# Remove colData from altExps for clarity
for (e in altExpNames(TSE)) {
  colData(altExp(TSE, e)) <- NULL
}

## --------- Cross-validations ----------------------------

# Define train/test split and store their binary indicators in colData
inds <- logical(ncol(TSE));
inds[seq(5000)] <- TRUE;
TSE$train.sample <- sample(inds)


# -------------- Add community types (CST) and enterosignatures  -------------------------------------

# Requires some manual analytics;
# adds CST in colData(TSE) in the end
# This indicates the determined community type for each sample
# source("communityping.R")
source("nmf.R") # "Enterosignatures" with NMF
# source("dmm.R")

# ------------- Save ------------------------------

# Save the joint data object
saveRDS(TSE, file="../data/TSE.rds")

# -------------------------------------------------------------

# R-4.2.3 installation
# --with-readline=no --with-x=no

#rm(list  = ls())

