
## Stuff to add taxids to blast results based on accessions.
### https://bioinf.shenwei.me/csvtk/
### https://bioinf.shenwei.me/taxonkit/tutorial/#add-taxonomy-information-to-blast-result
#### http://www.phytools.org/eqg/Exercise_3.2/
library(ape)
library(ggtree)
library(tidyverse)
library(TreeSummarizedExperiment)
library(mia)

# Read in Blast output of Resinfer database ARGs against the NCBI nucleotide database
# The blast output only includes rows with ARGs which were found in the FINRISK cohort samples,
# not the whole Resfinder database

resfinder_species <-
  read.csv(
    "../RESULTS/final_output_files/resfinder_results_withspecies_only_found_genes.out",
    sep = "\t"
  )

# Blast done with outfmt 6 and custom output, last one is staxids, 
# which has the ncbi tax id. 
# Species column was added using the staxid and following these tutorials
### https://bioinf.shenwei.me/csvtk/
### https://bioinf.shenwei.me/taxonkit/tutorial/#add-taxonomy-information-to-blast-result

# Add column names based on blast oufmt 6 custom output
colnames(resfinder_species) <-
  c(
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "staxids",
    "Species"
  )

# Remove spaces in the 'Species' column
resfinder_species$Species <- gsub(" ", "_", resfinder_species$Species) 

# Extract Genus information
resfinder_species$Genus <- gsub("_.*","", resfinder_species$Species)
resfinder_species$Genus <- gsub("[", "", resfinder_species$Genus, fixed = TRUE)
resfinder_species$Genus <- gsub("]", "", resfinder_species$Genus, fixed = TRUE)

# Read in Metaphlan's phylogenetic tree
tree <-
  read.tree("../RESULTS/final_output_files/mpa_v30_CHOCOPhlAn_201901_species_tree.nwk")

# Clean up taxa names in the tree
tree$tip.label <- gsub(".__", "", tree$tip.label)
tree$tip.label 

# Create a taxa table from the tree tip labels
metaphlan_tree_tax <-
  read.table(text = gsub("|", ",", tree$tip.label, fixed = TRUE),
             sep = ",")

# Add column names to the taxa table
taxtab <- metaphlan_tree_tax
colnames(taxtab) <-
  c("ID" ,
    "Domain",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species")

# Set tree tip labels to Genus for better matching
tree$tip.label <- taxtab$Genus

################
# Define specific genera for genes of interest

# CfxA6 Prevotella, parabacteroides, Bacteroides
cfxa6_genus<- taxtab[taxtab$Genus %in% resfinder_species[resfinder_species$qseqid %in% ("cfxA6_1_GQ342996"),]$Genus,]
cfxa6_genus$Gene <- "cfxA6"
# tetW - Faecalibacterium_prausnitzii, C. difficile, Bifidobacterium, Lactobacillus
tetW_genus<-taxtab[taxtab$Genus %in% resfinder_species[resfinder_species$qseqid %in% ("tet(W)_5_AJ427422"),]$Genus,]
tetW_genus$Gene <- "tet(W)"
# tetq - Bacteroides dominant, Prevotella, Parabacteroides
tetQ_genus<-taxtab[taxtab$Genus %in% resfinder_species[resfinder_species$qseqid %in% ("tet(Q)_1_L33696"),]$Genus,]
tetQ_genus$Gene <- "tet(Q)"
# tetO - Streptococcus, enterococcus, Faecalibacterium C. difficile
tetO_genus<-taxtab[taxtab$Genus %in% resfinder_species[resfinder_species$qseqid %in% ("tet(O)_3_Y07780"),]$Genus,]
tetO_genus$Gene <- "tet(O)"
# ErmB - Prevotella, Bacteroides, Salmonella, Parabacteroides
ermb_genus<-taxtab[taxtab$Genus %in% resfinder_species[resfinder_species$qseqid %in% ("erm(B)_26_AF080450"),]$Genus,]
ermb_genus$Gene <- "erm(B)"
############

#####Combine tables function to merge and remove duplicates#####
combine_tables <- function(tables_list) {
  # Combine the tables in the list
  combined_table <- do.call(rbind, tables_list)
  
  # Remove duplicate rows based on Genus and Gene columns
  unique_table <- unique(combined_table[, c("Genus", "Gene")])
  
  return(unique_table)
}
#############
# Example usage
tables_list <- list(tetO_genus, tetW_genus, tetQ_genus, ermb_genus, cfxa6_genus)

result_table <- combine_tables(tables_list)


####### Process and visualize the data #######

# Subset to only prevalent genera
# All prevalent genera
prevalent_genera <- mergeFeaturesByPrevalence(TSE[rowData(TSE)$Domain=="Bacteria"], rank="Genus", assay.type="relabundance", detection=0.1/100, prevalence=1/100) %>%
  rowData()

# prevalent_tree <- plot_tree_for_genus(prevalent_genera) 

# Subset tree based on prevalent genera
indices_selected_tips = prevalent_genera$Genus[prevalent_genera$Genus %in% tree$tip.label]
# Remove taxa that have unclassified in the name
indices_selected_tips <- indices_selected_tips[grep("unclassified", indices_selected_tips, invert = TRUE)]

# Then to prune or select the corresponding subtree:
pruned_tree = drop.tip(tree, tip = indices_selected_tips)
subtree = keep.tip(tree, tip = indices_selected_tips)

# Load color palette function
source("R_functions.R")
pal <- get.arg.gene.palette()

############### Create a modified data frame for visualization #############

# Combine genes for each genus
 aggregated_data_2 <- result_table %>%
  group_by(Genus) %>%
  summarize(Gene_List = paste(unique(Gene), collapse = ", ")) 

aggregated_data <- result_table %>%
  group_by(Genus) %>%
  summarize(First_gene = paste(unique(Gene))) 


# Add how many gene variants are found
aggreagted_all_species <- resfinder_species %>% 
   group_by(Genus) %>%
  summarize(Gene_number = length(unique(qseqid)))

aggregated_data <-merge(aggreagted_all_species, aggregated_data, by = "Genus", all.x = TRUE)
aggregated_data_2$Gene_List[is.na(aggregated_data_2$Gene_List)] <- "Others"
aggregated_data$First_gene[is.na(aggregated_data$First_gene)] <- "Others"

# Assuming aggregated_data is your data frame
unique_genes <- c("Others", "tet(O)", "tet(W)", "erm(B)", "tet(Q)", "cfxA6")

# Create columns for each unique gene
for (gene in unique_genes) {
  aggregated_data_2[[gene]] <- as.character(grepl(gene, aggregated_data_2$Gene_List, fixed = TRUE))
}

# Print the modified data frame
print(aggregated_data_2)


extra <- aggregated_data_2[,4:8] %>% data.frame
row.names(extra) <- aggregated_data_2$Genus

extra$tet.O. <- gsub("TRUE", "tet(O)", extra$tet.O., fixed =TRUE)
extra$tet.O. <- gsub("FALSE", NA, extra$tet.O.)
extra$tet.W. <- gsub("TRUE", "tet(W)", extra$tet.W., fixed =TRUE)
extra$tet.W. <- gsub("FALSE", NA, extra$tet.W.)
extra$erm.B. <- gsub("TRUE", "erm(B)", extra$erm.B., fixed =TRUE)
extra$erm.B. <-gsub("FALSE", NA, extra$erm.B. )
extra$cfxA6 <- gsub("TRUE", "cfxA6", extra$cfxA6, fixed = TRUE)
extra$cfxA6 <- gsub("FALSE", NA, extra$cfxA6)
extra$tet.Q. <- gsub("TRUE", "tet(Q)", extra$tet.Q., fixed =TRUE)
extra$tet.Q. <- gsub("FALSE", NA, extra$tet.Q.)


###########################################################

# Create and visualize the aggregated tree

aggregated_tree <- ggtree(subtree, layout = "circular") %<+% aggregated_data +
  geom_tiplab(aes(label = label), size = 6, offset = 0.1, align = TRUE,  linetype = NULL ) +
  geom_tippoint(aes(size=Gene_number)) +
  scale_color_manual(values = pal)+
  labs(fill = "Gene", size = "nr. of ARG variants")


# Create a heatmap to show gene presence for each genus

fig.tree <- gheatmap(aggregated_tree, extra, offset = 0,  
         colnames_position="top",
	 width = 0.14,
         colnames_angle=90, colnames_offset_y = 0.1, 
         hjust=0,
	 font.size=3,
	 colnames = FALSE) +
  scale_fill_manual(values = pal, na.value="white") +
  theme(legend.text = element_text(size=20), legend.title = element_text(size=20),
        legend.position = c(1,1), plot.margin = margin(3,3,3,3, "cm"))



