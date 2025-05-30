library(cowplot)
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

# First run to get family-ARG associations
Family_ARG_associations <- readRDS("Family_ARG_associations.rds")

# TSE <- readRDS("../data/TSE.rds")

source("R_functions.R")
source("FigureES.R") # Needed to calculate sigfams 
source("pcoas.R") # Metaphlan3 PcoA
source("enrichments.R") # Plots Metaphlan3 enrichment

# Metaphlan4 PCoA

# Calculate explained variance
e <- attr(reducedDim(altExp(TSE, "metaphlan4"), "PCoA-met4"), "eig")
rel_eig <- e / sum(e[e > 0])

colData(altExp(TSE, "metaphlan4")) <- colData(TSE)

colData(TSE)$Top_fam_meta4_clean <- gsub("^f__", "", colData(TSE)$Top_fam_meta4)
colData(altExp(TSE, "metaphlan4")) <- colData(TSE)

colData(TSE)$Top_fam_meta4_clean <- colData(TSE)$Top_fam_meta4_clean %>% as.factor()

# Extract the ordered family names from TSE
ordered_families <- sort(unique(colData(TSE)$Top_fam_meta4_clean))



###########
pal <- pal.family()

p <- scater::plotReducedDim(altExp(TSE, "metaphlan4"), "PCoA-met4",
                            colour_by = "Top_fam_meta4_clean",
                            point_alpha=0.5, # see ?"scater-plot-args"
                            point_size=3) +   # ?"scater-plot-args"
  scale_color_manual(values = pal, "Dominant Family") +
  # Add explained variance for each axis
  labs(x = paste("PCoA 1 (", round(100 * rel_eig[[1]], 1), "%", ")", sep = ""),
       y = paste("PCoA 2 (", round(100 * rel_eig[[2]], 1), "%", ")", sep = ""))

theme_set(theme_bw(20))
base.size <-20

metaphlan4.plot <- p  +
  theme_classic(base.size) +
  theme(legend.position=c(0.87, 0.17)) +
  coord_equal(clip = "on")  +
  theme(plot.margin  = unit(c(1, 1, 1, 1), "pt"),
        plot.title = element_text(size = base.size),
        axis.text.x = element_text(size = base.size),
        axis.text.y = element_text(size = base.size),
        axis.title.x = element_text(size = base.size),
        axis.title.y = element_text(size = base.size),
        legend.text = element_text(size=0.9*base.size),
        legend.title = element_text(size=1*base.size)				  
  ) +
  guides(color = guide_legend(override.aes = list(size=5))) 				    


# Calculate Top ARG enrichments by family abundance bin
nbreaks <- 5
sigfam <- setdiff(levels(TSE$Top_fam_meta4), "Other")

dfe <- calculateEnrichment(altExp(TSE, "Family-meta4", withColData=TRUE), varname="ARG_burd", assay.type="relabundance", features=sigfam, nbreaks=nbreaks) 

dfe$Level <-  gsub("^f__", "", dfe$Level)



# Set dfe$Level as a factor with levels matching Top_fam_meta4_clean order
dfe$Level <- factor(dfe$Level, levels = ordered_families, ordered = TRUE)


# Get the ordered unique levels of Top_fam_meta4_clean
ordered_families <- levels(dfe$Level)  # These should be the same as Top_fam_meta4_clean

dfe_clean <- dfe[dfe$Level%in%c("Bacteroidaceae","Lachnospiraceae", "Prevotellaceae", "Bifidobacteriaceae", "Enterobacteriaceae", "Rikenellaceae"),]



# High-ARG enrichments by bacterial family
fig3enr_meta4 <- ggplot(dfe_clean, aes(x=Bin, y=value, group=Level, color=Level)) +
  geom_line(size=3) +
  geom_point(size=5) +    
  scale_y_continuous(label=scales::percent, limits=range(dfe$value)) +
  geom_hline(yintercept=0) +
  scale_color_manual(values = pal) +
  labs(x="Abundance bin", y="ARG\nEnrichment (%)", title="High ARG load") +
  facet_wrap(Level ~ ., ncol=1) +
  theme(strip.text.x = element_text(size = 14,)) +
  theme(legend.position="none", plot.title = element_text(size = 20))


figS7_meta4 <- cowplot::plot_grid(metaphlan.plot, fig3enr, metaphlan4.plot, fig3enr_meta4, ncol=2, rel_widths=c(13,3), labels="auto", label_size=30)
s <- 800
library(Cairo)
CairoJPEG("../RESULTS/R_results/FigureS7_meta4.jpg", width=2*s, height=2.8*s)
print(figS7_meta4)
dev.off()

