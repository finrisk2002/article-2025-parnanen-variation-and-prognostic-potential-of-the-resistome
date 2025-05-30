# Load libraries
library(viridis)
library(mia)
library(miaViz)
library(RColorBrewer)

source("R_functions.R")

# Read data
TSE <- readRDS("../data/TSE.rds")

# Generate a bar plot of ARG class abundance
# (the necessary data groupings are now created in Carpentry already)
# tse <- altExp(TSE, "Class_top", withColData=TRUE)
tse <- altExp(TSE, "ClassPrevalent", withColData=TRUE)

# Define some name conversions
# TODO simplify..
arg.palette <- get.arg.palette()
names(dic) <- dic <- names(arg.palette)
dic <- gsub("Macrolide, Lincosamide, Streptogramin B", "MLSB", dic)
dic <- gsub("Macrolide, Aminoglycoside, Tetracycline, Quinolone, Amphenicol, Rifamycin", "MATQAR", dic)
arg.palette2 <- arg.palette
names(arg.palette2) <- dic[names(arg.palette2)]

base.size <- 30
theme_set(theme_bw(base.size))
ARG.barplot <- plotAbundance(tse,
                   assay.type="relabundance", 
                   order_sample_by = "Tetracycline",
		   order_rank_by="abund") +
       guides(fill = guide_legend(nrow = 7, title="ARG class")) +
       labs(x = paste0("Samples (n=", ncol(TSE), ")"),
            y = "Relative abundance (%)",
            caption  = "")  +
       scale_y_continuous(labels=scales::percent) + 
       theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
	     axis.text.y = element_text(size=0.8*base.size),
	     axis.title.y = element_text(size=base.size),
	     axis.title.x = element_text(size=base.size),
	     legend.text = element_text(size=0.8*base.size),
	     legend.title = element_text(size=base.size),	     	     	     	     
             legend.position = "right"
	     )  +
       scale_fill_manual("",
         values =  arg.palette[rownames(assay(tse, "relabundance"))],
         labels = dic[rownames(assay(tse, "relabundance"))],
	 name="ARG"	
	 ) +
	 theme(text = element_text(size=15))	 

pdf("../RESULTS/R_results/FigS4_ARGbarplot.pdf", width=15, height=8)
print(ARG.barplot)
dev.off()

# ===========================================================================================

# Define a different custom color palette and generate a bar plot of family abundance
pal <- c("#663333","#CC9999", "#CCCFFF", "#CC9966", "#330000", "#6699CC", "#666666", "#CCFFFF", "#000000")
tse <- altExp(TSE, "Family_aggregated", withColData=TRUE)

# TODO: miaViz could pick this from rows directly, too
# tse$Bacteroidaceae <- assay(tse, "relabundance")["Bacteroidaceae", ]
family.barplot <- plotAbundance(tse,
	           rank="Family",
                   assay.type="relabundance", 
                   order_sample_by = "Bacteroidaceae",
		   order_rank_by="abund"
		   ) +
       guides(fill = guide_legend(ncol  = 3)) +
       labs(x = "Samples",
         y = "Relative abundance (%)",
         subtitle = "Bacterial Family",
         caption = "")  +
       scale_y_continuous(labels=scales::percent) + 
       theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
	     legend.text = element_text(face = "italic"),
	     legend.position = "bottom")  +
       scale_fill_manual("", values = pal)

# family.barplot

# ===========================================================================================