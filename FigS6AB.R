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
library(circlize)
library(RColorBrewer)
theme_set(theme_bw(20))

# Load functions
source("R_functions.R")

# Read the data (ordination pre-calculated)
# First run to get family-ARG associations
TSE <- readRDS(file="../data/TSE.rds")

# PCoA with ARG profiles; coloring by dominant ARG 
tse <- altExp(TSE, "arg", withColData=TRUE)
tse$Top_ARGgene <- factor(tse$Top_ARGgene, levels=c("tet(O)", "tet(Q)", "tet(W)", "cfxA6", "erm(B)", "Other"))
tse$Top_ARGgene <- factor(tse$Top_ARGgene, levels=names(sort(table(tse$Top_ARGgene))))

# Swap PC2
# reducedDim(tse, "PCoA_ARG")[,2] <- -reducedDim(tse, "PCoA_ARG")[,2]

p <- scater::plotReducedDim(tse, "PCoA_ARG",
			    point_size=2,    # ?"scater-plot-args"
			    point_alpha=0.4,
                            colour_by = "Top_ARGgene")			    

# Calculate explained variance
e <- attr(reducedDim(tse, "PCoA_ARG"), "eig")
rel_eig <- e / sum(e[e > 0])
# Add explained variance for each axis
p <- p + labs(x = paste("PCoA 1 (", round(100 * rel_eig[[1]], 1), "%", ")", sep = ""),
              y = paste("PCoA 2 (", round(100 * rel_eig[[2]], 1), "%", ")", sep = ""))

arg.palette <- get.arg.palette()

names(dic) <- dic <- names(arg.palette)
dic <- gsub("Macrolide, Lincosamide, Streptogramin B", "MLSB", dic)
dic <- gsub("Macrolide, Aminoglycoside, Tetracycline, Quinolone, Amphenicol, Rifamycin", "MATQAR", dic)
arg.palette2 <- arg.palette
names(arg.palette2) <- dic[names(arg.palette2)]
pal <- get.arg.gene.palette()

dic.pal <- c("tet(O)"="Tetracyline: tet(O)",
         "tet(Q)"="Tetracyline: tet(Q)",
	 "tet(W)"="Tetracyline: tet(W)",
	 "cfxA6"="Beta-lactam: cfxA6",
	 "erm(B)"="MLSB: erm(B)",
	 "Other"="Other")

base.size <- 20

# -----------------------------------------------------------------------

ARG.plot <- p  + 
  guides(color = guide_legend(ncol = 1)) + 
  scale_color_manual(values = pal,
                     labels=dic.pal[names(pal)],
                     name="Dominant ARG") +
  theme_classic(20) +
  # labs(title="Resistome: population variation") +
  coord_equal()  +
          theme(
	     plot.margin  = unit(c(1, 1, 1, 1), "pt"),
	     legend.position=c(0.85, 0.85),
	     # axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
	     axis.text.y = element_text(size=base.size),
	     axis.text.x = element_text(size=base.size),	     
	     axis.title.y = element_text(size=base.size),
	     axis.title.x = element_text(size=base.size),
	     legend.text = element_text(size=0.9*base.size),
	     legend.title = element_text(size=base.size)
	     ) +
	     guides(color = guide_legend(override.aes = list(size=5))) 

tse <- altExp(TSE, "Class_top", withColData=TRUE)
tse <- altExp(tse, "ClassPrevalent")

# -----------------------------------------------------------------------

dfm <- melt(assay(tse, "relabundance"))
dfm$Var1 <- dic[as.character(dfm$Var1)]
dfm$Var1 <- factor(dfm$Var1, level=names(sort(sapply(split(dfm$value, dfm$Var1), mean))))
dfm <- dfm[!is.na(dfm$Var1),]
ARG.jitterplot <- ggplot(dfm, aes(x=Var1, y=value, color=Var1)) +
                    geom_jitter(size=1, alpha=0.5) +
		    scale_y_continuous(label=scales::percent) +
		    scale_color_manual(
                      values = arg.palette2[levels(dfm$Var1)]
                      #labels = dic[match(levels(dfm$Var1), dic)],
		      ) +
		    coord_flip() +
		    labs(x="", y="Relative abundance (%)") +
             theme(# axis.text.x = element_blank(),
	           #axis.ticks.x = element_blank(),
	           text = element_text(size=15), 
		   legend.position="none",
	     	   axis.text.y = element_text(size=base.size),
	     	   axis.text.x = element_text(size=base.size),		   
	     	   axis.title.y = element_text(size=base.size),
	     	   axis.title.x = element_text(size=base.size)
	     )  

# -----------------------------------------------------------------------

labsize <- 30
fig <-
  cowplot::plot_grid(
    ARG.plot,
    ARG.jitterplot,        
    nrow = 1,
    rel_widths=c(2,1),
    labels = "auto",
    label_size=labsize
    ) 
pdf("../RESULTS/R_results/FigS6AB.pdf", width=24, height=10)
print(fig)
dev.off()

