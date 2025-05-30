library(mia)
library(miaViz)
library(ggpubr)
library(patchwork)
library(grid)
library(microViz)
library(viridis)
library(multcomp)
library(phyloseq)

theme_set(theme_bw(20))

# Load functions
source("R_functions.R")

# Read the data (ordination pre-calculated)
TSE <- readRDS("../data/TSE.rds") 

# If needed, run first:
# Figure3.R

# Figure 3B: associations between ARG and driver taxa

# df_SP_axes <- cbind(reducedDim(TSE), colData(TSE)) %>%
df_SP_axes <- cbind(reducedDim(TSE), DF) %>%
                as.data.frame %>% 
                dplyr::rename(ARG_load = SUM_norm)

lm(log(ARG_load)~PC1+PC2+PC3, data = df_SP_axes) %>% summary
lm(log(ARG_load)~Bacteroidaceae+Enterobacteriaceae+Bifidobacteriaceae+Prevotellaceae+Lachnospiraceae, data = df_SP_axes) %>% summary

tops <- c("Bacteroidaceae","Enterobacteriaceae","Prevotellaceae","Bifidobacteriaceae","Lachnospiraceae")
# tops <- c("Bacteroidaceae","Enterobacteriaceae",""Prevotellaceae")
plots <- lapply(tops, function (i) plotfun(i, df_SP_axes)); names(plots) <- tops
# This was previously added on Prevotellaceae; not sure if really needed
# scale_color_manual(values = pal, "") + scale_fill_manual(values = pal, "") 

fig3B <- cowplot::plot_grid(
                   align = "h", 
                   plots[["Bacteroidaceae"]],
                   plots[["Enterobacteriaceae"]],
                   plots[["Prevotellaceae"]],
                   nrow=1, 
		   labels = c("c", "d", "e"))

pdf("../RESULTS/R_results/FigED1B.pdf", width=20, height=10)
print(fig3B)
dev.off()


# =============================================================================================

# Figure 3C: RDA plot

# ------------------------------------------------------------------

# This ordination plot is with microViz, which only supports phyloseq
# -> convert to phyloseq and operate with that.
# FIXME: if this is standard dbRDA, then use examples from OMA to rewrite.
ps <- makePhyloseqFromTreeSE(altExp(TSE, "Family", withColData=TRUE))
dd <- as.data.frame(colData(TSE)) %>%
        mutate(MEN=as.numeric(MEN))
sample_data(ps) <- sample_data(dd)

pal <-
  c(
    "#663333",
    "#6699CC",
    "khaki1",
    "#CCCFFF",
    "#330000",
    "#CC9966",
    "#CC9999",
    "khaki4",
    "lightcyan4"
  )


ordplot <- ps %>% 
  ord_calc(
    constraints = c(
      "BMI",
      "BL_AGE",
      "SUM_norm",
      "MEN"
    )
  ) %>%
  microViz::ord_plot(colour = "Top_fam",
                     plot_taxa = 1:8,
                     size = 1,
                     auto_caption = NA) +
  scale_color_manual(values = pal, "")  +
  # scale_fill_manual(values = pal, "") +
  theme(legend.position = c(c(0.1, 0.85)))



# ------------------------------------------------------------------

# Doesnt work for ARG data
# This ordination plot is with microViz, which only supports phyloseq
# -> convert to phyloseq and operate with that.
# FIXME: if this is standard dbRDA, then use examples from OMA to rewrite.
library(phyloseq)
library(microViz)

pal <-
  c(
    "#663333",
    "#6699CC",
    "khaki1",
    "#CCCFFF",
    "#330000",
    "#CC9966",
    "#CC9999",
    "khaki4",
    "lightcyan4"
  )


tse <- altExp(TSE, "arg_class", withColData=TRUE)
colnames(rowData(tse)) <- c("Family", "Class", "OTU") # FAKE; function requires "Rank" names
ps <- makePhyloseqFromTreeSE(tse, assay.type="relabundance")
dd <- as.data.frame(colData(TSE)) %>%
        mutate(MEN=as.numeric(MEN))
sample_data(ps) <- sample_data(dd)
ordplot2 <- ps %>% 
  ord_calc(
    constraints = c(
      "BMI",
      "BL_AGE",
      "SUM_norm",
      "MEN"
    )
  ) %>%
  microViz::ord_plot(colour = "Top_ARGgene",
                     plot_taxa = 1:5,
                     size = 1,
                     auto_caption = NA) +
  scale_color_manual(values = pal, "")  +
  theme(legend.position = c(c(0.1, 0.85)))

fig3C <-
  cowplot::plot_grid(
    ordplot,
    ordplot2,
    ncol = 2,
    labels = "auto"
  )

pdf("../RESULTS/R_results/FigED1C.pdf", width=20, height=10)
print(fig3C)
dev.off()
