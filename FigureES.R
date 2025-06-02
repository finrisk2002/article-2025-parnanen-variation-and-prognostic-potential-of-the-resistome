library(reshape2)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)
library(mia)

dic.nmf <- get.dic.nmf()

select <- dplyr::select

# The most significant families
sigfam <- Family_ARG_associations %>% filter(Variable=="ARG load" & FDR<0.05) %>%
                                      select(Family) %>%
				      unlist() %>%
				      as.character()
				      
# Dominant families
doms <- addPerSampleDominantFeatures(altExp(TSE, "Family"), name="top")$top

# Combine dominant families that are not significant or have small sample size
TSE$Signif_fam <- doms
TSE$Signif_fam[!TSE$Signif_fam %in% sigfam | TSE$Signif_fam %in% names(which(sort(table(doms)) < 50))] <- "Other"
TSE$Signif_fam <- droplevels(factor(TSE$Signif_fam))

# Order levels manually so that they match the Enterosignatures roughly
TSE$Signif_fam <- factor(TSE$Signif_fam, levels=c("Bacteroidaceae", "Lachnospiraceae", "Prevotellaceae", "Bifidobacteriaceae", "Enterobacteriaceae", "Rikenellaceae", "Other"))



# Loadings on NMF components (i.e. enterosignature profiles)

# Pick the NMF loadings and scores
H <- metadata(TSE)$NMF_loadings
W <- as.data.frame(colData(TSE)) %>% select(matches(c("nmf1", "nmf2", "nmf3", "nmf4", "nmf5")))
rownames(W) <- colnames(TSE)

# Store ES to altExp slot
m <- as.matrix(t(W))
rownames(m) <- dic.nmf[rownames(m)]
altExp(TSE, "ES") <- TreeSummarizedExperiment(assay=SimpleList(signal=m))
TSE$TopES <- factor(unname(dic.nmf[colnames(W)[apply(W, 1, which.max)]]), levels=c("ES-Bact", "ES-Firm", "ES-Prev", "ES-Bifi", "ES-Esch"))

# --------------------------------------------------------------

# Very slow!
selected.vars <- c("SUM_norm", "ARG_div", "Species_diversity")
library(stringr)

  # ARG classes vs. enterosignatures
  cors3 <- cor(t(assay(altExp(TSE, "Class_top"), "counts")), W, method="kendall")

  # Calculate all NMF associations at one go to deal with multiple testing properly
  colData(TSE) <- cbind(colData(TSE), t(assay(altExp(TSE, "Class_top"), "counts")))
  
  top.classes <- setdiff(unique(rownames(altExp(TSE, "Class_top"))), "Other")

  # Quick hack: this is needed for getExperimentCrossAssociation to work
  # TODO: open issue 
  colnames(TSE) <- colData(TSE)$Row.names
  
  #res3 <- getExperimentCrossAssociation(TSE,
  #                                   # Top ARG classes, count assay
  #                                   colData_variable1 = c(top.classes, selected.vars),#
  #				     # ... correlated with colData columns starting with "nmf"
  #                                     colData_variable2 = stringr::str_subset(names(colData(TSE)), "^nmf"),# 
  #				     # use Kendall's tau for robustness
  #                                    method = "kendall",
  #				     mode = "table",
  #				     p_adj_method="fdr",
  #				     test_significance = TRUE)

  res3 <- getCrossAssociation(TSE,
                                     # Top ARG classes, count assay
                                     col.var1 = c(top.classes, selected.vars),
				     # ... correlated with colData columns starting with "nmf"
                                     col.var2 = stringr::str_subset(names(colData(TSE)), "^nmf[1-5]$"),
				     # use Kendall's tau for robustness
                                     method = "kendall",
				     mode = "table",
				     p.adj.method="fdr",
				     test.signif = TRUE)

  colnames(res3) <- str_replace(colnames(res3), "Var1", "Variable")
  colnames(res3) <- str_replace(colnames(res3), "Var2", "Component")
  colnames(res3) <- str_replace(colnames(res3), "cor", "Tau")
  colnames(res3) <- str_replace(colnames(res3), "pval", "P")
  colnames(res3) <- str_replace(colnames(res3), "p_adj", "Padj")

  res <- res3 %>% filter(!Variable %in% selected.vars)
  res2 <- res3 %>% filter(Variable %in% selected.vars)


# ES profile across bacteria (normalized to 1 per ES)
x <- t(apply(H, 1, function (x) {x/sum(x)}))
dfm <- melt(x)
colnames(dfm) <- c("ES", "Genus", "score")
dfm <- dfm %>% filter(!Genus=="Other")
m <-mapTaxonomy(TSE, taxa = as.character(dfm$Genus), to="Phylum")
dfm$Phylum <- unname(m)
# Rearrange rows
dfm <- dfm %>% arrange(Phylum) %>% mutate(Genus=factor(Genus, levels=unique(Genus)))

# Alternative ES visualization
p1b <- ggplot(dfm, aes(x=Genus, y=score, fill=Phylum)) +
       geom_bar(stat="identity", position="dodge") +
       geom_text(data=(dfm%>%filter(score>0.05)), aes(label=Genus),
         check_overlap=TRUE, angle=45, vjust=0.5, y=0.05, hjust=0.05) +
       theme(axis.text.x = element_text(angle = 0)) +
       scale_y_continuous(limits=c(0,max(dfm$score)), label=scales::percent) + 
       facet_grid(. ~ ES) 

# ------------------------------------------------------------

library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
breaks <- seq(-0.15, 0.15, 0.03)
# col_fun <- colorRamp2(c(min(breaks), 0, max(breaks)), c("blue", "white", "red"))
col_fun <- colorRamp2(breaks, rev(brewer.pal(length(breaks), "RdBu")))

# ES profile across bacteria (normalized to 1 per ES)
# Remove genera whose relative importance in all ES remains under <x%
# ESprofiles <- Hscaled[, !colMeans(Hscaled<5/100)==1]
weights <- t(apply(H, 2, function (x) {x/sum(x)}))
# scores <- t(apply(W, 2, function (x) {x/sum(x)}))

rd <- rowData(altExp(TSE, "GenusPrevalent"))[rownames(weights), ]

# Add top Family info from TSE main
rd$Signif_fam <- rd$Family
# Indicate families that match with the main enterosignatures
rd$Signif_fam[!rd$Family %in% setdiff(unique(TSE$Signif_fam), "Rikenellaceae")] <- "Other"
fam2 <- tidyr::replace_na(as.character(rd$Signif_fam), "Other")
rd$Signif_fam <- factor(rd$Signif_fam)
fam <- as.character(rd$Signif_fam)
fam[is.na(fam)] <- "Other"

# Define row annotation
pal <- pal.family()
row_ha <- rowAnnotation(Family=fam, na_col="gray", col=list(Family=pal[unique(fam)]), show_annotation_name=TRUE, show_legend=TRUE)

p1 <- Heatmap(weights,
  name="Weight",
  row_names_gp = gpar(fontsize = 12),
  row_title_gp = gpar(fontsize = 12),  
  col = rev(brewer.pal(length(breaks), "Greys")),
  show_row_dend=FALSE,
  cluster_columns=FALSE,
  row_title="Genus",
  heatmap_height = unit(nrow(weights)/3, "cm"),
  show_column_names = FALSE, 
  bottom_annotation = HeatmapAnnotation(
          text=anno_text(dic.nmf[colnames(weights)], rot=90, location=0, just = "left", show_name=FALSE)
  	),
  right_annotation = row_ha
  )

# ------------------------------------------------------------

# Correlations and their p-values
C <- reshape::cast(res, Variable ~ Component, value="Tau")  %>% tibble::column_to_rownames("Variable")
P <- reshape::cast(res, Variable ~ Component, value="Padj") %>% tibble::column_to_rownames("Variable")

dic <- rownames(C); names(dic) <- dic
dic <- gsub("Macrolide..Lincosamide..Streptogramin.B", "MLSB", dic)
dic <- gsub("Beta.lactam", "Beta-lactam", dic)
dic <- gsub("Macrolide..Aminoglycoside..Tetracycline..Quinolone..Amphenicol..Rifamycin", "MATQAR", dic)
dic.args <- dic

Csort <- C[c(5,2,4,1,3),]
p2 <- Heatmap(Csort, # Sort as in Fig3b
  name="Tau",
  col = col_fun,
  cluster_columns=FALSE,
  cluster_rows=FALSE,      
  show_row_dend=FALSE,
  show_column_dend=FALSE,
  show_column_names = FALSE,   
  row_title="ARG class",
  row_labels=dic[rownames(Csort)],  
  heatmap_height = unit(0.5*nrow(Csort), "cm"),
    cell_fun = function(j, i, x, y, w, h, fill) {
      gb = textGrob("*")
      gb_w = convertWidth(grobWidth(gb), "mm")
      gb_h = convertHeight(grobHeight(gb), "mm")
      if(P[i, j] < 0.05) {    
        grid.text("*", x, y - gb_h*0.8 + gb_w*0, gp=gpar(fontsize=20, col="black"))	
      } 
    }
)

Csort <- C[c(5,2,4,1,3),]
colnames(Csort) <- dic.nmf[colnames(Csort)]
p2b <- Heatmap(Csort, # Sort as in Fig3b
  name="Tau",
  col = col_fun,
  cluster_columns=FALSE,
  cluster_rows=FALSE,      
  show_row_dend=FALSE,
  show_column_dend=FALSE,
  show_column_names = TRUE,   
  row_title="ARG class",
  row_labels=dic[rownames(Csort)],
  #col_labels=unname(),    
  # heatmap_height = unit(0.5*nrow(Csort), "cm"),
  heatmap_height = unit(0.9*nrow(Csort), "cm"),  
    cell_fun = function(j, i, x, y, w, h, fill) {
      gb = textGrob("*")
      gb_w = convertWidth(grobWidth(gb), "mm")
      gb_h = convertHeight(grobHeight(gb), "mm")
      if(P[i, j] < 0.05) {    
        grid.text("*", x, y - gb_h*0.8 + gb_w*0, gp=gpar(fontsize=20, col="black"))	
      } 
    }
)
# p2b

# -------------------------------------------------------------

# Correlations and their p-values
C2 <- reshape::cast(res2, Variable ~ Component, value="Tau") %>% tibble::column_to_rownames("Variable")
P2 <- reshape::cast(res2, Variable ~ Component, value="Padj") %>% tibble::column_to_rownames("Variable")
dic <- c(SUM_norm="ARG load", ARG_div="ARG diversity", Species_diversity="Species diversity")

C3 <- C2[c("ARG_div", "SUM_norm"), ]
P3 <- P2[c("ARG_div", "SUM_norm"), ]
p3 <- Heatmap(C3,
  name="Tau",
  col = col_fun,  
  cluster_rows=FALSE,
  cluster_columns=FALSE,    
  show_row_dend=FALSE,
  show_column_dend=FALSE,
  show_column_names = FALSE,   
  # row_title="ARG burden",
  row_labels=dic[rownames(C3)],
  heatmap_height = unit(0.5*nrow(C3), "cm"),
    cell_fun = function(j, i, x, y, w, h, fill) {
      gb = textGrob("*")
      gb_w = convertWidth(grobWidth(gb), "mm")
      gb_h = convertHeight(grobHeight(gb), "mm")
      if(P3[i, j] < 0.05) {
        grid.text("*", x, y - gb_h*0.8 + gb_w*0, gp=gpar(fontsize=20, col="darkgray"))	
      } 
    }
)

C4 <- C2["Species_diversity",]
P4 <- P2["Species_diversity",]
p4 <- Heatmap(C4,
  name="Tau",
  col = col_fun,  
  cluster_rows=FALSE,
  cluster_columns=FALSE,  
  show_row_dend=FALSE,
  show_column_dend=FALSE,
  show_column_names = FALSE,   
  row_title="",
  row_labels=dic[rownames(C4)],
  heatmap_height = unit(0.5*nrow(C4), "cm"),
    cell_fun = function(j, i, x, y, w, h, fill) {
      gb = textGrob("*")
      gb_w = convertWidth(grobWidth(gb), "mm")
      gb_h = convertHeight(grobHeight(gb), "mm")
      if(P4[i, j] < 0.05) {
        grid.text("*", x, y - gb_h*0.8 + gb_w*0, gp=gpar(fontsize=20, col="darkgray"))	
      } 
    }
)

# -------------------------------------------

# Save the Excel that is shown on the heatmap

tab <- res %>% select(-P)
tab$Variable <- dic.args[str_replace_all(as.character(tab$Variable), "[-\\,]", ".")]
tab$Component <- dic.nmf[tab$Component]
tab <- tab %>% rename(Enterosignature=Component) %>%
               mutate(Tau=round(Tau, 2)) %>%
               mutate(Padj=round(Padj, 16))

tab2 <- res2 %>% select(-P)
dic.variables <- as.character(unique(tab2$Variable)); names(dic.variables) <- as.character(dic.variables)
dic.variables <- gsub("SUM_norm", "ARG_load", dic.variables)
dic.variables <- gsub("ARG_div", "ARG_diversity", dic.variables)
tab2$Variable <- dic.variables[as.character(tab2$Variable)]
tab2$Component <- dic.nmf[tab2$Component]
tab2 <- tab2 %>% rename(Enterosignature=Component) %>%
                 mutate(Tau=round(Tau, 2)) %>%
                 mutate(Padj=round(Padj, 16))

tabs <- rbind(tab, tab2)
tabs$Padj <- str_replace(as.character(tabs$Padj), "^0$", "<1e-16")
tabs$Tau <- as.character(tabs$Tau)
tabs <- tabs %>% rename(FDR=Padj)

library(writexl)
tabs$Variable <- str_replace(tabs$Variable, "_", " ")
ES_ARG_associations <- tabs

# -------------------------------------------

# ES scores for all samples / ggplot
set.seed(22254)
scores <- t(apply(W, 1, function (x) {x/sum(x)}))
cap <- quantile(unlist(scores), 0.9)
scores[scores>cap] <- cap # cap on colors
colnames(scores) <- dic.nmf[colnames(scores)]
dd <- vegan::vegdist(scores, "euclidean")
hc <- hclust(dd, method="complete")
d <- melt(scores[hc$order,])
d$Var1 <- factor(d$Var1, levels=unique(d$Var1))
d$Var2 <- factor(d$Var2, levels=rev(unique(d$Var2)))
theme_set(theme_bw(30))
p <- ggplot(d, aes(x=Var1, y=Var2, fill=value)) +
       geom_tile() +
       scale_fill_gradient(low = "black", high = "white") +
       #labs(title="Samples (N=7095)", y="", x="") +
       labs(title="", y="", x="") +       
       theme(axis.text.x=element_blank(),
             axis.title.x=element_blank(),
             axis.ticks.x=element_blank()
	     ) +
       guides(fill = guide_legend(title="Weight", reverse=TRUE))
p.es.scores <- p


breaks <- seq(-0.15, 0.15, 0.05)
cap <- quantile(unlist(weights), 0.9)
weights[weights>cap] <- cap # cap on colors
p1horizontal <- Heatmap(t(weights),
  name="Weight",
  column_names_gp = gpar(fontsize = 10),
  column_title_gp = gpar(fontsize = 10),  
  col = rev(brewer.pal(length(breaks), "Greys")),
  show_column_dend=FALSE,
  cluster_rows=FALSE,
  column_title="",
  heatmap_width = unit(nrow(weights)/3, "cm"),
  show_row_names = FALSE, 
  left_annotation = rowAnnotation(
          text=anno_text(dic.nmf[colnames(weights)], rot=0, location=0, just = "left", show_name=FALSE)
  	),
  bottom_annotation = HeatmapAnnotation(Family=fam, na_col="gray", col=list(Family=pal[unique(fam)]), show_annotation_name=TRUE, show_legend=TRUE)
  )

#figES  <- p4 %v% p3 %v% p2 %v% p1
#figES2 <- p4 %v% p3 %v% p2b



