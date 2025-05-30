## ----setup--------------------------------------------------------------------
TSE <- readRDS(file="../data/TSE.rds")
library(reshape2)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)
library(mia)


## -----------------------------------------------------------------------------
# Pick the NMF loadings and scores
H <- metadata(TSE)$NMF_loadings
W <- as.data.frame(colData(TSE)) %>% select(starts_with("nmf"))


## ----esprofiles---------------------------------------------------------------
# ES profile across bacteria (normalized to 1 per ES)
Hscaled <- t(apply(H, 1, function (x) {x/sum(x)}))
Hscaled2 <- apply(H, 2, function (x) {x/sum(x)})

# Remove genera whose relative importance in all ES remains under <x%
# ESprofiles <- Hscaled[, !colMeans(Hscaled<5/100)==1]
x <- Hscaled2
breaks <- seq(0, round(max(x), 1), length.out = 9)
pheatmap(t(x),
  color = rev(brewer.pal(length(breaks), "Blues")),
  breaks = breaks,
  show_colnames = TRUE,
  cluster_cols = F)

x <- Hscaled
dfm <- melt(x)
colnames(dfm) <- c("ES", "Genus", "score")
dfm <- dfm %>% filter(!Genus=="Other")
m <-mapTaxonomy(TSE, taxa = as.character(dfm$Genus), to="Phylum")
dfm$Phylum <- unname(m)
# Rearrange rows
dfm <- dfm %>% arrange(Phylum) %>% mutate(Genus=factor(Genus, levels=unique(Genus)))

p <- ggplot(dfm, aes(x=Genus, y=score, fill=Phylum)) +
       geom_bar(stat="identity", position="dodge") +
       geom_text(data=(dfm%>%filter(score>0.05)), aes(label=Genus), check_overlap=TRUE, angle=45, vjust=0.5, y=0.1, hjust=0.05) +
       theme(axis.text.x = element_text(angle = 90)) +
       scale_y_continuous(limits=c(0,1), label=scales::percent) + 
       facet_grid(ES ~ .) 
print(p)


## ----eval=FALSE, echo=FALSE---------------------------------------------------
## # Same as barplot
## #dfm <- melt(H)
## ## df <- dfm %>% arrange(Var1, value) %>% mutate(Var2=factor(Var2, levels=unique(Var2)))
## #p <- ggplot(dfm, aes(y=Var2,x=value)) +
## #       geom_bar(position="dodge", stat="identity") +
## #       facet_grid(. ~ Var1)
## #print(p)


## ----argloadcor, fig.width=10, fig.height=2-----------------------------------
cors <- cor(W, TSE$SUM_norm, method="kendall")
colnames(cors) <- "Tau"
cors <- as.data.frame(cors)
cors$Component <- rownames(cors)

# Heatmap
breaks <- seq(-0.15, 0.15, 0.03)
d <- t(cors[, "Tau"])
colnames(d) <- colnames(W)
pheatmap(d,
  color = rev(brewer.pal(length(breaks), "RdBu")),
  breaks = breaks, cluster_cols=FALSE)


# Same as Barplot
#p <- ggplot(cors, aes(y=Tau, x=Component)) +
#  geom_bar(stat="identity") +
#  labs(y="NMF component")
#print(p)
# Just checking:
# Same magnitude when associating individual genera
# cors2 <- cor(t(assay(altExp(TSE, "Top_genus"), "clr")), TSE$SUM_norm, method="kendall") # For comparison with individual genera
#range(cors2[,1])


## ----csts, fig.width=15, fig.height=5-----------------------------------------
d <- cbind(as.data.frame(colData(TSE)), W)
dfm <- melt(d[, c("CST", colnames(W))])
p <- ggplot(dfm, aes(y=variable, x=value)) +
  geom_jitter(alpha=0.2, width=0.2) +
  facet_grid(. ~ CST)
print(p)


## -----------------------------------------------------------------------------
cors3 <- cor(t(assay(altExp(TSE, "Class_top"), "counts")), W, method="kendall")
# Heatmap
breaks <- seq(-0.15, 0.15, 0.03)
pheatmap(cors3,
  color = rev(brewer.pal(length(breaks), "RdBu")),
  breaks = breaks)

