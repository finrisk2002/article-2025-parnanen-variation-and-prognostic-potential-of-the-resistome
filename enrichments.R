
# Calculate Top ARG enrichments by family abundance bin
nbreaks <- 5
theme_set(theme_bw(20))
sigfam <- setdiff(levels(TSE$Signif_fam), "Other")
dfe <- calculateEnrichment(altExp(TSE, "Family", withColData=TRUE), varname="ARG_burd", assay.type="relabundance", features=sigfam, nbreaks=nbreaks) 

# ARG classes vs. enterosignatures
# corsf <- cor(t(assay(altExp(TSE, "Family_aggregated"), "counts")), TSE$SUM_norm, method="kendall")
# Calculate all NMF associations at one go to deal with multiple testing properly
F <- t(assay(altExp(TSE, "FamilyPrevalent"), "counts"))
colData(TSE) <- cbind(colData(TSE), F)
selected.vars <- "SUM_norm"

# Combine family abundances and colData
DF <- cbind(as.data.frame(colData(TSE)), as.data.frame(t(assay(altExp(TSE, "Family"), "relabundance"))))


pics <- pics2 <- list()
g <- NULL
gs <- levels(TSE$ARG_burd)
for (g in gs){
  print(g)
  high <- NULL
  for (tax in sigfam) {
    A <- assay(altExp(TSE, "Family"), "relabundance")[tax,]
    br <- unique(c(min(A), quantile(A[A>min(A)], probs=seq(0, 1, length=nbreaks))))
    cuts <- cut(A, breaks=br, include.lowest=TRUE, right=FALSE)
    # Observed / Expected gives the enrichment
    ptab <- table(cuts, TSE$ARG_burd)[, tolower(g)] / (mean(TSE$ARG_burd==tolower(g)) * table(cuts))  
    high <- rbind(high, ptab)
  }
  rownames(high) <- sigfam
  colnames(high) <- paste0("Q", 1:5)

}

# ----------------------------------------------------------------------------------------

theme_set(theme_bw(20))
pics <- list()

# High-ARG enrichments by bacterial family
fig3enr <- ggplot(dfe, aes(x=Bin, y=value, group=Level, color=Level)) +
    geom_line(size=3) +
    geom_point(size=5) +    
    scale_y_continuous(label=scales::percent, limits=range(dfe$value)) +
    geom_hline(yintercept=0) +
    scale_color_manual(values = pal.family()[as.character(levels(dfe$Level))], name="") +    
    labs(x="Abundance bin", y="ARG\nEnrichment (%)", title="High ARG load") +
    facet_wrap(Level ~ ., ncol=1) +
    theme(legend.position="none", plot.title = element_text(size = 20))

# High-ARG enrichments by ES
dic.nmf <- get.dic.nmf()
df2 <- calculateEnrichment(altExp(TSE, "ES", withColData=TRUE), varname="ARG_burd", assay.type="signal", features=NULL, nbreaks=nbreaks)

# With probabilistic credible intervals
# TODO finalize later:
# - add colors
# - also update the supplementary bacterial family figure
# - update corresponding tables
# - polish code
pics.brms <- calculateEnrichment.brms(altExp(TSE, "ES", withColData=TRUE), varname="ARG_burd", assay.type="signal", features=NULL, nbreaks=nbreaks)
pic <- plot_grid(pics.brms[[1]],# + scale_color_manual(values=pal.es()[names(pics.brms)[[1]]]),
                 pics.brms[[2]],
                 pics.brms[[3]],
                 pics.brms[[4]],
                 pics.brms[[5]],
		 nrow=1
                 )

# Library sizes do vary significantly between bins
#df2b <- calculateEnrichment2(altExp(TSE, "ES", withColData=TRUE), varname="ARG_burd", assay.type="signal", features=NULL, nbreaks=nbreaks) %>%
#        as.data.frame
#df2b$librarysize <- as.numeric(colData(TSE)[rownames(df2b), "V2"])

# Median library sizes per cut
libsizes <- calculateEnrichment3(altExp(TSE, "ES", withColData=TRUE), varname="ARG_burd", assay.type="signal", features=NULL, nbreaks=nbreaks) 

# Variations in library sizes are not associated with the enrichment values
cor.pval <- c()
for (feat in rownames(libsizes)) {
  if (!all(colnames(libsizes)==as.character(subset(df2, Level==feat)$Bin))) {stop("error")}
  cor.pval[[feat]] <- cor.test(as.numeric(libsizes[feat, ]), subset(df2, Level==feat)$value, method="kendall")$p.value
}
print(p.adjust(cor.pval))

#Library sizes do vary significantly between bins (Kruskal-Wallis
# test; p<0.05 for all enterosignatures); however, variations in
# library sizes are not associated with the enrichment values
# (Kendall's tau).
df2$Level <- factor(dic.nmf[df2$Level], levels=c("ES-Bact", "ES-Firm", "ES-Prev", "ES-Bifi", "ES-Esch"))

#---------------------------------------------------------------------------------------------------------------

fig3abswap <- ggplot(df2, aes(x=Bin, y=value, group=Level, color=Level)) +
    geom_line(size=1) +
    geom_point(size=3) +
    scale_y_continuous(label=scales::percent, limits=range(df2$value), breaks=seq(-0.4, 0.8, 0.4)) +
    geom_hline(yintercept=0) +    
    labs(x="Abundance bin", y="ARG\nEnrichment (%)") +
    facet_wrap(Level ~ ., ncol=5) +    
    scale_color_manual(values = pal.es(), name="") 

Ctot <- as.data.frame(rbind(C, C2))
Ctot$varname <- rownames(Ctot)
rownames(Ctot) <- NULL
colnames(Ctot) <- dic.nmf[colnames(Ctot)]
df <- melt(Ctot)

Ptot <- as.data.frame(rbind(P, P2))
Ptot$varname <- rownames(Ptot)
rownames(Ptot) <- NULL
colnames(Ptot) <- dic.nmf[colnames(Ptot)]
dfp <- melt(Ptot)

df$FDR <- factor(ifelse(dfp$value < 0.05, "p<0.05", "NS"))
colnames(df) <- c("variable", "ES", "value", "FDR")
dic <- c(dic.args, c("SUM_norm"="ARG load", "ARG_div"="ARG diversity", "Species_diversity"="Species diversity"))
df$variable <- factor(dic[df$variable], levels=rev(unname(dic[c(rev(rownames(C2)), rownames(C))])))

p <- ggplot(df, aes(x=variable, y=value, fill=FDR)) +
       geom_bar(stat="identity", position="dodge", color="black") +
       labs(y="Kendall's tau", x="") +       
       facet_grid(. ~ ES) +
       coord_flip() +
       scale_fill_manual(values=rev(c("darkgray", "white")))
# print(p)
pic.es.associations <- p


base.size <- 20
figa <- fig3abswap + theme(legend.position="none", plot.title = element_text(size = base.size),
                                  axis.text.x = element_text(size = 0.8*base.size),
                                  axis.text.y = element_text(size = base.size),
                                  axis.title.x = element_text(size = base.size),
                                  axis.title.y = element_text(size = base.size)	  
    				  )