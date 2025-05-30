# Resistome profiles for all samples
set.seed(22254)
resistome <- assay(altExp(TSE, "arg_class"), "relabundance")

# Check top-n most abundant ARG classes
topvar2 <- names(head(rev(sort(apply(resistome, 1, function (x) {mean(x)}))), 10))

pseudo <- sort(unique(as.vector(resistome)))[[2]]/2

# Show log10 relative abundance
scores <- t(log10(resistome[topvar2, ] + pseudo))

arg.palette <- get.arg.palette()
names(dic) <- dic <- names(arg.palette)
dic <- gsub("Macrolide, Lincosamide, Streptogramin B", "MLSB", dic)
dic <- gsub("Macrolide, Aminoglycoside, Tetracycline, Quinolone, Amphenicol, Rifamycin", "MATQAR", dic)
dic2 <- setdiff(colnames(scores), names(dic)); names(dic2) <- dic2
dic <- c(dic, dic2)
colnames(scores) <- dic[colnames(scores)]

cap <- quantile(unlist(scores), 0.9)
scores[scores>cap] <- cap # cap on colors
# scores[scores<(-cap)] <- -cap # cap on colors
dd <- vegan::vegdist(scores, "euclidean")
#hc <- hclust(dd, method="complete")
hc <- hclust(dd, method="ward.D2")
d <- melt(scores[hc$order,])
# d <- melt(scores[order(TSE$SUM_norm),])
d$Var1 <- factor(d$Var1, levels=unique(d$Var1))
d$Var2 <- factor(d$Var2, levels=rev(unique(d$Var2)))
p <- ggplot(d, aes(x=Var1, y=Var2, fill=value)) +
       geom_tile() +
       scale_fill_gradient(low = "lightgray", high = "darkred", 
                           breaks=log10(10^c(seq(-3, -1, 1), 100)),			   
                           labels=paste(c(100 * 10^seq(-3, -1, 1), 100), "%")
			   ) +              
       labs(x="Samples (n=7,095)", y="") +
       theme(axis.text.x=element_blank(),
             axis.title.x=element_blank(),
             axis.ticks.x=element_blank()
	     ) +
       guides(fill = guide_legend(title="Abundance", reverse=TRUE)) +
       labs(title=paste0("Samples (n=", ncol(resistome), ")"))
p.resistome.scores <- p
# print(p)


