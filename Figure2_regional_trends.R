# Load required libraries
library(mboost)
library(caret)
library(pdp)
library(mgcv)
library(mia)
library(dplyr)
library(cowplot)
set.seed(498)
source("R_functions.R")
TSE <- readRDS("../data/TSE.rds") 

selected.vars <- c(
      "PREVAL_",
      "BL_",
      "KY100_",
      "vaesto",
      "^MEN",
      "Species",
      "BL_AGE",
      "TULOT",
      "HFC_score",
      "KOULGR",
      "BMI",
      "KOL",
      "SUM_norm",
      "EAST",
      "ARG",
      "V2" # Library size
    )

FR02names <-
  read.table(
    "../RESULTS/final_output_files/FR02_pheno_names.txt",
    sep = "\t",
    fill = TRUE,
    header = TRUE
  )

FR02names <- rbind(FR02names, c(VARIABLE = "Population density (log10)", LONGNAME = "Population density (log10)"))

results <- NULL
# Variables to test
vars <- c("KY100_14", "KY100_21", "TULOT", "vaesto", "BL_AGE")
for (alue in unique(TSE$ALUE)) {
  for (men in unique(TSE$MEN)) {
  
    print(paste(alue, men))
    
    # Unit harmonization same for all analyses (TODO: move this to Carpentry)
    DF <- prepare.data.for.survival.analysis(TSE[, TSE$ALUE==alue & TSE$MEN==men]) 

    source("check.significance.R")
    df$alue <- rep(alue, nrow(df))    
    results <- rbind(results, df)
    
  }
}

# Also add East/West analysis
for (alue in unique(TSE$EAST)) {
  for (men in unique(TSE$MEN)) {

    print(paste(alue, men))
    
    # Unit harmonization same for all analyses (TODO: move this to Carpentry)
    DF <- prepare.data.for.survival.analysis(TSE[, TSE$EAST==alue & TSE$MEN==men])
    
    source("check.significance.R")

    if (alue==0) {
      df$alue <- rep("West", nrow(df))
    } else if (alue==1) {
      df$alue <- rep("East", nrow(df))
    }
    results <- rbind(results, df)

  }
}

# ---------------------------------------------------------------------------------------

df <- results %>% bind_rows()
df$alue <- factor(df$alue, levels=c("Helsinki", "Turku", "Oulu", "Kuopio", "Karelia", "Lapland", "East", "West"))
df$cov.estimate <- as.numeric(df$cov.estimate)
df$CI2.5 <- as.numeric(df$CI2.5)
df$CI97.5 <- as.numeric(df$CI97.5)
df$gender <- factor(ifelse(df$gender==1, "Men", "Women"), levels=c("Men", "Women"))

# Print the table
dfw <- df
names(dfw) <- gsub("VARIABLE", "Variable", names(dfw))
names(dfw) <- gsub("cov.estimate", "Effect size", names(dfw))
names(dfw) <- gsub("padj", "FDR", names(dfw))
names(dfw) <- gsub("alue", "Region", names(dfw))
names(dfw) <- gsub("gender", "Gender", names(dfw))
dfw[, "Effect size"] <- round(dfw[, "Effect size"], 3)
dfw[, "FDR"] <- round(dfw[, "FDR"], 3)
library(writexl)

write_xlsx(dfw[, c("Variable", "Region", "Gender", "Effect size", "FDR")],
  path = "../RESULTS/R_results/Table_S3.xlsx")

lo <- 0.9
hi <- 1.3    
pics <- list()
for (varname in vars) {

  fp <-
    ggplot(data = subset(df, variable==varname),
         aes(
           x = alue,
           y = cov.estimate,
	   color = gender,
           ymin = CI2.5,
           ymax = CI97.5
         ))  +
  geom_linerange(linewidth = 0.5,
                 #colour = "grey",
		 position=position_dodge(.4))  +
  geom_point(
    shape = 18,
    size = 5,
    # colour  = "black",
    stroke  = 1, position=position_dodge(.4)
  )  +
  labs(#title="Pairwise LMs: Antibiotic use adjusted",
       title=dic.varnames()[varname],
       x="",
       y="Coefficient, 95% CI") +
  geom_hline(yintercept = 1, lty = 2)  +  # add a dotted line at x = 1 after flip
  coord_flip()  +  # flip coordinates (puts labels on y axis)
  #scale_y_continuous(limits  = c(0.95, 1.25), breaks = seq(0.95, 1.2, 0.05)) +
  scale_y_continuous(limits  = c(lo, hi), breaks = seq(lo, hi, 0.1), trans="log2") +    
  theme_classic(16) +
  scale_color_manual(
      values = rev(c("chocolate4", "cornflowerblue")),
      name  = "", labels = levels(df$gender)
  ) +
  theme(
    panel.border  = element_blank(),
    panel.background  = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line  = element_line(colour  = "black"),
    axis.text.y  = element_text(size  = 16, colour  = "black"),
    axis.text.x.bottom  = element_text(size  = 16, colour  = "black"),
    axis.title.x  = element_text(size  = 16, colour  = "black"),
    title  = element_text(size  = 16, colour  = "black")
  )

  pics[[varname]] <- fp

}

library(cowplot)
leg <- get_legend(pics[[1]])
p <- plot_grid(pics[[1]] + theme(legend.position="none"),
               pics[[2]] + theme(legend.position="none"),
	       pics[[3]] + theme(legend.position="none"),
	       pics[[4]] + theme(legend.position="none"),
	       pics[[5]] + theme(legend.position="none"),
	       leg,
	       nrow=2)

# Slope values with credible intervals
pdf("../RESULTS/R_results/Figure3_regional_trends.pdf", width=20, height=10)
print(p)
dev.off()

# Per region and variable
pics <- list()

df <- colData(TSE) %>% as.data.frame
df$vaesto_quantile <- cut(df$vaesto, breaks=c(0, 10, 100, 1000, 10000, Inf), labels=c("1", "10", "100", "1000", "10000"), include.lowest=TRUE, right=FALSE)

# Sort regions by population density
df$Region <- factor(df$Region, levels=as.vector((df %>% group_by(Region) %>% summarise(pop=mean(vaesto, na.rm=TRUE)) %>% arrange(pop))$Region))

# Transform population density
df$vaesto <- log10(df$vaesto)

vars <- c(vegetables="KY100_14", poultry="KY100_21", income="TULOT", popdens="vaesto_quantile", Age="age_class")

for (myvar in vars) {

  # Duplicate to add another layer on East/West
  dfa <- df
  dfb <- df; dfb$ALUE <- ifelse(dfb$EAST==1, "East", "West")
  df2 <- rbind(dfa, dfb)
  df2$ALUE <- factor(df2$ALUE)

  df2$varname <- as.numeric(df2[[myvar]])
  
  df2$Region <- df2$ALUE
  
  df2 <- df2 %>% filter(!is.na(SUM_norm) & !is.na(varname) & !is.na(Region))  
  df2 <- (df2 %>% unite("group", MEN, Region) %>% mutate(group=str_replace(group, "0_", "Women/")) %>%
                                                  mutate(group=str_replace(group, "1_", "Men/")))

  ds <- NULL
  for (g in unique(df2$group)) {
    d0 <- (df2 %>% filter(group==g))[, c("SUM_norm", "varname")]
    res <- lm(log(SUM_norm) ~ varname, data = d0)
    d <- unique(cbind(Variable=dic.varnames()[myvar],
             Region=unlist(strsplit(g, "/"))[[2]],
             Gender=unlist(strsplit(g, "/"))[[1]],	     
             Level=d0$varname,
	     Prediction=predict(res)))
    d <- as.data.frame(d)	     
    # d$Level <- factor(d$Level)
    d$Prediction <- as.numeric(d$Prediction)
    ds <- rbind(ds, d)
  }

  ds$Gender <- factor(ds$Gender, levels=rev(c("Women", "Men")))
  ds$Region <- factor(ds$Region, levels=c("Helsinki","Turku","Oulu","Kuopio","Karelia","Lapland", "East", "West"))
  if (is.numeric(df2[,myvar])) {
    ds$Level <- factor(ds$Level, levels=sort(unique(ds$Level)))
  } else if (is.factor(df2[,myvar])) {
    ds$Level <- factor(levels(df2[, myvar])[as.numeric(ds$Level)], levels=levels(df2[, myvar]))
  }
  
  v <- seq(180,300,20)
  p <- ggplot(ds, aes(x=Level, y=exp(Prediction), group=Region, shape=Region, color=Gender)) +
       geom_point(size=3) +
       geom_line() +
       #scale_y_log10(breaks=v,labels=v) +
       scale_y_continuous(breaks=v,labels=v) +        
       facet_grid(. ~ Gender) +
       theme(legend.title=element_blank()) +
       guides(colour = guide_legend(reverse=TRUE)) + 
       scale_color_manual(
         values = rev(c("chocolate4", "cornflowerblue")),	 
         name  = "") +
       scale_shape_manual(
         values = 1:8,	 
         name  = "") +	 
       labs(y="ARG load (RPKM)", title=dic.varnames()[myvar]) 
             
  pics[[myvar]] <- p

}

leg <- get_legend(pics[[1]])
p <- plot_grid(pics[[1]] + theme(legend.position="none"),
               pics[[2]] + theme(legend.position="none"),
	       pics[[3]] + theme(legend.position="none"),
	       pics[[4]] + theme(legend.position="none"),
	       pics[[5]] + theme(legend.position="none"),
	       leg,
	       ncol=3
	       )
print(p)


pdf("../RESULTS/R_results/FigureS3_regional_trends_fitted.pdf", width=22, height=9)
print(p)
dev.off()


