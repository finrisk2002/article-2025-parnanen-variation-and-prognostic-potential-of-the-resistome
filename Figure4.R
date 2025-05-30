# Read rds and run boosted GLM
library("tidyverse")
library("survival")
library("mia")
library("tidybayes")
library("survminer")
library(stringr)

print("Load TSE")
TSE <- readRDS("../data/TSE.rds")
source("Dictionary.R")
source("R_functions.R")

# Prepare data frame for the survival analysis
# Note: logging and scaling for some variables!
DF <- prepare.data.for.survival.analysis(TSE)

# Event occurrence & time
DF$TIMEDIFF <- DF$DEATH_AGEDIFF
DF$EVENT <- DF$DEATH

# Run the survival model (Cox) with Stan using brms
# Variables to use and control
vars <- c("EVENT", "TIMEDIFF", "ARGload_log10", get.covariates())

# -------------------------------------------------------------------------

# Remove cases with negative f/u, only keep complete cases, scale variables to zero mean unit variance for comparability
dff <- DF %>% filter(!TIMEDIFF < 0) 
binaries <- c("MEN", "CURR_SMOKE", "BL_USE_RX_L", "BL_USE_RX_J01", "PREVAL_DIAB", "BP_TREAT") # Binary variables
dff <- cbind(dff[, c("EVENT", "TIMEDIFF")],
             dff[, binaries],
             dff[, setdiff(vars, c(binaries, "EVENT", "TIMEDIFF"))]
	     )
dff <- dff[complete.cases(dff),] # Important to have this *after* selecting the variables
	     
# N=6849

# Fit the model
set.seed(25624)
fit <- get.fitted.model(dff)
# Gather the same as in varnames, except the time and event variables
dfplot <- get.dfplot(fit)
# Only keep the ones where credible interval does not contain 1
dfplot <- dfplot[sign(dfplot$q5-1) == sign(dfplot$q95-1),]
# Sort by effects and clean up variable names
mytab <- dfplot %>% select(variable_print, median, q5, q95) %>% arrange(desc(median))

# ----------------------------------------------------

# Same, but continuous covariates scaled by sd
dff.scaled <- DF %>% filter(!TIMEDIFF < 0) 
dff.scaled <- cbind(dff.scaled[, c("EVENT", "TIMEDIFF")],
             dff.scaled[, binaries],
             scale(dff.scaled[, setdiff(vars, c(binaries, "EVENT", "TIMEDIFF"))]) 
	     )
dff.scaled <- dff.scaled[complete.cases(dff.scaled),]  # Important to have this *after* selecting the variables

# Fit the model
set.seed(25624)
fit.scaled <- get.fitted.model(dff.scaled)
# Gather the same as in varnames, except the time and event variables
dfplot.scaled <- get.dfplot(fit.scaled)
# Only keep the ones where credible interval does not contain 1
dfplot.scaled <- dfplot.scaled[sign(dfplot.scaled$q5-1) == sign(dfplot.scaled$q95-1),]
# Sort by effects and clean up variable names
mytab.scaled <- dfplot.scaled %>% select(variable_print, median, q5, q95) %>% arrange(desc(median))

# -------------------------------------------------------------------------

theme_set(theme_bw(20))
v <- c(0.7, 1, 1.4, 2, 3)
fp.totalmortality <- ggplot(dfplot, aes(y = variable_print, x = median, xmin = q5, xmax = q95)) + 
   geom_pointinterval() +
   labs(x="Hazard ratio", y="", title="All-cause mortality risk") +
   geom_vline(xintercept=1, color="darkgray", linetype=2) #+
   #scale_x_continuous(breaks=v, labels=v, trans="log2", limits=c(0.86, 3.3)) 
fp.totalmortality

# Extra check: control each food component as a covariate, one-by-one
skip <- TRUE
if (!skip) {
fits <- NULL
foods <- names(colData(TSE))[grepl("KY100", names(colData(TSE)))]
for (varname in foods) {
  print(varname)
  vars <- c("EVENT", "TIMEDIFF", "ARGload_log10", get.covariates())
  vars <- str_replace(vars, "KY100_14", varname)
  dff <- DF[, vars] %>% filter(!TIMEDIFF < 0); dff <- dff[complete.cases(dff),]
  colnames(dff) <- str_replace(colnames(dff), varname, "food")
  fit <- get.fitted.model(dff)
  # Gather the same as in varnames, except the time and event variables
  dfplot <- get.dfplot3(fit)
  # Only keep the ones where credible interval does not contain 1
  tmp <- filter(dfplot, variable=="ARGload_log10")
  tmp$food <- varname
  fits <- rbind(fits, tmp)  
}
}

# -------------------------------------------------------------------------------

# Discrete curve

# List the "signif" covariates with ARGload
# covariates <- setdiff(as.character(dfplot$variable), "ARG load (log10 RPKM)")

# replace continuous ARG load with discrete ARG burden for the visualization
DFdisc <- dff
DFdisc$ARG_burd <- DF[rownames(DFdisc), "ARG_burd"]
DFdisc$ARGload_log10 <- NULL

set.seed(45421)

base.size <- 22
pics <- list()
for (men in c(0, 1)) {

    lab <- ifelse(men==0, "Women", "Men")

    d <- DFdisc %>% filter(MEN==men);
    f <- survfit(Surv(TIMEDIFF, EVENT) ~ ARG_burd, data = d)

    p <- ggsurvplot(f,
            data=d,
  	    conf.int=TRUE,
  	    pval=TRUE,
  	    fun="event"
   	    )$plot +
  	  scale_y_continuous(limits=c(0, 0.35), label=scales::percent) +
  	  scale_x_continuous(limits=c(0, 18)) +	  
  	  scale_color_manual(values=c("blue", "darkred"), labels=c("Conventional", "High")) +
  	  scale_fill_manual(values=c("blue", "darkred"), labels=c("Conventional", "High")) +  
  	  labs(title=lab, x="Time (y)", y="Cumulative mortality (%)") +
  	  guides(fill=guide_legend(nrow=1,byrow=TRUE, reverse=TRUE,title="ARG burden", title.position="left"),
                color=guide_legend(nrow=1,byrow=TRUE, reverse=TRUE, title="ARG burden", title.position="left")) +  
  	  theme_classic(base.size) +
  	  theme(plot.title = element_text(size = base.size),
                                  axis.text.x = element_text(size = base.size),
                                  axis.text.y = element_text(size = base.size),
                                  axis.title.x = element_text(size = base.size),
                                  axis.title.y = element_text(size = base.size),
	                          legend.text = element_text(size=0.9*base.size),
	     			  legend.title = element_text(size=base.size),
                                  legend.position.inside=c(0.4, 0.9)
				  )
				  
    pics[[lab]] <- p
    

}

library(cowplot)
leg <- get_legend(pics[[1]])
#figB <- plot_grid(pics[[1]] + theme(legend.position=c(0.4, 0.8)),
figB0 <- plot_grid(pics[[1]] + theme(legend.position="none"),
		  pics[[2]] + theme(legend.position="none"),
		  nrow=1)
figB <- plot_grid(figB0, leg, rel_heights=c(0.9, 0.1), nrow=2)



figA <- fp.totalmortality + 
     theme_classic(base.size) +
  theme(plot.title = element_text(size = base.size),
                                  axis.text.x = element_text(size = base.size),
                                  axis.text.y = element_text(size = base.size),
                                  axis.title.x = element_text(size = base.size),
                                  axis.title.y = element_text(size = base.size),
	                          legend.text = element_text(size=0.9*base.size),
	     			  legend.title = element_text(size=base.size))

fig <- cowplot::plot_grid(
        figA,
	figB,
	nrow=2,
	rel_heights=c(0.5, 0.5),
	labels=c("a", "b"),
	# align="v",
	label_size = 30
	)
	
# Main Figure: Total mortality and ARG load
pdf("../RESULTS/R_results/Fig4_2row.pdf", width=12, height=12)
print(fig)
dev.off()


figBv2 <- plot_grid(pics[[1]] +
            guides(fill =guide_legend(nrow=2, byrow=TRUE, reverse=TRUE, title="ARG burden", title.position="top"),
                   color=guide_legend(nrow=2, byrow=TRUE, reverse=TRUE, title="ARG burden", title.position="top")) +
            theme(legend.position=c(0.4, 0.8)),
          pics[[2]] + theme(legend.position="none"))
	
fig2 <- cowplot::plot_grid(
        figA,
        figBv2,
	nrow=1,
	rel_heights=c(0.5, 0.5),
	labels=c("a", "b"),
	# align="v",
	label_size = 30
	)
	
# Main Figure: Total mortality and ARG load
pdf("../RESULTS/R_results/Fig4.pdf", width=20, height=5)
print(fig2)
dev.off()


# Store the results
li <- list(figure=fig,
	   table=mytab,
	   table.scaled=mytab.scaled
	   
	   )
saveRDS(li, file="totalmortality.rds")


# --------------------------------------------------------------

# Additional manual check: 
# Controlled HR per each (prevalent) family
# source("survival_totalmortality_byfamily.R") 






