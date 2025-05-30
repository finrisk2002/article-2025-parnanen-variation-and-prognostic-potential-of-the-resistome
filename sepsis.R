# Read rds and run boosted GLM
library("tidyverse")
library("survival")
library("mia")
library("tidybayes")
library("survminer")
library(stringr)

TSE <- readRDS("../data/TSE.rds") 
# source("GLMs_TreeSE.R")
source("Dictionary.R")
source("R_functions.R")

# Prepare data frame for the survival analysis
# Note: logging and scaling for some variables!
DF <- prepare.data.for.survival.analysis(TSE)

# Sepsis occurrence time
DF$TIMEDIFF <- DF$AB1_SEPSIS_BACT_YEAR+DF$AB1_SEPSIS_BACT_AGEDIFF
DF$EVENT <- DF$AB1_SEPSIS_BACT_AGE

# Run the survival model (Cox) with Stan using brms
# Variables to use and control
vars <- c("EVENT", "TIMEDIFF", "ARGload_log10", get.covariates())

# ----------------------------------------------------------

dff <- DF[, vars] %>% filter(!TIMEDIFF < 0);
dff <- dff[complete.cases(dff),]
binaries <- c("MEN", "CURR_SMOKE", "BL_USE_RX_L", "BL_USE_RX_J01", "PREVAL_DIAB", "BP_TREAT") # Binary variables
continuous <- setdiff(vars, c(binaries, "EVENT", "TIMEDIFF"))
dff <- cbind(dff[, c("EVENT", "TIMEDIFF")],
             dff[, binaries],
             dff[, continuous]
	     )

# Fit the model
fit <- get.fitted.model(dff)

# Gather the same as in varnames, except the time and event variables
dfplot <- get.dfplot(fit)
# Only keep the ones where credible interval does not contain 1
dfplot <- dfplot[sign(dfplot$q5-1) == sign(dfplot$q95-1),]
# survival_totalmortality_byfamily.R koodi sisältää tietoa miten vaihtaa muihin luottovaleihin tarvittaessa
    
# Sort by effects and clean up variable names
mytab <- dfplot %>% select(variable_print, median, q5, q95) %>% arrange(desc(median))
print(c(min(mytab$q5), max(mytab$q95)))

# -------------------------------------------------------------------------------

theme_set(theme_bw(20))
v <- c(0.5, 1, 2, 4)
fp.sepsis <- ggplot(dfplot, aes(y = variable_print, x = median, xmin = q5, xmax = q95)) + 
   geom_pointinterval() +
   labs(x="Hazard ratio", y="", title="Sepsis") +
   geom_vline(xintercept=1, color="darkgray", linetype=2) +
   scale_x_continuous(breaks=v, labels=v, trans="log2", limits=c(0.84, 5.55))   
fp.sepsis


# ----------------------------------------------------------

# With added scaling / LL
dff.scaled <- DF[, vars] %>% filter(!TIMEDIFF < 0);
dff.scaled <- dff.scaled[complete.cases(dff.scaled),]
dff.scaled <- cbind(dff.scaled[, c("EVENT", "TIMEDIFF")],
             dff.scaled[, binaries],
             scale(dff.scaled[, continuous])
	     )
fit.scaled <- get.fitted.model(dff.scaled)
dfplot.scaled <- get.dfplot(fit.scaled)
dfplot.scaled <- dfplot.scaled[sign(dfplot.scaled$q5-1) == sign(dfplot.scaled$q95-1),]
# Sort by effects and clean up variable names
mytab.scaled <- dfplot.scaled %>% select(variable_print, median, q5, q95) %>% arrange(desc(median))

# ----------------------------------------------------------

li <- list(table=mytab, table.scaled=mytab.scaled, figure=fp.sepsis)
saveRDS(li, file="sepsis.rds")
