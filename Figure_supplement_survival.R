# Figure_supplement_survival.R should be run after Figure4.R (uses saved data from others)

#Plots forest plot and saves tsv file of Cox ph with mortality and ARG burden
# Cox survival analysis
library("tidyverse")
library("survival")
library("mia")
source("R_functions.R")

# Generate the RDS files (slow)
#source("causes.R") # table.cause / 2-3 hours / # Cause-specific mortality and sepsis
#source("sepsis.R") # fp.sepsis, cox_model_final_print_sepsis

li1 <- readRDS(file="totalmortality.rds")
li2 <- readRDS(file="sepsis.rds")
li3 <- readRDS(file="causes.rds")

# Write results in a table
li <- list(Totalmortality=li1$table,
           Sepsis=li2$table,
	   Causes=li3$table
           )
writexl::write_xlsx(li, path = "../RESULTS/R_results/Extended_data_table_SURVIVAL.xlsx")

# ----------------------------------------------------------

# Supplementary Figure: Sepsis & Cause-specific effect sizes for ARG load
base.size <- 20
theme_set(theme_bw(base.size))

th <- theme(plot.margin  = unit(c(1, 1, 1, 1), "pt"),
	plot.title = element_text(size = base.size),  
        legend.position=c(0.8, 0.8),
        legend.text = element_text(size=0.7*base.size),
	legend.title = element_text(size=base.size),
                                  axis.text.x = element_text(size = base.size),
                                  axis.text.y = element_text(size = base.size),
                                  axis.title.x = element_text(size = base.size),
                                  axis.title.y = element_text(size = base.size)	  
				  ) 

figS <- cowplot::plot_grid(li2$figure + th,
                           li3$figure + th,
			   labels="auto",
			   align="h",
			   label_size = 20,
			   rel_widths=c(0.57, 0.43)) 

pdf("../RESULTS/R_results/FigS8.pdf", width=16, height=5)
print(figS)
dev.off()

# Compare to traditional Cox model
# Run Cox model and drop out non-significant variables (in step())
#cox_model<-
#  coxph(
#  Surv(DEATH_AGEDIFF, DEATH)  ~ . ,
#  data  = df_full[, varnames]
#) %>% step()
## Fit the final model with the selected variables
#cox_model_final <- cox_model %>% finalfit::fit2df(condense = FALSE) 
## effect sizes compare well between freq and probabilistic estimates 
#hr <- cox_model$coefficients
#hrf <- cox_model_final$HR; names(hrf) <- cox_model_final$explanatory
#hrstan <- exp(fixef(fit)[, "Estimate"])
#plot(hrf, hrstan[names(hrf)], xlab="Frequentist", ylab="Stan", type="n", xlim=c(0, 3), ylim=c(0, 3))
#text(hrf, hrstan[names(hrf)], names(hrf)); abline(0, 1)

