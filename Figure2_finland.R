#print("Load functions and data object")
source("R_functions.R")
#print("Load TSE")

library(writexl)

TSE <- readRDS("../data/TSE.rds")

source("GLM_plots_TreeSE.R")  # Must be run after GLM_TreeSE.R; only small part needed here to initialize
source("Linear_models_ARG_load_TreeSE.R")
age.labs <- age.labels()
# source("Boxplot_figures_TreeSE.R")
source("Figure2_map.R")

# Plot
df <- colData(TSE) %>% as.data.frame

# -----------------------------------------------------

theme_set(theme_classic(20))
lims <- c(200, 320)
labsize <- 30

# Whole country; bin-specific distributions (not lm)
source("Figure2_brms_lineplots.R") # lineplots
# Same but per region; unstable and very slow; generates the figure "Figure5_regional.pdf"
# source("Figure2_brms_lineplots_with_regions.R") # lineplots; sync with Linear_models_ARG_load_TreeSE.R & covariates

# Save the Excel
# Figure 3c ARG load vs. covariate assosiations
write_xlsx(xs, "../RESULTS/R_results/DataS3.xlsx") 

# Collect Figure 2:
source("Fig2panel.R")

pdf("../RESULTS/R_results/Fig2.pdf", width=24, height=15)
print(fig3)
dev.off()

# ------------------------------------------------------------------

# Ordinary linear models by region and variable, for comparison
source("Figure2_regional_trends.R")
# -> "Figure5_regional_lm.xlsx"
# -> "Figure5_regional_trends_fitted.pdf"

# stop("Run region-specific experiments manually from here onwards")
# Region-specific analyses; slow; run manually
# This is better in principle, with more rigorous log-linear model
# but obtaining level-specific predictions for simple visualization
# is in process
# source("lm_brms_regional.R")

