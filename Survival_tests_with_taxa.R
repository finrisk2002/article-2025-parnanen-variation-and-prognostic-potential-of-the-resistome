library(survival)

tse <- transformAssay(altExp(TSE, "FamilyPrevalent"), method = "log10", pseudocount = 0.00001)

df2 <-
  cbind(as.data.frame(colData(TSE)), as.data.frame(t(assay(
    (tse), "log10"
  ))))
df2$Age <- df2$BL_AGE/10

df_tmp<- df2 %>% dplyr::select(ends_with(c(
  taxa,
  c(
    "SUM_norm",
    "SysBP",
    "Age",
    "DEATH",
    "DEATH_AGEDIFF",
    "MEN",
    "CURR_SMOKE",
    "KY100_14",
    "TULOT",
    "BMI",
    "BL_USE_RX_L",
    "BP_TREAT",
    "PREVAL_DIAB",
    "BL_USE_RX_J01"
  ))
)) 

df_taxa <- df_tmp[complete.cases(df_tmp),]
df_taxa$ARG_load <- log10(df_taxa$SUM_norm)
df_taxa <- df_taxa %>% dplyr::select(-ends_with(c("RX_J01_AGE", "DEATH_AGE", "AB1_SEPSIS_BACT_AGE", "BL_AGE")))

# Get the list of taxa names
taxa <- as.data.frame(t(assay((tse), "log10"))) %>% names()

# Initialize variables to store results
all <- data.frame(matrix(NA, # Create empty data frame
                         nrow  = 4,
                         ncol  = length(taxa)))
all_ARG <- data.frame(matrix(NA, # Create empty data frame
                             nrow  = 4,
                             ncol  = length(taxa)))

p <- c(1:length(taxa))
p_ARG <- c(1:length(taxa))


# Covariates in the model
covariates <- c("DEATH", "DEATH_AGEDIFF", "ARG_load", "SysBP", "Age", "MEN", "CURR_SMOKE", "KY100_14", "TULOT", "BMI", "BL_USE_RX_L", "BP_TREAT", "PREVAL_DIAB")

# Loop through each taxa
for (i in 1:length(taxa)) {
  cox_model <- coxph(
    Surv(DEATH_AGEDIFF, DEATH) ~ .,
    data  = df_taxa %>% dplyr::select(ends_with(c(covariates)), taxa[i])
  )
  
  # Extract confidence intervals, p-values, and events
  x <- summary(cox_model)$conf.int[taxa[i],] %>% data.frame()
  x_ARG <- summary(cox_model)$conf.int["ARG_load",] %>% data.frame()
  y <- summary(cox_model)$coefficients[taxa[i], "Pr(>|z|)"]
  y_ARG <- summary(cox_model)$coefficients["ARG_load", "Pr(>|z|)"]
  # Set column names
  colnames(x) <- taxa[i]
  colnames(x_ARG) <- taxa[i]
  
  
  # Store results
  all <- cbind(all, x)
  all_ARG <- cbind(all_ARG, x_ARG)
  p[i] <- y
  p_ARG[i] <- y_ARG
  
  # Garbage collection
  gc()
}

# Create a dataframe with results

table_death_fam <- data.frame(t(all %>% dplyr::select(-starts_with("X"))))

table_death_fam$p <- p
table_death_fam$padj <- p.adjust(table_death_fam$p, method = "fdr")
table_death_fam <- round(table_death_fam, digits = 4)

table_death_ARG_w_fams <- data.frame(t(all_ARG %>% dplyr::select(-starts_with("X"))))
table_death_ARG_w_fams$p <- p_ARG
table_death_ARG_w_fams$padj <- p.adjust(table_death_ARG_w_fams$p, method = "fdr")
table_death_ARG_w_fams <- round(table_death_ARG_w_fams, digits = 4)

