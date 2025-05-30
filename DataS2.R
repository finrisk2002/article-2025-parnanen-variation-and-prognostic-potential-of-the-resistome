# Diet and taxa

TSE <- readRDS("../data/TSE.rds")
altExp(TSE, "FamilyPrevalent") <- mergeFeaturesByPrevalence(TSE[rowData(TSE)$Domain=="Bacteria"], rank="Family", assay.type="relabundance", detection=0.1/100, prevalence=1/100)

library(mia)
library(tidyverse)
TSE_family <- altExp(TSE, "FamilyPrevalent", withColData=TRUE) 
##########

# Pick abundance assay
df <- t(assay(TSE_family, "relabundance"))
row.names(df) <- TSE_family$Row.names

# Pick metadata
df1 <- data.frame(colData(TSE_family))
row.names(df1) <- TSE_family$Row.names
df1 <- df1 %>% dplyr::select((-ends_with(".1")))

# Merge
lm_df <- merge(df, df1, by = 0)
row.names(lm_df) <- TSE_family$Row.names

# Log-transform bacterial families
df <- lm_df %>% 
  mutate_at(vars(matches(c("ceae", "deae", "classified"))),
            ~ log10(. + min(.[. > 0]) / 2)) %>% 
  dplyr::select(-where(~ any(is.infinite(.))))

# Extract bacterial families
families <- df %>% dplyr::select(matches(c("ceae", "deae", "classified"))) %>% names

# Extract NMF components
nmfs <- df %>% dplyr::select(matches(c("nmf"))) %>% names

# List diet components
KY100s <- df %>% dplyr::select(matches(c("KY100", "HFC_score", "BMI", "KOL"))) %>% names

# Function to run GLM for a given response variable (family or NMF)
run_glm_analysis <- function(response_vars, dataset_name) {
  results_list <- list()
  
  
  for (response_var in response_vars) {
    cat("Processing:", response_var, "\n")
    
    # Initialize an empty data frame for storing results
    all_results <- data.frame()
    
    for (ky in KY100s) {
      print(paste("Running GLM for", response_var, " ~ ", ky))
      
      glm_model <- glm(as.formula(paste(response_var, "~", ky)), data = df)
      
      # Extract coefficients and p-values
      summary_model <- summary(glm_model)$coefficients
      est <- summary_model[2, 1]  # Coefficient estimate
      se <- summary_model[2, 2]   # Standard error
      t_val <- summary_model[2, 3] # T-value
      p_val <- summary_model[2, 4] # P-value
      
      # Store in a new row
      all_results <- rbind(all_results, data.frame(
        ResponseVariable = response_var,
        Covariate = ky,
        Estimate = est,
        StdError = se,
        Tvalue = t_val,
        p = p_val
      ))
    }
    
    results_list[[response_var]] <- all_results
  }
  
  # Combine results into one table
  final_table <- do.call(rbind, results_list)
  
  # Adjust p-values within each dataset separately
  final_table$FDR <- p.adjust(final_table$p, method = "fdr")
  
  # Save each table separately
  assign(paste0("final_table_", dataset_name), final_table, envir = .GlobalEnv)
}

# Run for bacterial families
run_glm_analysis(families, "families")
final_table_families <- final_table_families %>% arrange(FDR) %>% filter(FDR < 0.05) 
final_table_families <- replace_variable_with_dictionary(final_table_families, "Covariate", substitutions)

# Run for NMF components
run_glm_analysis(nmfs, "nmfs")

final_table_nmfs <- replace_variable_with_dictionary(final_table_nmfs, "Covariate", substitutions) 

final_table_nmfs <- replace_variable_with_dictionary(final_table_nmfs, "ResponseVariable", substitutions)

final_table_nmfs <- final_table_nmfs %>% arrange(FDR) %>% filter(FDR < 0.05) 

# Correlation matrix for diet

# Select only KY100s variables
ky100_df <- df %>% dplyr::select(all_of(c(KY100s, "PREVAL_RX_J01_NEVT")))

# Compute correlation matrix
cor_matrix <- cor(ky100_df, use = "pairwise.complete.obs", method = "pearson")

# Convert to a data frame for easier viewing
cor_matrix_df <- as.data.frame(cor_matrix)
cor_matrix_df$Covariate <- rownames(cor_matrix_df)
cor_matrix_df <- replace_variable_with_dictionary(cor_matrix_df, "Covariate", substitutions )


library(writexl)

# Convert row names to a column in the Correlation Matrix
cor_matrix_df_clean <- as.data.frame(cor_matrix_df)
cor_matrix_df_clean[is.na(cor_matrix_df_clean)] <- 0  # Replace NAs with 0
cor_matrix_df_clean$RowNames <- rownames(cor_matrix_df_clean)

# Convert row names to a column in the Families table
final_table_families$RowNames <- rownames(final_table_families)

# Convert row names to a column in the NMFs table
final_table_nmfs$RowNames <- rownames(final_table_nmfs)


# Combine all tables into a list
sheets <- list(
  "Significant Families" = final_table_families,
  "Significant NMFs" = final_table_nmfs,
  "Correlation Matrix" = cor_matrix_df_clean
)

# Save the list of data frames to an Excel file
write_xlsx(sheets, "../RESULTS/R_results/DataS2.xlsx")

