# Create Extended Data Tables
# Load data
source("loadData.R")
library(mia)
library(tidyverse)
## Antibiotic reimburements
df_NEVT <-
  colData(TSE) %>% data.frame() %>%
  dplyr::select(contains(c("J01"))) %>% 
  dplyr::select(contains(c("NEVT"))) %>%
  colSums() %>% data.frame()

colnames(df_NEVT) <- c("Reimbursements (nr of events)")
df_NEVT$Antibiotic <- rownames(df_NEVT)
df_NEVT <-
  df_NEVT %>% replace_variable_with_dictionary("Antibiotic", substitutions)
df_PREVAL <-
  colData(TSE) %>% data.frame() %>%
  dplyr::select(contains(c("J01"))) %>% 
  dplyr::select(-contains(c("NEVT", "AGE", "INCIDENT", "BL_USE"))) %>%  colSums() %>% data.frame
df_PREVAL$Antibiotic <- rownames(df_PREVAL)
df_PREVAL <-
  df_PREVAL %>% replace_variable_with_dictionary("Antibiotic", substitutions)
colnames(df_PREVAL) <-
  c("Prior reimbursements (nr of participants)")


library(psych)

antibiotic_use <- TSE %>% colData() %>% data.frame() %>%
  dplyr::select(contains("J01") & ends_with("NEVT")) %>%
  psych::describe() %>% data.frame %>%
  round(digits = 2) %>%  mutate(median_min_max = sprintf("%d (%d-%d)", median, min, max))

df_NEVT$"Prior reimbursements (nr of participants)" <- df_PREVAL[,1]

df_NEVT$Mean <- antibiotic_use$mean
df_NEVT$"Standard deviation" <- antibiotic_use$sd
df_NEVT$"Median (min-max)" <- antibiotic_use$median_min_max
  
## Correlations in the general stats, report p values
library(Hmisc)
cor.mats <- TSE %>% colData() %>% data.frame() %>%
  dplyr::select(c("vaesto", "KY100_14", "TULOT", "BL_AGE", "PREVAL_RX_J01_NEVT", "BMI", "KOL")) %>% as.matrix() %>% rcorr()
cor.p <- cor.mats$P %>% data.frame

cor.r <- cor.mats$r %>% data.frame

cor.n <- cor.mats$n %>% data.frame

# Create a function to format the entries with coefficient and p-value
format_entry <- function(r, p, n) {
  sprintf("%.2f (%.3f) [%.4f]", r, p, n)
}

# Create an empty matrix to store the formatted entries
combined_matrix <- matrix(nrow = nrow(cor.r), ncol = ncol(cor.r), dimnames = dimnames(cor.r))

# Fill the combined matrix with formatted entries
for (i in 1:nrow(combined_matrix)) {
  for (j in 1:ncol(combined_matrix)) {
    combined_matrix[i, j] <- format_entry(cor.r[i, j], cor.p[i, j], cor.n[i, j])
  }
}

# Convert to data frame for better viewing in R
combined_matrix_df <- as.data.frame(combined_matrix)

combined_matrix_df$rownames <- rownames(combined_matrix_df)
combined_matrix_df <-
  replace_variable_with_dictionary(combined_matrix_df, "rownames", substitutions)


## Antibiotic resistance gene classes
AB_class_table <- TSE %>% colData() %>% data.frame %>% dplyr::select("Top_class") %>% table %>% data.frame()

library(writexl)
#xs <- list() # Do not write unless code is modified, small edits made in Excel manually
#xs[["S1"]] <- combined_matrix_df # ED 1
# xs[["S2"]] <- df_NEVT #ED 2
#xs[["S3"]] <- AB_class_table # ED 6
#write_xlsx(xs, path = "../RESULTS/R_results/Extended_data_tables1-2,6.xlsx")

# Write in separate tables
write_xlsx(df_NEVT, path = "../RESULTS/R_results/Table_S1.xlsx")
write_xlsx(AB_class_table, path = "../RESULTS/R_results/AB_class_table.xlsx")
write_xlsx(combined_matrix_df, path = "../RESULTS/R_results/Table_S9.xlsx")


