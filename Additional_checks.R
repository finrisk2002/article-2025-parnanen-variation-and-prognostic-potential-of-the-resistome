# Extra checks
TSE <- readRDS("../data/TSE.rds") 
load("../data/FR_metagenomes.RData")

library(dplyr)
library(mia)
library(ggplot2)
library(jtools)

##########
# Women antibiotics
df <- colData(TSE) %>% as.data.frame()
lm((PREVAL_RX_J01_NEVT)~MEN+BL_AGE+TULOT, data=df) %>% summ( confint = TRUE)
lm((df$KY100_14)~MEN+BL_AGE+TULOT, data=df) %>% summ( confint = TRUE)

###### # Detected ARGs in samples versus BLANKs
# Checks the median numer of ARGs detected in samples versus negative controls
ARG_bt_res %>% dplyr::select(starts_with("820")) %>% colSums() %>% median() # Samples
ARG_bt_res %>% dplyr::select(matches("BLANK")) %>% colSums() %>% median() # Negative controls

# Prepare data
ARG_bt_res$GENE
tmp <- read.csv("../RESULTS/final_output_files/resfinder_phenotypes.txt", header = TRUE, sep = "\t", fill  = TRUE)
tmp$GENE <- tmp$Gene_accession.no

ARG_blanks <- merge(ARG_bt_res, tmp, by = "GENE") %>% dplyr::select(matches(c("Gene_accession", "BLANK", "Blank")))

metaphlan_table <- read.table("../RESULTS/final_output_files/metaphlan3_joined_bugs_list_species.tsv",
                              header = TRUE, row.names  = 1, check.names  = FALSE, sep = "\t")

metaphlan_blanks <- metaphlan_table %>% dplyr::select(matches(c("BLANK", "Blank")))
metaphlan_blanks$taxa <- rownames(metaphlan_blanks)


library(ggplot2)
library(dplyr)
library(writexl)

xs <- list()
xs[["MetaPhlAn taxa in neg. controls"]] <- metaphlan_blanks

xs[["ARGs in neg. controls"]] <- ARG_blanks

write_xlsx(xs, path = "../RESULTS/R_results/DataS4.xlsx")



##### Antibiotic use versus antibiotic class correlations ####
# Sanity checks for whether antibiotic purchases 
# correlate with antibiotic resistance of the respective classes

ARG <- altExp(TSE, "arg", withColData = TRUE)

# Add resistance class for prelavent ARGs

altExp(ARG, "ClassPrevalent", withColData = TRUE) <-
  agglomerateByPrevalence(
    ARG,
    rank = "Class",
    assay.type = "counts",
    detection = 0.1 / 100,
    prevalence = 5 / 100
  )
altExp(ARG, "ClassPrevalent")

# Extract abundance data (assuming "relabundance" contains the abundance values)
arg_abundance <-
  as.data.frame(assay(
    altExp(ARG, "ClassPrevalent", withColData = TRUE),
    "counts"
  ))

# Transpose so samples are rows and ARG classes are columns
arg_abundance <- t(arg_abundance)

# Convert to DataFrame (compatible with colData)
arg_abundance_df <- DataFrame(arg_abundance)

# Assign to colData
colData(TSE) <- cbind(colData(TSE), arg_abundance_df)

df <- colData(TSE) %>% as.data.frame()

# Libraries
library(broom)
library(purrr)
library(writexl)

# Add pseudo-count (0.5 * minimum nonzero value) to each antibiotic resistance variable
df$Tetracycline <- df$Tetracycline + 0.5 * min(df$Tetracycline[df$Tetracycline > 0], na.rm = TRUE)
df$Beta.lactam <- df$Beta.lactam + 0.5 * min(df$Beta.lactam[df$Beta.lactam > 0], na.rm = TRUE)
df$Folate.pathway.antagonist <- df$Folate.pathway.antagonist + 
  0.5 * min(df$Folate.pathway.antagonist[df$Folate.pathway.antagonist > 0], na.rm = TRUE)
df$Macrolide..Lincosamide..Streptogramin.B <- df$Macrolide..Lincosamide..Streptogramin.B + 
  0.5 * min(df$Macrolide..Lincosamide..Streptogramin.B[df$Macrolide..Lincosamide..Streptogramin.B > 0], na.rm = TRUE)

# Fit models using log10 transformation
models <- list(
  "Tetracycline" = lm(log10(Tetracycline) ~ PREVAL_RX_J01A_NEVT, data = df),
  "Beta-lactam_J01C" = lm(log10(Beta.lactam) ~ PREVAL_RX_J01C_NEVT, data = df),
  "Beta-lactam_J01D" = lm(log10(Beta.lactam) ~ PREVAL_RX_J01D_NEVT, data = df),
  "Sulphonamide" = lm(log10(Folate.pathway.antagonist) ~ PREVAL_RX_J01E_NEVT, data = df),
  "Macrolide" = lm(log10(Macrolide..Lincosamide..Streptogramin.B) ~ PREVAL_RX_J01F_NEVT, data = df)
)

# Extract estimates, confidence intervals, and p-values
results <- purrr::map_df(models, ~ tidy(.x, conf.int = TRUE), .id = "Antibiotic")

# Adjust p-values using FDR correction and compute exponentiated estimates and CIs
results <- results %>%
  filter(term != "(Intercept)") %>%
  mutate(
    p_adj = p.adjust(p.value, method = "fdr"),
    exp_estimate = 10^estimate,  # Exponentiate estimate (inverse of log10)
    percent_change = (exp_estimate - 1) * 100,  # Convert to % change
    exp_conf_low = 10^conf.low,  # Exponentiate lower CI bound
    exp_conf_high = 10^conf.high  # Exponentiate upper CI bound
  )

# Select relevant columns
results <- results %>%
  dplyr::select(
    Antibiotic, term, estimate, conf.low, conf.high,  # Estimate & 95% CI
    exp_estimate, exp_conf_low, exp_conf_high,  # Exponentiated values & 95% CI
    percent_change, std.error, p.value, p_adj  # Additional model details
  )
# Write results to an Excel file
write_xlsx(results, "../RESULTS/R_results/TableS12.xlsx")


###### Linear models just women #######
# This runs the models just for women, the tables integrated with DataS1.xlsx manually

df_train_num <-
  df_train_num[train$MEN == 0, ] %>% dplyr::select(-matches("MEN"))

result <-
  do.call(rbind, lapply(1:ncol(df_train_num), function(x) {
    fit <-
      lm(
        log10(SUM_norm) ~ df_train_num[, x] +
          PREVAL_RX_J01_NEVT + PREVAL_RX_J01A_NEVT +
          PREVAL_RX_J01A + PREVAL_RX_J01F + BL_USE_RX_J01,
        data = df_train_num
      )
    lm_with_cov.result <- fit %>% summary()
    CI <- fit %>% confint()
    CI_low <- 10 ^ CI[2, 1]
    CI_up <- 10 ^ CI[2, 2]
    cov.estimate <- lm_with_cov.result$coefficients[2, 1]
    cov.pvalue <- lm_with_cov.result$coefficients[2, 4]
    return(
      data.frame(
        pvalue = cov.pvalue,
        Estimate = cov.estimate,
        CI2.5 = CI[2, 1],
        CI97.5 = CI[2, 2],
        exp_Estimate = 10 ^ cov.estimate,
        exp_CI2.5 = CI_low,
        exp_CI97.5 = CI_up
      )
    )
  }))
rownames(result) <- colnames(df_train_num)
result0 <- result

result <-
  result %>% t() %>%
  data.frame() %>%
  dplyr::select(-contains(
    c("sum", "SUM_norm", "scaled", "ARG", "_KAIKKI", "GAMMA", "RX_J")
  )) %>%
  t() %>%
  data.frame()
result$padj <-
  p.adjust(result$pvalue, method = "fdr")

result$VARIABLE <- rownames(result)
result <- result %>% mutate_if(is.numeric, signif, digits = 3)

df <-
  data.frame(VARIABLE = FR02names$VARIABLE, LONGNAME = FR02names$LONGNAME)

result <- dplyr::right_join(df, result, by = join_by(VARIABLE))

sigs <- result[result$padj < 0.05, ]
df <-
  sigs %>% replace_variable_with_dictionary("VARIABLE", substitutions) %>%
  dplyr::arrange(Estimate) %>%
  mutate(VARIABLE = factor(VARIABLE, levels = VARIABLE))

## Without adjusting for antibiotic use

result_unadj <-
  do.call(rbind, lapply(1:ncol(df_train_num), function(x) {
    fit <-
      lm(log10(SUM_norm) ~ df_train_num[, x],
         data = df_train_num)
    lm_without_cov.result <- fit %>% summary()
    CI <- fit %>% confint()
    CI_low <- 10 ^ CI[2, 1]
    CI_up <- 10 ^ CI[2, 2]
    cov.estimate <- lm_without_cov.result$coefficients[2, 1]
    cov.pvalue <- lm_without_cov.result$coefficients[2, 4]
    return(
      data.frame(
        pvalue = cov.pvalue,
        Estimate = cov.estimate,
        CI2.5 = CI[2, 1],
        CI97.5 = CI[2, 2],
        exp_Estimate = 10 ^ cov.estimate,
        exp_CI2.5 = CI_low,
        exp_CI97.5 = CI_up
      )
    )
  }))


rownames(result_unadj) <- colnames(df_train_num)

result_unadj <-
  result_unadj %>% t() %>%
  data.frame() %>%
  dplyr::select(-contains(c(
    "sum", "SUM_norm", "ARG", "_KAIKKI", "GAMMA", "scaled"
  ))) %>%
  t() %>%
  data.frame()
result_unadj$padj <-
  p.adjust(result_unadj$pvalue, method = "fdr")


result_unadj <-
  result_unadj %>% mutate_if(is.numeric, signif, digits = 3)

result_unadj$VARIABLE <- rownames(result_unadj)

df <-
  data.frame(VARIABLE = FR02names$VARIABLE, LONGNAME = FR02names$LONGNAME)

result_unadj <-
  dplyr::right_join(df, result_unadj, by = join_by(VARIABLE))


sigs <- result_unadj[result_unadj$padj < 0.05, ]  %>%
  dplyr::select(-contains(c("NEVT")))

df <-
  sigs %>% replace_variable_with_dictionary("VARIABLE", substitutions) %>%
  dplyr::arrange(Estimate) %>%
  mutate(VARIABLE = factor(VARIABLE, levels = VARIABLE))


## ARG Diversity

result_div <-
  do.call(rbind, lapply(1:ncol(df_train_num), function(x) {
    fit <-
      lm(ARG_div ~ df_train_num[, x],
         data = df_train_num)
    lm_with_cov.result <- fit %>% summary()
    CI <- fit %>% confint()
    CI_low <- CI[2, 1]
    CI_up <- CI[2, 2]
    cov.estimate <- lm_with_cov.result$coefficients[2, 1]
    # pvalue <- lm.result$coefficients[2,4]
    cov.pvalue <- lm_with_cov.result$coefficients[2, 4]
    # estimate <- lm.result$coefficients[2,1]
    return(
      data.frame(
        pvalue = cov.pvalue,
        Estimate = cov.estimate,
        CI2.5 = CI_low,
        CI97.5 = CI_up
      )
    )
  }))

rownames(result_div) <- colnames(df_train_num)



result_div <- result_div %>% t() %>% data.frame() %>%
  dplyr::select(-ends_with(c(
    "sum", "SUM_norm", "_KAIKKI", "GAMMA", "ARG_div", "ARG_obs"
  ))) %>%
  dplyr::select(-contains(c("HELIUS", "eae"))) %>% t()  %>% data.frame()

result_div$padj_cov <-
  p.adjust(result_div$pvalue, method = "fdr")

result_div$VARIABLE <- rownames(result_div)


result_div <-
  result_div %>% mutate_if(is.numeric, signif, digits = 3)

df <-
  data.frame(VARIABLE = FR02names$VARIABLE, LONGNAME = FR02names$LONGNAME)

result_div <-
  dplyr::right_join(df, result_div, by = join_by(VARIABLE))

library(writexl)

xs <- list()
x <- result[order(result$padj), ]
rownames(x) <- NULL
xs[["Women ARG load, controlled for AB"]] <- x

x <- result_unadj[order(result_unadj$padj), ]
rownames(x) <- NULL
xs[["Women ARG load, not controlled for AB"]] <- x

x <- result_div[order(result_div$padj), ]
rownames(x) <- NULL
xs[["Women ARG diversity, not controlled for AB"]] <- x

# To do multiple sheets just give the named list of data frames here (xs)
write_xlsx(xs, path = "../RESULTS/R_results/DATAS1_women.xlsx") # These sheets are copied to the DataS1 manually

####### Critical genes ####
# This code add gene level info to colData so that manual check for critical genes can be done
# for example df might contain blandm, blakpc, blaoxa48 or mcr which are critical ARGs.

taxonomyRanks((ARG)) # check the ranks contain gene

altExp(ARG, "GenePrevalent", withColData = TRUE) <-
  agglomerateByPrevalence(
    ARG,
    rank = "Gene",
    assay.type = "relabundance",
    detection = 0.00001 / 100,
    prevalence = 0.01 / 100
  )
altExp(ARG, "GenePrevalent")

# Extract abundance data (assuming "relabundance" contains the abundance values)
arg_abundance <-
  as.data.frame(assay(
    altExp(ARG, "GenePrevalent", withColData = TRUE),
    "relabundance"
  ))

# Transpose so samples are rows and ARG classes are columns
arg_abundance <- t(arg_abundance)

# Convert to DataFrame (compatible with colData)
arg_abundance_df <- DataFrame(arg_abundance)

# Assign to colData
colData(TSE) <- cbind(colData(TSE), arg_abundance_df)

df <- colData(TSE) %>% as.data.frame()

df %>% dplyr::select()


### Checks with using a subset of samples with more than 0.2M reads

source("GLMs_TreeSE.R")

#ARG load
# fit model with more than 0.2 M reads
fit_glm_subset <-
  glm(log10(SUM_norm) ~ .,
      data = df_train_glm[df_train_glm$V2>0.2,],
      family = "gaussian")
summary_glm_subset <- summary(fit_glm_subset, correlation=TRUE)

# Apply substitutions to rownames of plot_glm_table
plot_glm_table_subset <- fit_glm_subset %>%
  jtools::summ(confint = TRUE, exp = TRUE) %>%
  .$coeftable %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  replace_variable_with_dictionary("rowname", substitutions) %>%
  column_to_rownames(var = "rowname")

# Add column named Covariate
plot_glm_table_subset$Covariate <- rownames(plot_glm_table_subset)

# remove Intercept row
df <- plot_glm_table_subset[2:length(plot_glm_table_subset$Covariate), ]

# Sort
df <- df %>% arrange(`exp(Est.)`) %>%
  mutate(Covariate=factor(Covariate, levels=unique(Covariate)))


# Plot forest plot
forestplot.without.bacfam_subet <- ggplot(df, aes(y  = Covariate, x  = `exp(Est.)`))  +
  geom_point(shape = 18, size = 5)  +
  geom_errorbarh(aes(xmin  = `2.5%`, xmax  = `97.5%`), height  = 0.25, color="black")  +
  geom_vline(
    xintercept  = 1,
    color  = "darkgray",
    linetype  = "dashed",
    cex  = 2
  )  +
  # labs(x="Coefficient (95% CI)", y="", title="Drivers of ARG load") +
  theme_classic()  +
  theme(
    panel.border  = element_blank(),
    panel.background  = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line  = element_line(colour  = "black"),
    axis.text.y  = element_text(size = 0.9*base.size, colour  = "black"),
    axis.text.x.bottom  = element_text(size  = base.size, colour  = "black"),
    axis.title.x  = element_text(size  = base.size, colour  = "black"),
    title  = element_text(size  = base.size, colour  = "black")
  )

print(forestplot.without.bacfam_subet)


# ARG div
# fit model with more than 0.2 M reads
fit_glm_subset_div <-
  glm(ARG_div ~ .,
      data = df_train_glm_div[df_train_glm_div$V2>0.2,],
      family = "gaussian")
summary_glm_subset <- summary(fit_glm_subset_div, correlation=TRUE)

# Apply substitutions to rownames of plot_glm_table
plot_glm_table_subset_div <- fit_glm_subset_div %>%
  jtools::summ(confint = TRUE, exp = TRUE) %>%
  .$coeftable %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  replace_variable_with_dictionary("rowname", substitutions) %>%
  column_to_rownames(var = "rowname")

# Add column named Covariate
plot_glm_table_subset_div$Covariate <- rownames(plot_glm_table_subset_div)

# remove Intercept row
df <- plot_glm_table_subset_div[2:length(plot_glm_table_subset_div$Covariate), ]

# Sort
df <- df %>% arrange(`exp(Est.)`) %>%
  mutate(Covariate=factor(Covariate, levels=unique(Covariate)))


# Plot forest plot
forestplot.without.bacfam_subet_div <- ggplot(df, aes(y  = Covariate, x  = `exp(Est.)`))  +
  geom_point(shape = 18, size = 5)  +
  geom_errorbarh(aes(xmin  = `2.5%`, xmax  = `97.5%`), height  = 0.25, color="black")  +
  geom_vline(
    xintercept  = 1,
    color  = "darkgray",
    linetype  = "dashed",
    cex  = 2
  )  +
  # labs(x="Coefficient (95% CI)", y="", title="Drivers of ARG load") +
  theme_classic()  +
  theme(
    panel.border  = element_blank(),
    panel.background  = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line  = element_line(colour  = "black"),
    axis.text.y  = element_text(size = 0.9*base.size, colour  = "black"),
    axis.text.x.bottom  = element_text(size  = base.size, colour  = "black"),
    axis.title.x  = element_text(size  = base.size, colour  = "black"),
    title  = element_text(size  = base.size, colour  = "black")
  )

print(forestplot.without.bacfam_subet_div)
cowplot::plot_grid(forestplot.without.bacfam_subet, 
                   forestplot.without.bacfam_subet_div, labels="auto",
                   ncol=2, label_size=40)



pdf("../RESULTS/R_results/FigS10.pdf", width=30, height=12)
cowplot::plot_grid(forestplot.without.bacfam_subet, 
                   forestplot.without.bacfam_subet_div, labels="auto",
                   ncol=2, label_size=40)

dev.off()

s <- 350
library(Cairo)
CairoJPEG("../RESULTS/R_results/FigureS10.jpg", width=4*s, height=1.7*s)
cowplot::plot_grid(forestplot.without.bacfam_subet, 
                   forestplot.without.bacfam_subet_div, labels="auto",
                   ncol=2)
dev.off()

blaOXA.486