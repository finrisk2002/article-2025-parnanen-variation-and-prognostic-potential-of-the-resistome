source("loadData.R")
source("GLMs_TreeSE.R")
source("R_functions.R")

FR02names <-
  read.table(
    "../RESULTS/final_output_files/FR02_pheno_names.txt",
    sep = "\t",
    fill = TRUE,
    header = TRUE
  )

# train and test are created in GLMs_TreeSE.R
# we use the same here for consistency and to avoid
# data leaking in CV

# Linear regression for ARG load and metadata
result <-
  do.call(rbind, lapply(1:ncol(df_train_num), function(x) {
    fit <-
      lm(log10(SUM_norm) ~ df_train_num[, x]+
           PREVAL_RX_J01_NEVT + PREVAL_RX_J01A_NEVT +
           PREVAL_RX_J01A + PREVAL_RX_J01F + BL_USE_RX_J01,
         data = df_train_num)
    lm_with_cov.result <- fit %>% summary()
    CI <- fit %>% confint()
    CI_low <- 10^CI[2, 1]
    CI_up <- 10^CI[2, 2]
    cov.estimate <- lm_with_cov.result$coefficients[2, 1]
    cov.pvalue <- lm_with_cov.result$coefficients[2, 4]
    return(
      data.frame(
        pvalue = cov.pvalue,
        Estimate = cov.estimate,
        CI2.5 = CI[2, 1],
        CI97.5 = CI[2, 2],
        exp_Estimate = 10^cov.estimate,
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
  dplyr::select(-contains(c("sum", "SUM_norm", "scaled", "ARG", "_KAIKKI", "GAMMA", "RX_J"))) %>%
  t() %>%
  data.frame()
result$padj <-
  p.adjust(result$pvalue, method = "fdr")

result$VARIABLE <- rownames(result)
result <- result %>% mutate_if(is.numeric, signif, digits =3)

df <-
  data.frame(VARIABLE = FR02names$VARIABLE, LONGNAME = FR02names$LONGNAME)

result <- dplyr::right_join(df, result, by = join_by(VARIABLE))

sigs <- result[result$padj<0.05,] 
df <- sigs %>% replace_variable_with_dictionary("VARIABLE", substitutions) %>% 
  dplyr::arrange(Estimate) %>% 
  mutate(VARIABLE=factor(VARIABLE, levels=VARIABLE)) 
  
## Without adjusting for antibiotic use

result_unadj <-
  do.call(rbind, lapply(1:ncol(df_train_num), function(x) {
    fit <-
      lm(log10(SUM_norm) ~ df_train_num[, x],
         data = df_train_num)
    lm_without_cov.result <- fit %>% summary()
    CI <- fit %>% confint()
    CI_low <- 10^CI[2, 1]
    CI_up <- 10^CI[2, 2]
    cov.estimate <- lm_without_cov.result$coefficients[2, 1]
    cov.pvalue <- lm_without_cov.result$coefficients[2, 4]
    return(
      data.frame(
        pvalue = cov.pvalue,
        Estimate = cov.estimate,
        CI2.5 = CI[2, 1],
        CI97.5 = CI[2, 2],
        exp_Estimate = 10^cov.estimate,
        exp_CI2.5 = CI_low,
        exp_CI97.5 = CI_up
      )
    )
  }))


rownames(result_unadj) <- colnames(df_train_num)

result_unadj <-
  result_unadj %>% t() %>%
  data.frame() %>%
  dplyr::select(-contains(c("sum", "SUM_norm", "ARG","_KAIKKI", "GAMMA", "scaled"))) %>%
  t() %>%
  data.frame()
result_unadj$padj <-
  p.adjust(result_unadj$pvalue, method = "fdr")


result_unadj <- result_unadj %>% mutate_if(is.numeric, signif, digits =3)

result_unadj$VARIABLE <- rownames(result_unadj)

df <-
  data.frame(VARIABLE = FR02names$VARIABLE, LONGNAME = FR02names$LONGNAME)

result_unadj <-
  dplyr::right_join(df, result_unadj, by = join_by(VARIABLE))


sigs <- result_unadj[result_unadj$padj<0.05,]  %>%
  dplyr::select(-contains(c("NEVT")))
  
df <- sigs %>% replace_variable_with_dictionary("VARIABLE", substitutions) %>% 
  dplyr::arrange(Estimate) %>% 
  mutate(VARIABLE=factor(VARIABLE, levels=VARIABLE)) 


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


result_div <- result_div %>% mutate_if(is.numeric, signif, digits =3)

df <-
  data.frame(VARIABLE = FR02names$VARIABLE, LONGNAME = FR02names$LONGNAME)

result_div <-
  dplyr::right_join(df, result_div, by = join_by(VARIABLE))


# ------------------------------------------------------------------------

# TODO: would be better to move to Carpentry to avoid recalculations?
# Generates the table for Family - ARG load associations
xxx <- getCrossAssociation(altExp(TSE, "FamilyPrevalent", withColData=TRUE),
                           assay.type1="counts",
		           col.var2=c("SUM_norm", "ARG_div", "Species_diversity"), 
                           method="kendall",
			   mode="table",
			   test.signif=TRUE,
			   p.adj.method="fdr"
			   )
xxx <- xxx[, !colnames(xxx)=="pval"]
colnames(xxx) <- c("Family", "Variable", "Tau", "FDR")
xxx$Tau <- round(xxx$Tau, 3)
xxx$FDR <- round(xxx$FDR, 3)
xxx <- xxx[!xxx$Family=="Other",]
xxx$Variable <- str_replace(xxx$Variable, "ARG_div", "ARG diversity")
xxx$Variable <- str_replace(xxx$Variable, "SUM_norm", "ARG load")
xxx$Variable <- str_replace(xxx$Variable, "Species_diversity", "Species diversity")
xxx <- xxx[, c(2, 1, 3, 4)] %>% arrange(Variable, Family)
Family_ARG_associations <- xxx
saveRDS(Family_ARG_associations, file="Family_ARG_associations.rds")

# ------------------------------------------------------------------------

# Let us use writexl library as it is free of Java dependencies and works better across multiple systems
# TODO: manuscript Tables should be numbered instead of naming sheets

library(writexl)

xs <- list()
x <- result[order(result$padj),]; rownames(x) <- NULL
xs[["ARG load, controlled for AB"]] <- x

x <- result_unadj[order(result_unadj$padj),]; rownames(x) <- NULL
xs[["ARG load, not controlled for AB"]] <- x

x <- result_div[order(result_div$padj),]; rownames(x) <- NULL
xs[["ARG diversity, not controlled for AB"]] <- x

x <- Family_ARG_associations; rownames(x) <- NULL
xs[["Family associations with ARG"]] <- x

# To do multiple sheets just give the named list of data frames here (xs)
write_xlsx(xs, path = "../RESULTS/R_results/DataS1.xlsx")



