# analyse hospital esposure in the last year, manually add the result for K27_bin to text.


df2 <- cbind(as.data.frame(colData(TSE)), as.data.frame(t(assay(altExp(TSE, "FamilyPrevalent"), "relabundance"))))

train <- df2[df2$train.sample==TRUE,] %>%  mutate_at(vars(matches(c("ceae", "deae", "classified"))),
                                                     ~ log10(. + min(.[. > 0]) / 2)) %>% 
  dplyr::select(-where(~ any(is.infinite(.))))

test <- df2[df2$train.sample==FALSE,] %>%  mutate_at(vars(matches(c("ceae", "deae", "classified"))),
                                                     ~ log10(. + min(.[. > 0]) / 2)) %>% 
  dplyr::select(-where(~ any(is.infinite(.))))

# Create a character vector of ATC codes for later use
ATC_1level <- paste("RX_", LETTERS, sep = "")

selected.vars <- c(
  "PREVAL_",
  "BL_",
  "KY100_",
  "vaesto",
  "Species_diversity",
  "Enterobacteriaceae",
  "^MEN",
  "family",
  "BL_AGE",
  "TULOT",
  "BMI",
  "KOL",
  "SUM_norm",
  "KOULGR",
  "HFC_SCORE",
  "EAST",
  "ARG",
  "V2" # Library size
)

# New
df_train_num <- train   %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(matches(selected.vars)) %>%
  dplyr::select(-matches(c("scaled", "^MENTALDIS")))
df_train_num$K27_bin <- train$K27_bin %>% as.numeric()
row.names(df_train_num) <- train$Row.names

df_train_num$vaesto <- log10(df_train_num$vaesto)
df_train_num$V2 <- df_train_num$V2/1000000 # Library size - transform to Million reads

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

