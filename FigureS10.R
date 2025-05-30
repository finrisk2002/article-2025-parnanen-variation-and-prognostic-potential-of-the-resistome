# Pick the ARG data for training and test set
set <- altExp(TSE, "arg", withColData = TRUE)
train <- set[, TSE$train.sample]
test <- set[, !TSE$train.sample]

fitControl <- trainControl(method = "cv")

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

row.names(df_train_num) <- train$Row.names

df_train_num$vaesto <- log10(df_train_num$vaesto)
df_train_num$V2 <- df_train_num$V2/1000000 # Library size - transform to Million reads

df_test_num <- test   %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(matches(selected.vars)) %>%
  dplyr::select(-matches(c("scaled")))

row.names(df_test_num) <- test$Row.names

df_test_num$vaesto <- log10(df_test_num$vaesto)
df_test_num$V2 <- df_test_num$V2/1000000 # Library size - transform to Million reads

# GLM boost on train set

# Set row names and remove redundant column
row.names(df_train_num) <- df_train_num$Row.names
df_train_num <-
  df_train_num %>% dplyr::select(-matches("Row.names"))

row.names(df_test_num) <- df_test_num$Row.names
df_test_num <- df_test_num %>% dplyr::select(-matches("Row.names"))

# Train GLMBoost model for SUM_norm
glmBoostModel <-
  caret::train(
    log10(SUM_norm)  ~  .,
    data = df_train_num %>% dplyr::select(-ends_with(c(
      "ARG_div", "ARG_obs", "KOULGR", "PREVAL_HIBP", "PREVAL_IHD"
    ))) ,
    method  = "glmboost",
    tuneLength = 5,
    center = TRUE,
    na.action  = na.omit
  )

# Train GLMBoost model for ARG_div
glmBoostModel_div <-
  caret::train(
    ARG_div  ~  .,
    data = df_train_num %>% dplyr::select(-matches(c(
      "SUM_norm", "ARG_obs", "KOULGR", "PREVAL_HIBP", "PREVAL_IHD"
    ))),
    method  = "glmboost",
    tuneLength = 5,
    center = TRUE,
    na.action  = na.omit
  )

# Model summary
glmboostsummary <- summary(glmBoostModel)
glmboostsummary

glmboostsummary_div <- summary(glmBoostModel_div)
glmboostsummary_div

# Merge tables for coefficients
cf <-
  summary(glmBoostModel)$object %>%
  coef(off2int = TRUE) %>%
  as.list %>%
  as.matrix %>%
  as.data.frame %>%
  dplyr::select(-(matches("Intercept")))
selprob <- summary(glmBoostModel)$selprob  %>% as.data.frame
coefs_table <- merge(cf, selprob, by = 0)

cf_div <-
  summary(glmBoostModel_div)$object %>% coef(off2int = TRUE)  %>% as.list %>% as.matrix %>% as.data.frame %>% dplyr::select(-(matches("Intercept")))
selprob <- summary(glmBoostModel_div)$selprob  %>% as.data.frame

coefs_table_div <- merge(cf_div, selprob, by = 0)

colnames(coefs_table) <-
  c("Variable", "Coefficient", "Selection probablity")
colnames(coefs_table_div) <-
  c("Variable", "Coefficient", "Selection probablity")
coefs_table
coefs_table_div


# Regular GLM with glmboost covariates

covariates <-
  c(coefs_table$Variable, "SUM_norm")
covariates_div <-
  c(coefs_table_div$Variable, "ARG_div")

df_train_glm <- df_train_num %>%
  dplyr::select(all_of(covariates)) %>% as.data.frame()
fit_glm_200 <-
  glm(log10(SUM_norm) ~ .,
      data = df_train_glm[df_train_glm$V2>0.2,],
      family = "gaussian")
summary_glm <- summary(fit_glm, correlation=TRUE)
library(jtools)
summ(fit_glm, digits = 3, confint = TRUE)

df_train_glm_div <- df_train_num %>%
  dplyr::select(all_of(covariates_div)) %>% as.data.frame()



df_test_glm <- df_test_num %>%
  dplyr::select(all_of(covariates)) %>% as.matrix %>% as.data.frame



fit_glm_div_200 <-
  glm(ARG_div ~ .,
      data = df_train_glm_div[df_train_glm_div$V2>0.2,],
      family = "gaussian")
summary_glm <- summary(fit_glm_div)
library(jtools)
summ(fit_glm_div, digits = 3, confint = TRUE)
