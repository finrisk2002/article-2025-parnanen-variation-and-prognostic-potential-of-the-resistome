# glmboost
# Load required libraries
library(mboost)
library(caret)
library(pdp)
library(mgcv)
library(mia)
library(dplyr)
library(jtools)
set.seed(481311)

# TSE <- readRDS("../data/TSE.rds")

# Pick the ARG data for training and test set
set <- altExp(TSE, "arg", withColData = TRUE)
train <- set[, TSE$train.sample]
test <- set[,!TSE$train.sample]

fitControl <- trainControl(method = "cv")

df2 <-
  cbind(as.data.frame(colData(TSE)), as.data.frame(t(assay(
    altExp(TSE, "FamilyPrevalent"), "relabundance"
  ))))

train <-
  df2[df2$train.sample == TRUE, ] %>%  mutate_at(vars(matches(c(
    "ceae", "deae", "classified"
  ))),
  ~ log10(. + min(.[. > 0]) / 2)) %>%
  dplyr::select(-where( ~ any(is.infinite(.))))

test <-
  df2[df2$train.sample == FALSE, ] %>%  mutate_at(vars(matches(c(
    "ceae", "deae", "classified"
  ))),
  ~ log10(. + min(.[. > 0]) / 2)) %>%
  dplyr::select(-where( ~ any(is.infinite(.))))

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
df_train_num$V2 <-
  df_train_num$V2 / 1000000 # Library size - transform to Million reads

df_test_num <- test   %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(matches(selected.vars)) %>%
  dplyr::select(-matches(c("scaled")))

row.names(df_test_num) <- test$Row.names

df_test_num$vaesto <- log10(df_test_num$vaesto)
df_test_num$V2 <-
  df_test_num$V2 / 1000000 # Library size - transform to Million reads

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
    data = df_train_num %>% dplyr::select(-ends_with(
      c("ARG_div", "ARG_obs", "KOULGR", "PREVAL_HIBP", "PREVAL_IHD")
    )) ,
    method  = "glmboost",
    tuneLength = 5,
    center = TRUE,
    na.action  = na.omit
  )

# Train GLMBoost model for ARG_div
glmBoostModel_div <-
  caret::train(
    ARG_div  ~  .,
    data = df_train_num %>% dplyr::select(-matches(
      c("SUM_norm", "ARG_obs", "KOULGR", "PREVAL_HIBP", "PREVAL_IHD")
    )),
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
fit_glm <-
  glm(log10(SUM_norm) ~ .,
      data = df_train_glm,
      family = "gaussian")
summary_glm <- summary(fit_glm, correlation = TRUE)
library(jtools)
summ(fit_glm, digits = 3, confint = TRUE)


df_train_glm <- df_train_num %>%
  dplyr::select(all_of(covariates)) %>% as.data.frame() %>% dplyr::select(-matches(c("KY100_16", "KY100_9"))) # Remove foods which correlate


df_train_glm_div <- df_train_num %>%
  dplyr::select(all_of(covariates_div)) %>% as.data.frame()



df_test_glm <- df_test_num %>%
  dplyr::select(all_of(covariates)) %>% as.matrix %>% as.data.frame



fit_glm_div <-
  glm(ARG_div ~ .,
      data = df_train_glm_div,
      family = "gaussian")
summary_glm <- summary(fit_glm_div)
library(jtools)
summ(fit_glm_div, digits = 3, confint = TRUE)

##### R2 calculations ####
test <- test
#############
# R2 ARG div +family div
fit_glm_div_family_only <-
  glm(ARG_div ~ Species_diversity,
      data = df_train_glm_div,
      family = "gaussian")

# Predicted vs observed plot
predictions <-
  predict(fit_glm_div_family_only , test) %>% data.frame()

rownames(test) <- test$Row.names
predictions <-
  predict(fit_glm_div_family_only , test[test %>% dplyr::select(matches(c("ARG_div", "Species_diversity"))) %>% complete.cases(), ]) %>% data.frame
test.target <- subset(test, select = ARG_div)
test.target <- test.target %>%
  t %>% subset(select = rownames(predictions)) %>% t

# RMSE
sqrt(mean((test.target[, 1] - predictions[, 1]) ^ 2))
# R2
cor(test.target, predictions) ^ 2


# ARG load Enterosignatures#

fit_glm_nmf <-
  glm((SUM_norm) ~ nmf1 + nmf2 + nmf3 + nmf4 + nmf5,
      data = train ,
      family = "gaussian")
summary_glm <- summary(fit_glm_nmf)
library(jtools)
summ(fit_glm_nmf, digits = 3, confint = TRUE)

# Predicted vs observed plot
predictions <- predict(fit_glm_nmf, test) %>% data.frame()

test.target <- subset(test, select = SUM_norm)
test.target <- test.target %>%
  t %>% subset(select = rownames(predictions)) %>% t

# RMSE
sqrt(mean((test.target[, 1] - predictions[, 1]) ^ 2))
# R2
cor(test.target, predictions) ^ 2

######

# ARG load family diversity
fit_glm_div_family_only <-
  glm(SUM_norm ~ Species_diversity,
      data = train ,
      family = "gaussian")

# Predicted vs observed plot
predictions <-
  predict(fit_glm_div_family_only , test) %>% data.frame()

rownames(test) <- test$Row.names
predictions <-
  predict(fit_glm_div_family_only , test[test %>% dplyr::select(matches(c("SUM_norm", "Species_diversity"))) %>% complete.cases(), ]) %>% data.frame
test.target <- subset(test, select = SUM_norm)
test.target <- test.target %>%
  t %>% subset(select = rownames(predictions)) %>% t

# RMSE
sqrt(mean((test.target[, 1] - predictions[, 1]) ^ 2))
# R2
cor(test.target, predictions) ^ 2



########
# Enterosignatures diversity#

fit_glm_nmf <-
  glm((ARG_div) ~ nmf1 + nmf2 + nmf3 + nmf4 + nmf5,
      data = train ,
      family = "gaussian")
summary_glm <- summary(fit_glm_nmf)
library(jtools)
summ(fit_glm_nmf, digits = 3, confint = TRUE)

# Predicted vs observed plot
predictions <- predict(fit_glm_nmf, test) %>% data.frame()

test.target <- subset(test, select = SUM_norm)
test.target <- test.target %>%
  t %>% subset(select = rownames(predictions)) %>% t

# RMSE
sqrt(mean((test.target[, 1] - predictions[, 1]) ^ 2))
# R2
cor(test.target, predictions) ^ 2
#

## East/West ##

fit_glm_east <-
  glm(log10(SUM_norm) ~ EAST,
      data = train ,
      family = "gaussian")
summary(fit_glm_east)
library(jtools)
summ(fit_glm_east, digits = 3, confint = TRUE)

predictions <- predict(fit_glm_east, test) %>% data.frame()

test.target <- subset(test, select = SUM_norm)
test.target <- test.target %>%
  t %>% subset(select = rownames(predictions)) %>% t

# RMSE
sqrt(mean((test.target[, 1] - predictions[, 1]) ^ 2))
# R2
cor(test.target, predictions) ^ 2


### Diet ##

diet_train <-
  train  %>% dplyr::select(matches(c("KY100", "SUM_norm")))

# Fit GLM

fit_glm_diet <-
  glm(log10(SUM_norm) ~ .,
      data = diet_train,
      family = "gaussian")


# Predicted vs observed plot
diet_test <- test %>% dplyr::select(matches(c("KY100", "SUM_norm")))
diet_test <- diet_test[complete.cases(diet_test), ]

predictions <- predict(fit_glm_diet, diet_test) %>% data.frame()

test.target <- subset(diet_test, select = SUM_norm)
test.target <- test.target %>%
  t %>% subset(select = rownames(predictions)) %>% t

# RMSE
sqrt(mean((test.target[, 1] - predictions[, 1]) ^ 2))
# R2
cor(test.target, predictions) ^ 2

####

# Only antibiotics
fit_glm_ab <-
  glm(log10(SUM_norm) ~ .,
      data = train  %>% dplyr::select(matches(c(
        "SUM_norm", "PREVAL_RX_J01"
      ))),
      family = "gaussian")
summary_glm <- summary(fit_glm_ab)
library(jtools)
summ(fit_glm_ab, digits = 3, confint = TRUE)

# Predicted vs observed plot
predictions <- predict(fit_glm_ab, test) %>% data.frame()

test.target <- subset(test, select = SUM_norm)
test.target <- test.target %>%
  t %>% subset(select = rownames(predictions)) %>% t

# RMSE
sqrt(mean((test.target[, 1] - predictions[, 1]) ^ 2))
# R2
cor(test.target, predictions) ^ 2


################ GLM with just bacterial family

df2 <-
  cbind(as.data.frame(colData(TSE)), as.data.frame(t(assay(
    altExp(TSE, "FamilyPrevalent"), "relabundance"
  ))))

train <- df2[df2$train.sample == TRUE, ]
test <- df2[df2$train.sample == FALSE, ]

df_train_num_taxa <- train %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(matches(c("SUM_norm", "ae", "classified"))) %>%
  dplyr::select(-matches(c("vaesto"))) %>%
  mutate_at(vars(matches(c(
    "ceae", "deae", "classified"
  ))),
  ~ log10(. + min(.[. > 0]) / 2)) %>%
  dplyr::select(-where( ~ any(is.infinite(.))))

fit_glm_family <-
  glm(log10(SUM_norm) ~ .,
      data = df_train_num_taxa,
      family = "gaussian")
summary_glm <- summary(fit_glm_family)
library(jtools)
summ(fit_glm_family, digits = 3, confint = TRUE)

### Predicted vs observed for only bacterial family

# Predicted vs observed plot
predictions <-
  predict(fit_glm_family, test, na.action = na.omit) %>% data.frame()

test.target <- subset(test, select = SUM_norm)
test.target <- test.target %>%
  t %>% subset(select = rownames(predictions)) %>% t

# RMSE
sqrt(mean((test.target[, 1] - predictions[, 1]) ^ 2))
# R2
cor(test.target, predictions) ^ 2

###########################################

####### GLM boost with all other covariates AND bacterial family

df2 <-
  cbind(as.data.frame(colData(TSE)), as.data.frame(t(assay(
    altExp(TSE, "FamilyPrevalent"), "relabundance"
  ))))
train <- df2[df2$train.sample == TRUE, ]
test <- df2[df2$train.sample == FALSE, ]

df2

df_train_num_taxa <- train %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(matches(c(selected.vars, "ae", "classified"))) %>%
  dplyr::select(-matches(c(
    "scaled", "^MENTALDIS", "Helicobacteraceae", "KY100_"
  ))) %>%
  mutate_at(vars(matches(c(
    "ceae", "deae", "classified"
  ))),
  ~ log10(. + min(.[. > 0]) / 2)) %>%
  dplyr::select(-where( ~ any(is.infinite(.))))

df_train_num_taxa$KY100_14 <-
  train$KY100_14 # Add back poultry and fresh veg
df_train_num_taxa$KY100_21 <- train$KY100_21

df_test_num_taxa <- test %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(matches(c(selected.vars, "ae", "classified"))) %>%
  dplyr::select(-matches(c(
    "scaled", "^MENTALDIS", "Helicobacteraceae", "KY100_"
  ))) %>%
  mutate_at(vars(matches(c(
    "ceae", "deae", "classified"
  ))),
  ~ log10(. + min(.[. > 0]) / 2)) %>%
  dplyr::select(-where( ~ any(is.infinite(.))))

df_test_num_taxa$KY100_14 <-
  test$KY100_14 # Add back poultry and fresh veg
df_test_num_taxa$KY100_21 <- test$KY100_21

# Set row names and remove redundant column
row.names(df_train_num_taxa) <- df_train_num_taxa$Row.names
df_train_num_taxa <-
  df_train_num_taxa %>% dplyr::select(-matches("Row.names"))

row.names(df_test_num_taxa) <- df_test_num_taxa$Row.names
df_test_num_taxa <-
  df_test_num_taxa %>% dplyr::select(-matches("Row.names"))

# Train GLMBoost model for SUM_norm
glmBoostModel_family <-
  caret::train(
    log10(SUM_norm)  ~  .,
    data = df_train_num_taxa %>% dplyr::select(-matches(
      c("ARG_div", "ARG_obs", "PREVAL_HIBP", "PREVAL_IHD")
    )),
    method  = "glmboost",
    tuneLength = 5,
    center = TRUE,
    na.action  = na.omit
  )

# Model summary
glmboostsummary_family <- summary(glmBoostModel_family)
glmboostsummary_family


# Merge tables for coefficients
cf <-
  summary(glmBoostModel_family)$object %>%
  coef(off2int = TRUE) %>%
  as.list %>%
  as.matrix %>%
  as.data.frame %>%
  dplyr::select(-(matches("Intercept")))
selprob <- summary(glmBoostModel_family)$selprob  %>% as.data.frame
coefs_table <- merge(cf, selprob, by = 0)


colnames(coefs_table) <-
  c("Variable", "Coefficient", "Selection probablity")

coefs_table

# Regular GLM with glmboost covariates

covariates_family <-
  c(coefs_table$Variable, "SUM_norm")

df_train_glm_family <- df_train_num_taxa %>%
  dplyr::select(all_of(covariates_family)) %>% as.data.frame()


fit_glm_family <-
  glm(log10(SUM_norm) ~ .,
      data = df_train_glm_family,
      family = "gaussian")
summary_glm <- summary(fit_glm_family)
library(jtools)
summ(fit_glm_family, digits = 3, confint = TRUE)


#########################

# Bacterial family tests

# Diet
df_fam_test <- df_train_num_taxa %>%
  dplyr::select(-matches(c(
    "scaled", "MENTALDIS", "Helicobacteraceae"
  ))) %>%
  dplyr::select(matches(
    c(
      "ae",
      "classified",
      "KY100_21",
      "KY100_14",
      "SUM_norm",
      "PREVAL_RX_J01"
    )
  )) %>% na.omit()

fit_glm_family_test_diet <-
  glm(log10(SUM_norm) ~ .,
      data = df_fam_test,
      family = "gaussian") %>% step
summary_glm <- summary(fit_glm_family_test_diet, correlation = TRUE)
library(jtools)
summ(fit_glm_family_test_diet,
     digits = 3,
     confint = TRUE)

# Geo and demography
df_fam_test <- df_train_num_taxa %>%
  dplyr::select(-matches(c(
    "scaled", "MENTALDIS", "Helicobacteraceae"
  ))) %>%
  dplyr::select(matches(
    c(
      "ae",
      "classified",
      "MEN",
      "TULOT",
      "EAST",
      "SUM_norm",
      "PREVAL_RX_J01"
    )
  )) %>% na.omit()

fit_glm_family_test_geo_demo <-
  glm(log10(SUM_norm) ~ .,
      data = df_fam_test,
      family = "gaussian") %>% step
summary_glm <-
  summary(fit_glm_family_test_geo_demo, correlation = TRUE)
library(jtools)
summ(fit_glm_family_test_geo_demo,
     digits = 3,
     confint = TRUE)


############## GLM boost with meta4 families ########
altExp(TSE, "metaphlan4FamilyPrevalent") <-
  mergeFeaturesByPrevalence(
    altExp(TSE, "metaphlan4")[rowData(altExp(TSE, "metaphlan4"))$kingdom == "k__Bacteria"],
    rank =
      "family",
    assay.type = "relabundance",
    detection = 0.1 / 100,
    prevalence = 1 / 100
  )


df2 <-
  cbind(as.data.frame(colData(TSE)), as.data.frame(t(assay(
    altExp(TSE, "metaphlan4FamilyPrevalent"), "relabundance"
  ))))
train <- df2[df2$train.sample == TRUE, ]
test <- df2[df2$train.sample == FALSE, ]

df2

df_train_num_taxa_meta4 <- train %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(matches(c(selected.vars, "ae", "classified"))) %>%
  dplyr::select(-matches(c(
    "scaled", "^MENTALDIS", "Helicobacteraceae", "KY100_"
  ))) %>%
  mutate_at(vars(matches(c(
    "ceae", "deae", "classified"
  ))),
  ~ log10(. + min(.[. > 0]) / 2)) %>%
  dplyr::select(-where( ~ any(is.infinite(.))))

df_train_num_taxa_meta4$KY100_14 <-
  train$KY100_14 # Add back poultry and fresh veg
df_train_num_taxa_meta4$KY100_21 <- train$KY100_21

df_test_num_taxa_meta4 <- test %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(matches(c(selected.vars, "ae", "classified"))) %>%
  dplyr::select(-matches(c(
    "scaled", "^MENTALDIS", "Helicobacteraceae", "KY100_"
  ))) %>%
  mutate_at(vars(matches(c(
    "ceae", "deae", "classified"
  ))),
  ~ log10(. + min(.[. > 0]) / 2)) %>%
  dplyr::select(-where( ~ any(is.infinite(.))))

df_test_num_taxa_meta4$KY100_14 <-
  test$KY100_14 # Add back poultry and fresh veg
df_test_num_taxa_meta4$KY100_21 <- test$KY100_21

# Set row names and remove redundant column
row.names(df_train_num_taxa_meta4) <-
  df_train_num_taxa_meta4$Row.names
df_train_num_taxa_meta4 <-
  df_train_num_taxa_meta4 %>% dplyr::select(-matches("Row.names"))

row.names(df_test_num_taxa_meta4) <-
  df_test_num_taxa_meta4$Row.names
df_test_num_taxa_meta4 <-
  df_test_num_taxa_meta4 %>% dplyr::select(-matches("Row.names"))

# Train GLMBoost model for SUM_norm
glmBoostModel_family_meta4 <-
  caret::train(
    log10(SUM_norm)  ~  .,
    data = df_train_num_taxa_meta4 %>% dplyr::select(-matches(
      c("ARG_div", "ARG_obs", "PREVAL_HIBP", "PREVAL_IHD")
    )),
    method  = "glmboost",
    tuneLength = 5,
    center = TRUE,
    na.action  = na.omit
  )

# Model summary
glmboostsummary_family_meta4 <- summary(glmBoostModel_family_meta4)
glmboostsummary_family_meta4


# Merge tables for coefficients
cf <-
  summary(glmBoostModel_family_meta4)$object %>%
  coef(off2int = TRUE) %>%
  as.list %>%
  as.matrix %>%
  as.data.frame %>%
  dplyr::select(-(matches("Intercept")))
selprob <-
  summary(glmBoostModel_family_meta4)$selprob  %>% as.data.frame
coefs_table <- merge(cf, selprob, by = 0)


colnames(coefs_table) <-
  c("Variable", "Coefficient", "Selection probablity")

coefs_table

# Regular GLM with glmboost covariates

covariates_family_meta4 <-
  c(coefs_table$Variable, "SUM_norm")

df_train_glm_family_meta4 <- df_train_num_taxa_meta4 %>%
  dplyr::select(all_of(covariates_family_meta4)) %>% as.data.frame()

# Clean names
colnames(df_train_glm_family_meta4) <-
  gsub("^f__", "", colnames(df_train_glm_family_meta4))


fit_glm_family_meta4 <-
  glm(log10(SUM_norm) ~ .,
      data = df_train_glm_family_meta4,
      family = "gaussian")
summary_glm <- summary(fit_glm_family_meta4)
library(jtools)

library(writexl)

# Extract GLM summary as a data frame
glm_summary_df <-
  summ(
    fit_glm_family_meta4,
    digits = 3,
    confint = TRUE,
    output = "data.frame"
  )$coeftable %>% as.data.frame()
glm_summary_df$Variable <- rownames(glm_summary_df)
coeftable <-
  replace_variable_with_dictionary(glm_summary_df, "Variable", substitutions)
# Save to an Excel file
write_xlsx(coeftable, "../RESULTS/R_results/TableS11.xlsx")

# Extract GLM summary as a data frame
glm_summary_df <- summ(
  fit_glm_family,
  digits = 3,
  confint = TRUE,
  output = "data.frame"
)$coeftable %>% as.data.frame()
glm_summary_df$Variable <- rownames(glm_summary_df)
coeftable <-
  replace_variable_with_dictionary(glm_summary_df, "Variable", substitutions)
# Save to an Excel file
write_xlsx(coeftable, "../RESULTS/R_results/TableS10.xlsx")

# Check hospital exposure
source("Hospital_exposure.R")
