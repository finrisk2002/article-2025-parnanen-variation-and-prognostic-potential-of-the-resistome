library("tidyverse")
library("survival")
library("mia")
library("tidybayes")

# Read rds and run boosted GLM
TSE <- readRDS("../data/TSE.rds") 
# source("GLMs_TreeSE.R")
source("Dictionary.R")
source("R_functions.R")

# Prepare data frame for the survival analysis
# Note: logging and scaling for some variables!

# Causes of death
DF <- prepare.data.for.survival.analysis(TSE) 


df <- DF %>% 
  mutate(
    Infectious  = case_when(
      K_TPKS_level1 == "A" ~ 1,
      K_VKS_level1 == "A" ~ 1,!K_TPKS_level1 == "A" ~ 0,!K_VKS_level1 == "A" ~ 0
    )
  ) %>%
  mutate(
    Cardiovascular  = case_when(
      K_TPKS_level1 == "I" ~ 1,
      K_VKS_level1 == "I" ~ 1,!K_TPKS_level1 == "I" ~ 0,!K_VKS_level1 == "I" ~ 0
    )
  ) %>%
  mutate(
    Respiratory  = case_when(
      K_TPKS_level1 == "J" ~ 1,
      K_VKS_level1 == "J" ~ 1,!K_TPKS_level1 == "J" ~ 0,!K_VKS_level1 == "J" ~ 0
    ) ) %>%
  mutate(
    Pneumonia = case_when(
      K_VKS == "J189" ~ 1,
      K_VKS == "J159" ~ 1,
      K_VKS == "J181" ~ 1,
      K_VKS == "J188" ~ 1,
      K_VKS == "J180" ~ 1, !K_VKS == "J189" ~ 0,  !K_VKS == "J159" ~ 0,  !K_VKS == "J181" ~ 1, !K_VKS == "J188" ~ 1,
    )
  ) %>%
  mutate(
    Gastrointestinal  = case_when(
      K_TPKS_level1 == "K" ~ 1,
      K_VKS_level1 == "K" ~ 1,!K_TPKS_level1 == "K" ~ 0,!K_VKS_level1 == "K" ~ 0
    )
  ) %>%
  mutate(
    Cancer  = case_when(
      K_TPKS_level1 == "C" ~ 1,
      K_VKS_level1 == "C" ~ 1,!K_TPKS_level1 == "C" ~ 0,!K_VKS_level1 == "C" ~ 0
    )
  ) %>%
  mutate(
    Neurological  = case_when(
      K_TPKS_level1 == "G" ~ 1,
      K_VKS_level1 == "G" ~ 1,!K_TPKS_level1 == "G" ~ 0,!K_VKS_level1 == "G" ~ 0
    )
  ) %>%
  mutate(
    Physical_trauma  = case_when(
      K_TPKS_level1 == "V" ~ 1, # 
      K_TPKS_level1 == "W" ~ 1,
      K_TPKS_level1 == "X" ~ 1,
      K_TPKS_level1 == "Y" ~ 1,
      !K_TPKS_level1 == "V" ~ 0,!K_VKS_level1 == "W" ~ 0,!K_VKS_level1 == "X" ~ 0,!K_VKS_level1 == "Y" ~ 0
    )
  )
df$All <- df$DEATH
df$AB_K_TPKS <-
  df$Respiratory + df$Cancer + df$Infectious + df$Gastrointestinal > 0
df$AB_K_TPKS <- as.numeric(df$AB_K_TPKS)

df_full <- df

causes <-
  c(
    "Cardiovascular",
    "Respiratory",
    "Cancer",
    "Gastrointestinal",
    "Neurological",
    "Physical_trauma",
 #   "Infectious",
    "All"
  )


library(brms)
set.seed(6322)
table.cause <- NULL

for (cause in causes) {

  print(cause)

  # Only keep complete cases for this cause
  # Same as in all-cause model
  sel <- c(cause, "DEATH_AGEDIFF", "ARGload_log10", get.covariates())  

  dff <- df_full[, sel]
  dff$EVENT <- dff[, cause]; dff[, cause] <- NULL
  dff$TIMEDIFF <- dff$DEATH_AGEDIFF; dff$DEATH_AGEDIFF <- NULL  
  dff <- dff[complete.cases(dff), ] %>%
               filter(TIMEDIFF>0) # Also exclude past events


  # Added scaling / LL
  binaries <- c("MEN", "CURR_SMOKE", "BL_USE_RX_L", "BL_USE_RX_J01", "PREVAL_DIAB", "BP_TREAT") # Binary variables
  continuous <- setdiff(sel, c(binaries, "EVENT", "TIMEDIFF", "DEATH_AGEDIFF", cause))
  dff <- cbind(dff[, c("EVENT", "TIMEDIFF")],
               dff[, binaries],
               scale(dff[, continuous])
	      )

  # Fit the model
  fit <- get.fitted.model(dff)
  # summary(fit) # cox_model

  library(tidybayes)
  # get_variables(fit)
  # Gather the same as in varnames, except the time and event variables
  dfc <- fit %>% gather_draws(b_ARGload_log10) %>%
    summarise_draws() %>%
    arrange(median) %>%
    select(-mean) %>%
    # Convert to the original domain 
    mutate(median=exp(median)) %>%
    # survival_totalmortality_byfamily.R koodi sisältää tietoa miten vaihtaa muihin luottovaleihin tarvittaessa
    mutate(q5=exp(q5)) %>%
    mutate(q95=exp(q95)) 
  dfc$cause <- rep(cause, nrow(dfc))
  dfc$N_events <- sum(dff[, "EVENT"])
  dfc$N_total <- nrow(dff)  
  table.cause <- rbind(table.cause, dfc)
  gc()
}


# Sort
table.cause0 <- table.cause # Store original

library(tidyverse)
library(stringr)
library(dplyr)
table.cause <- table.cause0

table.cause$cause <- apply(table.cause, 1, function (x) {gsub("= ", "=", paste0(x[["cause"]], " (N=", x[["N_events"]], ")"))})
table.cause <- table.cause %>%
               mutate(cause = str_replace(cause, "_", " ")) %>%
               arrange(median) %>%		 
	       as.data.frame() %>%
	       mutate(cause=factor(cause, levels=unique(cause))) %>%
	       select(cause, median, q5, q95, N_events, N_total)

theme_set(theme_bw(20))
v <- c(0.7, 1, 1.4, 2, 3)
fp.causes <- ggplot(table.cause, aes(y = cause, x = median, xmin = q5, xmax = q95)) + 
   geom_pointinterval() +
   labs(x="Hazard ratio", y="", title="Cause specific mortality") +
   geom_vline(xintercept=1, color="darkgray", linetype=2) +
   scale_x_continuous(trans="log2")
   #scale_x_continuous(breaks=v, labels=v, trans="log2", limits=c(0.59, 3.1))      
   #  geom_vline(xintercept=1, color="darkgray", linetype=2) +
   #  scale_x_continuous(breaks=v, labels=v, trans="log2", limits=c(0.8, 1.8))   
#theme_set(theme_bw(20))
#v <- c(0.7, 1, 1.4, 2, 3)
#fp.causes <- ggplot(table.cause, aes(y = Cause, x = `exp(coef)`, xmin = `lower .95`, xmax = `upper .95`)) + 

# fp.causes

li <- list(table=table.cause[rev(seq(nrow(table.cause))), ], figure=fp.causes)
saveRDS(li, file="causes.rds")	

