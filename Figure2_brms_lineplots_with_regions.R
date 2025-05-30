library(dplyr)
library(tidyr)
library(tidybayes)
library(ggplot2)
library(cowplot)
library(rstan)
library(brms)
library(ggrepel)
library(RColorBrewer)
library(stringr)
library(posterior)
theme_set(theme_tidybayes() + panel_border())
# These options help Stan run faster:
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# -----------------------------------------

#  df %>% group_by(Region) %>% summarise(ave=mean(vaesto, na.rm=TRUE)) %>% arrange(ave)
## A tibble: 6 × 2
#  Region     ave
#  <fct>    <dbl>
#1 Lapland   718.
#2 Oulu      862.
#3 Karelia   927.
#4 Kuopio   1067.
#5 Turku    2504.
#6 Helsinki 4603.

df <- colData(TSE) %>% as.data.frame
df$vaesto_quantile <- cut(df$vaesto, breaks=c(0, 10, 100, 1000, 10000, Inf), labels=c("1", "10", "100", "1000", "10000"), include.lowest=TRUE, right=FALSE)
vars <- c(vegetables="KY100_14", poultry="KY100_21", income="TULOT", popdens="vaesto_quantile", Age="age_class")

# Sort regions by population density
df$Region <- factor(df$Region, levels=as.vector((df %>% group_by(Region) %>% summarise(pop=mean(vaesto, na.rm=TRUE)) %>% arrange(pop))$Region))

library(tidyverse)
library(mia)
library(brms)

models <- list()
# models <- readRDS(file=paste0("models.Rds"))
set.seed(3436436)

# Transform population density
df$vaesto <- log10(df$vaesto)

for (myvar in vars) {

    print(myvar)

    # Pick complete cases
    df2 <- df
    df2$varname <- df2[[myvar]]
    df2 <- df2 %>% filter(!is.na(SUM_norm) & !is.na(varname) & !is.na(ALUE))
    df2 <- (df2 %>% unite("sexregion", MEN, ALUE) %>% mutate(sexregion=str_replace(sexregion, "0_", "Women/")) %>%
                                                      mutate(sexregion=str_replace(sexregion, "1_", "Men/")))
    
    # Run model for each sex and region separately for now (later, a proper hierarchical model)
    # m <- lapply(unique(df2$sexregion), function (g) {brm(SUM_norm ~ factor(varname)-1, family = lognormal(), data = (df2 %>% filter(sexregion==g)))})

    # Control same covariates than in the linear model
    m <- lapply(unique(df2$sexregion), function (g) {brm(SUM_norm ~ factor(varname)-1 + PREVAL_RX_J01_NEVT + PREVAL_RX_J01A_NEVT + PREVAL_RX_J01A + PREVAL_RX_J01F + BL_USE_RX_J01,
                                                         family = lognormal(), data = (df2 %>% filter(sexregion==g)))})
							 
    names(m) <- unique(df2$sexregion)
    models[[myvar]] <- m
    rm("m")

}

#saveRDS(models, file=paste0("models.Rds"))
models <- readRDS(paste0("models.Rds"))

# log(SUM_norm) ~ df_train_num[, x] + PREVAL_RX_J01_NEVT + PREVAL_RX_J01A_NEVT + PREVAL_RX_J01A + PREVAL_RX_J01F + BL_USE_RX_J01,
  
# ---------------------------------------------------------------------------------------

posteriors <- list()
for (myvar in names(models)) {
  m <- models[[myvar]]
  # Extract posterior draws and summaries for all variables
  post <- do.call("rbind", lapply(names(m), function (nam) {exp(t(apply(as_draws_array(m[[nam]], variable = "^b_", regex = TRUE), 3, function (x) keyvars(x, iq=0.1, oq=0.1)))) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="varname") %>%
    mutate(varname=str_replace(varname, "b_factorvarname", "")) %>%
    mutate(group=rep(nam, n()))
    })) %>%
    separate("group", c("Sex", "Region")) %>%
    mutate(Sex=factor(Sex)) %>%
    mutate(Region=factor(Region, levels=c("Helsinki","Turku","Oulu","Kuopio","Karelia","Lapland")))  
  posteriors[[myvar]] <- post
}

posteriors[["KY100_14"]]$varname <- as.numeric(posteriors[["KY100_14"]]$varname)
posteriors[["KY100_21"]]$varname <- as.numeric(posteriors[["KY100_21"]]$varname)
posteriors[["TULOT"]]$varname <- factor(posteriors[["TULOT"]]$varname)
posteriors[["vaesto_quantile"]]$varname <- as.numeric(posteriors[["vaesto_quantile"]]$varname)
posteriors[["age_class"]]$varname <- factor(posteriors[["age_class"]]$varname)

# -----------------------------------------

labsize <- 20
age.labs <- age.labels()
source("R_functions.R")
leg <- get_legend(brplot2(posteriors[[1]]))
li <- lapply(names(posteriors), function (i) {brplot2(posteriors[[i]]) +
             labs(title=names(vars)[vars==i]) +
             theme(legend.position="none")})

lineplots <- plot_grid(
          li[[1]] + labs(title="Fresh vegetables"),
          li[[2]] + labs(title="Poultry consumption", y=""),
	  li[[3]] + labs(title="Household income", y=""),
	  li[[4]] + labs(title="Population density", x="Level (n/km2)", y="") +
	            scale_x_log10(labels=scales::trans_format('log10',scales::math_format(10^.x))),
	  li[[5]] + labs(title="Age", x="Age decade", y="") +
	            scale_x_discrete(labels=age.labs),
	  nrow=5,
	  labels=c("", "", "", ""),
	  label_size=labsize)

# lineplots

# Combine GLM forestplot and prediction plot
lineplots <- lineplots + annotation_custom(
                   ggplotGrob(ggpubr::as_ggplot(leg)),
                       xmin = 0.9, xmax=0.95,
                       ymin = 0.11, ymax = 0.16)

pdf("../RESULTS/R_results/Figure5_regional.pdf", width=20, height=30)
plot(lineplots)
dev.off()


library(writexl)
xs <- list()

# Write the posterior results to a file
for (case in names(posteriors)) {

  x0 <- posteriors[[case]][, c("varname", "Sex", "Region", "mean", "outer.quantile.low", "outer.quantile.high")]

  x <- x0 %>% 
           mutate(mean=round(mean, 0)) %>%
           mutate(low=round(outer.quantile.low, 0)) %>%
           mutate(high=round(outer.quantile.high, 0)) %>%
           mutate("ARG load"=paste0(mean, " (", low, "-", high , ")")) %>%
	   dplyr::rename(Level=varname) %>%
	   mutate(Level=as.character(Level)) %>%
	   dplyr::select(Region, Sex, Level, "ARG load") %>%
	   reshape::cast(Region + Sex ~ Level)

  # Add a new sheet for each iteration
  xs[[paste0(dic.varnames()[case], "_regional")]] <- x

}

xs.regional <- xs

tmp <- write_xlsx(xs, path="../RESULTS/R_results/Supplementary_TableFig5_regional.xlsx")
