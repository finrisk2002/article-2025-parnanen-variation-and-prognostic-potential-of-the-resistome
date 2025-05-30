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
library(mia)
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

#df$vaesto_quantile <- cut(df$vaesto, quantile(df$vaesto, probs=seq(0, 1, length=5), na.rm=TRUE))
#df$vaesto_quantile <- cut(df$vaesto, breaks=c(0, 10^seq(2, 4), Inf))
# df$vaesto_quantile <- cut(df$vaesto, breaks=c(0, 200, 2000, Inf), labels=c("<200", "<2000", "<20000"), include.lowest=TRUE, right=FALSE)
# df$vaesto_quantile <- cut(df$vaesto, breaks=c(0, 10, 100, 1000, 10000, Inf), labels=c("<10", "<100", "<1000", "<10000", "<20000"), include.lowest=TRUE, right=FALSE)
#df$vaesto_quantile <- cut(df$vaesto, breaks=c(0, 10, 100, 1000, 10000, Inf), labels=c("<10", "<100", "<1e3", "<1e4", "<2e4"), include.lowest=TRUE, right=FALSE)
df$vaesto_quantile <- cut(df$vaesto, breaks=c(0, 10, 100, 1000, 10000, Inf), labels=c("1", "10", "100", "1000", "10000"), include.lowest=TRUE, right=FALSE)

#vars <- c(vegetables="KY100_14", poultry="KY100_21", income="TULOT", populationdensity="vaesto")

vars <- c(vegetables="KY100_14", poultry="KY100_21", income="TULOT", popdens="vaesto_quantile", Age="age_class")

# Sort regions by population density
df$Region <- factor(df$Region, levels=as.vector((df %>% group_by(Region) %>%
               summarise(pop=mean(vaesto, na.rm=TRUE)) %>%
	       arrange(pop))$Region))

models <- list()
for (myvar in vars) {

  print(myvar)

  # Pick complete cases
  df2 <- df
  df2$varname <- df2[[myvar]]
  df2 <- df2 %>% filter(!is.na(SUM_norm) & !is.na(varname))

  # Run model for each sex separately for now
  m <- lapply(unique(df2$MEN), function (g) {brm(SUM_norm ~ factor(varname)-1, family = lognormal(), data = (df2 %>% filter(MEN==g)))})
  names(m) <- unique(df2$MEN)
  models[[myvar]] <- m

}

posteriors <- list()
for (myvar in vars) {
  m <- models[[myvar]]
  # Extract posterior draws and summaries for all variables
  post <- do.call("rbind", lapply(names(m), function (nam) {exp(t(apply(as_draws_array(m[[nam]], variable = "^b_", regex = TRUE), 3, function (x) keyvars(x, iq=0.1, oq=0.1)))) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="varname") %>%
    mutate(varname=str_replace(varname, "b_factorvarname", "")) %>%
    mutate(sex=rep(nam, n()))}))
  posteriors[[myvar]] <- post
}
# posteriors[["Region"]]$varname <- factor(posteriors[["Region"]]$varname, levels=levels(df$Region))
posteriors[["vaesto_quantile"]]$varname <- as.numeric(posteriors[["vaesto_quantile"]]$varname)

# -----------------------------------------

leg <- get_legend(brplot(posteriors[[1]]))
li <- lapply(names(posteriors), function (i) {brplot(posteriors[[i]]) +
             labs(title=names(vars)[vars==i]) +
             theme(legend.position="none")})

lineplots <- plot_grid(
          li[[1]] + labs(title="Fresh vegetable\nconsumption"),
          li[[2]] + labs(title="Poultry\nconsumption", y=""),
	  li[[3]] + labs(title="Household\nincome", y=""),
	  li[[4]] + labs(title="Population\ndensity", x="Level (n/km2)", y="") +
	            scale_x_log10(labels=scales::trans_format('log10',scales::math_format(10^.x))),
	  li[[5]] + labs(title="Age", x="Age decade", y="") + scale_x_discrete(labels=age.labs),
	  nrow=1,
	  labels=c("c", "", "", ""),
	  label_size=labsize)


# Combine GLM forestplot and prediction plot
lineplots <- lineplots + annotation_custom(
                   ggplotGrob(ggpubr::as_ggplot(leg)),
                       xmin = 0.08, xmax=0.13,
                       ymin = 0.6, ymax = 0.8)


library(writexl)
xs <- list()

# Write the posterior results to a file
for (case in names(posteriors)) {
  x0 <- posteriors[[case]][, c("varname", "sex", "mean", "outer.quantile.low", "outer.quantile.high")]

  # Sorting of the levels
  ord <- as.character((subset(x0, sex==0) %>% arrange(mean))$varname)

  x <- x0 %>% mutate(sex=case_when(sex==0 ~ "Female", sex==1 ~ "Male")) %>%
           mutate(mean=round(mean, 0)) %>%
           mutate(low=round(outer.quantile.low, 0)) %>%
           mutate(high=round(outer.quantile.high, 0)) %>%
           mutate("ARG load"=paste0(mean, " (", low, "-", high , ")")) %>%
	   dplyr::rename(Level=varname) %>%
	   mutate(Level=as.character(Level)) %>%
	   dplyr::rename(Sex=sex) %>%
	   dplyr::select(Level, Sex, "ARG load") %>%
	   reshape::cast(Level ~ Sex)

  # Sort
  x <- x[match(ord, x$Level),]
  x <- x %>% arrange(Level)

  # Add a new sheet for each iteration
  xs[[dic.varnames()[case]]] <- x

}


