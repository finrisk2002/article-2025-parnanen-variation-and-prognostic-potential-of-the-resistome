library(tidyverse)
library(brms)

# models <- readRDS(file=paste0("models.Rds"))
set.seed(3436436)

df <- colData(TSE) %>% as.data.frame
df$vaesto_quantile <- cut(df$vaesto, breaks=c(0, 10, 100, 1000, 10000, Inf), labels=c("1", "10", "100", "1000", "10000"), include.lowest=TRUE, right=FALSE)

# Sort regions by population density
df$Region <- factor(df$Region, levels=as.vector((df %>% group_by(Region) %>% summarise(pop=mean(vaesto, na.rm=TRUE)) %>% arrange(pop))$Region))

# Transform population density
df$vaesto <- log10(df$vaesto)

vars <- c(vegetables="KY100_14", poultry="KY100_21", income="TULOT", popdens="vaesto_quantile", Age="age_class")


models <- list()
posteriors <- list()
dfs.all <- NULL
for (myvar in vars) {

  print(myvar)

  # Pick complete cases
  df2 <- df
  df2$varname <- as.numeric(df2[[myvar]])

  # All regions combined
  df2$Region <- factor(rep("All", nrow(df2)))
  df2 <- df2 %>% filter(!is.na(SUM_norm) & !is.na(varname) & !is.na(Region))  
  df2 <- (df2 %>% unite("group", MEN, Region) %>% mutate(group=str_replace(group, "0_", "Women/")) %>%
                                                  mutate(group=str_replace(group, "1_", "Men/")))
  
  res <- get.post(df2)
  # models[[paste0(myvar, "/All")]] <- res$m0
  posteriors[[paste0(myvar, "/All")]] <- res$post

  # Per East/West
  df2 <- df
  df2$varname <- as.numeric(df2[[myvar]])  
  df2$Region <- df2$EAST
  df2 <- df2 %>% filter(!is.na(SUM_norm) & !is.na(varname) & !is.na(Region))  
  df2 <- (df2 %>% unite("group", MEN, Region) %>% mutate(group=str_replace(group, "0_", "Women/")) %>%
                                                  mutate(group=str_replace(group, "1_", "Men/")))
  
  res <- get.post(df2)
  # models[[paste0(myvar, "/Regional")]]    <- res$m0
  posteriors[[paste0(myvar, "/Regional1")]] <- res$post

  # Per Region
  df2 <- df
  df2$varname <- as.numeric(df2[[myvar]])  
  df2$Region <- df2$ALUE
  df2 <- df2 %>% filter(!is.na(SUM_norm) & !is.na(varname) & !is.na(Region))  
  df2 <- (df2 %>% unite("group", MEN, Region) %>% mutate(group=str_replace(group, "0_", "Women/")) %>%
                                                  mutate(group=str_replace(group, "1_", "Men/")))
  
  res <- get.post(df2)
  # models[[paste0(myvar, "/Regional")]] <- res$m0
  posteriors[[paste0(myvar, "/Regional2")]] <- res$post

  dfs <- NULL
  fits <- res$m0
  for (k in 1:length(fits)) {
    nam <- names(fits)[[k]]
    e <- conditional_effects(fits[[nam]], effects="varname")[[1]]
    d <- e[, c("varname", "estimate__", "lower__", "upper__")]
    d$group <- rep(nam, nrow(e))
    d$variable <- rep(myvar, nrow(e))
    dfs <- rbind(dfs, d)
  }
  dfs.all <- rbind(dfs.all, dfs)

}

ddd <- separate(dfs.all, "group", c("Sex", "Region")) %>%
         mutate(Sex=factor(Sex)) %>%
         mutate(Region=factor(Region, levels=c("Helsinki", "Turku", "Oulu", "Kuopio", "Karelia", "Lapland"))) %>%
	 mutate(variable=factor(dic.varnames()[variable], levels=c("Fresh vegetables", "Poultry consumption", "Household income", "Population density", "Age")))

theme_set(theme_bw(20))
# Identify posterior predictions that match the levels
inds <- which((ddd$varname - round(ddd$varname))==0)
#p <- ggplot(ddd[seq(1, 6000, 20),], aes(x=varname, y=exp(estimate__), color=Sex, shape=Region)) +
p <- ggplot(ddd[inds,], aes(x=varname, y=exp(estimate__), color=Sex, shape=Region)) +
       #geom_smooth() +
       geom_line() +
       geom_point(size=2) +       
       labs(title=paste(myvar, nam, sep = "/")) +
       facet_grid(variable ~ Sex) +
       scale_color_manual(
          values = rev(c("chocolate4", "cornflowerblue")),
          name  = ""
        ) +
	labs(x="Level", y="ARG load (RPKM)", title="") +
	guides(colour = guide_legend(reverse=TRUE))

print(p)
fig.lm.regional <- p

# --------------------------------------------

# Plot

# LM for a brief comparison
mx <- lapply(unique(df2$group), function (g) {
  lm(log(SUM_norm) ~ varname + PREVAL_RX_J01_NEVT + PREVAL_RX_J01A_NEVT + PREVAL_RX_J01A + PREVAL_RX_J01F + BL_USE_RX_J01, data = (df2 %>% filter(group==g)))
})






# Correspondence with the probabilistic model:
# cbind(subset(post, Sex=="Women")$mean, exp(mx[[1]]$coef)) # lognormal
# print(cbind(res$post$mean, mx[[1]]$coef)) # gaussian   
plot(cbind(res$post$mean, mx[[1]]$coef)); abline(0,1) # gaussian

# ----------------

#df <- bind_rows(posteriors)
#df$variable <- factor(df$variable)
#nam <- names(posteriors)[[1]]
#post <- posteriors[[nam]] %>% arrange(mean) %>% mutate(varname=factor(varname, levels=unique(varname)))
#brplot3(post %>% filter(!varname=="b_Intercept")) + labs(title=dic.varnames()[nam])

# Include only selected regional splits
pp <-  posteriors# [!grepl("Regional1", names(posteriors))]
dfp <- lapply(pp, function (post) {subset(post, varname=="b_varname")}) %>%
        bind_rows %>%
	mutate(Region=str_replace(Region, "0", "West")) %>%
        mutate(Region=str_replace(Region, "1", "East")) %>%
	mutate(Region=factor(Region, levels=c("Helsinki", "Turku", "Oulu", "Kuopio", "Karelia", "Lapland", "West", "East", "All"))) %>%
	mutate(variable=dic.varnames()[variable])

dfp <- dfp %>%
	mutate(variable=factor(variable, levels=unique(dfp$variable)))

# TODO add the table
# write_excel()

source("R_functions.R")
p <- brplot4(dfp)
# leg <- get_legend(p)
print(p)
fig.lm.regional.efs <- p + guides(color="none")

# plot_grid(fig.lm.regional, fig.lm.regional.efs, rel_widths=c(2,1))

pdf(file="../RESULTS/R_results/FigED5regionalA.pdf", width=12, height=12)
print(p)
dev.off()

saveRDS(posteriors, file="regional.models.Rds")


