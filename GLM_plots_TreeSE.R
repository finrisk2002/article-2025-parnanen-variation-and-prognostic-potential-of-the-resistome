# Plots
set.seed(25421)

library(tidyverse)
source("GLMs_TreeSE.R")
source("R_functions.R")
base.size <- 22
theme_set(theme_bw())

# Apply substitutions to rownames of plot_glm_table
plot_glm_table <- fit_glm %>%
  jtools::summ(confint = TRUE, exp = TRUE) %>%
  .$coeftable %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  replace_variable_with_dictionary("rowname", substitutions) %>%
  column_to_rownames(var = "rowname")

# Add column named Covariate
plot_glm_table$Covariate <- rownames(plot_glm_table)

# remove Intercept row
df <- plot_glm_table[2:length(plot_glm_table$Covariate), ]

# Sort
df <- df %>% arrange(`exp(Est.)`) %>%
             mutate(Covariate=factor(Covariate, levels=unique(Covariate)))

write_xlsx(df, "../RESULTS/R_results/Supplementary_TableFig3b.xlsx")

# Plot forest plot
forestplot.without.bacfam <- ggplot(df, aes(y  = Covariate, x  = `exp(Est.)`))  +
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

print(forestplot.without.bacfam)

#######################################################################################

# Predicted vs observed plot
predictions <- predict(glmBoostModel, newdata = df_test_num)
test.target <- subset(df_test_num, select = SUM_norm)
test.target <- test.target %>%
  t %>% subset(select = rownames(predictions)) %>% t

# RMSE
sqrt(mean((test.target - predictions) ^ 2))
# R2
cor(test.target, predictions) ^ 2

# Plot predicted vs observed
df <- data.frame(observed  = log10(as.vector(test.target)),
                 predicted = unname(unlist(predictions)))
		 
ranx <- range(df$predicted)
rany <- range(df$observed)
prediction.scatterplot.argload <- ggplot(df,
                         aes(x = predicted,
                             y = observed))  +
  geom_point(size = 1, alpha = 0.7, color = "darkgrey") +
  annotate(
    "text",
    label  = paste0("R2=", round(cor(df$predicted, df$observed)^2, digits = 2)),    
    x  = 2.3,
    y  = 1.2,
    size  = unit(8, "pt")
  )  +

  theme_bw()  +
  theme(
    panel.border  = element_blank(),
    panel.background  = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line  = element_line(colour  = "black"),
    axis.text.y  = element_text(size  = base.size, colour  = "black"),
    axis.text.x.bottom  = element_text(size  = base.size, colour  = "black"),
    axis.title.x  = element_text(size  = base.size, colour  = "black")
  ) 

######################################################################################

# Diversity predicted vs observed plot
predictions <- predict(glmBoostModel_div, newdata = df_test_num)
test.target <- subset(df_test_num, select = ARG_div)
test.target <- test.target %>%
  t %>% subset(select = rownames(predictions)) %>% t

# RMSE
sqrt(mean((test.target - predictions) ^ 2))
# R2
cor(test.target, predictions) ^ 2

# Plot predicted vs observed
df <- data.frame(observed = test.target, predicted = predictions)

prediction.scatterplot.diversity <- ggplot(df,
    aes(x  = predicted, y = ARG_div))  +
  geom_point(size = 1, alpha = 1, color = "grey")  +
  annotate(
    "text",
    label  = paste("R2 = ", round(cor(test.target, predictions) ^ 2, digits = 2)),
    x  = 2.4,
    y  = 1.2,    
    size  = unit(8, "pt")
  )  +
  geom_abline(
    intercept  = 0,
    slope  = 1,
    color  = "black",
    linewidth  = 1
  )  +
  labs(x="Predicted (Shannon)",
       y="Observed (Shannon)")  +
  theme_bw()  +
  theme(
    panel.border  = element_blank(),
    panel.background  = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line  = element_line(colour  = "black"),
    axis.text.y  = element_text(size  = base.size, colour  = "black"),
    axis.text.x.bottom  = element_text(size = base.size, colour  = "black"),
    axis.title.x  = element_text(size = base.size, colour  = "black"),
    axis.title.y  = element_text(size = base.size, colour  = "black")    
  )

######################################################################################

# Apply substitutions to rownames of plot_glm_table
plot_glm_table <- fit_glm_div %>%
  jtools::summ(confint = TRUE) %>%
  .$coeftable %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  replace_variable_with_dictionary("rowname", substitutions) %>%
  column_to_rownames(var = "rowname")

# Add column named Covariate
plot_glm_table$Covariate <- rownames(plot_glm_table)

# remove Intercept row
df <- plot_glm_table[2:length(plot_glm_table$Covariate), ]

# Sort
df <- df %>% arrange(`Est.`) %>% mutate(Covariate=factor(Covariate, levels=unique(Covariate)))

# Plot forest plot
forestplot.diversity <- ggplot(df, aes(y  = Covariate, x  = `Est.`))  +
  geom_point(shape  = 18, size  = 5)  +
  geom_errorbarh(aes(xmin  = `2.5%`, xmax  = `97.5%`), height  = 0.25)  +
  geom_vline(
    xintercept  = 0,
    color  = "darkgray",
    linetype  = "dashed",
    cex  = 2,
    alpha  = 0.5
  )  +
  xlab("Coefficient (95% CI)")  +
  ylab(" ")  +
  theme_bw()  +
  theme(
    panel.border  = element_blank(),
    panel.background  = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line  = element_line(colour  = "black"),
    axis.text.y  = element_text(size  = base.size, colour  = "black"),
    axis.text.x.bottom  = element_text(size  = base.size, colour  = "black"),
    axis.title.x  = element_text(size  = base.size, colour  = "black")
  )



#### With family

# Apply substitutions to rownames of plot_glm_table
plot_glm_table <- fit_glm_family %>%
  jtools::summ(confint = TRUE, exp = TRUE) %>%
  .$coeftable %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  replace_variable_with_dictionary("rowname", substitutions) %>%
  column_to_rownames(var = "rowname")

# Add column named Covariate
plot_glm_table$Covariate <- rownames(plot_glm_table)

# remove Intercept row
df <- plot_glm_table[2:length(plot_glm_table$Covariate), ]

# Sort
df <- df %>% arrange(`exp(Est.)`) %>%
   mutate(Covariate=factor(Covariate, levels=unique(Covariate)))

write_xlsx(df, "../RESULTS/R_results/Supplementary_TableFigS2b.xlsx")

# Plot forest plot
forestplot.bacfam <- ggplot(df, aes(y  = Covariate, x  = `exp(Est.)`))  +
  geom_point(shape = 18, size = 5)  +
  geom_errorbarh(aes(xmin  = `2.5%`, xmax  = `97.5%`), height  = 0.25, color="black")  +
  geom_vline(
    xintercept  = 1,
    color  = "darkgray",
    linetype  = "dashed",
    cex  = 2,
    alpha=0.5
  )  +
  labs(x="Coefficient (95% CI)", y="") +
  theme_bw()  +
  theme(
    panel.border  = element_blank(),
    panel.background  = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line  = element_line(colour  = "black"),
    axis.text.y  = element_text(size  = base.size, colour  = "black"),
    axis.text.x.bottom  = element_text(size  = base.size, colour  = "black"),
    axis.title.x  = element_text(size  = base.size, colour  = "black")
  )
  
print(forestplot.bacfam)

# Predicted vs observed plot
predictions <- predict(glmBoostModel_family, newdata = df_test_num_taxa)
test.target <- subset(df_test_num_taxa, select = SUM_norm)

test.target <- test.target %>%
  t %>% subset(select = rownames(predictions)) %>% t

# RMSE
sqrt(mean((test.target - predictions) ^ 2))
# R2
cor(test.target, predictions) ^ 2

# Plot predicted vs observed
df <- data.frame(observed = test.target, predicted = predictions)
ranx <- range(10^df$predicted)
rany <- range(df$SUM_norm)
prediction.scatterplot.family <- ggplot(df, 
  aes(x  = predicted,
      y  = log10(SUM_norm)))  +
  geom_point(size = 1, alpha = 0.7, color = "darkgrey") +
  annotate(
    "text",
    label  = paste0("R2=", round(cor(df$predicted, log10(df$SUM_norm))^2, digits = 2)),
    x  = 2.5,
    y  = 1.2,
    size  = unit(8, "pt")
  ) +				  
  geom_abline(
    intercept  = 0,
    slope  = 1,
    color  = "black",
    linewidth  = 1
  )  +
  labs(x="Predicted (log10 RPKM)",
       y="Observed (log10 RPKM)")  +
       
  theme_bw()  +
  theme(
    panel.border  = element_blank(),
    panel.background  = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line  = element_line(colour  = "black"),
    axis.text.y  = element_text(size = base.size, colour  = "black"),
    axis.text.x.bottom  = element_text(size = base.size, colour  = "black"),
    axis.title.x  = element_text(size = base.size, colour  = "black"),
    axis.title.y  = element_text(size = base.size, colour  = "black")    
  )

pbacfam <- forestplot.bacfam + annotation_custom(
                   ggplotGrob(prediction.scatterplot.family),
                       xmin = 1.02, xmax = 1.10, 
                       ymin = 2, ymax = 16)

# Without bacterial families:
# forestplot.without.bacfam

# Combine GLM forestplot and prediction plot
pargdiv <- forestplot.diversity + annotation_custom(
                   ggplotGrob(prediction.scatterplot.diversity),
                       xmin = 0.1, xmax = 0.3, 
                       ymin = 2.4, ymax = 16)

# Final figure
figS2 <- cowplot::plot_grid(pargdiv, pbacfam, labels="auto", ncol=2, label_size=40)
pdf("../RESULTS/R_results/FigS2.pdf", width=30, height=12)
print(figS2)
dev.off()

### Test excels #####

# Apply substitutions to rownames of table
plot_glm_table <- fit_glm_family_test_geo_demo  %>%
  jtools::summ(confint = TRUE, exp = TRUE) %>%
  .$coeftable %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  replace_variable_with_dictionary("rowname", substitutions) %>%
  column_to_rownames(var = "rowname")

# Add column named Covariate
plot_glm_table$Covariate <- rownames(plot_glm_table)

# remove Intercept row
df <- plot_glm_table[2:length(plot_glm_table$Covariate), ]

# Sort
df <- df %>% arrange(`exp(Est.)`) %>%
  mutate(Covariate=factor(Covariate, levels=unique(Covariate)))


write_xlsx(df, "../RESULTS/R_results/Supplementary_TableFamilyControls_geo.xlsx")

# Apply substitutions to rownames of table
plot_glm_table <- fit_glm_family_test_diet   %>%
  jtools::summ(confint = TRUE, exp = TRUE) %>%
  .$coeftable %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  replace_variable_with_dictionary("rowname", substitutions) %>%
  column_to_rownames(var = "rowname")

# Add column named Covariate
plot_glm_table$Covariate <- rownames(plot_glm_table)

# remove Intercept row
df <- plot_glm_table[2:length(plot_glm_table$Covariate), ]

# Sort
df <- df %>% arrange(`exp(Est.)`) %>%
  mutate(Covariate=factor(Covariate, levels=unique(Covariate)))

write_xlsx(df, "../RESULTS/R_results/Supplementary_TableFamilyControls_diet.xlsx")