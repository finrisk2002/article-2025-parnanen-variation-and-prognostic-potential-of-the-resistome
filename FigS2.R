#source("GLMs_TreeSE.R")

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


predicted_div.plot <- ggplot(df,   # Draw plot using ggplot2 package
                             aes(x  = predicted,
                                 y  = ARG_div))  +
  geom_point(size = 1, alpha = 1, color = "grey")  +
  annotate(
    "text",
    label  = paste0("R2=", round(cor(test.target, predictions) ^ 2, digits = 2)),
    x  = 3,
    y  = 0.8,
    size  = unit(8, "pt")
  )  +
  geom_abline(
    intercept  = 0,
    slope  = 1,
    color  = "black",
    linewidth  = 1
  )  +
  ylab("Observed ARG diversity")  +
  xlab("Predicted ARG diversity") +
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
plot2 <- ggplot(df, aes(y  = Covariate, x  = `Est.`))  +
  geom_point(shape  = 18, size  = 5)  +
  geom_errorbarh(aes(xmin  = `2.5%`, xmax  = `97.5%`), height  = 0.25)  +
  geom_vline(
    xintercept  = 0,
    color  = "darkgrey",
    linetype  = "dashed",
    cex  = 1,
    alpha  = 0.5
  )  +
  xlab("Coefficient (95% CI)")  +
  ylab(" ") +
  theme_classic()  +
  theme(
    panel.border  = element_blank(),
    panel.background  = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line  = element_line(colour  = "black"),
    axis.text.y  = element_text(size = 0.9*base.size, colour  = "black"),
    axis.text.x.bottom  = element_text(size  = 1.2*base.size, colour  = "black"),
    axis.title.x  = element_text(size  = 1.2*base.size, colour  = "black"),
    title  = element_text(size  = base.size, colour  = "black")
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

# Plot forest plot
plot1_family <- ggplot(df, aes(y  = Covariate, x  = `exp(Est.)`))  +
  geom_point(shape = 18, size = 5)  +
  geom_errorbarh(aes(xmin  = `2.5%`, xmax  = `97.5%`), height  = 0.25, color="black")  +
  geom_vline(
    xintercept  = 1,
    color  = "darkgray",
    linetype  = "dashed",
    cex  = 2
  )  +
  labs(x="Coefficient (95% CI)", y="") +
  theme_classic()  +
  theme(
    panel.border  = element_blank(),
    panel.background  = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line  = element_line(colour  = "black"),
    axis.text.y  = element_text(size = 0.9*base.size, colour  = "black"),
    axis.text.x.bottom  = element_text(size  = 1.2*base.size, colour  = "black"),
    axis.title.x  = element_text(size  = 1.2*base.size, colour  = "black"),
    title  = element_text(size  = base.size, colour  = "black")
  )

print(plot1_family)

#pdf("../RESULTS/R_results/Fig4A.pdf", width=20, height=10)
#print(plot1)
#dev.off()

# Predicted vs observed plot
#predictions <- predict(glmBoostModel_family, newdata = df_test_num_taxa)
#test.target <- subset(df_test_num, select = SUM_norm)
#predictions <- predict(glmBoostModel_family, newdata = df_test_num)
#test.target <- subset(df_test_num, select = SUM_norm)
predictions <- predict(glmBoostModel_family, newdata = df_test_num_taxa)
test.target <- subset(df_test_num_taxa, select = SUM_norm)

test.target <- test.target %>%
  t %>% subset(select = rownames(predictions)) %>% t

# RMSE
sqrt(mean((test.target - predictions) ^ 2))
# R2
cor(test.target, predictions) ^ 2


# Plot predicted vs observed
base.size <- 20
df <- data.frame(observed = test.target, predicted = predictions)
ranx <- range(10^df$predicted)
rany <- range(df$SUM_norm)
predicted.plot_family <- ggplot(df,   # Draw plot using ggplot2 package
                                 aes(x  = predicted,
                                     y  = log10(SUM_norm)))  +
  geom_point(size = 1, alpha = 0.7, color = "darkgrey") +
  annotate(
    "text",
    # label  = paste0("R2=", round(cor(test.target, predictions)^2, digits = 3)),
    label  = paste0("R2=", round(cor(df$predicted, log10(df$SUM_norm))^2, digits = 2)),
    x  = 2.5,
    y  = 1.2,
    size  = unit(8, "pt")
  ) +
  theme(
    plot.title = element_text(size = 0.5*base.size),
    title = element_text(size = 0.5*base.size),
    axis.text.x = element_text(size = 0.7*base.size),				  
    axis.text.y = element_text(size = 0.7*base.size),
    axis.title.x = element_text(size = 0.6*base.size),
    axis.title.y = element_text(size = 0.6*base.size)
  ) 

predicted.plot_family <- predicted.plot_family +  
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
    axis.text.y  = element_text(size  = base.size, colour  = "black"),
    axis.text.x.bottom  = element_text(size  = base.size, colour  = "black"),
    axis.title.x  = element_text(size  = base.size, colour  = "black")
  ) 



# Diversity plot for supplement

supplement_div_plot <- plot2 + annotation_custom(
  ggplotGrob(predicted_div.plot),
  xmin = 0.07, xmax = 0.42, 
  ymin = 2, ymax = 20) 

supplement_fam_plot <- plot1_family + annotation_custom(
  ggplotGrob(predicted.plot_family),
  xmin = 1.02, xmax = 1.16,
  ymin = 2, ymax= 21) 

figS1 <- cowplot::plot_grid(supplement_div_plot, supplement_fam_plot, labels="auto")
pdf("../RESULTS/R_results/FigS1.pdf", width=24, height=11)
print(figS1)
dev.off()

