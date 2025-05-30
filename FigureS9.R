df <- TSE %>% colData() %>% as.data.frame()
df$SUM_norm 
df$ARG_obs
df$V2 <- as.numeric(df$V2)
# Load required package
library(ggplot2)
library(cowplot)


fit <- lm(ARG_div~V2, data=df) %>% summary()
plot_df <- fit$coefficients %>% as.data.frame()

# Plot cumulative Shannon diversity vs. cumulative RPKM
b <- ggplot(df, aes(x = (V2), y = ARG_div)) +
  geom_line(color = "grey") +
  geom_smooth(method = "loess", color = "black", se = FALSE) +
  labs(
    x = "Sequencing depth (reads)",
    y = "Shannon Diversity",
    title = "Saturation of ARG Diversity"
  ) +
  annotate(
    "text",
    label  = paste0("Estimate=", round(plot_df$Estimate[2], digits = 8), ", p=", round(plot_df$`Pr(>|t|)`[2], digits=3)),    
    y  = 4,
    x  = 3000000,
    size  = unit(6, "pt")
  ) +
  
  xlim(0, 5000000) +
  theme_minimal(20)

fit <- lm(SUM_norm~V2, data=df) %>% summary()
plot_df <- fit$coefficients %>% as.data.frame()

a <- ggplot(df, aes(x = (V2), y = SUM_norm)) +
  geom_line(color = "grey") +
  geom_smooth(method = "loess", color = "black", se = FALSE) +
  labs(
    x = "Sequencing depth (reads)",
    y = "ARG load (RPKM)",
    title = "Saturation of ARG load"
  ) +
  annotate(
    "text",
    label  = paste0("Estimate=", round(plot_df$Estimate[2], digits = 8), ", p=", round(plot_df$`Pr(>|t|)`[2], digits=3)),    
    y  = 1000,
    x  = 3000000,
    size  = unit(6, "pt")
  ) +
  xlim(0, 5000000) +
  theme_minimal(20)

fit <- lm(ARG_obs~V2, data=df) %>% summary()
plot_df <- fit$coefficients %>% as.data.frame()

c <- ggplot(df, aes(x = (V2), y = ARG_obs)) +
  geom_line(color = "grey") +
  geom_smooth(method = "loess", color = "black", se = FALSE) +
  labs(
    x = "Sequencing depth (reads)",
    y = "Observed ARGs",
    title = "Saturation of observed ARGs"
  ) +
  annotate(
    "text",
    label  = paste0("Estimate=", round(plot_df$Estimate[2], digits = 8), ", p=", round(plot_df$`Pr(>|t|)`[2], digits=3)),    
    y  = 175,
    x  = 3000000,
    size  = unit(6, "pt")
  ) +
  xlim(0, 5000000) +
  theme_minimal(20)

df <- df %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Sample") %>%
  dplyr::select(Sample, ReadCount = V2) %>%
  mutate(ReadCount = as.numeric(ReadCount)) %>%
  arrange(ReadCount) %>%
  mutate(Sample = factor(Sample, levels = Sample))  # order factor levels by ReadCount

d <- ggplot(df, aes(x = Sample, y = (ReadCount))) +
  geom_point(stat = "identity",color = "black", size = 0.000001) +
  theme_classic(20) +
  labs(x = "Sample", y = "log10(Read count)",  title = "Sequencing depth") +
  ylim(0, 5000000) +
  theme(axis.text.x = element_blank())


  

plot_grid(a,b,c,d, labels = "auto")

s <- 200
library(Cairo)
CairoJPEG("../RESULTS/R_results/FigureS9.jpg", width=4*s, height=4*s)
plot_grid(a,b,c,d, labels = "auto", nrow=2)
dev.off()


