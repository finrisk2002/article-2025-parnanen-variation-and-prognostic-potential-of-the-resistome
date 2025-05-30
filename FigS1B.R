# Descriptive figures of dataset
library(dplyr)
library(mia)
library(ggplot2)

print("Load functions and data object")
source("R_functions.R")  # load the functions defined in the file "R_functions.R
#TSE <- readRDS("../data/TSE.rds")

# Plot cohort participants by age and sex
df <- colData(TSE) %>% as.data.frame

age.labs <- age.labels()

df <-
  mutate(df,
         age_class = cut(
           BL_AGE,
           breaks = seq(20, 80, by = 10),
	   right=FALSE,
           labels = age.labs,
           include.lowest = TRUE
         ))

# Add in TSE
TSE$age_class <- df$age_class

Age_sex_barplot <-
  ggplot(df, aes(x = age_class, fill = as.factor(MEN)))  +
   labs(x = "Age (y)", y = "Frequency (n)") +
   geom_bar(alpha = 0.5,
           position = "dodge2",
           color = "black") +
  scale_fill_manual(
    values = c("chocolate4", "cornflowerblue"),
    name  = "",
    labels = c("Women", "Men")
  )  +
  theme(
    panel.border  = element_blank(),
    panel.background  = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line  = element_line(colour  = "black"),
    axis.text = element_text(size = 12)
  )

df_compl <- df %>% dplyr::select(ends_with(c("MEN", "TULOT")))  
df_complete <- df_compl[complete.cases(df_compl),]


lm_summary <- lm(vaesto~MEN, data=df) %>% summary()
lm_summary$coefficients[2,1]
lm_summary$coefficients[2,4]



Income_class_barplot <-
  ggplot(df_complete, aes(x = as.factor(TULOT), fill = as.factor(MEN)))  +
  labs(x="Household income level", y = "Frequency (n)")  +
  geom_bar(alpha = 0.5,
           position = "dodge2",
           color = "black") +
  scale_fill_manual(
    values = c("chocolate4", "cornflowerblue"),
    name  = "Sex",
    labels = c("Female", "Male")
  )  +
  theme(
    panel.border  = element_blank(),
    panel.background  = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line  = element_line(colour  = "black"),
    legend.position = "none"
  )

df <-
  mutate(df,
         AB_class = cut(
           PREVAL_RX_J01_NEVT,
           breaks = c(0, 0.1, 1, 5, 10, 20, 85),
           labels = c("0", "1", "2-\n5", "6-\n10",
                      "11-\n20", ">20"),
           include.lowest = TRUE
         ))

Antibiotic_use_barplot <-
  ggplot(df, aes(x = as.factor(AB_class), fill = as.factor(MEN)))  +
  labs(x="Abx reimbursements (n)", y = "Frequency (n)") +
  geom_bar(alpha = 0.5,
           position = "dodge2",
           color = "black") +
  scale_fill_manual(
    values = c("chocolate4", "cornflowerblue"),
    name  = "Sex",
    labels = c("Female", "Male")
  )  +
  theme_bw()  +
  theme(
    panel.border  = element_blank(),
    panel.background  = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line  = element_line(colour  = "black"),
    legend.position = "none"
  )

lm_summary <- lm(PREVAL_RX_J01_NEVT~TULOT+MEN, data=df) %>% summary
lm_summary$coefficients[2,1]
lm_summary$coefficients[2,4]



lm_summary <- lm(PREVAL_RX_J01_NEVT~TULOT+MEN, data=df) %>% summary
lm_summary$coefficients[2,1]
lm_summary$coefficients[2,4]
summ(lm(PREVAL_RX_J01_NEVT~TULOT+MEN, data=df), confint = TRUE)

ab_income <-
  ggplot(df, aes(x = TULOT,
                 y = PREVAL_RX_J01_NEVT,
                 fill = as.factor(MEN)))  +
  labs(x="Household income level", y = "Antibiotic\nreimbursements") +
  geom_smooth(method = "lm", color = "black") +
  scale_fill_manual(
    values = c("chocolate4", "cornflowerblue"),
    name  = "Sex",
    labels = c("Female", "Male")
  ) +
  scale_x_continuous(limits = c(1, 9),
                     breaks = seq(1, 9, 1)) +
  ylim(c(0,5))+
  annotate(
    "text",
    label  = paste0("Estimate: ", round(lm_summary$coefficients[2,1], digits=3), ", p=", round(lm_summary$coefficients[2,4], digits=3)),
    x  = 4,
    y  = 5,
    size  = unit(5, "pt")
  ) +
  theme(legend.position = "none")

lm_summary <- lm(vaesto~TULOT+MEN, data=df) %>% summary
lm_summary$coefficients[2,1]
lm_summary$coefficients[2,4]
summ(lm(vaesto~TULOT+MEN, data=df), confint = TRUE)
income_pop_dens <-
  ggplot(df, aes(y = vaesto,
                 x = TULOT,
                 fill = as.factor(MEN)))  +
  labs(x="Household income level", y="Population density") +
  geom_smooth(method = "lm", color = "black") +

  scale_fill_manual(
    values = c("chocolate4", "cornflowerblue"),
    name  = "Sex",
    labels = c("Female", "Male")
  ) +
  annotate(
    "text",
    label  = paste0("Estimate: ", round(lm_summary$coefficients[2,1]), ", p<0.001"),
    x  = 3.5,
    y  = 2200,
    size  = unit(5, "pt")
  )   +
  scale_x_continuous(limits = c(1, 9),
                     breaks = seq(1, 9, 1)) +
  theme(legend.position = "none")



lm_summary <- lm(TULOT~KY100_14+MEN, data=df) %>% summary
summ(lm(TULOT~KY100_14+MEN, data=df), confint = TRUE)
lm_summary$coefficients[2,1]
lm_summary$coefficients[2,4]

veg_income <-
  ggplot(df, aes(y = KY100_14,
                 x = TULOT,
                 fill = as.factor(MEN)))  +
  labs(x="Household income level", y="Fresh veg.\nconsumption") +
  scale_y_continuous(limits = c(1, 5),
                     breaks = seq(1, 5, 1)) +
  scale_x_continuous(limits = c(1, 9),
                     breaks = seq(1, 9, 1)) +
  geom_smooth(method = "lm", color = "black") +
  scale_fill_manual(
    values = c("chocolate4", "cornflowerblue"),
    name  = "Sex",
    labels = c("Female", "Male")
  ) +
  annotate(
    "text",
    label  = paste0("Estimate: ", round(lm_summary$coefficients[2,1], digits=3), ", p<0.001"),
    y  = 4.5,
    x  = 4,
    size  = unit(5, "pt")
  )   +
  theme(legend.position = "none")

lm_summary <- lm(TULOT~KY100_21+MEN, data=df) %>% summary
lm_summary$coefficients[2,1]
lm_summary$coefficients[2,4]
summ(lm(TULOT~KY100_21+MEN, data=df), confint = TRUE)
poultry_income <-
  ggplot(df, aes(y = KY100_21,
                 x = TULOT,
                 fill = as.factor(MEN)))  +
  labs(x="Household income level", y="Poultry\nconsumption") +
  scale_y_continuous(limits = c(1, 5),
                     breaks = seq(1, 5, 1)) +
  scale_x_continuous(limits = c(1, 9),
                     breaks = seq(1, 9, 1)) +
  geom_smooth(method = "lm", color = "black") +
  scale_fill_manual(
    values = c("chocolate4", "cornflowerblue"),
    name  = "Sex",
    labels = c("Female", "Male")
  ) +
  annotate(
    "text",
    label  = paste0("Estimate: ", round(lm_summary$coefficients[2,1], digits=3), ", p<0.001"),
    y  = 4.5,
    x  = 4,
    size  = unit(5, "pt")
  )   +
  theme(legend.position = "none")





general_stats_plot <-
  cowplot::plot_grid(
    labels = c("b","c","d","e","f","g"),
    nrow = 3,
    ncol = 2,
    align = "hv",
    Age_sex_barplot +
      ylim(c(0, 2000)) +
      theme(axis.title = element_text(size = 16), legend.position = c(0.3,0.8)),
    income_pop_dens +
      theme_classic() +
      theme(legend.position = "none", axis.title = element_text(size = 16)),
    Income_class_barplot  +
      theme_classic() +
      theme(
        legend.position = "none",
        axis.title = element_text(size = 16)
      ),
    ab_income +
      theme_classic() +
      theme(legend.position = "none", axis.title = element_text(size = 16)),
    poultry_income +
      theme_classic()+
      theme(legend.position = "none", axis.title = element_text(size = 16)),
    veg_income  +
      theme_classic() +
      theme(legend.position = "none", axis.title = element_text(size = 16))
    
  )

theme_set(theme_bw(25))
general_stats_plot

png(filename = "../RESULTS/R_results/General_stats_TreeSE.png",
    width = 1000,
    height = 1500)
print(general_stats_plot)
dev.off()

pdf("../RESULTS/R_results/Fig2B.pdf", width=10, height=14)
print(general_stats_plot)
dev.off()

#### Check gender vs abx ###

df <- colData(TSE) %>% as.data.frame()
library(jtools)
glm(PREVAL_RX_J01_NEVT~BL_AGE+MEN+TULOT, data=df) %>% summ(confint=TRUE)
glm(KY100_14~BL_AGE+MEN+TULOT, data=df) %>% summ(confint=TRUE, digits=6)
glm(vaesto~MEN+BL_AGE+TULOT, data=df) %>% summ(confint=TRUE, digits=6)
