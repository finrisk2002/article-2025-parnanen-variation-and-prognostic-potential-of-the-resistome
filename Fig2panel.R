labsize <- 30
base.size <- 20

predicted.burden.plot <- prediction.scatterplot.argload +  
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
  ) +
  theme_set(theme_bw(20))

# Combine GLM forestplot and prediction plot
pglm <- forestplot.without.bacfam + annotation_custom(
                   ggplotGrob(predicted.burden.plot),
                       xmin = 1.02, xmax = 1.12, 
                       ymin = 2, ymax = 20)

ptop <- cowplot::plot_grid(pmap, pglm,
                           rel_widths=c(2,3),
			   nrow=1,
			   labels=c("a", "b"),
			   label_size=labsize)

fig3 <- cowplot::plot_grid(
  ptop,
  lineplots,
  nrow=2,
  rel_heights=c(8,4),
  label_size=labsize)

