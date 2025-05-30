# Plot Metaphlan
pal <- pal.family()

# Plot the PCoA reduced dimension, which has been stored in the TSE object

# Calculate explained variance
e <- attr(reducedDim(TSE, "PCoA"), "eig")
rel_eig <- e / sum(e[e > 0])

p <- scater::plotReducedDim(TSE, "PCoA",
                            colour_by = "Signif_fam",
			    point_alpha=0.5, # see ?"scater-plot-args"
			    point_size=3) +   # ?"scater-plot-args"
                 # Add explained variance for each axis
                 labs(x = paste("PCoA 1 (", round(100 * rel_eig[[1]], 1), "%", ")", sep = ""),
                      y = paste("PCoA 2 (", round(100 * rel_eig[[2]], 1), "%", ")", sep = ""))

p2 <- scater::plotReducedDim(TSE, "PCoA",
                            colour_by = "TopES",
			    point_alpha=0.5, # see ?"scater-plot-args"
			    point_size=3)  +    # ?"scater-plot-args"
                  # Add explained variance for each axis
                  labs(x = paste("PCoA 1 (", round(100 * rel_eig[[1]], 1), "%", ")", sep = ""),
                       y = paste("PCoA 2 (", round(100 * rel_eig[[2]], 1), "%", ")", sep = ""))

base.size <- 40
# Possibly to supplements (genus level PCoA)
metaphlan.plot <- p  +
  scale_color_manual(values = pal.family(), name="Dominant Family") +
  theme_classic(base.size) +
  theme(legend.position=c(0.87, 0.17)) +
  coord_equal(clip = "on")  +
  theme(plot.margin  = unit(c(1, 1, 1, 1), "pt"),
        plot.title = element_text(size = base.size),
                                  axis.text.x = element_text(size = base.size),
                                  axis.text.y = element_text(size = base.size),
                                  axis.title.x = element_text(size = base.size),
                                  axis.title.y = element_text(size = base.size),
	                          legend.text = element_text(size=0.9*base.size),
	     			  legend.title = element_text(size=1*base.size)				  
				  ) +
	     guides(color = guide_legend(override.aes = list(size=5))) 				    


