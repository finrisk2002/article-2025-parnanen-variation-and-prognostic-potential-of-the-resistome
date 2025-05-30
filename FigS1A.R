library(geofi)
library(sf)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(svglite)
## Cannot be run, contains sensitive location data ##
theme_set(theme_tufte(base_family = "sans", base_size = 18) + 
theme(panel.border = element_rect(colour = "black", fill = NA), panel.background = element_rect(fill = "white"), plot.background = element_rect(fill="white"), axis.text = element_text(colour = "black", size = 18)))

spatial_file_path = "spatial_file.txt"
spatial_data <- read.csv(spatial_file_path, sep="\t", header = T)
spatial_data <- spatial_data[complete.cases(spatial_data),]
pop_grid_2005 <- readRDS("pop_grid_2005.rds")

lonlat_to_pop_dens <- function(pointsDF,
                            grid = pop_grid_2005,
                            name_cols = c("vaesto", "miehet", "naiset", "ika_0_14", "ika_15_64", "ika_65_")) {
    ## Convert points data.frame to an sf POINTS object
    pts <- st_as_sf(pointsDF, coords = c("EASTING", "NORTHING"), crs = 3067)
    ## Find population densities of grids intersected by each point
    pop_dens <- data.frame(grid)[name_cols]
    ii <- as.integer(st_intersects(pts, grid))
    pop_dens <- pop_dens[ii,]
	#approximate all NAs to pop density = 1 (since the data point exists in the data)
	pop_dens[is.na(pop_dens$vaesto),"vaesto"] <- 1
	pop_dens[is.na(pop_dens)] <- -1
	pop_dens
}

pop_dens_df <- lonlat_to_pop_dens(spatial_data)
pop_dens_df$FID <- spatial_data$FID
write.csv(pop_dens_df, "FR02_pop_dens_df.csv", row.names=FALSE)
spatial_data$vaesto <- log10(pop_dens_df$vaesto)

jitter_spatial <- spatial_data[which(pop_dens_df$vaesto > 10),]
set.seed(42)
jitter_spatial$EASTING <- jitter(jitter_spatial$EASTING, amount = 10000)
jitter_spatial$NORTHING <- jitter(jitter_spatial$NORTHING, amount = 10000)
jitter_spatial <- st_as_sf(jitter_spatial, coords = c("EASTING", "NORTHING"), crs = 3067)

#Finland's national borders downloaded from https://geodata.lib.utexas.edu/catalog/stanford-mb800bj2863
fin_borders <- st_read("mb800bj2863.shp")
fin_borders <- st_transform(fin_borders, crs = 3067)

p1 <- ggplot() + geom_sf(data = fin_borders, alpha = 0) + 
geom_sf(data = jitter_spatial, aes(colour=vaesto)) + 
scale_colour_viridis_c(name = bquote("Population\ndensity per"~km^2), limits = c(0,5), breaks = c(0,1,2,3,4), labels = c(0, 10, 100, 1000, 10000), alpha = 0.7) +
theme(panel.border = element_blank(), axis.text=element_text(size=12))
ggsave("vaesto_ggplot.pdf", p1)