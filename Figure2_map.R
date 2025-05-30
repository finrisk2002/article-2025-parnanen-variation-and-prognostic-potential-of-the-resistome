library(mia)
library(dplyr)
library(ggplot2)
library(geofi)
library(sf)

# PK R code for maps
# Figure5_map_2.R
# yearly_map_in_sorvi.R

# ------------------------------------

library(dplyr)
library(mia)
library(reshape2)
library(ggplot2)

dfprint <- colData(TSE) %>%
  as.data.frame %>%    
  group_by(Region) %>%
  summarise(N = n(),
            freq=mean(ARG_burd=="high"),
	    mean.argdiversity=mean(ARG_div, na.rm=TRUE),
	    median.argdiversity=median(ARG_div, na.rm=TRUE), 	    
            HighARGPrevalence=100 * mean(ARG_burd=="high"),	    
	    mean.rpkm=exp(mean(log(SUM_norm), na.rm=TRUE)), # Use geometric mean
	    mean.rpkm.arithmetic=mean(SUM_norm, na.rm=TRUE), # Use arithmetic mean
	    median.rpkm=median(SUM_norm, na.rm=TRUE), # Use median (without log)	    	    
	    qlow.rpkm=exp(quantile(log(SUM_norm), na.rm=TRUE, 0.05)), 
	    qhigh.rpkm=exp(quantile(log(SUM_norm), na.rm=TRUE, 0.95)),

	    qlow.rpkm.diversity =exp(quantile(ARG_div, na.rm=TRUE, 0.05)), 
	    qhigh.rpkm.diversity=exp(quantile(ARG_div, na.rm=TRUE, 0.95))

 	    ) %>% 
  arrange(HighARGPrevalence) %>% 
  mutate(Region=factor(Region, levels=unique(Region))) %>%
  # Contrast the highARG prevalence and mean RPKM to Lapland for all other regions  
  mutate(HighARGRatio=HighARGPrevalence/HighARGPrevalence[Region=="Lapland"]) %>%
  mutate(MedianARGRatio=median.rpkm/median.rpkm[Region=="Lapland"]) %>%    
  mutate_if(is.numeric, function (x) round(x, 2)) %>%
  mutate(MedianARG=paste0(round(median.rpkm), " (", round(qlow.rpkm), "-", round(qhigh.rpkm), ")", sep = "")) %>%

  mutate(MeanARG=paste0(round(mean.rpkm, 1), " (", round(qlow.rpkm, 1), "-", round(qhigh.rpkm, 1), ")", sep = "")) %>%

  mutate(MeanARGdiversity=paste0(round(mean.argdiversity, 1), " (", round(qlow.rpkm.diversity, 1), "-", round(qhigh.rpkm.diversity, 1), ")", sep = "")) %>%

  mutate(MedianARGdiversity=paste0(round(median.argdiversity, 2), " (", round(qlow.rpkm.diversity, 2), "-", round(qhigh.rpkm.diversity, 2), ")", sep = "")) %>%
  
  mutate(Region=dic.region[as.character(Region)]) %>%
  mutate(HighARGPrevalence=round(HighARGPrevalence, 1)) %>%
  mutate(HighARGRatio=round(HighARGRatio, 2)) %>%
  mutate(MedianARGRatio=round(MedianARGRatio, 2)) %>%    
  select(Region, N, MedianARGdiversity, MedianARG, HighARGPrevalence, MedianARGRatio, HighARGRatio) 

# Region specific ARG loads
writexl::write_xlsx(dfprint, path = "../RESULTS/R_results/Table_S2.xlsx")

# -------------------------------------------------------------------------------------------

d1 <- dfprint %>% dplyr::select(Region, HighARGRatio) %>%
                     dplyr::rename(value="HighARGRatio") %>%
                     arrange(value)
d1$variable <- rep("High ARG", nrow(d1))

d2 <- dfprint %>% dplyr::select(Region, MedianARGRatio) %>%
              dplyr::rename(value="MedianARGRatio");
d2$variable <- rep("Median ARG", nrow(d2))

d <- rbind(d1, d2)
d <- d %>% mutate(Region=factor(Region, levels=d1$Region)) %>%
           mutate(variable=factor(variable, levels=rev(unique(variable)))) %>%
	   mutate(Region=dic.region[as.character(d$Region)])

theme_set(theme_bw(20))
p.ratios <- ggplot(d, aes(x=Region, y=value-1, fill=variable)) +
       geom_bar(position="dodge", color="black", stat="identity") +
       scale_fill_manual(values=c("lightgray", "darkgray"), guide=guide_legend(reverse=TRUE)) +
       # Must set values manually due to some geom_bar limitations on automation
       scale_y_continuous(trans="identity", breaks=seq(0,1,0.2), labels=seq(1,2,0.2)) +              
       coord_flip() +
       labs(fill="", x="", y="Ratio") +
       theme(legend.position.inside = c(0.8, 0.15))

pdf("../RESULTS/R_results/Figure3_ratios.pdf", width=7, height=7)
print(p.ratios)
dev.off()

# -------------------------------------------------------------------------------------------

d <- colData(TSE) %>%
  as.data.frame %>%    
  mutate(ARGload=cut(SUM_norm, breaks=quantile(SUM_norm, seq(0, 1, 0.2)), labels=paste0("Q", seq(5)))) %>%
  group_by(Region, ARGload) %>% 
  summarize(n=n()) %>%
  filter(!is.na(ARGload)) %>%
  mutate(Region=dic.region[as.character(Region)]) %>%  
  mutate(Region=factor(Region, levels=unname(d1$Region))) %>%
  mutate(f=n/sum(n)) 

p.ratios1 <- ggplot(d) +
  geom_line(aes(x=ARGload, y=f, group=Region, color=Region)) +
  geom_point(aes(x=ARGload, y=f, group=Region, color=Region)) +  
  facet_wrap(Region ~ .) +
  labs(x="Quantile", y="Fraction (%)") +
  scale_y_continuous(labels=scales::percent)

p.ratios2 <- ggplot(d) +
  geom_bar(aes(x=Region, y=f, fill=ARGload), stat="identity", position="dodge", color="black") +
  #geom_bar(aes(x=ARGload, y=f, fill=Region), stat="identity", position="dodge", color="black") +  
  #geom_bar(aes(x=Region, y=f, fill=ARGload), stat="identity", position="dodge", color="black") +  
  # geom_bar(aes(x=ARGload, y=f, group=Region), stat="identity", position="dodge") +  
  # facet_wrap(Region ~ .) +
  scale_fill_manual(values=rev(gray.colors(5)), guide=guide_legend(reverse=TRUE)) + 
  labs(x="", y="Fraction (%)", fill="Quantile") +
  scale_y_continuous(labels=scales::percent)

# ------------------------------------

#### General code ####

dfr <- dfprint
dfr$Region <- as.character(dfr$Region)
dfr$municipality <- c("Rovaniemi", "Kontiolahti", "Kuopio", "Utajärvi", "Aura", "Helsinki")

# mun <- geofi::get_municipalities(year = 2013)
mun <- geofi::get_municipality_pop(year = 2013)

#### Add empty columns
mun$freq.relative <- NA
mun$freq <- NA
# mun$mean.rpkm <- NA
mun$share_of_cases <- NA

# Update "mun"
clean_and_join_regions("Lapland", 1, "HighARGRatio", "HighARGPrevalence")
clean_and_join_regions("North Karelia", 2, "HighARGRatio", "HighARGPrevalence")
clean_and_join_regions("Pohjois-Savo", 3, "HighARGRatio", "HighARGPrevalence")
clean_and_join_regions(c("North Ostrobothnia", "Kainuu"), 4, "HighARGRatio", "HighARGPrevalence")
clean_and_join_municipalities(c("Turku", "Loimaa", "Ypäjä", "Aura", "Pöytyä", "Oripää", "Punkalaidun"), 5, "HighARGRatio", "HighARGPrevalence")
clean_and_join_municipalities(c("Helsinki", "Vantaa"), 6, "HighARGRatio", "HighARGPrevalence")
mun$HighARGPrevalence <- mun$freq
mun$HighARGRatio <- mun$freq.relative

# Kontiolahti is the central location for North Karelia | Rovaniemi or Sodankylä for Lapland
centrals <- geofi::municipality_central_localities

# Add together populations for South Karelia and North Karelia
dfr$Region <- as.character(dfr$Region)
dfr$municipality <- NA
dfr$municipality <- c("Rovaniemi", "Kontiolahti", "Kuopio", "Utajärvi", "Aura", "Helsinki")
dfr$municipality <- as.factor(dfr$municipality)

# Combine point locations to dfr
central_points <- dplyr::left_join(dfr, centrals, by = c("municipality" = "municipality_name_fi"))
breaks2 <- seq(0.06, 0.14, 0.02)

base.size <- 30
# theme_set(theme_bw(base.size))
theme_set(theme_classic(base.size)) 
# TODO pick colors from probabilistic estimates
pmap <- ggplot2::ggplot() +
  geom_sf(data = mun, aes(fill = HighARGPrevalence/100), color="lightgray", alpha = 1) +
  scale_fill_gradient2(low = "darkblue", mid="gray", high = "darkred",
                       midpoint=0.1, na.value = "white", breaks = breaks2,
		       limits=range(breaks2), label=scales::percent) +
  # Other points than Helsinki might be redundant
  geom_sf(data = central_points, aes(size = 1,
                                     color = HighARGPrevalence/100,
                                     geometry = geometry),
          show.legend = FALSE,
          alpha = 0.6) +
  scale_colour_gradient2(low = "darkblue",
                         mid = "gray",
                         midpoint = 0.1,
                         high = "darkred") +
  guides(fill = guide_legend(title = "")) +
  # See this to fix the map font issue https://r-graphics.org/recipe-output-fonts-pdf
  theme(
        plot.title = element_text(size = base.size),
                                  axis.text.x = element_text(size = base.size, family="serif"),
                                  axis.text.y = element_text(size = base.size, family="serif"),
                                  axis.title.x = element_text(size = base.size),
                                  axis.title.y = element_text(size = base.size),
	                          legend.text = element_text(size=0.9*base.size)
				  ) +
  labs(x="", y="") +
  geom_segment(aes(x=140000, xend=620000, y=7050000, yend=6650000), linetype=2, color="darkgray", size=2) +
  guides(fill = guide_legend(reverse=TRUE, title="High ARG prevalence"))

# -----------------------------------------------------------------------------

##### Points map ####
#
#points_lapland <- st_sample(x = mun[which(mun$maakunta_name_en == "Lapland"),],
#                            size = mun$share_of_cases[which(mun$maakunta_name_en == "Lapland")])
#points_karelia <- st_sample(x = mun[which(mun$maakunta_name_en == "North Karelia"),],
#                            size = mun$share_of_cases[which(mun$maakunta_name_en == "North Karelia")])
#points_savo <- st_sample(x = mun[which(mun$maakunta_name_en == "Pohjois-Savo"),],
#                         size = mun$share_of_cases[which(mun$maakunta_name_en == "Pohjois-Savo")])
#points_oulu <- st_sample(x = mun[which(mun$maakunta_name_en %in% c("North Ostrobothnia", "Kainuu")),],
#                         size = mun$share_of_cases[which(mun$maakunta_name_en %in% c("North Ostrobothnia", "Kainuu"))])
#points_turku <- st_sample(x = mun[which(mun$municipality_name_fi %in% c("Turku", "Loimaa", "Ypäjä", "Aura", "Pöytyä", "Oripää", "Punkalaidun")),],
#                          size = mun$share_of_cases[which(mun$municipality_name_fi %in% c("Turku", "Loimaa", "Ypäjä", "Aura", "Pöytyä", "Oripää", "Punkalaidun"))])
#points_helsinki <- st_sample(x = mun[which(mun$municipality_name_fi %in% c("Helsinki", "Vantaa")),],
#                             size = mun$share_of_cases[which(mun$municipality_name_fi %in% c("Helsinki", "Vantaa"))])
#plot(st_geometry(mun))
#plot(points_lapland,  col = "green",   add = TRUE, pch = 16, cex = 0.5)
#plot(points_karelia,  col = "yellow",  add = TRUE, pch = 16, cex = 0.5)
#plot(points_savo,     col = "orange",  add = TRUE, pch = 16, cex = 0.5)
#plot(points_oulu,     col = "blue",    add = TRUE, pch = 16, cex = 0.5)
#plot(points_turku,    col = "magenta", add = TRUE, pch = 16, cex = 0.5)
#plot(points_helsinki, col = "grey",    add = TRUE, pch = 16, cex = 0.5)

# (Geometric) mean burden
#> exp(mean(log(df$SUM_norm)))
#[1] 267.5819
		    
# print(pmap)

# Frequency per region plot; now replaced by lineplots in Figure5_brms_lineplots.R
# source("Figure5_brms.R")
#pdf("../RESULTS/R_results/Fig5map.pdf", width=15, height=10)
#cowplot::plot_grid(pfreq, pmap, nrow=1, labels="auto", label_size=22)
#dev.off()
# 
# cowplot::plot_grid(fig5, lm_plots, nrow=2, labels=c("", "e"), rel_heights = c(0.5, 0.5), label_size=labsize)

#### Bubble map ####
# Calculate populations for different regions and municipalities
#dfr$pop <- NA
#for (i in seq_along(dfr)) {
#  if (is.na(dfr$Region[i])) {
#    dfr$pop[i] <- pop$vaesto[which(pop$municipality_name_en == dfr$municipality#[i])]
#  } else {
#    # Region population is a sum of populations of municipalities in the region
#    dfr$pop[i] <- sum(pop$vaesto[which(pop$maakunta_name_en == dfr$Region[i])])
#  }
#}

#### Choropleth map ####

# Modify and compare to Salosensaari 2021 map which has broader regions
# Examples geofi for 3 different things:
# 1) Salosensaari scatterplot thingy
# 2) Regions with different solid colors thing
# 3) Regions with heatmap gradient colors
# -> Including municipalities and larger regions like in FINRISK
#-> See Figure5_map

# Eastern Finland: North Karelia and Northern Savo 
# Southwestern Finland: Turku and Loimaa regions in 
# Capital region: cities of Helsinki and Vantaa 
# Northwestern Finland: province of Northern Ostrobothnia
# Northwestern Finland: province of Kainuu
# Province of Lapland
#
# mutate(across('ALUE', str_replace, '2', 'Oulu')) %>%
#  mutate(across('ALUE', str_replace, '3', 'Kuopio')) %>%    
#  mutate(across('ALUE', str_replace, '4', 'Turku')) %>%
#  mutate(across('ALUE', str_replace, '5', 'Helsinki')) %>%
#  mutate(across('ALUE', str_replace, '6', 'Karelia')) %>%  
#  mutate(across('ALUE', str_replace, '7', 'Lapland')) %>%