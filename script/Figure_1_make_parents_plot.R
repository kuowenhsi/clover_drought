library(tidyverse)
library(sf)
library(sp)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(scatterpie)
# library(maps)
library(ggpubr)
library(raster)

setwd("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/Drought_F3_paper/Env")

parent_loc <- read_csv("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/GBS_wild_population/cyano_buf_climate_data_419.csv")%>%
  dplyr::select(genotype, Latitude, Longitude)%>%
  filter(genotype %in% c("DMN_010", "STL_0701", "GFL_007"))%>%
  st_as_sf(coords = c("Longitude", "Latitude"), agr = "constant",  crs = 4326)


usa_state <- ne_states(country = "United States of America", returnclass = "sf")%>%st_crop(xmin = -135, ymin = -55, xmax = -65, ymax = 60)%>%st_geometry()

canada_state <- ne_states(country = "canada", returnclass = "sf")%>%st_crop(xmin = -135, ymin = -55, xmax = -65, ymax = 60)%>%st_geometry()
mexico_state <- ne_states(country = "mexico", returnclass = "sf")%>%st_crop(xmin = -135, ymin = -55, xmax = -65, ymax = 60)%>%st_geometry()

p <- ggplot(data = parent_loc)+
  geom_sf(data = usa_state, fill = NA, color = "gray75")+
  geom_sf(data = canada_state, fill = NA, color = "gray75")+
  geom_sf(data = mexico_state, fill = NA, color = "gray75")+
  geom_sf(aes(color = genotype, shape = genotype), size = 2)+
  theme_bw()+
  theme(legend.position="bottom", legend.box="vertical", panel.grid = element_blank())+
  coord_sf(ylim = c(2e6,7e6),expand = FALSE, crs = 3857)

p

ggsave("just_parents_location.png", width = 8, height = 6, dpi = 600)

list.files("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/Drought_F3_paper/Env/wc2.1_2.5m_tavg")
tmp_layer <- stack(c(paste0("/Users/kuowenhsi/OneDrive\ -\ Washington\ University\ in\ St.\ Louis/Drought_F3_paper/Env/wc2.1_2.5m_tavg/", list.files("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/Drought_F3_paper/Env/wc2.1_2.5m_tavg"))))%>%raster::crop(extent(-135,-55,15,60))

plot(tmp_layer)


tmp_extracted <- raster::extract(tmp_layer, parent_loc, df = TRUE)%>%
  mutate(ID = parent_loc$genotype)%>%
  pivot_longer(cols = 2:13, names_to = "month", values_to = "avg_temp")%>%
  mutate(month = as.integer(str_remove(month, "wc2.1_2.5m_tavg_")))%>%
  group_by(ID)%>%
  mutate(GDD5 = (avg_temp - 5))%>%
  mutate(GDD5 = ifelse(GDD5 > 0, GDD5*30, 0))%>%
  mutate(GDD5 = cumsum(GDD5))%>%
  ungroup()

########

prep_layer <- stack(c(paste0("/Users/kuowenhsi/OneDrive\ -\ Washington\ University\ in\ St.\ Louis/Drought_F3_paper/Env/wc2.1_2.5m_prec/", list.files("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/Drought_F3_paper/Env/wc2.1_2.5m_prec"))))%>%raster::crop(extent(-135,-55,15,60))

plot(prep_layer)


prep_extracted <- raster::extract(prep_layer, parent_loc, df = TRUE)%>%
  mutate(ID = parent_loc$genotype)%>%
  pivot_longer(cols = 2:13, names_to = "month", values_to = "avg_prep")%>%
  mutate(month = as.integer(str_remove(month, "wc2.1_2.5m_prec_")))
  
comb_extracted <- left_join(tmp_extracted, prep_extracted, by = c("ID", "month"))%>%
  mutate(avg_prep = avg_prep/10, GDD5 = GDD5/100)%>%
  pivot_longer(cols = 3:5, names_to = "env_varibles", values_to = "values")


pGDD <- ggplot(data = filter(comb_extracted, env_varibles == "GDD5"), aes (x = month, y = values*100))+
  geom_hline(yintercept = 1000, linetype = "dashed")+
  geom_point(aes(color = ID), size = 2)+
  geom_line(aes(color = ID), alpha = 1)+
  scale_x_continuous(breaks = 1:12)+
  scale_y_continuous(expression("Cumulative Growing Degree Days (> 5"*degree*C*")"), n.breaks = 10)+
  scale_color_manual("Locations", values = c("blue", "red2", "green4"), labels = c(DMN_010 = "Duluth, MN", GFL_007 = "Gainsville, FL", STL_0701 = "St. Louis, MO"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank())

pGDD

ggsave("GDD5_monthly_20240530.png", width = 8, height = 4, dpi = 600)

p_legend <- ggplot(data = filter(comb_extracted, ID == "DMN_010"), aes(x = month, y = values))+
  geom_point(aes(color = env_varibles), size = 1)+
  geom_line(aes(color = env_varibles), alpha = 0.7)+
  scale_y_continuous("Temperature (°C)",sec.axis = sec_axis(~ . * 10, name = "Precipitation (mm)"), limits = c(-15,58))+
  scale_x_continuous(breaks = c(1:6)*2)+
  scale_color_manual(name = "", values = c("blue", "red2", "green4"),labels = c("avg. prep.", "avg. temp.", "GDD5"))+
  labs(title = "DMN_010")+
  theme_bw()+
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "#F8766D60"), plot.title = element_text(color = "#F8766D", hjust = 0.5, face = "bold",size = 8), axis.title = element_text(size = 8), axis.text = element_text(size = 6),plot.background = element_rect(fill = "#FFFFFF90", color = NA), legend.background = element_rect(fill = NA, color = NA))


p_DMN010 <- ggplot(data = filter(comb_extracted, ID == "DMN_010"), aes(x = month, y = values))+
  geom_point(aes(color = env_varibles), size = 1)+
  geom_line(aes(color = env_varibles), alpha = 0.7)+
  scale_y_continuous("Temperature (°C)",sec.axis = sec_axis(~ . * 10, name = "Precipitation (mm)"), limits = c(-15,58))+
  scale_x_continuous(breaks = c(1:6)*2)+
  scale_color_manual(name = "", values = c("blue", "red2", "green4"),labels = c("avg. prep.", "avg. temp.", "GDD5"))+
  labs(title = "DMN_010")+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = "none", panel.background = element_rect(fill = "#F8766D60"), plot.title = element_text(color = "#F8766D", hjust = 0.5, face = "bold",size = 8), axis.title = element_text(size = 8), axis.text = element_text(size = 6),plot.background = element_rect(fill = "#FFFFFF90", color = NA))

p_DMN010


p_STL0701 <- ggplot(data = filter(comb_extracted, ID == "STL_0701"), aes(x = month, y = values))+
  geom_point(aes(color = env_varibles))+
  geom_line(aes(color = env_varibles), alpha = 0.7)+
  scale_y_continuous("Temperature (°C)",sec.axis = sec_axis(~ . * 10, name = "Precipitation (mm)"), limits = c(-15,58))+
  scale_x_continuous(breaks = c(1:6)*2)+
  scale_color_manual(name = "", values = c("blue", "red2", "green4"),labels = c("avg. prep.", "avg. temp.", "GDD5"))+
  labs(title = "STL_0701")+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = "none", panel.background = element_rect(fill = "#619CFF60"), plot.title = element_text(color = "#619CFF", hjust = 0.5, face = "bold", size = 8), axis.title = element_text(size = 8), axis.text = element_text(size = 6),plot.background = element_rect(fill = "#FFFFFF90", color = NA))

p_STL0701

p_GFL007 <- ggplot(data = filter(comb_extracted, ID == "GFL_007"), aes(x = month, y = values))+
  geom_point(aes(color = env_varibles))+
  geom_line(aes(color = env_varibles), alpha = 0.7)+
  scale_y_continuous("Temperature (°C)",sec.axis = sec_axis(~ . * 10, name = "Precipitation (mm)"), limits = c(-15,58))+
  scale_x_continuous(breaks = c(1:6)*2)+
  scale_color_manual(name = "", values = c("blue", "red2", "green4"),labels = c("avg. prep.", "avg. temp.", "GDD5"))+
  labs(title = "GFL_007")+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = "none", panel.background = element_rect(fill = "#00BA3860"), plot.title = element_text(color = "#00BA38", hjust = 0.5, face = "bold", size = 8), axis.title = element_text(size = 8), axis.text = element_text(size = 6),plot.background = element_rect(fill = "#FFFFFF90", color = NA))

p_GFL007


p + annotation_custom(ggplotGrob(p_DMN010), xmin = -10e6, xmax = -7.7e6,ymin = 5e6,ymax = 7e6)+
  annotation_custom(ggplotGrob(p_STL0701), xmin = -12.8e6, xmax = -10.5e6,ymin = 3.8e6,ymax = 5.8e6)+
  annotation_custom(ggplotGrob(p_GFL007), xmin = -11.5e6, xmax = -9.2e6,ymin = 2e6,ymax = 4e6)+
  annotation_custom(get_legend(p_legend), xmin = -16e6, xmax = -12.5e6,ymin = 6.4e6,ymax = 6.7e6)


#############
env_layer <- raster("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/Undergrad/Grace/Env_varibles/current_30arcsec_growingDegDays5.tif")

#############

env_layer_croped <- crop(env_layer, extent(-135,-65,25,50))/10
# plot(env_layer_croped)

env_layer_df <- as.data.frame(env_layer_croped, xy=TRUE)%>% drop_na()
head(env_layer_df)
max(env_layer_croped)

p_env <- ggplot(data = parent_loc)+
  geom_raster(data = env_layer_df, aes(x=x, y=y,fill=current_30arcsec_growingDegDays5))+
  geom_sf(data = usa_state, fill = NA)+
  geom_sf(data = canada_state, fill = NA)+
  geom_sf(data = mexico_state, fill = NA)+
  geom_sf(aes(color = genotype, shape = genotype), size = 2, show.legend = FALSE)+
  scale_fill_viridis_c(option = "inferno", name = "")+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.key.height = unit(0.1, "in"), legend.key.width = unit(0.4, "in"), legend.margin = margin(t=-0.02, b=-0.02, unit = "in"), legend.position = "top", panel.background = element_rect(fill = "white"), axis.title = element_text(size = 8), legend.background = element_rect(fill = NA, color = NA), plot.background = element_rect(fill = NA, color = NA))+
  xlab(expression("Growing Degree Days 5 "*degree*C))+
  ylab("")+
  coord_sf(xlim = c(-130, -65), ylim = c(25, 50), expand = FALSE, label_axes = "")

# p_env

out_p <- p + annotation_custom(ggplotGrob(p_DMN010), xmin = -10e6, xmax = -7.7e6,ymin = 5e6,ymax = 7e6)+
  annotation_custom(ggplotGrob(p_STL0701), xmin = -12.8e6, xmax = -10.5e6,ymin = 3.8e6,ymax = 5.8e6)+
  annotation_custom(ggplotGrob(p_GFL007), xmin = -11.5e6, xmax = -9.2e6,ymin = 2e6,ymax = 4e6)+
  annotation_custom(get_legend(p_legend), xmin = -16e6, xmax = -12.5e6,ymin = 6.4e6,ymax = 6.7e6)+
  annotation_custom(ggplotGrob(p_env), xmin = -16e6, ymin = 2e6, xmax = -11e6, ymax = 4e6)


ggsave("clim_genotype_month.png", width = 8, height = 6, dpi = 600)
