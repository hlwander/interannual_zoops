# Redundancy analysis for all BVR zoop data 2014-2021

pacman::p_load(vegan, tidyr, data.table, lubridate, rLakeAnalyzer, car, ggnewscale)

#read in zoop data
all_zoops_dens <- read.csv("Output/all_zoops_dens.csv",header = TRUE)

#list of all taxa
taxa <- unique(all_zoops_dens$Taxon)

#taxa as cols, dates as rows 
all_zoops <- all_zoops_dens |> 
  select(DateTime, Taxon, dens) |> 
  filter(Taxon %in% taxa) |> 
  mutate(DateTime = as.Date(DateTime)) |> 
  pivot_wider(names_from = Taxon, values_from = dens) |> 
  mutate_all(~replace(., is.na(.), 0)) |> 
  ungroup() |> group_by(DateTime) |>
  summarise(Bosmina = mean(Bosmina),
            Ceriodaphnia = mean(Ceriodaphnia),
            Daphnia = mean(Daphnia),
            Cyclopoida = mean(Cyclopoida),
            Nauplii = mean(Nauplii),
            Ascomorpha = mean(Ascomorpha),
            Conochilus = mean(Conochilus),
            Keratella = mean(Keratella),
            Kellicottia = mean(Kellicottia),
            Polyarthra = mean(Polyarthra)) |> 
  ungroup() |>
  filter(!DateTime %in% c("2021-07-07", "2021-07-08", "2020-12-02", 
                          "2020-08-13", "2019-07-18", "2016-08-03",
                          "2016-05-20", "2014-04-30", "2014-04-17", 
                          "2014-04-04", "2021-06-16", "2019-07-25",
                          "2019-07-10")) #removing the dates that don't have corresponding ysi/ctd casts

#select only data cols
zoops_dens <- all_zoops |> select(Bosmina:Polyarthra)

#hellinger transform data
zoop_dens_trans <- labdsv::hellinger(zoops_dens)

#list of dates to match up with env data (n=68)
dates <- unique(all_zoops$DateTime)

#read in ctd
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/13/27ceda6bc7fdec2e7d79a6e4fe16ffdf" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

ctd <-read.csv(infile1,header=T) |> 
  mutate(DateTime = as.Date(DateTime)) |>
  filter(DateTime %in% dates & 
           Depth_m > 0 & Site ==50 & Reservoir %in% c("BVR"),
         !is.na(Temp_C)) |> 
  select(Reservoir, Site, Depth_m, DateTime,
         Temp_C, DO_mgL)

# read in ysi data
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/198/11/6e5a0344231de7fcebbe6dc2bed0a1c3" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(500, getOption("timeout"))))

ysi <- read.csv(infile1,header=T) |> 
  mutate(DateTime = as.Date(DateTime)) |>
  filter(DateTime %in% dates &
           !DateTime %in% as.Date(ctd$DateTime) &
           Depth_m > 0 &
           Site ==50 & Reservoir %in% c("BVR")) |> 
  select(DateTime, Depth_m, 
         Temp_C, DO_mgL)

#select every 0.5m from casts for thermocline calc below
ctd_final_temp_do <- ctd |>
  dplyr::mutate(rdepth = plyr::round_any(Depth_m, 0.5)) |>
  dplyr::group_by(DateTime, rdepth) |>
  dplyr::summarise(temp = mean(Temp_C, na.rm=T),
                   DO = mean(DO_mgL, na.rm=T)) |>
  dplyr::rename(depth = rdepth) 

#rename columns
ctd_final_temp_do <- ctd_final_temp_do |> 
  rename(Depth_m = depth,
         Temp_C = temp,
         DO_mgL = DO)

#clean up ysi
ysi_clean <- ysi |> 
  select(DateTime, Depth_m,Temp_C, DO_mgL) |> 
  filter(!is.na(Temp_C)) |>
  distinct(DateTime, Depth_m, .keep_all = TRUE) #one repeated cast so getting rid of the second occurrence of 2021-06-15

#combine ctd and ysi dfs 
temp_final <-ctd_final_temp_do |>   
  select(DateTime, Depth_m,Temp_C) |> 
  full_join(ysi_clean, multiple = 'all', 
            by = c('DateTime', 'Depth_m',"Temp_C")) |>
  arrange(DateTime, Depth_m) |> 
  select(-c(DO_mgL)) |>
  mutate(Depth_m = ifelse(Depth_m==0, 0.1, Depth_m))

do_final <- ctd_final_temp_do |> 
  select(DateTime, Depth_m, DO_mgL) |> 
  full_join(ysi_clean, multiple = 'all', 
            by = c('DateTime', 'Depth_m',"DO_mgL")) |>
  arrange(DateTime, Depth_m) |> 
  select(-c(Temp_C)) |>
  mutate(Depth_m = ifelse(Depth_m==0, 0.1, Depth_m))

#calculate thermocline depth using combined ctd and ysi df
ctd_thermo_depth <- temp_final |> 
  group_by(DateTime) |> 
  summarise(therm_depth = thermo.depth(Temp_C,Depth_m)) |>
  mutate(therm_depth = ifelse(is.nan(therm_depth), 0, therm_depth)) 
#replacing the three NAs with 0 becuase they are essentially isothermal

#calculate oxycline depth
ctd_oxy_depth <- do_final |> 
  group_by(DateTime) |> 
  summarise(oxy_depth = thermo.depth(DO_mgL,Depth_m)) |> 
  ungroup() |>
  group_by(DateTime) |> 
  summarise(oxy_depth = mean(oxy_depth, na.rm=T)) |>
  dplyr::filter(DateTime %in% dates)  |>
  mutate(oxy_depth = ifelse(is.nan(oxy_depth), 0, oxy_depth)) 

#round ctd to nearest m
ctd_final <- ctd |>
  mutate(rdepth = plyr::round_any(Depth_m, 1),
         month = month(DateTime),
         year = year(DateTime)) |>
  select(-Depth_m) |> 
  rename(Depth_m = rdepth) |> 
  left_join(ctd_thermo_depth, by = "DateTime") |> 
  group_by(DateTime) |> 
  summarise(
    Temp_C_epi = ifelse(
      first(therm_depth) == 0, 
      mean(Temp_C, na.rm = TRUE),
      mean(Temp_C[Depth_m < plyr::round_any(first(therm_depth), 1)], na.rm = TRUE)
    ),
    Temp_C_hypo = mean(Temp_C[Depth_m >= plyr::round_any(first(therm_depth), 1)], na.rm = TRUE),
    DO_mgL_epi = ifelse(
      first(therm_depth) == 0, 
      mean(DO_mgL, na.rm = TRUE),
      mean(DO_mgL[Depth_m < plyr::round_any(first(therm_depth), 1)], na.rm = TRUE)
    ),
    .groups = "drop"
  )

ysi_final <- ysi |>
  left_join(ctd_thermo_depth, by = "DateTime") |>
  group_by(DateTime) |>
  summarise(
    Temp_C_epi = ifelse(
      first(therm_depth) == 0,
      mean(Temp_C, na.rm = TRUE),
      mean(Temp_C[Depth_m < plyr::round_any(first(therm_depth), 1)], na.rm = TRUE)
    ),
    Temp_C_hypo = mean(
      Temp_C[Depth_m >= plyr::round_any(first(therm_depth), 1)],
      na.rm = TRUE
    ),
    DO_mgL_epi = ifelse(
      first(therm_depth) == 0,
      mean(DO_mgL, na.rm = TRUE),
      mean(DO_mgL[Depth_m < plyr::round_any(first(therm_depth), 1)], na.rm = TRUE)
    ),
    .groups = "drop"
  )

#combine ctd and ysi to account for missing ctd days
profiles <- bind_rows(ctd_final, ysi_final) |> arrange(DateTime) 

#download bathymetry
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/1254/1/f7fa2a06e1229ee75ea39eb586577184" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

bathymetry <- readr::read_csv(infile1, show_col_types = F)  |>
  dplyr::select(Reservoir, Depth_m, SA_m2) |>
  dplyr::filter(Reservoir == "BVR")

#read in water level 
water_level <- read.csv("./Output/BVR_WaterLevel_2014_2022_interp.csv") |> 
  select(Date, WaterLevel_m) |> 
  mutate(WaterLevel_m = as.numeric(WaterLevel_m)) |> 
  filter(Date %in% dates) |> 
  group_by(Date) |> 
  summarise(waterlevel = mean(WaterLevel_m,na.rm=T)) |> 
  arrange(Date) |>
  mutate(wl_cv = sd(waterlevel)/mean(waterlevel)) |> 
  ungroup() 

wl <- read.csv("./Output/BVR_WaterLevel_2014_2022_interp.csv") |> 
  select(Date, WaterLevel_m) |> 
  mutate(WaterLevel_m = as.numeric(WaterLevel_m)) |> 
  mutate(Date = as.Date(Date)) |> 
  filter(Date %in% dates)

#combine temp_final and wl
temp_final_wl <- temp_final |> 
  dplyr::rename(Date = DateTime) |> 
  dplyr::full_join(wl, multiple = 'all', by = 'Date') |> 
  dplyr::filter(!is.na(Temp_C)) |>
  mutate(Reservoir = "BVR")

#Create a dataframe with bathymetry at each date
flexible_bathy <- temp_final_wl |> # takes the depth at each day
  dplyr::ungroup() |>
  dplyr::select(-Depth_m) |> 
  dplyr::distinct(Date, WaterLevel_m, Reservoir) |>
  dplyr::full_join(bathymetry, multiple = 'all', by = 'Reservoir') |>
  dplyr::group_by(Date) |>
  dplyr::mutate(Depth_m = Depth_m - (max(Depth_m) - mean(unique(WaterLevel_m))),
                WaterLevel_m = mean(WaterLevel_m)) |>
  dplyr::filter(Depth_m>=0) |>
  dplyr::distinct()

#initialize SS df
schmidts <- data.frame("Date" = unique(temp_final_wl$Date),
                       "SS" = NA)

#for loop for physical metrics
for(i in 1:length(unique(flexible_bathy$Date))) {
  baths <- flexible_bathy |>
    dplyr::filter(Date==unique(flexible_bathy$Date)[i])
  
  temps <- temp_final_wl |>
    dplyr::filter(Date == unique(flexible_bathy$Date)[i],
                  # cannot have an observation at a depth shallower than the
                  # shallowest bathymetry (returns NA below) so these are filtered out
                  Depth_m >= min(baths$Depth_m))
  
  #calculate schmidt stability
  schmidts[i,2] <- rLakeAnalyzer::schmidt.stability(wtr = temps$Temp_C,
                                                    depths = temps$Depth_m,
                                                    bthA = baths$SA_m2,
                                                    bthD = baths$Depth_m,
                                                    sal = rep(0,length(temps$Temp_C)))
  
}

#now calculate month/yearly means
physics <- schmidts |> 
  group_by(Date) |> 
  dplyr::summarise(SS = mean(SS)) |>
  dplyr::filter(Date %in% dates)

#read in chem from edi
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/199/11/509f39850b6f95628d10889d66885b76" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

#note that for 3 dates the epi and hypo values are the same bc the wc is isothermal
chem <- read.csv(infile1, header = TRUE) |> 
  mutate(DateTime = as.Date(DateTime)) |>
  filter(DateTime %in% dates, Site == 50, Reservoir == "BVR") |>
  mutate(month = month(DateTime), year = year(DateTime)) |>
  select(Depth_m, Rep, TN_ugL, TP_ugL, DateTime) |>
  left_join(ctd_thermo_depth |> distinct(DateTime, .keep_all = TRUE),
            by = "DateTime") |>
  group_by(DateTime) |>
  summarise(TN_ugL_epi = ifelse(first(therm_depth) == 0,mean(TN_ugL, na.rm = TRUE),
      mean(TN_ugL[Depth_m > plyr::round_any(first(therm_depth), 1)], na.rm = TRUE)),
    TN_ugL_hypo = ifelse(first(therm_depth) == 0,mean(TN_ugL, na.rm = TRUE),
      mean(TN_ugL[Depth_m <= plyr::round_any(first(therm_depth), 1)], na.rm = TRUE)),
    TP_ugL_epi = ifelse(first(therm_depth) == 0,mean(TP_ugL, na.rm = TRUE),
      mean(TP_ugL[Depth_m > plyr::round_any(first(therm_depth), 1)], na.rm = TRUE)),
    TP_ugL_hypo = ifelse(first(therm_depth) == 0,mean(TP_ugL, na.rm = TRUE),
      mean(TP_ugL[Depth_m <= plyr::round_any(first(therm_depth), 1)], na.rm = TRUE)),
    .groups = "drop") |>
  arrange(DateTime)


#read in secchi data (not using secchi bc missing 7 obs on zoop collection days)
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/198/11/81f396b3e910d3359907b7264e689052" 
#infile1 <- tempfile()
#try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

#secchi <-read.csv(infile1) |> 
#  mutate(DateTime = as.Date(DateTime)) |> 
#  filter(Reservoir == "BVR" & Site == 50 &
#           DateTime %in% dates) |> 
#  distinct() |> 
#  mutate(year = format(DateTime, "%Y"),
#         month = format(DateTime, "%m")) |> 
#  group_by(DateTime) |>
#  summarise(secchi = mean(Secchi_m))

#read in nldas met data
nldas <- read.csv("./inputs/BVR_GLM_NLDAS_010113_123121_GMTadjusted.csv") |> 
  dplyr::mutate(time = date(time)) |> 
  dplyr::group_by(time) |>
  dplyr::summarise(AirTemp = mean(AirTemp),
                   Shortwave = mean(ShortWave),
                   Longwave = mean(LongWave),
                   RelHum = mean(RelHum),
                   WindSpeed = mean(WindSpeed),
                   Rain = mean(Rain)) |>
  rename(DateTime = time) |>
  dplyr::filter(DateTime %in% dates)

#read in fp data (not doing this either bc only n=51 obs)
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/272/8/0359840d24028e6522f8998bd41b544e" 
#infile1 <- tempfile()
#try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

#fp <- read.csv(infile1) |> 
#  dplyr::mutate(DateTime = as.Date(DateTime)) |>
#  dplyr::filter(DateTime %in% dates & 
#                  Site ==50 & Reservoir %in% c("BVR")) |> 
#  dplyr::mutate(rdepth = plyr::round_any(Depth_m, 1)) |>
#  dplyr::mutate(month = month(DateTime),
#                year = year(DateTime)) |>
#  dplyr::select(-Depth_m) |> 
#  dplyr::rename(Depth_m = rdepth) |> 
#  select(!c(CastID,Temp_C:Flag_RFU_470nm)) |> 
#  dplyr::left_join(ctd_thermo_depth, by = "DateTime") |> 
#  dplyr::group_by(DateTime) |> 
#  dplyr::summarise(Green_ugL = mean(GreenAlgae_ugL[Depth_m < plyr::round_any(therm_depth,1)], na.rm=T), 
#                   Bluegreen_ugL = mean(Bluegreens_ugL[Depth_m < plyr::round_any(therm_depth,1)], na.rm=T), 
#                   Brown_ugL = mean(BrownAlgae_ugL[Depth_m < plyr::round_any(therm_depth,1)], na.rm=T),
#                   Mixed_ugL = mean(MixedAlgae_ugL[Depth_m < plyr::round_any(therm_depth,1)], na.rm=T),
#                   Total_ugL = mean(TotalConc_ugL[Depth_m < plyr::round_any(therm_depth,1)], na.rm=T))

#calculate residence time (volume / total inflow)
#EDI bvr bathy = 1357140.6 m3 (vol at full pond)
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/1254/1/f7fa2a06e1229ee75ea39eb586577184" 
infile1 <- tempfile()
try(download.file(inUrl1, infile1, timeout = max(300, getOption("timeout"))))

vol <- readr::read_csv(infile1, show_col_types = FALSE) |>
  dplyr::select(Reservoir, Depth_m, Volume_below_L) |>
  dplyr::filter(Reservoir == "BVR") |>
  dplyr::mutate(Volume_m3 = Volume_below_L / 1000) |>
  dplyr::select(-Volume_below_L)

daily_wl <- read.csv("./Output/BVR_WaterLevel_2014_2022_interp.csv") |>
  dplyr::filter(Date <= "2021-12-31")

# Maximum depth with nonzero volume (e.g., 13 m in your case)
max_nonzero_depth <- max(vol$Depth_m[vol$Volume_m3 > 0], na.rm = TRUE)

# Calculate volume at each daily depth using 13 − depth mapping
daily_wl <- daily_wl %>%
  rowwise() %>%
  mutate(Volume_m3 = {
    # Map measured depth to depth in vol table
    lookup_depth <- max_nonzero_depth - WaterLevel_m
    
    # Clamp within valid range
    lookup_depth <- max(min(lookup_depth, max(vol$Depth_m)), min(vol$Depth_m))
    
    # Interpolate volume at that lookup depth
    approx(
      x = vol$Depth_m,
      y = vol$Volume_m3,
      xout = lookup_depth,
      rule = 2
    )$y
  }) %>%
  ungroup()

res_time <- read.csv("Output/BVR_flow_calcs_NLDAS_2014-2021.csv") |> 
  dplyr::filter(time %in% dates) |>  
  dplyr::left_join(daily_wl, by = c("time" = "Date")) |>  # join volumes
  dplyr::mutate(res_time_d = Volume_m3 / Q_m3pd) |> 
  dplyr::group_by(time) |> 
  dplyr::summarise(res_time = mean(res_time_d), .groups = "drop") |>
  dplyr::rename(DateTime = time)
# the 2021-12-06 res time is super high but I have no reason to think it's incorrect....

# create combined env df
all_drivers <- bind_cols(ctd_thermo_depth, ctd_oxy_depth[!colnames(ctd_oxy_depth) %in% c("DateTime")], 
                         profiles[!colnames(profiles) %in% c("DateTime", "Temp_C_epi")], #colinear
                         water_level[!colnames(water_level) %in% c("Date")],
                         physics[!colnames(physics) %in% c("Date")],
                         chem[!colnames(chem) %in% c("DateTime")],
                         nldas[!colnames(nldas) %in% c("DateTime", "AirTemp","Longwave")], #these vars are colinear so removing
                         res_time[!colnames(res_time) %in% c("DateTime")])
#write.csv(all_drivers, "./Output/RDA_env_drivers.csv", row.names=FALSE)
                         
#------------------------------------------------------------------------------#
# Make sure the rows match between predictors and responses
stopifnot(nrow(all_drivers) == nrow(all_zoops))

#check gradient lengths to determine whether to do RDA vs. CCA
dca <- decorana(zoop_dens_trans)   
# axis lengths are < 3, so RDA should be good

#drop datetime cols
all_drivers_num <- all_drivers[ , !names(all_drivers) %in% c("DateTime")]

# Run the RDA
rda_model <- rda(zoop_dens_trans ~ ., data = all_drivers_num)

#next see whether env vars are colinear (VIF>10)
vif.cca(rda_model) #remove airtemp, longwave, Temp_C_epi from above env df

# 39.3% of total zooplankton variation is explained by env variables; 60.6% is unexplained
summary(rda_model)

# ANOVA to test significance
anova(rda_model, by = "axis", permutations = 999) 
#RDA1 and RDA2 explain a significant fraction of the variance!

#term
anova(rda_model, by = "term", permutations = 999)
#SS, water level, hypo temp, and oxycline depth are significant drivers of zoop community structure

# Site (sample) scores
site_scores <- vegan::scores(rda_model, display = "sites", choices = c(1,2)) |>
  as.data.frame() |>
  mutate(year = format(all_zoops$DateTime, "%Y"))

# Species scores
species_scores <- vegan::scores(rda_model, display = "species", choices = c(1,2)) |>
  as.data.frame() |>
  mutate(Species = rownames(.))

# List of significant drivers from anova(rda_model, by="term")
sig_drivers <- c("oxy_depth", "Temp_C_hypo", "waterlevel", "SS")

# Environmental (biplot) scores
env_scores <- vegan::scores(rda_model, display = "bp", choices = c(1,2)) %>%
  as.data.frame() %>%
  mutate(Env = rownames(.),
         Significant = ifelse(Env %in% sig_drivers, "yes", "no"),
         RDA1_end = RDA1 * arrow_mult,
         RDA2_end = RDA2 * arrow_mult)

# Multiply arrows by a constant for visibility
arrow_mult <- 2
env_scores <- env_scores %>%
  mutate(RDA1_end = RDA1 * arrow_mult,
         RDA2_end = RDA2 * arrow_mult)

ggplot() +
  geom_point(data = site_scores, aes(x = RDA1, y = RDA2, color=year), size = 2) +
  #geom_text(data = species_scores, aes(x = RDA1, y = RDA2, label = Species),
  #          color = "black", size = 3, vjust = -0.5) +
  scale_color_viridis_d(name = "year") + 
  ggnewscale::new_scale_color() +
  geom_segment(data = env_scores, 
               aes(x = 0, y = 0, xend = RDA1_end, yend = RDA2_end, color = Significant),
               arrow = arrow(length = unit(0.3,"cm")), show.legend = FALSE) +
  #geom_text_repel(data = env_scores, 
  #                aes(x = RDA1_end, y = RDA2_end, label = Env),
  #               color = ifelse(env_scores$Significant == "yes", "red", "grey50"),
  #                size = 3) +
  geom_text(data = env_scores, 
            aes(x = RDA1_end, y = RDA2_end, label = Env, color = Significant),
            size = 2, vjust = -0.5, show.legend = FALSE) +
  scale_color_manual(values = c("yes" = "red", "no" = "grey50")) +
  xlab(paste0("RDA1 (", round(summary(rda_model)$cont$importance[2,1]*100,1), "%)")) +
  ylab(paste0("RDA2 (", round(summary(rda_model)$cont$importance[2,2]*100,1), "%)")) +
  theme_minimal() +
  theme(text = element_text(size=10), 
        legend.position = "right",
        axis.text = element_text(size=6, color="black"), 
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)) 
#ggsave("Figures/zoop_RDA.jpg", width=5, height=4) 

#now see how much variation the 4 sig env vars make up
rda_model_sig <- rda(zoop_dens_trans ~ oxy_depth + SS + Temp_C_hypo + waterlevel,
                 data = all_drivers_num)

# Overall RDA results
summary(rda_model_sig) #these vars explain 23.4% of the total variation (collectively)

#variation of each driver individually
varpart(zoop_dens_trans,
        all_drivers_num$oxy_depth,
        all_drivers_num$SS,
        all_drivers_num$Temp_C_hypo,
        all_drivers_num$waterlevel)
#oxycline depth = 46%; SS = 8.7%; hypo temp = 2.7%; water level = 4.2% 