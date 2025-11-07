#code to visualize seasonal and multi-annual zoop patterns

#read in packages, starting with the pacman package (un-comment and run 
#the next line if pacman is not already installed)
#install.packages('pacman')
pacman::p_load(tidyverse, NatParksPalettes, rLakeAnalyzer,
               EDIutils, xml2, gsheet, dplyr, zoo, lubridate)

#------------------------------------------------------------------------------#
#Pull in environmental data for 2014-2023

#read in zoop df to get dates
zoop <- read.csv("Output/all_zoops_dens.csv")
  zoop_dates <- unique(zoop$DateTime)

#dataframe with all zoop dates
all_dates_df <- data.frame(DateTime = as.Date(zoop_dates))
#just gonna keep the first MSN obs

#read in ctd (200.14)
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/14/0432a298a90b2b662f26c46071f66b8a" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

ctd <-read.csv(infile1,header=T) |> 
  mutate(DateTime = as.Date(DateTime)) |>
  filter(DateTime %in% all_dates_df$DateTime & 
           Depth_m > 0 & Site ==50 & Reservoir %in% c("BVR"),
         !is.na(Temp_C)) |> 
  select(Reservoir, Site, Depth_m, DateTime,
         Temp_C, DO_mgL)  #removing chl a bc fp and ctd aren't so comparable so can't fill in missing days

  
#sub missing ctd data with ysi data (198.12)
inUrl1  <-  "https://pasta.lternet.edu/package/data/eml/edi/198/12/e0181372c6d4cbc66eca4644f8385470" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(500, getOption("timeout"))))

ysi <- read.csv(infile1,header=T) |> 
  mutate(DateTime = as.Date(DateTime)) |>
  filter(!DateTime %in% ctd$DateTime & 
           Depth_m > 0 &
           Site ==50 & Reservoir %in% c("BVR")) |> 
  select(Reservoir, Site, DateTime, Depth_m, 
         Temp_C, DO_mgL) 

#select every 0.5m from casts for thermocline calc below
ctd_final_temp_do <- ctd |>
  dplyr::mutate(rdepth = plyr::round_any(Depth_m, 0.5)) |>
  dplyr::group_by(DateTime, rdepth, Reservoir, Site) |>
  dplyr::summarise(temp = mean(Temp_C, na.rm=T),
                   DO = mean(DO_mgL, na.rm=T)) |>
  dplyr::rename(depth = rdepth) |> ungroup()

#rename columns
ctd_final_temp_do <- ctd_final_temp_do |> 
  rename(Depth_m = depth,
         Temp_C = temp,
         DO_mgL = DO)

#clean up ysi
ysi_clean <- ysi |> 
  select(Reservoir,DateTime, Depth_m,Temp_C, DO_mgL) |> 
  filter(!is.na(Temp_C)) |>
  distinct(Reservoir, DateTime, Depth_m, .keep_all = TRUE)

#combine ctd and ysi dfs 
temp_final <- ctd_final_temp_do |> 
  select(DateTime, Depth_m,Temp_C) |> 
  full_join(ysi_clean, multiple = 'all', 
            by = c('DateTime', 'Depth_m',"Temp_C")) |>
  arrange(DateTime, Depth_m) |> 
  select(-c(Reservoir,DO_mgL))

do_final <- ctd_final_temp_do |> 
  select(DateTime, Depth_m, DO_mgL) |> 
  full_join(ysi_clean, multiple = 'all', 
            by = c('DateTime', 'Depth_m',"DO_mgL")) |>
  arrange(DateTime, Depth_m) |> 
  select(-c(Reservoir,Temp_C))

#calculate thermocline depth using combined ctd and ysi df
ctd_thermo_depth_all <- temp_final |> 
  group_by(DateTime) |> 
  summarise(therm_depth = thermo.depth(Temp_C,Depth_m))

ctd_thermo_depth <- ctd_thermo_depth_all |> 
  ungroup() |>
  group_by(DateTime) |>
  summarise(therm_depth = mean(therm_depth)) |>
  dplyr::filter(DateTime %in% all_dates_df$DateTime)

#calculate oxycline depth
ctd_oxy_depth <- do_final |> 
  group_by(DateTime) |> 
  summarise(oxy_depth = thermo.depth(DO_mgL,Depth_m)) |> 
  ungroup() |>
  group_by(DateTime) |> 
  summarise(oxy_depth = mean(oxy_depth, na.rm=T)) |>
  dplyr::filter(DateTime %in% all_dates_df$DateTime)

#round ctd to nearest m
ctd_final <- ctd |>
  dplyr::mutate(rdepth = plyr::round_any(Depth_m, 1)) |>
  dplyr::select(-Depth_m) |> 
  dplyr::rename(Depth_m = rdepth) |> 
  dplyr::left_join(ctd_thermo_depth, by = "DateTime") |> 
  group_by(DateTime) |> 
  dplyr::summarise(Temp_C_epi = mean(Temp_C[Depth_m < plyr::round_any(therm_depth,1)], na.rm=T),
                   Temp_C_hypo = mean(Temp_C[Depth_m >= plyr::round_any(therm_depth,1)], na.rm=T),
                   DO_mgL_epi = mean(DO_mgL[Depth_m < plyr::round_any(therm_depth,1)], na.rm=T))

ysi_final <- ysi |> 
  dplyr::left_join(ctd_thermo_depth, by = "DateTime") |> 
  group_by(DateTime) |> 
  summarise(Temp_C_epi = mean(Temp_C[Depth_m < plyr::round_any(therm_depth,1)], na.rm=T),
            Temp_C_hypo = mean(Temp_C[Depth_m >= plyr::round_any(therm_depth,1)], na.rm=T),
            DO_mgL_epi = mean(DO_mgL[Depth_m < plyr::round_any(therm_depth,1)], na.rm=T))
  
#combine ctd and ysi to account for missing ctd days
profiles <- bind_rows(ctd_final, ysi_final) |> arrange(DateTime) |>
  dplyr::filter(DateTime %in% all_dates_df$DateTime)

available_dates <- unique(profiles$DateTime)
missing_dates <- setdiff(all_dates_df$DateTime, available_dates)

nearest_numeric <- sapply(missing_dates, function(md) {
  available_dates[which.min(abs(available_dates - md))]})

# Convert numeric back to Date
nearest_dates <- as.Date(nearest_numeric, origin = "1970-01-01")

nearest_df <- data.frame(
  missing_date = missing_dates,
  nearest_available = nearest_dates) |>
  mutate(diff = abs(as.numeric(nearest_available - missing_date)),
         replacement_date = ifelse(diff <= 7, nearest_available, as.Date(NA)),
         replacement_date = as.Date(replacement_date))

#------------------------------------------------------------------------------#
#download bathymetry
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/1254/1/f7fa2a06e1229ee75ea39eb586577184" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

bathymetry <- readr::read_csv(infile1, show_col_types = F)  |>
  dplyr::select(Reservoir, Depth_m, SA_m2) |>
  dplyr::filter(Reservoir == "BVR")

#read in water level (from thermistors aka bvrplatform edi 725.4)
water_level <- read.csv("./Output/BVR_WaterLevel_2014_2023_interp.csv") |> 
  select(Date, WaterLevel_m) |> 
  mutate(WaterLevel_m = as.numeric(WaterLevel_m)) |> 
  filter(Date >= "2014-01-01" & Date < "2024-01-01") |> 
  group_by(Date) |> 
  summarise(waterlevel = mean(WaterLevel_m,na.rm=T)) |> 
  arrange(Date) |> ungroup() |> 
  mutate(Date = as.Date(Date)) |> 
  dplyr::filter(Date %in% all_dates_df$DateTime) |>
  rename(DateTime = Date)

wl <- read.csv("./Output/BVR_WaterLevel_2014_2023_interp.csv") |> 
  select(Date, WaterLevel_m) |> 
  mutate(WaterLevel_m = as.numeric(WaterLevel_m)) |> 
  mutate(Date = as.Date(Date)) 

#add res column to wl df
wl$Reservoir <- "BVR"

#only select wl on days where we have temp obs
wl <- wl |> filter(Date %in% unique(temp_final$DateTime))

#combine temp_final and wl
temp_final_wl <- temp_final |> 
  dplyr::rename(Date = DateTime) |> 
  dplyr::full_join(wl, multiple = 'all', by = 'Date') |> 
  dplyr::filter(!is.na(Temp_C))

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
                        "SS" = NA,
                        "BF" = NA)

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
                                                  sal = rep(0,length(temps$Temp_C)))}

#now calculate month/yearly means
physics <- schmidts |> 
  group_by(Date) |> 
  dplyr::summarise(SS = mean(SS)) |>
  dplyr::filter(Date %in% all_dates_df$DateTime) |>
  rename(DateTime = Date)
  
#------------------------------------------------------------------------------#
#read in chem from edi (199.12)
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/199/12/a33a5283120c56e90ea414e76d5b7ddb" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

chem <-read.csv(infile1,header=T) |> 
  dplyr::mutate(DateTime = as.Date(DateTime)) |>
  dplyr::filter(DateTime %in% all_dates_df$DateTime & 
          Site ==50 & Reservoir %in% c("BVR")) |> 
  dplyr::select(Reservoir, Depth_m, Rep,
         TN_ugL, TP_ugL, DateTime) |>
  dplyr::left_join(ctd_thermo_depth, by = "DateTime") |>
  dplyr::group_by(DateTime) |> 
  dplyr::summarise(TN_ugL_epi = mean(TN_ugL[Depth_m > plyr::round_any(therm_depth,1)], na.rm=T),
                   TN_ugL_hypo = mean(TN_ugL[Depth_m <= plyr::round_any(therm_depth,1)], na.rm=T),
                   TP_ugL_epi = mean(TP_ugL[Depth_m > plyr::round_any(therm_depth,1)], na.rm=T),
                   TP_ugL_hypo = mean(TP_ugL[Depth_m <= plyr::round_any(therm_depth,1)], na.rm=T)) |> 
  dplyr::arrange(DateTime) |>
  dplyr::filter(DateTime %in% all_dates_df$DateTime)

#-----------------------------------------------------------------------------#
#read in secchi data (198.12)
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/198/12/80bec97dc53d85b0298a72bb1a098442" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

secchi <-read.csv(infile1) |> 
  mutate(DateTime = as.Date(DateTime)) |> 
  filter(Reservoir == "BVR" & Site == 50 &
           DateTime %in% all_dates_df$DateTime) |> 
  distinct() |> 
  group_by(DateTime) |>
  summarise(secchi = mean(Secchi_m)) |>
  mutate(DateTime = as.Date(ifelse(DateTime == as.Date("2014-09-25"), 
                                   as.Date("2014-10-23"), DateTime))) |>
  dplyr::filter(DateTime %in% all_dates_df$DateTime)
#counting the second 2014 secchi in sep as an oct obs

available_dates <- unique(secchi$DateTime)
missing_dates <- setdiff(all_dates_df$DateTime, available_dates)

nearest_numeric <- sapply(missing_dates, function(md) {
  available_dates[which.min(abs(available_dates - md))]})

# Convert numeric back to Date
nearest_dates <- as.Date(nearest_numeric, origin = "1970-01-01")

nearest_df <- data.frame(
  missing_date = missing_dates,
  nearest_available = nearest_dates) |>
  mutate(diff = abs(as.numeric(nearest_available - missing_date)),
         replacement_date = ifelse(diff <= 7, nearest_available, as.Date(NA)),
         replacement_date = as.Date(replacement_date)) 

#can't replace with closest date bc these dates already correspond with zoop samples
#replacement_secchi <- nearest_df |>
#  filter(!is.na(replacement_date)) |>
#  left_join(secchi, by = c("replacement_date" = "DateTime")) |>
#  mutate(DateTime = missing_date) |> #overwrite datetime to match the zoop date
#  select(-replacement_date, -diff, -missing_date, -nearest_available)

#combine dfs
#secchi_filled <- bind_rows(secchi, replacement_secchi) |>
#  arrange(DateTime)

# Left join the combined data to include all dates
#secchi_filled <- all_dates_df |>
#  left_join(secchi_filled, by = c("DateTime")) |>
#  arrange(DateTime)

#------------------------------------------------------------------------------#
#read in nldas met data
nldas <- read.csv("./inputs/BVR_GLM_NLDAS_2013-2023_hourly.csv") |> 
  dplyr::mutate(time = date(time)) |> 
  dplyr::filter(time %in% all_dates_df$DateTime) |> 
  dplyr::group_by(time) |>
  dplyr::summarise(AirTemp = mean(AirTemp),
                   Shortwave = mean(ShortWave),
                   Longwave = mean(LongWave),
                   RelHum = mean(RelHum),
                   WindSpeed = mean(WindSpeed),
                   Rain = mean(Rain)) |>
  rename(DateTime = time) |>
  dplyr::filter(DateTime %in% all_dates_df$DateTime)

#------------------------------------------------------------------------------#
#read in fp data
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/272/8/0359840d24028e6522f8998bd41b544e" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

fp <- read.csv(infile1) |> 
  dplyr::mutate(DateTime = as.Date(DateTime)) |>
  dplyr::filter(DateTime %in% all_dates_df$DateTime & 
                  Site ==50 & Reservoir %in% c("BVR")) |> 
  dplyr::mutate(rdepth = plyr::round_any(Depth_m, 1)) |>
  dplyr::mutate(month = month(DateTime),
                year = year(DateTime)) |>
  dplyr::select(-Depth_m) |> 
  dplyr::rename(Depth_m = rdepth) |> 
  select(!c(CastID,Temp_C:Flag_RFU_470nm)) |> 
  dplyr::left_join(ctd_thermo_depth, by = "DateTime") |> 
  dplyr::group_by(DateTime) |> 
  dplyr::summarise(Green_ugL = mean(GreenAlgae_ugL[Depth_m < plyr::round_any(therm_depth,1)], na.rm=T), 
                   Bluegreen_ugL = mean(Bluegreens_ugL[Depth_m < plyr::round_any(therm_depth,1)], na.rm=T), 
                   Brown_ugL = mean(BrownAlgae_ugL[Depth_m < plyr::round_any(therm_depth,1)], na.rm=T),
                   Mixed_ugL = mean(MixedAlgae_ugL[Depth_m < plyr::round_any(therm_depth,1)], na.rm=T),
                   Total_ugL = mean(TotalConc_ugL[Depth_m < plyr::round_any(therm_depth,1)], na.rm=T))

available_dates <- unique(fp$DateTime)
missing_dates <- setdiff(all_dates_df$DateTime, available_dates)

nearest_numeric <- sapply(missing_dates, function(md) {
  available_dates[which.min(abs(available_dates - md))]})

# Convert numeric back to Date
nearest_dates <- as.Date(nearest_numeric, origin = "1970-01-01")

nearest_df <- data.frame(
  missing_date = missing_dates,
  nearest_available = nearest_dates) |>
  mutate(diff = abs(as.numeric(nearest_available - missing_date)),
         replacement_date = ifelse(diff <= 7, nearest_available, as.Date(NA)),
         replacement_date = as.Date(replacement_date))

#can't replace with closest date bc these dates already correspond with zoop samples
#replacement_fp <- nearest_df |>
#  filter(!is.na(replacement_date)) |>
#  left_join(fp, by = c("replacement_date" = "DateTime")) |>
#  mutate(DateTime = missing_date) |> #overwrite datetime to match the zoop date
#  select(-replacement_date, -diff, -missing_date, -nearest_available)

#combine dfs
#fp_filled <- bind_rows(fp, replacement_fp) |>
#  arrange(DateTime)

# Left join the combined data to include all dates
#fp_filled <- all_dates_df |>
#  left_join(fp_filled, by = c("DateTime")) |>
#  arrange(DateTime)

#------------------------------------------------------------------------------#
#now merge all vars
all_drivers <- all_dates_df |>
  left_join(chem, by = "DateTime") |>
  left_join(profiles, by = "DateTime") |>
  left_join(water_level, by = "DateTime") |>
  left_join(ctd_thermo_depth, by = "DateTime") |>
  left_join(ctd_oxy_depth, by = "DateTime") |>
  left_join(physics, by = "DateTime") |>
  left_join(fp, by = "DateTime") |>
  left_join(secchi, by = "DateTime") |>
  left_join(nldas, by = "DateTime") |>
  mutate(across(where(is.numeric), ~ifelse(is.nan(.x), NA, .x))) 
  
#----------------------------------------
#add doy cols
all_drivers <- all_drivers |>
  mutate(DOY = yday(DateTime), 
         Year = year(DateTime),
         Month = month(DateTime))

# remove non-numeric column before imputation
all_drivers_num <- all_drivers[, sapply(all_drivers, is.numeric)]

#interpolate within each year first to maintain seasonality
all_drivers_fill <- all_drivers_num |>
  group_by(Year) |> 
  arrange(DOY) |>
  mutate(across(where(is.numeric), 
                ~ na.approx(.x, x = DOY, na.rm = FALSE, maxgap = 15))) |>
  ungroup()

#fill remaining NAs using monthly averages
all_drivers_fill <- all_drivers_fill |>
  group_by(Month) |>
  mutate(across(where(is.numeric), ~ 
                  ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))) |>
  ungroup() |>
  arrange(Year, Month, DOY) 

#add in date col again
all_drivers_fill$DateTime <- all_drivers$DateTime 

#drop unused cols
all_drivers_fill <- all_drivers_fill |>
  filter(!Month %in% c(3,12)) |> #dropping 4 days that couldn't be imputed
  dplyr::select(-c(Year,Month,DOY)) 

#export csv; n = 83 days; secchi has a lot of NAs...
write.csv(all_drivers_fill, "./Output/all_drivers.csv", row.names=FALSE)
