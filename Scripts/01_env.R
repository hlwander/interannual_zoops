#code to visualize seasonal and multiannual zoop patterns

#read in packages, starting with the pacman package if needed (uncomment and run 
#the next line if pacman is not already installed)
#install.packages('pacman')
pacman::p_load(tidyverse, NatParksPalettes, rLakeAnalyzer,
               EDIutils, xml2, gsheet)

#------------------------------------------------------------------------------#
#Pull in environmental data for 2014-2023

#list of dates only between May-Oct 
dates_list <- c(seq(as.Date("2014-05-01"), as.Date("2014-10-31"), by="days"),
                seq(as.Date("2015-05-01"), as.Date("2015-10-31"), by="days"),
                seq(as.Date("2016-05-01"), as.Date("2016-10-31"), by="days"),
                seq(as.Date("2019-05-01"), as.Date("2019-10-31"), by="days"),
                seq(as.Date("2020-05-01"), as.Date("2020-10-31"), by="days"),
                seq(as.Date("2021-05-01"), as.Date("2021-10-31"), by="days"),
                seq(as.Date("2023-05-01"), as.Date("2023-10-31"), by="days"))

#read in zoop df to get dates
zoop <- read.csv("Output/all_zoops_dens.csv")
  zoop_dates <- unique(zoop$DateTime)

#read in ctd (200.14)
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/14/0432a298a90b2b662f26c46071f66b8a" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

ctd <-read.csv(infile1,header=T) |> 
  mutate(DateTime = as.Date(DateTime)) |>
  filter(DateTime %in% dates_list & 
           Depth_m > 0 & Site ==50 & Reservoir %in% c("BVR"),
         !is.na(Temp_C)) |> 
  select(Reservoir, Site, Depth_m, DateTime,
         Temp_C, DO_mgL) #removing chl a bc fp and ctd aren't so comparable so can't fill in missing days

unique(ctd$DateTime)

#missing Aug 2015, May 2019, May 2020, Aug + Sep + Oct 2021
missing_dates <- c(seq(as.Date("2015-08-01"), as.Date("2015-08-31"), by="days"),
                   seq(as.Date("2019-05-01"), as.Date("2019-05-31"), by="days"),
                   seq(as.Date("2020-05-01"), as.Date("2020-05-31"), by="days"),
                   seq(as.Date("2021-08-01"), as.Date("2021-10-31"), by="days"))

#sub missing ctd data with ysi data (198.12)
inUrl1  <-  "https://pasta.lternet.edu/package/data/eml/edi/198/12/e0181372c6d4cbc66eca4644f8385470" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(500, getOption("timeout"))))

ysi <- read.csv(infile1,header=T) |> 
  mutate(DateTime = as.Date(DateTime)) |>
  filter(DateTime %in% missing_dates & 
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
  dplyr::rename(depth = rdepth) 

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
temp_final <-ctd_final_temp_do |> 
  select(DateTime, Depth_m,Temp_C) |> 
  full_join(ysi_clean, multiple = 'all', 
            by = c('DateTime', 'Depth_m',"Temp_C")) |>
  arrange(DateTime, Depth_m) |> 
  select(-c(Reservoir.x,Reservoir.y,DO_mgL))

do_final <- ctd_final_temp_do |> 
  select(DateTime, Depth_m, DO_mgL) |> 
  full_join(ysi_clean, multiple = 'all', 
            by = c('DateTime', 'Depth_m',"DO_mgL")) |>
  arrange(DateTime, Depth_m) |> 
  select(-c(Reservoir.x,Reservoir.y,Temp_C))

#calculate thermocline depth using combined ctd and ysi df
ctd_thermo_depth_all <- temp_final |> 
  group_by(DateTime) |> 
  summarise(therm_depth = thermo.depth(Temp_C,Depth_m)) |> 
  mutate(month = month(DateTime),
         year = year(DateTime)) 

ctd_thermo_depth <- ctd_thermo_depth_all |> 
  ungroup() |>
  #group_by(month, year) |> 
  group_by(DateTime) |>
  summarise(therm_depth = mean(therm_depth)) |>
  dplyr::filter(DateTime %in% zoop_dates)

#calculate oxycline depth
ctd_oxy_depth <- do_final |> 
  group_by(DateTime) |> 
  summarise(oxy_depth = thermo.depth(DO_mgL,Depth_m)) |> 
  mutate(month = month(DateTime),
         year = year(DateTime)) |> 
  ungroup() |>
  group_by(DateTime) |> 
  #group_by(month, year) |> 
  summarise(oxy_depth = mean(oxy_depth, na.rm=T)) |>
  dplyr::filter(DateTime %in% zoop_dates)

#round ctd to nearest m
ctd_final <- ctd |>
  dplyr::mutate(rdepth = plyr::round_any(Depth_m, 1)) |>
  dplyr::mutate(month = month(DateTime),
                year = year(DateTime)) |>
  dplyr::select(-Depth_m) |> 
  dplyr::rename(Depth_m = rdepth) |> 
  dplyr::left_join(ctd_thermo_depth, by = "DateTime") |> #c("month","year")
  group_by(DateTime) |> 
  #dplyr::group_by(month, year) |> 
  dplyr::summarise(Temp_C_epi = mean(Temp_C[Depth_m < plyr::round_any(therm_depth,1)], na.rm=T),
                   Temp_C_hypo = mean(Temp_C[Depth_m >= plyr::round_any(therm_depth,1)], na.rm=T),
                   DO_mgL_epi = mean(DO_mgL[Depth_m < plyr::round_any(therm_depth,1)], na.rm=T))

ysi_final <- ysi |> 
  mutate(month = month(DateTime),
                year = year(DateTime)) |> 
  dplyr::left_join(ctd_thermo_depth, by = "DateTime") |> #c("month","year")
  #group_by(month, year) |>   
  group_by(DateTime) |> 
  summarise(Temp_C_epi = mean(Temp_C[Depth_m < plyr::round_any(therm_depth,1)], na.rm=T),
            Temp_C_hypo = mean(Temp_C[Depth_m >= plyr::round_any(therm_depth,1)], na.rm=T),
            DO_mgL_epi = mean(DO_mgL[Depth_m < plyr::round_any(therm_depth,1)], na.rm=T))
  
#combine ctd and ysi to account for missing ctd days
profiles <- bind_rows(ctd_final, ysi_final) |> arrange(DateTime) |>
  dplyr::filter(DateTime %in% zoop_dates)
  #arrange(month, year) 

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
  mutate(month = month(Date),
         year = year(Date)) |> 
  group_by(Date) |> 
  #group_by(month, year) |> 
  summarise(waterlevel = mean(WaterLevel_m,na.rm=T)) |> 
  #filter(month %in% c(5,6,7,8,9),
  #       year %in% c(2014:2016,2019:2021)) |> 
  filter(month(Date) %in% c(5,6,7,8,9),
         year(Date) %in% c(2014:2016,2019:2021)) |> 
  #arrange(month,year) |> 
  arrange(Date) |>
  ungroup() |> 
 # group_by(year) |> 
  group_by(Date) |> 
  mutate(wl_cv = sd(waterlevel)/mean(waterlevel)) |> 
  ungroup() |>
  dplyr::filter(Date %in% zoop_dates)

wl <- read.csv("./Output/BVR_WaterLevel_2014_2023_interp.csv") |> 
  select(Date, WaterLevel_m) |> 
  mutate(WaterLevel_m = as.numeric(WaterLevel_m)) |> 
  mutate(Date = as.Date(Date)) |> 
  filter(Date >= "2014-01-01" & Date < "2024-01-01" &
           month(Date) %in% c(5,6,7,8,9,10))

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
                                                  sal = rep(0,length(temps$Temp_C)))
  
}

#now calculate month/yearly means
physics <- schmidts |> 
  dplyr::mutate(month = month(Date),
                year = year(Date)) |> 
 # dplyr::group_by(month, year) |> 
  group_by(Date) |> 
  dplyr::summarise(SS = mean(SS)) |>
  dplyr::filter(Date %in% zoop_dates)
  
#read in chem from edi (199.12)
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/199/12/a33a5283120c56e90ea414e76d5b7ddb" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

chem <-read.csv(infile1,header=T) |> 
  dplyr::mutate(DateTime = as.Date(DateTime)) |>
  dplyr::filter(DateTime %in% dates_list & 
          Site ==50 & Reservoir %in% c("BVR")) |> 
  dplyr::mutate(month = month(DateTime),
                year = year(DateTime)) |> 
  dplyr::select(Reservoir, Depth_m, Rep,
         TN_ugL, TP_ugL, DateTime) |> #month, year
  dplyr::left_join(ctd_thermo_depth, by = "DateTime") |>  #c("month","year")
  #dplyr::group_by(month, year) |> 
  dplyr::group_by(DateTime) |> 
  dplyr::summarise(TN_ugL_epi = mean(TN_ugL[Depth_m > plyr::round_any(therm_depth,1)], na.rm=T),
                   TN_ugL_hypo = mean(TN_ugL[Depth_m <= plyr::round_any(therm_depth,1)], na.rm=T),
                   TP_ugL_epi = mean(TP_ugL[Depth_m > plyr::round_any(therm_depth,1)], na.rm=T),
                   TP_ugL_hypo = mean(TP_ugL[Depth_m <= plyr::round_any(therm_depth,1)], na.rm=T)) |> 
  dplyr::arrange(DateTime) |>
  dplyr::filter(DateTime %in% zoop_dates)
  #dplyr::arrange(month, year)

#read in secchi data (198.12)
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/198/12/80bec97dc53d85b0298a72bb1a098442" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

secchi <-read.csv(infile1) |> 
  mutate(DateTime = as.Date(DateTime)) |> 
  filter(Reservoir == "BVR" & Site == 50 &
           DateTime %in% dates_list) |> 
  distinct() |> 
  mutate(year = format(DateTime, "%Y"),
         month = format(DateTime, "%m")) |> 
  #group_by(year, month) |> 
  group_by(DateTime) |>
  summarise(secchi = mean(Secchi_m)) |>
  mutate(DateTime = as.Date(ifelse(DateTime == as.Date("2014-09-25"), 
                                   as.Date("2014-10-01"), DateTime))) 
#counting the second 2014 secchi in sep as an oct obs

#missing secchi obs (miss forest to fill later)
s_missing <- data.frame("year" = "2015",
                        "month" = "09",
                        "secchi" = NA)

#combine dfs
secchi_df <- bind_rows(secchi, s_missing) |> 
  arrange(DateTime) |>
  dplyr::filter(DateTime %in% zoop_dates)

#read in nldas met data
nldas <- read.csv("./inputs/BVR_GLM_NLDAS_010113_123121_GMTadjusted.csv") |> 
  dplyr::mutate(time = date(time)) |> 
  dplyr::filter(time %in% dates_list) |> 
  dplyr::mutate(month = month(time),
                year = year(time)) |> 
  dplyr::group_by(time) |>
  dplyr::summarise(AirTemp = mean(AirTemp),
                   Shortwave = mean(ShortWave),
                   Longwave = mean(LongWave),
                   RelHum = mean(RelHum),
                   WindSpeed = mean(WindSpeed),
                   Rain = mean(Rain)) |>
  rename(DateTime = time) |>
  dplyr::filter(DateTime %in% zoop_dates)

#------------------------------------------------------------------------------#
#read in fp data
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/272/8/0359840d24028e6522f8998bd41b544e" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

fp <- read.csv(infile1) |> 
  dplyr::mutate(DateTime = as.Date(DateTime)) |>
  dplyr::filter(DateTime %in% dates_list & 
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

#add in NAs for missing months
fp_missing <- data.frame("month" = c(9,10,5),
                         "year" = c(2015,2015,2020),
                         "Green_ugL" = c(NA,NA,NA),
                         "Bluegreen_ugL" = c(NA,NA,NA),
                         "Brown_ugL" = c(NA,NA,NA),
                         "Mixed_ugL" = c(NA,NA,NA),
                         "Total_ugL" = c(NA,NA,NA))

#combine dfs
fp_df <-bind_rows(fp, fp_missing) |> 
  arrange(DateTime) |>
  dplyr::filter(DateTime %in% zoop_dates)


#------------------------------------------------------------------------------#
#now average across years
all_drivers <- bind_cols(chem, profiles[!colnames(profiles) %in% c("DateTime")], #n=54
            water_level[!colnames(water_level) %in%                              #n=61
                          c("Date")],
            ctd_thermo_depth[!colnames(ctd_thermo_depth) %in%                    #n=50
                           c("DateTime")],
            ctd_oxy_depth[!colnames(ctd_oxy_depth) %in%                          #n=50
                           c("DateTime")],
            physics[!colnames(physics) %in%                                      #n=50
                           c("Date")],
            fp_df[!colnames(fp_df) %in%                                          #n=38
                           c("DateTime")],
            secchi_df[!colnames(secchi_df) %in%                                  #n=48
                           c("DateTime")],
            nldas[!colnames(nldas) %in%                                          #n=61
                          c("DateTime")],
            res_time[!colnames(res_time) %in%                                    #n=61 
                          c("DateTime")])

#------------------------------------------------------------------------------#
#miss forest to fill in missing env vars

library(missForest)

# remove non-numeric columns (e.g., Date) before imputation
all_drivers_num <- all_drivers[, sapply(all_drivers, is.numeric)]

# run MissForest
set.seed(42)
env_imputed <- missForest(all_drivers_num, maxiter = 10, ntree = 100)

# replace imputed columns back
env_sub_imputed <- all_drivers
env_sub_imputed[, names(all_drivers_num)] <- env_imputed$ximp




#export csv
write.csv(all_drivers, "./Output/all_drivers.csv", row.names=FALSE)
