#code to visualize seasonal and multiannual zoop patterns

#read in packages
pacman::p_load(tidyverse, NatParksPalettes,rLakeAnalyzer)

#------------------------------------------------------------------------------#
#Pull in environmental data for 2014-2021

#list of dates only between May-Sep 
dates_list <- c(seq(as.Date("2014-05-01"), as.Date("2014-09-30"), by="days"),
                seq(as.Date("2015-05-01"), as.Date("2015-09-30"), by="days"),
                seq(as.Date("2016-05-01"), as.Date("2016-09-30"), by="days"),
                seq(as.Date("2019-05-01"), as.Date("2019-09-30"), by="days"),
                seq(as.Date("2020-05-01"), as.Date("2020-09-30"), by="days"),
                seq(as.Date("2021-05-01"), as.Date("2021-09-30"), by="days"))


#read in ctd
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/13/27ceda6bc7fdec2e7d79a6e4fe16ffdf" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

ctd <-read.csv(infile1,header=T) |> 
  mutate(DateTime = as.Date(DateTime)) |>
  filter(DateTime %in% dates_list & 
           Depth_m > 0 & Site ==50 & Reservoir %in% c("BVR")) |> 
  select(Reservoir, Site, DateTime, Depth_m, 
         Temp_C, DO_mgL, DOsat_percent) #removing chl a bc fp and ctd aren't so comparable

#round ctd to nearest m (need to see if I can use ysi bc some missing months w/ ctd)
ctd_final <- ctd |>
  dplyr::mutate(rdepth = plyr::round_any(Depth_m, 1)) |>
  dplyr::mutate(month = month(DateTime),
                year = year(DateTime)) |>
  dplyr::select(-Depth_m) |> 
  dplyr::rename(Depth_m = rdepth) |> 
  dplyr::group_by(month, year) |> 
  dplyr::filter(Depth_m %in% c(1,last(Depth_m))) |> 
  dplyr::summarise(Temp_C_epi = mean(Temp_C[Depth_m==1], na.rm=T),
                   Temp_C_hypo = mean(Temp_C[Depth_m!=1], na.rm=T),
                   DO_mgL_epi = mean(DO_mgL[Depth_m==1], na.rm=T),
                   DO_mgL_hypo = mean(DO_mgL[Depth_m!=1], na.rm=T))


#missing Aug 2015, May 2019, May 2020, Aug + Sep 2021
missing_dates <- c(seq(as.Date("2015-08-01"), as.Date("2015-08-31"), by="days"),
                   seq(as.Date("2019-05-01"), as.Date("2019-05-31"), by="days"),
                   seq(as.Date("2020-05-01"), as.Date("2020-05-31"), by="days"),
                   seq(as.Date("2021-08-01"), as.Date("2021-09-30"), by="days"))

#see of we can sub missing ctd data with ysi data
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/198/11/6e5a0344231de7fcebbe6dc2bed0a1c3" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))

ysi <- read.csv(infile1,header=T) |> 
  mutate(DateTime = as.Date(DateTime)) |>
  filter(DateTime %in% missing_dates & 
           Depth_m %in%  c(1,11) & 
           Site ==50 & Reservoir %in% c("BVR")) |> 
  select(Reservoir, Site, DateTime, Depth_m, 
         Temp_C, DO_mgL) |> 
  mutate(month = month(DateTime),
                year = year(DateTime)) |> 
  group_by(month, year) |>          
  summarise(Temp_C_epi = mean(Temp_C[Depth_m==1], na.rm=T),
            Temp_C_hypo = mean(Temp_C[Depth_m==11], na.rm=T),
            DO_mgL_epi = mean(DO_mgL[Depth_m==1], na.rm=T),
            DO_mgL_hypo = mean(DO_mgL[Depth_m==11], na.rm=T)) 
  
#combine ctd and ysi
profiles <- bind_rows(ctd_final, ysi) |> arrange(month, year) 

#select every 0.5m from casts for thermocline calc below
ctd_final_temp <- ctd |>
  dplyr::mutate(rdepth = plyr::round_any(Depth_m, 0.5)) |>
  dplyr::group_by(DateTime, rdepth, Reservoir, Site) |>
  dplyr::summarise(value = mean(Temp_C)) |>
  dplyr::rename(depth = rdepth) 

#calculate thermocline depth (missing some profiles...)
ctd_thermo_depth <- ctd_final_temp |> 
  group_by(DateTime) |> 
  summarise(therm_depth = thermo.depth(value,depth)) |> 
  mutate(month = month(DateTime),
                 year = year(DateTime)) |> 
  ungroup() |>
  group_by(month, year) |> 
  summarise(therm_depth = mean(therm_depth))

#read in chem from edi
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/199/11/509f39850b6f95628d10889d66885b76" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

chem <-read.csv(infile1,header=T) |> 
  dplyr::mutate(DateTime = as.Date(DateTime)) |>
  dplyr::filter(DateTime %in% dates_list & 
          Site ==50 & Reservoir %in% c("BVR")) |> 
  dplyr::mutate(month = month(DateTime),
                year = year(DateTime)) |> 
  dplyr::select(Reservoir, Depth_m, Rep,
         TN_ugL, TP_ugL,NH4_ugL ,NO3NO2_ugL, SRP_ugL, month, year) |> 
  dplyr::group_by(month, year) |> 
  dplyr::filter(Depth_m <=0.1 | Depth_m > 7) |>   #drop data when only one random depth was analyzed
  dplyr::filter(Depth_m %in%  c(0.1,max(Depth_m))) |> 
  dplyr::summarise(TN_ugL_epi = mean(TN_ugL[Depth_m==0.1], na.rm=T),
                   TN_ugL_hypo = mean(TN_ugL[Depth_m!=0.1], na.rm=T),
                   TP_ugL_epi = mean(TP_ugL[Depth_m==0.1], na.rm=T),
                   TP_ugL_hypo = mean(TP_ugL[Depth_m!=0.1], na.rm=T),
                   NH4_ugL_epi = mean(NH4_ugL[Depth_m==0.1], na.rm=T),
                   NH4_ugL_hypo = mean(NH4_ugL[Depth_m!=0.1], na.rm=T),
                   NO3NO2_ugL_epi = mean(NO3NO2_ugL[Depth_m==0.1], na.rm=T),
                   NO3NO2_ugL_hypo = mean(NO3NO2_ugL[Depth_m!=0.1], na.rm=T),
                   SRP_ugL_epi = mean(SRP_ugL[Depth_m==0.1], na.rm=T),
                   SRP_ugL_hypo = mean(SRP_ugL[Depth_m!=0.1], na.rm=T)) |> 
  dplyr::arrange(month, year)

#average the 10m tp from earlier in aug with the 9m tp from a later aug date so no NA (aug 2021)
chem$TP_ugL_hypo[is.nan(chem$TP_ugL_hypo)] <-  mean(c(26.4, 23.5))

#read in met data
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/389/7/02d36541de9088f2dd99d79dc3a7a853" 
#infile1 <- tempfile()
#try(download.file(inUrl1,infile1,method="curl"))
#if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")

#avg obs for each 24 hr period - no met station until jul 2015...
#met <-read.csv(infile1,header=T) |> 
#  dplyr::mutate(DateTime = as.POSIXct(DateTime, format="%Y-%m-%d %H:%M:%S")) |> 
#  dplyr::filter(DateTime %in% dates_list) |> 
#  dplyr::select(DateTime, AirTemp_C_Average, WindSpeed_Average_m_s,
#                            Rain_Total_mm, PAR_umolm2s_Average) |> 
#  dplyr::mutate(month = month(DateTime),
#                year = year(DateTime)) |> 
#  dplyr::group_by(month, year) |> 
#  dplyr::summarise(AirTemp_C = mean(AirTemp_C_Average),
#                   WindSpeed = mean(WindSpeed_Average_m_s),
#                   Rain_mm = mean(Rain_Total_mm),
#                   PAR_umolm2s = mean(PAR_umolm2s_Average))

#------------------------------------------------------------------------------#
#make an environmental driver df
env_drivers <- bind_cols(chem, profiles[!colnames(profiles) %in% c("month", "year")])
  
#export env csv
write.csv(env_drivers, "./Output/env.csv", row.names=FALSE)

#------------------------------------------------------------------------------#
#quick plots to visualize data
ggplot(env_drivers, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"), Temp_C_epi)) +
  geom_point() + geom_line() + theme_bw() + xlab("Date")

ggplot(env_drivers, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"), Temp_C_hypo)) +
  geom_point() + geom_line() + theme_bw() + xlab("Date")

ggplot(env_drivers, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"), DO_mgL_epi)) +
  geom_point() + geom_line() + theme_bw() + xlab("Date")

ggplot(env_drivers, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"), DO_mgL_hypo)) +
  geom_point() + geom_line() + theme_bw() + xlab("Date")

ggplot(env_drivers, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"), TN_ugL_epi)) +
  geom_point() + geom_line() + theme_bw() + xlab("Date")

ggplot(env_drivers, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"), TN_ugL_hypo)) +
  geom_point() + geom_line() + theme_bw() + xlab("Date")

ggplot(env_drivers, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"), TP_ugL_epi)) +
  geom_point() + geom_line() + theme_bw() + xlab("Date")

ggplot(env_drivers, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"), TP_ugL_hypo)) +
  geom_point() + geom_line() + theme_bw() + xlab("Date")

ggplot(env_drivers, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"), NH4_ugL_epi)) +
  geom_point() + geom_line() + theme_bw() + xlab("Date")

ggplot(env_drivers, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"), NH4_ugL_hypo)) +
  geom_point() + geom_line() + theme_bw() + xlab("Date")

ggplot(env_drivers, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"), NO3NO2_ugL_epi)) +
  geom_point() + geom_line() + theme_bw() + xlab("Date")

ggplot(env_drivers, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"), NO3NO2_ugL_hypo)) +
  geom_point() + geom_line() + theme_bw() + xlab("Date")

ggplot(env_drivers, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"), SRP_ugL_epi)) +
  geom_point() + geom_line() + theme_bw() + xlab("Date")

ggplot(env_drivers, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"), SRP_ugL_hypo)) +
  geom_point() + geom_line() + theme_bw() + xlab("Date")
                          

