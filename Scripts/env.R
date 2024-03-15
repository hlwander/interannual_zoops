#code to visualize seasonal and multiannual zoop patterns

#read in packages
pacman::p_load(tidyverse, NatParksPalettes, rLakeAnalyzer,
               EDIutils, xml2, gsheet)

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
         Temp_C, DO_mgL) #removing chl a bc fp and ctd aren't so comparable so can't fill in missing days

#round ctd to nearest m (also use ysi bc some missing months w/ ctd)
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
                   DO_mgL_epi = mean(DO_mgL[Depth_m==1], na.rm=T))


#missing Aug 2015, May 2019, May 2020, Aug + Sep 2021
missing_dates <- c(seq(as.Date("2015-08-01"), as.Date("2015-08-31"), by="days"),
                   seq(as.Date("2019-05-01"), as.Date("2019-05-31"), by="days"),
                   seq(as.Date("2020-05-01"), as.Date("2020-05-31"), by="days"),
                   seq(as.Date("2021-08-01"), as.Date("2021-09-30"), by="days"))

#sub missing ctd data with ysi data
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/198/11/6e5a0344231de7fcebbe6dc2bed0a1c3" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))

ysi <- read.csv(infile1,header=T) |> 
  mutate(DateTime = as.Date(DateTime)) |>
  filter(DateTime %in% missing_dates & 
           Depth_m > 0 &
           Site ==50 & Reservoir %in% c("BVR")) |> 
  select(Reservoir, Site, DateTime, Depth_m, 
         Temp_C, DO_mgL)

ysi_final <- ysi |> 
  filter(Depth_m %in%  c(1,11)) |> 
  mutate(month = month(DateTime),
                year = year(DateTime)) |> 
  group_by(month, year) |>          
  summarise(Temp_C_epi = mean(Temp_C[Depth_m==1], na.rm=T),
            Temp_C_hypo = mean(Temp_C[Depth_m==11], na.rm=T),
            DO_mgL_epi = mean(DO_mgL[Depth_m == 1], na.rm=T)) 
  
#combine ctd and ysi
profiles <- bind_rows(ctd_final, ysi_final) |> arrange(month, year) 

#select every 0.5m from casts for thermocline calc below
ctd_final_temp_do <- ctd |>
  dplyr::mutate(rdepth = plyr::round_any(Depth_m, 0.5)) |>
  dplyr::group_by(DateTime, rdepth, Reservoir, Site) |>
  dplyr::summarise(temp = mean(Temp_C),
                   DO = mean(DO_mgL)) |>
  dplyr::rename(depth = rdepth) 

#rename columns
ctd_final_temp_do <- ctd_final_temp_do |> 
  rename(Depth_m = depth,
         Temp_C = temp,
         DO_mgL = DO)

#clean up ysi
ysi_clean <- ysi |> 
  select(Reservoir,DateTime, Depth_m,Temp_C, DO_mgL) |> 
  filter(!is.na(Temp_C))

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

#calculate thermocline depth (missing n=5 profiles...) - bring in ysi profiles too??
ctd_thermo_depth <- temp_final |> 
  group_by(DateTime) |> 
  summarise(therm_depth = thermo.depth(Temp_C,Depth_m)) |> 
  mutate(month = month(DateTime),
                 year = year(DateTime)) |> 
  ungroup() |>
  group_by(month, year) |> 
  summarise(therm_depth = mean(therm_depth))

#calculate oxycline depth
ctd_oxy_depth <- do_final |> 
  group_by(DateTime) |> 
  summarise(oxy_depth = thermo.depth(DO_mgL,Depth_m)) |> 
  mutate(month = month(DateTime),
         year = year(DateTime)) |> 
  ungroup() |>
  group_by(month, year) |> 
  summarise(oxy_depth = mean(oxy_depth))

#calculate anoxic depth
anoxic_depth <- do_final |> 
  group_by(DateTime) |> 
  mutate(AD = first(Depth_m[DO_mgL <= 2])) |> 
  mutate(month = month(DateTime),
         year = year(DateTime)) |> 
  group_by(month, year) |> 
  summarise(anoxic_depth = mean(AD,na.rm=T))

#plot anoxic depth
ggplot(anoxic_depth, aes(yday(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d")),
                         anoxic_depth, color=as.factor(year))) +
  geom_point() + geom_line() + theme_bw() + xlab("doy") +
  scale_x_continuous(labels = scales::date_format("%b",tz="EST5EDT")) +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("Acadia", 6))
ggsave("Figures/anoxic_depth_vs_doy.jpg", width=6, height=3) 

#download bathymetry
infile <- tempfile()
try(download.file("https://pasta.lternet.edu/package/data/eml/edi/1254/1/f7fa2a06e1229ee75ea39eb586577184",
                  infile, method="curl"))
if (is.na(file.size(infile))) download.file(historic_file,infile,method="auto")

bathymetry <- readr::read_csv(infile, show_col_types = F)  |>
  dplyr::select(Reservoir, Depth_m, SA_m2) |>
  dplyr::filter(Reservoir == "BVR")

#read in water level (THIS LINK WILL NEED TO BE UPDATED ONCE QAQC IS DONE)
#gsheet_url <- 'https://docs.google.com/spreadsheets/d/1DDF-KZPuGBOjO2rB-owdod14N9bRs00XSZmmaASO7g8/edit#gid=0'

#using daily wl file for now

water_level <- read.csv("./Output/BVR_WaterLevel_2014_2022_interp.csv") |> 
  select(Date, WaterLevel_m) |> 
  mutate(WaterLevel_m = as.numeric(WaterLevel_m)) |> 
  filter(Date >= "2014-01-01" & Date < "2022-01-01") |> 
  mutate(month = month(Date),
         year = year(Date)) |> 
  group_by(month, year) |> 
  summarise(waterlevel = mean(WaterLevel_m,na.rm=T)) |> 
  filter(month %in% c(5,6,7,8,9),
         year %in% c(2014:2016,2019:2021)) |> 
  arrange(month,year) |> ungroup() |> 
  group_by(year) |> 
  mutate(diff = max(waterlevel) - min(waterlevel)) |> 
  ungroup()

ggplot(water_level, aes(yday(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d")),
                        waterlevel, color=as.factor(year))) +
  geom_point() + geom_line() + theme_bw() + xlab("doy") +
  scale_x_continuous(labels = scales::date_format("%b",tz="EST5EDT")) +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("Glacier", 6))
ggsave("Figures/waterlevel_vs_doy.jpg", width=6, height=3) 

wl <- read.csv("./Output/BVR_WaterLevel_2014_2022_interp.csv") |> 
  select(Date, WaterLevel_m) |> 
  mutate(WaterLevel_m = as.numeric(WaterLevel_m)) |> 
  mutate(Date = as.Date(Date)) |> 
  filter(Date >= "2014-01-01" & Date < "2022-01-01" &
           month(Date) %in% c(5,6,7,8,9))

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
  dplyr::rename(depth = Depth_m) |> 
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

#for loop for physcial metrics
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

#new wide df of temp for all dates
temp_final_wl_wide <- temp_final_wl |>
  pivot_wider(names_from = Depth_m, values_from = Temp_C,
              names_glue = "wtr_{Depth_m}") |> 
  select(-c(WaterLevel_m, Reservoir))

#calculate max buoyancy frequency 
bf <- rLakeAnalyzer::ts.buoyancy.freq(wtr = temp_final_wl_wide, 
                                          at.thermo = T, na.rm=T) 

schmidts[,3] <- bf$n2

#now calculate month/yearly means
physics <- schmidts |> 
  dplyr::mutate(month = month(Date),
                year = year(Date)) |> 
  dplyr::group_by(month, year) |> 
  dplyr::summarise(SS = mean(SS),
                   BF = mean(BF, na.rm=T))
  
#plot metrics over time
ggplot(physics, aes(yday(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d")),
                        SS, color=as.factor(year))) +
  geom_point() + geom_line() + theme_bw() + xlab("doy") +
  scale_x_continuous(labels = scales::date_format("%b",tz="EST5EDT")) +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("DeathValley", 6))

ggplot(physics, aes(yday(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d")),
                    BF, color=as.factor(year))) +
  geom_point() + geom_line() + theme_bw() + xlab("doy") +
  scale_x_continuous(labels = scales::date_format("%b",tz="EST5EDT")) +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("WindCave", 6))

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
         TN_ugL, TP_ugL, month, year) |> #NH4_ugL ,NO3NO2_ugL, SRP_ugL,
  dplyr::group_by(month, year) |> 
  dplyr::filter(Depth_m <=0.1 | Depth_m > 7) |>   #drop data when only one random depth was analyzed
  dplyr::filter(Depth_m %in%  c(0.1,max(Depth_m))) |> 
  dplyr::summarise(TN_ugL_epi = mean(TN_ugL[Depth_m==0.1], na.rm=T),
                   TN_ugL_hypo = mean(TN_ugL[Depth_m!=0.1], na.rm=T),
                   TP_ugL_epi = mean(TP_ugL[Depth_m==0.1], na.rm=T),
                   TP_ugL_hypo = mean(TP_ugL[Depth_m!=0.1], na.rm=T)) |> #,
                   #NH4_ugL_epi = mean(NH4_ugL[Depth_m==0.1], na.rm=T),
                   #NH4_ugL_hypo = mean(NH4_ugL[Depth_m!=0.1], na.rm=T),
                   #NO3NO2_ugL_epi = mean(NO3NO2_ugL[Depth_m==0.1], na.rm=T),
                   #NO3NO2_ugL_hypo = mean(NO3NO2_ugL[Depth_m!=0.1], na.rm=T),
                   #SRP_ugL_epi = mean(SRP_ugL[Depth_m==0.1], na.rm=T),
                   #SRP_ugL_hypo = mean(SRP_ugL[Depth_m!=0.1], na.rm=T)) |> 
  dplyr::arrange(month, year)

#average the 10m tp from earlier in aug with the 9m tp from a later aug date so no NA (aug 2021)
chem$TP_ugL_hypo[is.nan(chem$TP_ugL_hypo)] <-  mean(c(26.4, 23.5))

#convert fp from wide to long
chem_long <- chem |>
  pivot_longer(cols = c(TN_ugL_epi,TP_ugL_epi), 
               names_to = "variable")  

#plot chem over time
ggplot(data=subset(chem_long, !variable %in% c("TN_ugL_epi")), 
       aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"), value,
                    color = variable)) +
  geom_vline(xintercept = as.Date("2014-01-01")) +
  geom_vline(xintercept = as.Date("2015-01-01")) +
  geom_vline(xintercept = as.Date("2016-01-01")) +
  geom_vline(xintercept = as.Date("2017-01-01")) +
  geom_vline(xintercept = as.Date("2019-01-01")) +
  geom_vline(xintercept = as.Date("2020-01-01")) +
  geom_vline(xintercept = as.Date("2021-01-01")) +
  geom_point() + geom_line() + theme_bw() + xlab("Date") +
  annotate("text", x=as.Date("2014-07-01"), y=23, label= "2014") +
  annotate("text", x=as.Date("2015-07-01"), y=23, label= "2015") +
  annotate("text", x=as.Date("2016-07-01"), y=23, label= "2016") +
  annotate("text", x=as.Date("2019-07-01"), y=23, label= "2019") +
  annotate("text", x=as.Date("2020-07-01"), y=23, label= "2020") +
  annotate("text", x=as.Date("2021-07-01"), y=23, label= "2021") +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("Banff", 4), 
                     labels=c("NH4_ugL","NO3NO2_ugL","SRP_ugL","TP_ugL"))
ggsave("Figures/chem_timeseries_notn.jpg", width=6, height=3) 

#and now tn
ggplot(data=subset(chem_long, variable %in% c("TN_ugL_epi")), 
       aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"), value,
           color = variable)) +
  geom_vline(xintercept = as.Date("2014-01-01")) +
  geom_vline(xintercept = as.Date("2015-01-01")) +
  geom_vline(xintercept = as.Date("2016-01-01")) +
  geom_vline(xintercept = as.Date("2017-01-01")) +
  geom_vline(xintercept = as.Date("2019-01-01")) +
  geom_vline(xintercept = as.Date("2020-01-01")) +
  geom_vline(xintercept = as.Date("2021-01-01")) +
  geom_point() + geom_line() + theme_bw() + xlab("Date") +
  annotate("text", x=as.Date("2014-07-01"), y=400, label= "2014") +
  annotate("text", x=as.Date("2015-07-01"), y=400, label= "2015") +
  annotate("text", x=as.Date("2016-07-01"), y=400, label= "2016") +
  annotate("text", x=as.Date("2019-07-01"), y=400, label= "2019") +
  annotate("text", x=as.Date("2020-07-01"), y=400, label= "2020") +
  annotate("text", x=as.Date("2021-07-01"), y=400, label= "2021") 
ggsave("Figures/chem_timeseries_tn.jpg", width=6, height=3) 

#read in secchi data
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/198/11/81f396b3e910d3359907b7264e689052" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))

secchi <-read.csv(infile1) |> 
  mutate(DateTime = as.Date(DateTime)) |> 
  filter(Reservoir == "BVR" & Site == 50 &
           DateTime %in% dates_list) |> 
  distinct() |> 
  mutate(year = format(DateTime, "%Y"),
         month = format(DateTime, "%m")) |> 
  group_by(year, month) |> 
  summarise(secchi = mean(Secchi_m))

ggplot(secchi, aes(as.Date("2019-12-31") + 
                     yday(as.Date(paste0(year,"-",month,"-01"), 
                                  "%Y-%m-%d")), secchi, color=year)) + 
  geom_line() + geom_point() + theme_bw() + xlab("") +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("Cuyahoga", 6)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") 
ggsave("Figures/secchi_vs_doy.jpg", width=6, height=4)
#there is a 23oct2015 3m secchi obs for bvr, but nothing in sep...

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
env_drivers <- bind_cols(chem, anoxic_depth[!colnames(anoxic_depth) %in% 
                                              c("month", "year")],
                         profiles[!colnames(profiles) %in% c("month", "year")],
                         water_level[!colnames(water_level) %in% 
                                       c("month", "year")],
                         ctd_thermo_depth[!colnames(ctd_thermo_depth) %in%
                                            c("month","year")],
                         ctd_oxy_depth[!colnames(ctd_oxy_depth) %in%
                                         c("month","year")],
                         physics[!colnames(physics) %in% 
                                   c("month", "year")])#,
                        # secchi[!colnames(secchi) %in% 
                        #          c("month", "year")])

#export env csv
write.csv(env_drivers, "./Output/env.csv", row.names=FALSE)

#------------------------------------------------------------------------------#
#read in fp data
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/272/8/0359840d24028e6522f8998bd41b544e" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))

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
  dplyr::group_by(month, year) |> 
  dplyr::filter(Depth_m %in% c(1)) |> #only doing epi bc hypo isn't so useful
  dplyr::summarise(Green_ugL = mean(GreenAlgae_ugL, na.rm=T),
                   Bluegreen_ugL = mean(Bluegreens_ugL, na.rm=T),
                   Brown_ugL = mean(BrownAlgae_ugL, na.rm=T),
                   Mixed_ugL = mean(MixedAlgae_ugL, na.rm=T),
                   Total_ugL = mean(TotalConc_ugL, na.rm=T))

#save phyto df
write.csv(fp, "Output/phytos.csv", row.names=FALSE)

#convert fp from wide to long
fp_long <- fp |>
  pivot_longer(cols = Green_ugL:Total_ugL, 
               names_to = "variable")  

#plot phytos over time
ggplot(fp_long, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"), value,
                    color = variable)) +
  geom_vline(xintercept = as.Date("2014-01-01")) +
  geom_vline(xintercept = as.Date("2015-01-01")) +
  geom_vline(xintercept = as.Date("2016-01-01")) +
  geom_vline(xintercept = as.Date("2017-01-01")) +
  geom_vline(xintercept = as.Date("2019-01-01")) +
  geom_vline(xintercept = as.Date("2020-01-01")) +
  geom_vline(xintercept = as.Date("2021-01-01")) +
  geom_point() + geom_line() + theme_bw() + xlab("Date") +
  annotate("text", x=as.Date("2014-07-01"), y=23, label= "2014") +
  annotate("text", x=as.Date("2015-07-01"), y=23, label= "2015") +
  annotate("text", x=as.Date("2016-07-01"), y=23, label= "2016") +
  annotate("text", x=as.Date("2019-07-01"), y=23, label= "2019") +
  annotate("text", x=as.Date("2020-07-01"), y=23, label= "2020") +
  annotate("text", x=as.Date("2021-07-01"), y=23, label= "2021") +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("SmokyMtns", 6), 
                     labels=c("Bluegreen_ugL","Brown_ugL","Green_ugL",
                              "Mixed_ugL","Total_ugL"))
ggsave("Figures/phyto_succession.jpg", width=6, height=3) 

#remove green, brown, and total to see what other groups are doing
ggplot(data=subset(fp_long, variable %in% c("Mixed_ugL","Bluegreen_ugL")),
       aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"), value,
                    color = variable)) +
  geom_vline(xintercept = as.Date("2014-01-01")) +
  geom_vline(xintercept = as.Date("2015-01-01")) +
  geom_vline(xintercept = as.Date("2016-01-01")) +
  geom_vline(xintercept = as.Date("2017-01-01")) +
  geom_vline(xintercept = as.Date("2019-01-01")) +
  geom_vline(xintercept = as.Date("2020-01-01")) +
  geom_vline(xintercept = as.Date("2021-01-01")) +
  geom_point() + geom_line() + theme_bw() + xlab("Date") +
  annotate("text", x=as.Date("2014-07-01"), y=6.1, label= "2014") +
  annotate("text", x=as.Date("2015-07-01"), y=6.1, label= "2015") +
  annotate("text", x=as.Date("2016-07-01"), y=6.1, label= "2016") +
  annotate("text", x=as.Date("2019-07-01"), y=6.1, label= "2019") +
  annotate("text", x=as.Date("2020-07-01"), y=6.1, label= "2020") +
  annotate("text", x=as.Date("2021-07-01"), y=6.1, label= "2021") +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("BryceCanyon",3), 
                     labels=c("Bluegreen_ugL","Mixed_ugL","Yellow_ugL"))
ggsave("Figures/phyto_succession_no_brown_or_greens.jpg", width=6, height=3) 
  
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
                          
#-----------------------------------------------------------------------------#
#standardize for each year and variable (n=30 1s)
phytos_std <- fp_long |> 
  filter(!variable %in% c("Total_ugL")) |> 
  group_by(year, month) |> 
  summarise(Green = mean(value[variable=="Green_ugL"]),
            Bluegreen = mean(value[variable=="Bluegreen_ugL"]),
            Brown = mean(value[variable=="Brown_ugL"]),
            Mixed = mean(value[variable=="Mixed_ugL"])) |> 
  pivot_longer(-c(year,month),
               names_to = c("variable"))  |> 
  ungroup() |> group_by(variable,year) |>
  mutate(min_val = min(value),
         max_val = max(value)) |> 
  mutate(standardized_abund = (value - min_val) / (max_val - min_val))

#playing around with order/layering of taxa
phytos_std$variable <- factor(phytos_std$variable,
                                     levels = c("Mixed","Bluegreen","Brown","Green"))
                                     #levels = c("Green","Brown","Bluegreen","Mixed","Yellow"))

#phyto succession shaded line plot
ggplot(phytos_std, 
       aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"),
           standardized_abund, color=variable)) +
  geom_area(aes(color = variable, fill = variable),
            position = "stack", stat="identity",
            alpha=0.7) +
  facet_wrap(~year, scales = "free")+
  #scale_color_manual(values = NatParksPalettes::natparks.pals("Volcanoes", 4))+
  #scale_fill_manual(values = NatParksPalettes::natparks.pals("Volcanoes", 4)) +
  scale_color_manual(values = c("#9C1D4B","#147B80","#8A5414","#6EB579"))+
  scale_fill_manual(values = c("#9C1D4B","#147B80","#8A5414","#6EB579")) +
                    #breaks = c("Cladocera","Copepoda","Rotifera"))+
  scale_x_date(expand = c(0,0),
               labels = scales::date_format("%b",tz="EST5EDT")) +
  scale_y_continuous(expand = c(0,0))+
  xlab("") + ylab("standardized density") +
  guides(color= "none",
         fill = guide_legend(ncol=1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        text = element_text(size=8), 
        axis.text.y = element_text(size = 8),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.x = element_text(face = "bold",hjust = 0),
        axis.text.x = element_text(angle=90),
        strip.background.x = element_blank(),
        axis.title.y = element_text(size = 9),
        plot.margin = unit(c(0, 1, 0, 0), "cm"),
        panel.background = element_rect(
          fill = "white"),
        panel.spacing = unit(0.5, "lines"))
ggsave("Figures/BVR_phyto_succession_stacked.jpg", width=6, height=4) 

#standardize for each year and variable (n=30 1s)
chem_std <- chem_long |> 
  group_by(year, month) |> 
  summarise(TN = mean(value[variable=="TN_ugL_epi"]),
            TP = mean(value[variable=="TP_ugL_epi"])) |> 
  pivot_longer(-c(year,month),
               names_to = c("variable"))  |> 
  ungroup() |> group_by(variable,year) |>
  mutate(min_val = min(value),
         max_val = max(value)) |> 
  mutate(standardized_value = (value - min_val) / (max_val - min_val))

#phyto succession shaded line plot
ggplot(chem_std, 
       aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"),
           standardized_value, color=variable)) +
  geom_area(aes(color = variable, fill = variable),
            position = "stack", stat="identity") +
  facet_wrap(~year, scales = "free_x")+
  scale_color_manual(values = NatParksPalettes::natparks.pals("DeathValley", 2))+
  scale_fill_manual(values = NatParksPalettes::natparks.pals("DeathValley", 2)) +
  #breaks = c("Cladocera","Copepoda","Rotifera"))+
  scale_x_date(expand = c(0,0),
               labels = scales::date_format("%b",tz="EST5EDT")) +
  scale_y_continuous(expand = c(0,0))+
  xlab("") + ylab("standardized density") +
  guides(color= "none",
         fill = guide_legend(ncol=1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        text = element_text(size=8), 
        axis.text.y = element_text(size = 8),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.x = element_text(face = "bold",hjust = 0),
        axis.text.x = element_text(angle=90),
        strip.background.x = element_blank(),
        axis.title.y = element_text(size = 9),
        plot.margin = unit(c(0, 1, 0, 0), "cm"),
        panel.background = element_rect(
          fill = "white"),
        panel.spacing = unit(0.5, "lines"))
ggsave("Figures/BVR_total_np_succession.jpg", width=6, height=4) 

