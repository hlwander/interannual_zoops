#code to visualize seasonal and multiannual zoop patterns

#read in packages
pacman::p_load(tidyverse, NatParksPalettes)

#read in zoop data from EDI
inUrl1  <- "https://pasta-s.lternet.edu/package/data/eml/edi/1090/14/c7a04035b0a99adc489f5b6daec1cd52" 
infile1 <-  tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")

zoops <- read.csv(infile1, header=T) |>
  filter(CollectionMethod=="Tow" & Reservoir %in% c("BVR", "FCR") &
           StartDepth_m > 7.1) |> 
  select(-c(Site,EndDepth_m,CollectionMethod))

#split data into pre 2019 and post
zoops_2016_2018 <- zoops[as.Date(zoops$DateTime)<"2019-01-01",]
zoops_2019_2021 <- zoops[as.Date(zoops$DateTime)>="2019-01-01",]

#add daphnia (D. catawba, D. ambigua), calanoida (diaptomus) for 2014-2018 data
zoops_pre <- zoops_2016_2018 |> 
  group_by(Reservoir, DateTime, StartDepth_m) |> 
  summarise(Daphnia = sum(Density_IndPerL[
    Taxon %in% c("D. catawba","D. ambigua")]),
    Calanoida = sum(Density_IndPerL[
      Taxon %in% c("Diaptomus")]),
    Cyclopoida = sum(Density_IndPerL[
      Taxon %in% c("Cyclopoids")]),
    Nauplii = sum(Density_IndPerL[
      Taxon %in% c("Nauplii")]),
    Bosmina = sum(Density_IndPerL[
      Taxon %in% c("Bosmina")]),
    Ceriodaphnia = sum(Density_IndPerL[
      Taxon %in% c("Ceriodaphnia")]),
    Conochilus = sum(Density_IndPerL[
      Taxon %in% c("Conochilus")]),
    Keratella = sum(Density_IndPerL[
      Taxon %in% c("Keratella")]),
    Trichocerca = sum(Density_IndPerL[
      Taxon %in% c("Trichocerca")]),
    Kellicottia = sum(Density_IndPerL[
      Taxon %in% c("Kellicottia")]),
    Lecane = sum(Density_IndPerL[
      Taxon %in% c("Lecane")]),
    Rotifera = sum(Density_IndPerL[
      Taxon %in% c("Total Rotifers")]),
    Cladocera = sum(Density_IndPerL[
      Taxon %in% c("Bosmina","D. catawba", "Chydorus","D. ambigua",
                   "Diaphanosoma","Ceriodaphnia")]),
    Copepoda = sum(Density_IndPerL[
      Taxon %in% c("Diaptomus","Nauplii", "Cyclopoids")]))

#convert back to long
zoops_final_pre <- zoops_pre |> 
  group_by(Reservoir, DateTime, StartDepth_m) |> 
  pivot_longer(cols=Daphnia:Copepoda,
               names_to = c("Taxon"),
               values_to = "Density_IndPerL") |> 
  filter(hour(DateTime) %in% c(9,10,11,12,13,14)) |> #drop nighttime samples
  mutate(DateTime = as.Date(DateTime)) |> 
  mutate(Density_IndPerL = Density_IndPerL * (1/0.031))  #10m bvr neteff from 2016 (n=2) - note that 7m neteff was also 0.31
  #avg from 2020 and 2021 is 0.021...
 
#list common taxa between pre and post
taxa <- c("Bosmina", "Daphnia", "Ceriodaphnia",
          "Cyclopoida","Calanoida", "nauplius", 
          "Conochilidae","Keratella", "Rotifera",
          "Trichocercidae","Kellicottia", "Lecane",
          "Cladocera", "Copepoda")

#average reps when appropriate
zoops_final_post <- zoops_2019_2021 |> 
  mutate(DateTime = as.POSIXct(DateTime, format="%Y-%m-%d %H:%M:%S", tz="UTC")) |> 
  filter(hour(DateTime) %in% c(9,10,11,12,13,14)) |> #drop nighttime samples
  filter(Taxon %in% c(taxa)) |> 
  mutate(DateTime = as.Date(DateTime)) |> 
  mutate(Taxon = ifelse(Taxon=="nauplius", "Nauplii",
                 ifelse(Taxon=="Trichocercidae", "Trichocerca", 
                 ifelse(Taxon=="Conochilidae", "Conochilus", Taxon)))) |> 
  group_by(Reservoir, DateTime, StartDepth_m, Taxon) |> 
  summarise(Density_IndPerL = mean(Density_IndPerL))

#combine all zoop data
all_zoops <- bind_rows(zoops_final_pre, zoops_final_post) |> 
  mutate_all(~replace(., is.nan(.), NA)) |>  #replace NAN with NA
  ungroup() |> select(-StartDepth_m)   #dropping, but note that depths range from 7.5-11.5m....

#new df with just total density, avg by month + calculate sd
bvr_total_zoops <- all_zoops |> group_by(Reservoir, DateTime) |> 
  summarise(Total = sum(Density_IndPerL[Taxon %in% c("Cladocera","Copeoda","Rotifera")])) |> 
  ungroup() |> filter(Reservoir=="BVR") |> 
  mutate(year = format(DateTime, "%Y"),
         month = format(DateTime, "%m")) |> 
  group_by(year, month) |> 
  summarise(Total_avg = mean(Total),
            Total_sd = sd(Total)) |> ungroup() |> 
  filter(!year %in% "2022") #drop 2022

#add doy and year column
all_zoops$doy <- yday(all_zoops$DateTime)
all_zoops$year <- year(all_zoops$DateTime)

bvr_total_zoops$doy <- unique(all_zoops$DateTime[all_zoops$Reservoir=="BVR"])

#look at doy on x and year by color
ggplot(bvr_total_zoops, aes(as.Date(paste0("2023-",month,"-01")),
                            Total_avg, color=as.factor(year))) + 
  theme_bw() + xlab("Month") + ylab ("Zooplankton (#/L)") +
  theme(text = element_text(size=14), 
        axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.direction = "vertical",
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = c(0.88,0.9), legend.spacing = unit(-0.5, 'cm'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_point() +
  geom_line(size=1.2) + scale_x_date(date_breaks="2 month", date_labels="%b") +
  geom_errorbar(aes(ymin = Total_avg -Total_sd, ymax = Total_avg+Total_sd)) +
  scale_color_manual("",values=natparks.pals("Banff", 6)) 


#look at doy on x and year by color
ggplot(all_zoops, aes(doy, Density_IndPerL, color=as.factor(year))) + 
  geom_point() + theme_bw() + geom_line() +
  facet_wrap(~Taxon+Reservoir, scales="free_y", nrow=3) 

ggplot(data=subset(all_zoops, Reservoir=="FCR" & Taxon=="Cladocera"),
                   aes(doy, Density_IndPerL, color=as.factor(year))) + 
  geom_point() + theme_bw() + geom_line() +
  facet_wrap(~Taxon+Reservoir, scales="free_y", nrow=3) 

#not sure if I should look at all 11 taxa or just cladocerans, copepods, and rotifers

ggplot(data=subset(all_zoops, Taxon %in% c("Cladocera", "Copepoda", "Rotifera")),
                   aes(doy, Density_IndPerL, color=as.factor(year))) + 
  geom_point() + theme_bw() + geom_line() +
  facet_wrap(~Taxon+Reservoir, scales="free_y", nrow=3) 

#3 taxa full timeseries
ggplot(data=subset(all_zoops, Taxon %in% c("Cladocera", "Copepoda", "Rotifera") & Reservoir=="BVR"),
       aes(DateTime, Density_IndPerL, color=as.factor(Taxon))) + 
  geom_point() + theme_bw() + geom_line() +
  facet_wrap(~Taxon, scales="free_y", nrow=3) 




#------------------------------------------------------------------------------#
#Pull in environmental data to make sure we have data associated with each zoop sample

dates <- unique(all_zoops$DateTime)

#read in ctd
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/13/27ceda6bc7fdec2e7d79a6e4fe16ffdf" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

ctd <-read.csv(infile1,header=T) |> 
  mutate(DateTime = as.Date(DateTime)) |>
  filter(DateTime >= "2014-04-04" & DateTime <= "2022-07-01" & 
           Depth_m > 0 & Site ==50 & Reservoir %in% c("FCR", "BVR")) |> 
  select(Reservoir, Site, DateTime, Depth_m, 
         Temp_C, DO_mgL, DOsat_percent, Chla_ugL)

#round ctd to nearest m
ctd_final <- ctd |>
  dplyr::mutate(rdepth = plyr::round_any(Depth_m, 1)) |>
  dplyr::group_by(Reservoir, DateTime, rdepth) |>
  dplyr::summarise(Temp_C = mean(Temp_C),
                   DO_mgL = mean(DO_mgL),
                   DOsat_percent = mean(DOsat_percent),
                   Chla_ugL = mean(Chla_ugL)) |>
  dplyr::rename(Depth_m = rdepth) |> 
  dplyr::filter(Depth_m %in% ifelse(Reservoir=="BVR", c(1,10),
                                    c(1,9)))

#read in chem from edi
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/199/11/509f39850b6f95628d10889d66885b76" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

chem <-read.csv(infile1,header=T) |> 
  dplyr::mutate(DateTime = as.Date(DateTime)) |>
  dplyr::filter(DateTime >= "2014-04-04" & DateTime <= "2022-07-01" & 
          Site ==50 & Reservoir %in% c("FCR", "BVR")) |> 
  dplyr::select(Reservoir, DateTime, Depth_m, Rep,
         TN_ugL, TP_ugL,NH4_ugL ,NO3NO2_ugL, SRP_ugL) |> 
  dplyr::group_by(DateTime, Reservoir) |> 
  dplyr::filter(Depth_m %in% ifelse(Reservoir=="BVR", c(0.1,last(Depth_m[Reservoir=="BVR"])),
                                    c(0.1,last(Depth_m[Reservoir=="FCR"])))) |> 
  dplyr::filter(Depth_m <=0.1 | Depth_m > 7) #drop data when only one random depth was analyzed


#plot surface (1m) and bottom (9 or 10m) variables
ggplot(ctd_final, aes(DateTime, Temp_C, color=as.factor(Depth_m))) +
  geom_point() + geom_line() + theme_bw() + facet_wrap(~Reservoir, nrow=2)

ggplot(ctd_final, aes(DateTime, DO_mgL, color=as.factor(Depth_m))) +
  geom_point() + geom_line() + theme_bw() + facet_wrap(~Reservoir, nrow=2)

ggplot(ctd_final, aes(DateTime, DOsat_percent, color=as.factor(Depth_m))) +
  geom_point() + geom_line() + theme_bw() + facet_wrap(~Reservoir, nrow=2)

ggplot(ctd_final, aes(DateTime, Chla_ugL, color=as.factor(Depth_m))) +
  geom_point() + geom_line() + theme_bw() + 
  facet_wrap(~Reservoir, nrow=2, scales="free_y")

ggplot(chem, aes(DateTime, TN_ugL, color=as.factor(Depth_m))) +
  geom_point() + geom_line() + theme_bw() + 
  facet_wrap(~Reservoir, nrow=2, scales="free_y")

ggplot(chem, aes(DateTime, TP_ugL, color=as.factor(Depth_m))) +
  geom_point() + geom_line() + theme_bw() + 
  facet_wrap(~Reservoir, nrow=2, scales="free_y")

