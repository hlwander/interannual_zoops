#zoop biomass succession figs

#read in packages
pacman::p_load(zoo, dplR, dplyr, tidyverse, ggplot2, ggpubr, sf, lubridate)

#read in zoop data from EDI
inUrl1  <-  "https://pasta.lternet.edu/package/data/eml/edi/197/3/9eb6db370194bd3b2824726d89a008a6" 
infile1 <-  tempfile()
download.file(inUrl1,infile1,method="curl")

zoops <- read.csv(infile1, header=T) |>
  filter(CollectionMethod=="Tow" & Reservoir %in% c("BVR") &
           StartDepth_m > 7.1) |> 
  select(-c(Site,EndDepth_m,CollectionMethod))

#split data into pre 2019 and post
zoops_2016_2018 <- zoops[as.Date(zoops$DateTime)<"2019-01-01",]
zoops_2019_2021 <- zoops[as.Date(zoops$DateTime)>="2019-01-01",]

#add daphnia (D. catawba, D. ambigua), calanoida (diaptomus) for 2014-2018 data
zoops_pre <- zoops_2016_2018 |> 
  group_by(Reservoir, DateTime, StartDepth_m) |> 
  summarise(Daphnia_biom = sum(Biomass_ugL[
    Taxon %in% c("D. catawba","D. ambigua")],na.rm=T),
    Daphnia_sd = sd(Biomass_ugL[
      Taxon %in% c("D. catawba","D. ambigua")],na.rm=T),
    Calanoida_biom = sum(Biomass_ugL[
      Taxon %in% c("Diaptomus")],na.rm=T),
    Calanoida_sd = sd(Biomass_ugL[
      Taxon %in% c("Diaptomus")],na.rm=T),
    Cyclopoida_biom = sum(Biomass_ugL[
      Taxon %in% c("Cyclopoids")],na.rm=T),
    Cyclopoida_sd = sd(Biomass_ugL[
      Taxon %in% c("Cyclopoids")],na.rm=T),
    Nauplii_biom = sum(Biomass_ugL[
      Taxon %in% c("Nauplii")],na.rm=T),
    Nauplii_sd = sd(Biomass_ugL[
      Taxon %in% c("Nauplii")],na.rm=T),
    Bosmina_biom = sum(Biomass_ugL[
      Taxon %in% c("Bosmina")],na.rm=T),
    Bosmina_sd = sd(Biomass_ugL[
      Taxon %in% c("Bosmina")],na.rm=T),
    #Chydorus = sum(Biomass_ugL[
    #  Taxon %in% c("Chydorus")]),
    Ceriodaphnia_biom = sum(Biomass_ugL[
      Taxon %in% c("Ceriodaphnia")],na.rm=T),
    Ceriodaphnia_sd = sd(Biomass_ugL[
      Taxon %in% c("Ceriodaphnia")],na.rm=T),
    #Diaphanosoma = sum(Biomass_ugL[
    #  Taxon %in% c("Diaphanosoma")]),
    Ascomorpha_biom = sum(Biomass_ugL[
      Taxon %in% c("Ascomorpha")],na.rm=T),
    Ascomorpha_sd = sd(Biomass_ugL[
      Taxon %in% c("Ascomorpha")],na.rm=T),
    #Asplanchna = sum(Biomass_ugL[
    #  Taxon %in% c("Asplanchna")]),
    Conochilus_biom = sum(Biomass_ugL[
      Taxon %in% c("Conochilus")],na.rm=T),
    Conochilus_sd = sd(Biomass_ugL[
      Taxon %in% c("Conochilus")],na.rm=T),
    Keratella_biom = sum(Biomass_ugL[
      Taxon %in% c("Keratella")],na.rm=T),
    Keratella_sd = sd(Biomass_ugL[
      Taxon %in% c("Keratella")],na.rm=T),
    Trichocerca_biom = sum(Biomass_ugL[
      Taxon %in% c("Trichocerca")],na.rm=T),
    Trichocerca_sd = sd(Biomass_ugL[
      Taxon %in% c("Trichocerca")],na.rm=T),
    Kellicottia_biom = sum(Biomass_ugL[
      Taxon %in% c("Kellicottia")],na.rm=T),
    Kellicottia_sd = sd(Biomass_ugL[
      Taxon %in% c("Kellicottia")],na.rm=T),
    # Lecane = sum(Biomass_ugL[
    #   Taxon %in% c("Lecane")]),
    Polyarthra_biom = sum(Biomass_ugL[
      Taxon %in% c("Polyarthra")],na.rm=T),
    Polyarthra_sd = sd(Biomass_ugL[
      Taxon %in% c("Polyarthra")],na.rm=T),
    Rotifera_biom = sum(Biomass_ugL[
      Taxon %in% c("Total rotifers")],na.rm=T),
    Rotifera_sd = sd(Biomass_ugL[
      Taxon %in% c("Total rotifers")],na.rm=T),
    Cladocera_biom = sum(Biomass_ugL[
      Taxon %in% c("Bosmina","D. catawba", "Chydorus","D. ambigua",
                   "Diaphanosoma","Ceriodaphnia")],na.rm=T),
    Cladocera_sd = sd(Biomass_ugL[
      Taxon %in% c("Bosmina","D. catawba", "Chydorus","D. ambigua",
                   "Diaphanosoma","Ceriodaphnia")],na.rm=T),
    Copepoda_biom = sum(Biomass_ugL[
      Taxon %in% c("Diaptomus","Nauplii", "Cyclopoids")],na.rm=T),
    Copepoda_sd = sd(Biomass_ugL[
      Taxon %in% c("Diaptomus","Nauplii", "Cyclopoids")],na.rm=T)) |> 
  pivot_longer(-c(Reservoir,DateTime,StartDepth_m),
               names_to = c("Taxon", ".value"),
               names_sep="_" )  |> 
  filter(hour(DateTime) %in% c(9,10,11,12,13,14)) |> #drop nighttime samples
  mutate(DateTime = as.Date(DateTime)) |> 
  mutate(biom = biom * (1/0.031))  #10m bvr neteff from 2016 (n=2) - note that 7m neteff was also 0.031
#avg from 2020 and 2021 is 0.021...

#list common taxa between pre and post
taxa <- c("Bosmina", "Daphnia", "Ceriodaphnia",
          "Cyclopoida","Calanoida", "Nauplius", 
          "Conochilus","Keratella", "Rotifera",
          "Trichocerca","Kellicottia","Ascomorpha",
          "Polyarthra","Cladocera", "Copepoda") 
#Diaphanosoma, Lecane, Asplanchna, Chydorus

#average reps when appropriate
zoops_post <- zoops_2019_2021 |> 
  mutate(DateTime = as.POSIXct(DateTime, format="%Y-%m-%d %H:%M:%S", tz="UTC")) |> 
  filter(hour(DateTime) %in% c(9,10,11,12,13,14) &  #drop nighttime samples
           !year(DateTime) %in% c(2022)) |> 
  filter(Taxon %in% c(taxa)) |> 
  mutate(DateTime = as.Date(DateTime)) |> 
  mutate(Taxon = ifelse(Taxon == "Nauplius", "Nauplii", Taxon)) |> 
  group_by(Reservoir, DateTime, StartDepth_m, Taxon) |> 
  summarise(biom = mean(Biomass_ugL,na.rm=T),
            sd = sd(Biomass_ugL,na.rm=T))

#combine all zoop data
all_zoops <- bind_rows(zoops_pre, zoops_post) |> 
  mutate_all(~replace(., is.nan(.), NA)) |>  #replace NAN with NA
  ungroup() |> select(-StartDepth_m) #dropping, but note that depths range from 7.5-11.5m....

#write all_zoops
write.csv(all_zoops, paste0("Output/all_zoops_biom.csv"),row.names = FALSE)

#add column for pre vs post
all_zoops$data <- ifelse(all_zoops$DateTime<="2019-01-01","pre","post")

#------------------------------------------------------------------------------#
# new df for total zoop biom
zoops_total <- all_zoops |> 
  group_by(DateTime) |> 
  summarise(Total = sum(biom[Taxon %in% c("Cladocera","Copeoda","Rotifera")],na.rm=T),
            sd = mean(sd[Taxon %in% c("Cladocera","Copeoda","Rotifera")],na.rm=T)) |> 
  ungroup() |> 
  mutate(year = format(DateTime, "%Y"),
         month = format(DateTime, "%m")) |> 
  group_by(year, month) |>
  summarise(Total_avg = mean(Total,na.rm=T),
            Total_sd = mean(sd,na.rm=T)) |> 
  ungroup() |> 
  mutate(Biom_std =(Total_avg - min(Total_avg,na.rm=T)) / 
           (max(Total_avg,na.rm=T) - min(Total_avg,na.rm=T)))

# 3 group zoop biom
zoops_3_groups <- all_zoops |> group_by(Taxon, DateTime) |> 
  filter(Taxon %in% c("Cladocera","Copepoda","Rotifera")) |> 
  mutate(year = format(DateTime, "%Y"),
         month = format(DateTime, "%m")) |> 
  ungroup() |> group_by(year, month) |> 
  summarise(Cladocera_avg = mean(biom[Taxon=="Cladocera"],na.rm=T),
            Cladocera_sd = mean(sd[Taxon=="Cladocera"],na.rm=T),
            Copepoda_avg = mean(biom[Taxon=="Copepoda"],na.rm=T),
            Copepoda_sd = mean(sd[Taxon=="Copepoda"],na.rm=T),
            Rotifera_avg = mean(biom[Taxon=="Rotifera"],na.rm=T),
            Rotifera_sd = mean(sd[Taxon=="Rotifera"],na.rm=T)) |> 
  pivot_longer(-c(year,month),
               names_to = c("Taxon", ".value"),
               names_sep="_" )  |> 
  ungroup() |> group_by(Taxon) |>
  mutate(min_biom = min(avg,na.rm=T),
         max_biom = max(avg,na.rm=T)) |> 
  mutate(standardized_biom = (avg - min_biom) / (max_biom - min_biom))

# 12 group zoop biom
zoops_12_groups <- all_zoops |> group_by(Taxon, DateTime) |> 
  filter(!Taxon %in% c("Cladocera","Copepoda","Rotifera")) |> 
  mutate(year = format(DateTime, "%Y"),
         month = format(DateTime, "%m")) |> 
  group_by(year, month) |> 
  summarise(Daphnia_avg = mean(biom[Taxon=="Daphnia"],na.rm=T),
            Daphnia_sd = mean(sd[Taxon=="Daphnia"],na.rm=T),
            Calanoida_avg = mean(biom[Taxon=="Calanoida"],na.rm=T),
            Calanoida_sd = mean(sd[Taxon=="Calanoida"],na.rm=T),
            Cyclopoida_avg = mean(biom[Taxon=="Cyclopoida"],na.rm=T),
            Cyclopoida_sd = mean(sd[Taxon=="Cyclopoida"],na.rm=T),
            Nauplii_avg = mean(biom[Taxon=="Nauplii"],na.rm=T),
            Nauplii_sd = mean(sd[Taxon=="Nauplii"],na.rm=T),
            Bosmina_avg = mean(biom[Taxon=="Bosmina"],na.rm=T),
            Bosmina_sd = mean(sd[Taxon=="Bosmina"],na.rm=T),
            Ceriodaphnia_avg = mean(biom[Taxon=="Ceriodaphnia"],na.rm=T),
            Ceriodaphnia_sd = mean(sd[Taxon=="Ceriodaphnia"],na.rm=T),
            Ascomorpha_avg = mean(biom[Taxon=="Ascomorpha"],na.rm=T),
            Ascomorpha_sd = mean(sd[Taxon=="Ascomorpha"],na.rm=T),
            Conochilus_avg = mean(biom[Taxon=="Conochilus"],na.rm=T),
            Conochilus_sd = mean(sd[Taxon=="Conochilus"],na.rm=T),
            Keratella_avg = mean(biom[Taxon=="Keratella"],na.rm=T),
            Keratella_sd = mean(sd[Taxon=="Keratella"],na.rm=T),
            Trichocerca_avg = mean(biom[Taxon=="Trichocerca"],na.rm=T),
            Trichocerca_sd = mean(sd[Taxon=="Trichocerca"],na.rm=T),
            Kellicottia_avg = mean(biom[Taxon=="Kellicottia"],na.rm=T),
            Kellicottia_sd = mean(sd[Taxon=="Kellicottia"],na.rm=T),
            Polyarthra_avg = mean(biom[Taxon=="Polyarthra"],na.rm=T),
            Polyarthra_sd = mean(sd[Taxon=="Polyarthra"],na.rm=T)) |> 
  pivot_longer(-c(year,month),
               names_to = c("Taxon", ".value"),
               names_sep="_" )  |> 
  ungroup() |> group_by(Taxon) |>
  mutate(min_biom = min(avg,na.rm=T),
         max_biom = max(avg,na.rm=T)) |> 
  mutate(standardized_biom = (avg - min_biom) / (max_biom - min_biom))

#------------------------------------------------------------------------------#
#plots
ggplot(zoops_total, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"),Total_avg)) +
  geom_vline(xintercept = as.Date("2014-01-01")) +
  geom_vline(xintercept = as.Date("2015-01-01")) +
  geom_vline(xintercept = as.Date("2016-01-01")) +
  geom_vline(xintercept = as.Date("2017-01-01")) +
  geom_vline(xintercept = as.Date("2019-01-01")) +
  geom_vline(xintercept = as.Date("2020-01-01")) +
  geom_vline(xintercept = as.Date("2021-01-01")) +
  geom_point() + geom_line() + theme_bw() + xlab("Date") +
  annotate("text", x=as.Date("2014-07-01"), y=2000, label= "2014") +
  annotate("text", x=as.Date("2015-07-01"), y=2000, label= "2015") +
  annotate("text", x=as.Date("2016-07-01"), y=2000, label= "2016") +
  annotate("text", x=as.Date("2019-07-01"), y=2000, label= "2019") +
  annotate("text", x=as.Date("2020-07-01"), y=2000, label= "2020") +
  annotate("text", x=as.Date("2021-07-01"), y=2000, label= "2021") +
  geom_errorbar(aes(ymin = Total_avg -Total_sd, ymax = Total_avg+Total_sd), 
                color="#009999") 
ggsave("Figures/zoop_total_biom.jpg", width=6, height=3) 

ggplot(zoops_3_groups, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"),
                           avg, color=Taxon)) +
  geom_vline(xintercept = as.Date("2014-01-01")) +
  geom_vline(xintercept = as.Date("2015-01-01")) +
  geom_vline(xintercept = as.Date("2016-01-01")) +
  geom_vline(xintercept = as.Date("2017-01-01")) +
  geom_vline(xintercept = as.Date("2019-01-01")) +
  geom_vline(xintercept = as.Date("2020-01-01")) +
  geom_vline(xintercept = as.Date("2021-01-01")) +
  geom_point() + geom_line() + theme_bw() + xlab("Date") +
  annotate("text", x=as.Date("2014-07-01"), y=5000, label= "2014") +
  annotate("text", x=as.Date("2015-07-01"), y=5000, label= "2015") +
  annotate("text", x=as.Date("2016-07-01"), y=5000, label= "2016") +
  annotate("text", x=as.Date("2019-07-01"), y=5000, label= "2019") +
  annotate("text", x=as.Date("2020-07-01"), y=5000, label= "2020") +
  annotate("text", x=as.Date("2021-07-01"), y=5000, label= "2021") +
  geom_errorbar(aes(ymin = avg -sd, ymax = avg+sd)) +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("Saguaro", 3))
ggsave("Figures/zoop_3taxa_biom.jpg", width=6, height=3) 

ggplot(zoops_12_groups, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"),
                            avg, color=Taxon)) +
  geom_vline(xintercept = as.Date("2014-01-01")) +
  geom_vline(xintercept = as.Date("2015-01-01")) +
  geom_vline(xintercept = as.Date("2016-01-01")) +
  geom_vline(xintercept = as.Date("2017-01-01")) +
  geom_vline(xintercept = as.Date("2019-01-01")) +
  geom_vline(xintercept = as.Date("2020-01-01")) +
  geom_vline(xintercept = as.Date("2021-01-01")) +
  geom_point() + geom_line() + theme_bw() + xlab("Date") +
  annotate("text", x=as.Date("2014-07-01"), y=5000, label= "2014") +
  annotate("text", x=as.Date("2015-07-01"), y=5000, label= "2015") +
  annotate("text", x=as.Date("2016-07-01"), y=5000, label= "2016") +
  annotate("text", x=as.Date("2019-07-01"), y=5000, label= "2019") +
  annotate("text", x=as.Date("2020-07-01"), y=5000, label= "2020") +
  annotate("text", x=as.Date("2021-07-01"), y=5000, label= "2021") +
  geom_errorbar(aes(ymin = avg -sd, ymax = avg+sd)) +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("Torres", 12))
ggsave("Figures/zoop_12taxa_biom.jpg", width=6, height=3) 

#------------------------------------------------------------------------------#
#standardized biomass
ggplot(zoops_total, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"),Biom_std)) +
  geom_vline(xintercept = as.Date("2014-01-01")) +
  geom_vline(xintercept = as.Date("2015-01-01")) +
  geom_vline(xintercept = as.Date("2016-01-01")) +
  geom_vline(xintercept = as.Date("2017-01-01")) +
  geom_vline(xintercept = as.Date("2019-01-01")) +
  geom_vline(xintercept = as.Date("2020-01-01")) +
  geom_vline(xintercept = as.Date("2021-01-01")) +
  geom_point(color="darkblue") + geom_line(color="darkblue") + theme_bw() + xlab("Date") +
  annotate("text", x=as.Date("2014-07-01"), y=1.1, label= "2014") +
  annotate("text", x=as.Date("2015-07-01"), y=1.1, label= "2015") +
  annotate("text", x=as.Date("2016-07-01"), y=1.1, label= "2016") +
  annotate("text", x=as.Date("2019-07-01"), y=1.1, label= "2019") +
  annotate("text", x=as.Date("2020-07-01"), y=1.1, label= "2020") +
  annotate("text", x=as.Date("2021-07-01"), y=1.1, label= "2021") 
ggsave("Figures/zoop_total_std_biom.jpg", width=6, height=3) 

ggplot(zoops_3_groups, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"),
                           standardized_biom, color=Taxon)) +
  geom_vline(xintercept = as.Date("2014-01-01")) +
  geom_vline(xintercept = as.Date("2015-01-01")) +
  geom_vline(xintercept = as.Date("2016-01-01")) +
  geom_vline(xintercept = as.Date("2017-01-01")) +
  geom_vline(xintercept = as.Date("2019-01-01")) +
  geom_vline(xintercept = as.Date("2020-01-01")) +
  geom_vline(xintercept = as.Date("2021-01-01")) +
  geom_point() + geom_line() + theme_bw() + xlab("Date") +
  annotate("text", x=as.Date("2014-07-01"), y=1.1, label= "2014") +
  annotate("text", x=as.Date("2015-07-01"), y=1.1, label= "2015") +
  annotate("text", x=as.Date("2016-07-01"), y=1.1, label= "2016") +
  annotate("text", x=as.Date("2019-07-01"), y=1.1, label= "2019") +
  annotate("text", x=as.Date("2020-07-01"), y=1.1, label= "2020") +
  annotate("text", x=as.Date("2021-07-01"), y=1.1, label= "2021") +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("Saguaro", 3))
ggsave("Figures/zoop_3taxa_std_biom.jpg", width=6, height=3) 

ggplot(zoops_12_groups, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"),
                            standardized_biom, color=Taxon)) +
  geom_vline(xintercept = as.Date("2014-01-01")) +
  geom_vline(xintercept = as.Date("2015-01-01")) +
  geom_vline(xintercept = as.Date("2016-01-01")) +
  geom_vline(xintercept = as.Date("2017-01-01")) +
  geom_vline(xintercept = as.Date("2019-01-01")) +
  geom_vline(xintercept = as.Date("2020-01-01")) +
  geom_vline(xintercept = as.Date("2021-01-01")) +
  geom_point() + geom_line() + theme_bw() + xlab("Date") +
  annotate("text", x=as.Date("2014-07-01"), y=1.1, label= "2014") +
  annotate("text", x=as.Date("2015-07-01"), y=1.1, label= "2015") +
  annotate("text", x=as.Date("2016-07-01"), y=1.1, label= "2016") +
  annotate("text", x=as.Date("2019-07-01"), y=1.1, label= "2019") +
  annotate("text", x=as.Date("2020-07-01"), y=1.1, label= "2020") +
  annotate("text", x=as.Date("2021-07-01"), y=1.1, label= "2021") +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("Torres", 12))
ggsave("Figures/zoop_12taxa_std_biom.jpg", width=6, height=3) 

#------------------------------------------------------------------------------#
#biomass vs. doy
ggplot(zoops_total, aes(yday(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d")),
                        Biom_std, color=year)) +
  geom_point() + geom_line() + theme_bw() + xlab("doy") +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("KingsCanyon", 6))
ggsave("Figures/zoop_total_std_biom_vs_doy.jpg", width=6, height=3) 

ggplot(zoops_3_groups, aes(yday(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d")),
                           standardized_biom, color=year)) +
  facet_wrap(~Taxon, nrow=3)+
  geom_point() + geom_line() + theme_bw() + xlab("doy") +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("KingsCanyon", 6))
ggsave("Figures/zoop_3taxa_std_biom_vs_doy.jpg", width=6, height=3) 

ggplot(zoops_12_groups, aes(yday(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d")),
                            standardized_biom, color=year)) +
  facet_wrap(~Taxon, nrow=4)+
  geom_point() + geom_line() + theme_bw() + xlab("doy") +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("KingsCanyon", 6))
ggsave("Figures/zoop_12taxa_std_biom_vs_doy.jpg", width=6, height=3) 
