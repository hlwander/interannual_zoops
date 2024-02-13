#zoop size succession figs

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
  summarise(Daphnia_size = mean(MeanLength_mm[
    Taxon %in% c("D. catawba","D. ambigua")],na.rm=T),
    Daphnia_sd = sd(MeanLength_mm[
      Taxon %in% c("D. catawba","D. ambigua")],na.rm=T),
    Calanoida_size = mean(MeanLength_mm[
      Taxon %in% c("Diaptomus")],na.rm=T),
    Calanoida_sd = sd(MeanLength_mm[
      Taxon %in% c("Diaptomus")],na.rm=T),
    Cyclopoida_size = mean(MeanLength_mm[
      Taxon %in% c("Cyclopoids")],na.rm=T),
    Cyclopoida_sd = sd(MeanLength_mm[
      Taxon %in% c("Cyclopoids")],na.rm=T),
    Nauplii_size = mean(MeanLength_mm[
      Taxon %in% c("Nauplii")],na.rm=T),
    Nauplii_sd = sd(MeanLength_mm[
      Taxon %in% c("Nauplii")],na.rm=T),
    Bosmina_size = mean(MeanLength_mm[
      Taxon %in% c("Bosmina")],na.rm=T),
    Bosmina_sd = sd(MeanLength_mm[
      Taxon %in% c("Bosmina")],na.rm=T),
    #Chydorus = mean(MeanLength_mm[
    #  Taxon %in% c("Chydorus")]),
    Ceriodaphnia_size = mean(MeanLength_mm[
      Taxon %in% c("Ceriodaphnia")],na.rm=T),
    Ceriodaphnia_sd = sd(MeanLength_mm[
      Taxon %in% c("Ceriodaphnia")],na.rm=T),
    #Diaphanosoma = mean(MeanLength_mm[
    #  Taxon %in% c("Diaphanosoma")]),
    Ascomorpha_size = mean(MeanLength_mm[
      Taxon %in% c("Ascomorpha")],na.rm=T),
    Ascomorpha_sd = sd(MeanLength_mm[
      Taxon %in% c("Ascomorpha")],na.rm=T),
    #Asplanchna = mean(MeanLength_mm[
    #  Taxon %in% c("Asplanchna")]),
    Conochilus_size = mean(MeanLength_mm[
      Taxon %in% c("Conochilus")],na.rm=T),
    Conochilus_sd = sd(MeanLength_mm[
      Taxon %in% c("Conochilus")],na.rm=T),
    Keratella_size = mean(MeanLength_mm[
      Taxon %in% c("Keratella")],na.rm=T),
    Keratella_sd = sd(MeanLength_mm[
      Taxon %in% c("Keratella")],na.rm=T),
    Trichocerca_size = mean(MeanLength_mm[
      Taxon %in% c("Trichocerca")],na.rm=T),
    Trichocerca_sd = sd(MeanLength_mm[
      Taxon %in% c("Trichocerca")],na.rm=T),
    Kellicottia_size = mean(MeanLength_mm[
      Taxon %in% c("Kellicottia")],na.rm=T),
    Kellicottia_sd = sd(MeanLength_mm[
      Taxon %in% c("Kellicottia")],na.rm=T),
    # Lecane = mean(MeanLength_mm[
    #   Taxon %in% c("Lecane")]),
    Polyarthra_size = mean(MeanLength_mm[
      Taxon %in% c("Polyarthra")],na.rm=T),
    Polyarthra_sd = sd(MeanLength_mm[
      Taxon %in% c("Polyarthra")],na.rm=T),
    Rotifera_size = mean(MeanLength_mm[
      Taxon %in% c("Total rotifers")],na.rm=T),
    Rotifera_sd = sd(MeanLength_mm[
      Taxon %in% c("Total rotifers")],na.rm=T),
    Cladocera_size = mean(MeanLength_mm[
      Taxon %in% c("Bosmina","D. catawba", "Chydorus","D. ambigua",
                   "Diaphanosoma","Ceriodaphnia")],na.rm=T),
    Cladocera_sd = sd(MeanLength_mm[
      Taxon %in% c("Bosmina","D. catawba", "Chydorus","D. ambigua",
                   "Diaphanosoma","Ceriodaphnia")],na.rm=T),
    Copepoda_size = mean(MeanLength_mm[
      Taxon %in% c("Diaptomus","Nauplii", "Cyclopoids")],na.rm=T),
    Copepoda_sd = sd(MeanLength_mm[
      Taxon %in% c("Diaptomus","Nauplii", "Cyclopoids")],na.rm=T)) |> 
  pivot_longer(-c(Reservoir,DateTime,StartDepth_m),
               names_to = c("Taxon", ".value"),
               names_sep="_" )  |> 
  filter(hour(DateTime) %in% c(9,10,11,12,13,14)) |> #drop nighttime samples
  mutate(DateTime = as.Date(DateTime)) 

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
  summarise(size = mean(MeanLength_mm,na.rm=T),
            sd = sd(MeanLength_mm,na.rm=T))

#combine all zoop data
all_zoops <- bind_rows(zoops_pre, zoops_post) |> 
  mutate_all(~replace(., is.nan(.), NA)) |>  #replace NAN with NA
  ungroup() |> select(-StartDepth_m) #dropping, but note that depths range from 7.5-11.5m....

#write all_zoops
write.csv(all_zoops, paste0("Output/all_zoops_size.csv"),row.names = FALSE)

#add column for pre vs post
all_zoops$data <- ifelse(all_zoops$DateTime<="2019-01-01","pre","post")

#------------------------------------------------------------------------------#
# new df for total zoop size
zoops_total <- all_zoops |> 
  group_by(DateTime) |> 
  summarise(Total = mean(size[Taxon %in% c("Cladocera","Copeoda","Rotifera")],na.rm=T),
            sd = mean(sd[Taxon %in% c("Cladocera","Copeoda","Rotifera")],na.rm=T)) |> 
  ungroup() |> 
  mutate(year = format(DateTime, "%Y"),
         month = format(DateTime, "%m")) |> 
  group_by(year, month) |>
  summarise(Total_avg = mean(Total,na.rm=T),
            Total_sd = mean(sd,na.rm=T)) 

# 3 group zoop size
zoops_3_groups <- all_zoops |> group_by(Taxon, DateTime) |> 
  filter(Taxon %in% c("Cladocera","Copepoda","Rotifera")) |> 
  mutate(year = format(DateTime, "%Y"),
         month = format(DateTime, "%m")) |> 
  ungroup() |> group_by(year, month) |> 
  summarise(Cladocera_avg = mean(size[Taxon=="Cladocera"],na.rm=T),
            Cladocera_sd = mean(sd[Taxon=="Cladocera"],na.rm=T),
            Copepoda_avg = mean(size[Taxon=="Copepoda"],na.rm=T),
            Copepoda_sd = mean(sd[Taxon=="Copepoda"],na.rm=T),
            Rotifera_avg = mean(size[Taxon=="Rotifera"],na.rm=T),
            Rotifera_sd = mean(sd[Taxon=="Rotifera"],na.rm=T)) |> 
  pivot_longer(-c(year,month),
               names_to = c("Taxon", ".value"),
               names_sep="_" ) 

# 12 group zoop size
zoops_12_groups <- all_zoops |> group_by(Taxon, DateTime) |> 
  filter(!Taxon %in% c("Cladocera","Copepoda","Rotifera")) |> 
  mutate(year = format(DateTime, "%Y"),
         month = format(DateTime, "%m")) |> 
  group_by(year, month) |> 
  summarise(Daphnia_avg = mean(size[Taxon=="Daphnia"],na.rm=T),
            Daphnia_sd = mean(sd[Taxon=="Daphnia"],na.rm=T),
            Calanoida_avg = mean(size[Taxon=="Calanoida"],na.rm=T),
            Calanoida_sd = mean(sd[Taxon=="Calanoida"],na.rm=T),
            Cyclopoida_avg = mean(size[Taxon=="Cyclopoida"],na.rm=T),
            Cyclopoida_sd = mean(sd[Taxon=="Cyclopoida"],na.rm=T),
            Nauplii_avg = mean(size[Taxon=="Nauplii"],na.rm=T),
            Nauplii_sd = mean(sd[Taxon=="Nauplii"],na.rm=T),
            Bosmina_avg = mean(size[Taxon=="Bosmina"],na.rm=T),
            Bosmina_sd = mean(sd[Taxon=="Bosmina"],na.rm=T),
            Ceriodaphnia_avg = mean(size[Taxon=="Ceriodaphnia"],na.rm=T),
            Ceriodaphnia_sd = mean(sd[Taxon=="Ceriodaphnia"],na.rm=T),
            Ascomorpha_avg = mean(size[Taxon=="Ascomorpha"],na.rm=T),
            Ascomorpha_sd = mean(sd[Taxon=="Ascomorpha"],na.rm=T),
            Conochilus_avg = mean(size[Taxon=="Conochilus"],na.rm=T),
            Conochilus_sd = mean(sd[Taxon=="Conochilus"],na.rm=T),
            Keratella_avg = mean(size[Taxon=="Keratella"],na.rm=T),
            Keratella_sd = mean(sd[Taxon=="Keratella"],na.rm=T),
            Trichocerca_avg = mean(size[Taxon=="Trichocerca"],na.rm=T),
            Trichocerca_sd = mean(sd[Taxon=="Trichocerca"],na.rm=T),
            Kellicottia_avg = mean(size[Taxon=="Kellicottia"],na.rm=T),
            Kellicottia_sd = mean(sd[Taxon=="Kellicottia"],na.rm=T),
            Polyarthra_avg = mean(size[Taxon=="Polyarthra"],na.rm=T),
            Polyarthra_sd = mean(sd[Taxon=="Polyarthra"],na.rm=T)) |> 
  pivot_longer(-c(year,month),
               names_to = c("Taxon", ".value"),
               names_sep="_" )

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
  annotate("text", x=as.Date("2014-07-01"), y=0.8, label= "2014") +
  annotate("text", x=as.Date("2015-07-01"), y=0.8, label= "2015") +
  annotate("text", x=as.Date("2016-07-01"), y=0.8, label= "2016") +
  annotate("text", x=as.Date("2019-07-01"), y=0.8, label= "2019") +
  annotate("text", x=as.Date("2020-07-01"), y=0.8, label= "2020") +
  annotate("text", x=as.Date("2021-07-01"), y=0.8, label= "2021") +
  geom_errorbar(aes(ymin = Total_avg -Total_sd, ymax = Total_avg+Total_sd), 
                color="#009999") 
ggsave("Figures/zoop_total_size.jpg", width=6, height=3) 

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
  annotate("text", x=as.Date("2014-07-01"), y=1, label= "2014") +
  annotate("text", x=as.Date("2015-07-01"), y=1, label= "2015") +
  annotate("text", x=as.Date("2016-07-01"), y=1, label= "2016") +
  annotate("text", x=as.Date("2019-07-01"), y=1, label= "2019") +
  annotate("text", x=as.Date("2020-07-01"), y=1, label= "2020") +
  annotate("text", x=as.Date("2021-07-01"), y=1, label= "2021") +
  geom_errorbar(aes(ymin = avg -sd, ymax = avg+sd)) +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("Saguaro", 3))
ggsave("Figures/zoop_3taxa_size.jpg", width=6, height=3) 

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
  annotate("text", x=as.Date("2014-07-01"), y=1, label= "2014") +
  annotate("text", x=as.Date("2015-07-01"), y=1, label= "2015") +
  annotate("text", x=as.Date("2016-07-01"), y=1, label= "2016") +
  annotate("text", x=as.Date("2019-07-01"), y=1, label= "2019") +
  annotate("text", x=as.Date("2020-07-01"), y=1, label= "2020") +
  annotate("text", x=as.Date("2021-07-01"), y=1, label= "2021") +
  geom_errorbar(aes(ymin = avg -sd, ymax = avg+sd)) +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("Torres", 12))
ggsave("Figures/zoop_12taxa_size.jpg", width=6, height=3) 

#------------------------------------------------------------------------------#
#size vs. doy
ggplot(zoops_total, aes(yday(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d")),
                        Total_avg, color=year)) +
  geom_point() + geom_line() + theme_bw() + xlab("doy") +
  geom_errorbar(aes(ymin = Total_avg -Total_sd, ymax = Total_avg+Total_sd)) +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("KingsCanyon", 6))
ggsave("Figures/zoop_size_vs_doy.jpg", width=6, height=3) 

ggplot(zoops_3_groups, aes(yday(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d")),
                           avg, color=year)) +
  facet_wrap(~Taxon, nrow=3)+
  geom_point() + geom_line() + theme_bw() + xlab("doy") +
  geom_errorbar(aes(ymin = avg -sd, ymax = avg+sd)) +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("KingsCanyon", 6))
ggsave("Figures/zoop_3taxa_size_vs_doy.jpg", width=6, height=3) 

ggplot(zoops_12_groups, aes(yday(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d")),
                            avg, color=year)) +
  facet_wrap(~Taxon, nrow=4)+
  geom_point() + geom_line() + theme_bw() + xlab("doy") +
  geom_errorbar(aes(ymin = avg -sd, ymax = avg+sd)) +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("KingsCanyon", 6))
ggsave("Figures/zoop_12taxa_size_vs_doy.jpg", width=6, height=3) 

