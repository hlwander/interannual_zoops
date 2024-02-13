#zoop density succession figs

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
  summarise(Daphnia_dens = sum(Density_IndPerL[
    Taxon %in% c("D. catawba","D. ambigua")]),
    Daphnia_sd = sd(Density_IndPerL[
      Taxon %in% c("D. catawba","D. ambigua")]),
    Calanoida_dens = sum(Density_IndPerL[
      Taxon %in% c("Diaptomus")]),
    Calanoida_sd = sd(Density_IndPerL[
      Taxon %in% c("Diaptomus")]),
    Cyclopoida_dens = sum(Density_IndPerL[
      Taxon %in% c("Cyclopoids")]),
    Cyclopoida_sd = sd(Density_IndPerL[
      Taxon %in% c("Cyclopoids")]),
    Nauplii_dens = sum(Density_IndPerL[
      Taxon %in% c("Nauplii")]),
    Nauplii_sd = sd(Density_IndPerL[
      Taxon %in% c("Nauplii")]),
    Bosmina_dens = sum(Density_IndPerL[
      Taxon %in% c("Bosmina")]),
    Bosmina_sd = sd(Density_IndPerL[
      Taxon %in% c("Bosmina")]),
    #Chydorus = sum(Density_IndPerL[
    #  Taxon %in% c("Chydorus")]),
    Ceriodaphnia_dens = sum(Density_IndPerL[
      Taxon %in% c("Ceriodaphnia")]),
    Ceriodaphnia_sd = sd(Density_IndPerL[
      Taxon %in% c("Ceriodaphnia")]),
    #Diaphanosoma = sum(Density_IndPerL[
    #  Taxon %in% c("Diaphanosoma")]),
    Ascomorpha_dens = sum(Density_IndPerL[
      Taxon %in% c("Ascomorpha")]),
    Ascomorpha_sd = sd(Density_IndPerL[
      Taxon %in% c("Ascomorpha")]),
    #Asplanchna = sum(Density_IndPerL[
    #  Taxon %in% c("Asplanchna")]),
    Conochilus_dens = sum(Density_IndPerL[
      Taxon %in% c("Conochilus")]),
    Conochilus_sd = sd(Density_IndPerL[
      Taxon %in% c("Conochilus")]),
    Keratella_dens = sum(Density_IndPerL[
      Taxon %in% c("Keratella")]),
    Keratella_sd = sd(Density_IndPerL[
      Taxon %in% c("Keratella")]),
    Trichocerca_dens = sum(Density_IndPerL[
      Taxon %in% c("Trichocerca")]),
    Trichocerca_sd = sd(Density_IndPerL[
      Taxon %in% c("Trichocerca")]),
    Kellicottia_dens = sum(Density_IndPerL[
      Taxon %in% c("Kellicottia")]),
    Kellicottia_sd = sd(Density_IndPerL[
      Taxon %in% c("Kellicottia")]),
   # Lecane = sum(Density_IndPerL[
   #   Taxon %in% c("Lecane")]),
    Polyarthra_dens = sum(Density_IndPerL[
      Taxon %in% c("Polyarthra")]),
    Polyarthra_sd = sd(Density_IndPerL[
     Taxon %in% c("Polyarthra")]),
    Rotifera_dens = sum(Density_IndPerL[
      Taxon %in% c("Total rotifers")]),
    Rotifera_sd = sd(Density_IndPerL[
     Taxon %in% c("Total rotifers")]),
    Cladocera_dens = sum(Density_IndPerL[
      Taxon %in% c("Bosmina","D. catawba", "Chydorus","D. ambigua",
                   "Diaphanosoma","Ceriodaphnia")]),
    Cladocera_sd = sd(Density_IndPerL[
     Taxon %in% c("Bosmina","D. catawba", "Chydorus","D. ambigua",
                  "Diaphanosoma","Ceriodaphnia")]),
    Copepoda_dens = sum(Density_IndPerL[
      Taxon %in% c("Diaptomus","Nauplii", "Cyclopoids")]),
    Copepoda_sd = sd(Density_IndPerL[
     Taxon %in% c("Diaptomus","Nauplii", "Cyclopoids")])) |> 
    pivot_longer(-c(Reservoir,DateTime,StartDepth_m),
               names_to = c("Taxon", ".value"),
               names_sep="_" )  |> 
  filter(hour(DateTime) %in% c(9,10,11,12,13,14)) |> #drop nighttime samples
  mutate(DateTime = as.Date(DateTime)) |> 
  mutate(dens = dens * (1/0.031))  #10m bvr neteff from 2016 (n=2) - note that 7m neteff was also 0.031
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
  summarise(dens = mean(Density_IndPerL),
            sd = sd(Density_IndPerL))

#combine all zoop data
all_zoops <- bind_rows(zoops_pre, zoops_post) |> 
  mutate_all(~replace(., is.nan(.), NA)) |>  #replace NAN with NA
  ungroup() |> select(-StartDepth_m) #dropping, but note that depths range from 7.5-11.5m....

#write all_zoops
write.csv(all_zoops, paste0("Output/all_zoops_dens.csv"),row.names = FALSE)

#add column for pre vs post
all_zoops$data <- ifelse(all_zoops$DateTime<="2019-01-01","pre","post")

#calculate proportion of total density for each taxa by day
all_zoops <- all_zoops |> group_by(DateTime) |> 
  mutate(n = sum(dens)) |> ungroup() |> 
           group_by(DateTime,Taxon) |> 
    mutate(prop = dens / sum(n))

#order data levels
all_zoops$data <- factor(all_zoops$data, levels=c("pre", "post"))

#create bar plots comparing taxa densities pre vs post
ggplot(data=subset(all_zoops, Taxon %in% 
                     c("Cladocera","Copepoda","Rotifera")), 
              aes(x=data, fill=Taxon ))+
  geom_bar(width = 1, stat="identity", aes(y=prop))
  #coord_polar("y")

#clads
ggplot(data=subset(all_zoops, Taxon %in% 
                     c("Daphnia","Ceriodaphnia","Bosmina")), 
       aes(x=data, fill=Taxon ))+theme_bw() +
  geom_bar(width = 1, stat="identity", aes(y=prop))

#copes
ggplot(data=subset(all_zoops, Taxon %in% 
                     c("Calanoida","Cyclopoida","Nauplii")), 
       aes(x=data, fill=Taxon ))+ theme_bw() +
  geom_bar(width = 1, stat="identity", aes(y=prop))

#rots
ggplot(data=subset(all_zoops, Taxon %in% 
                     c("Conochilus","Keratella","Kellicottia",
                       "Lecan","Trichocerca")), 
       aes(x=data, fill=Taxon ))+theme_bw() +
  geom_bar(width = 1, stat="identity", aes(y=prop))

#------------------------------------------------------------------------------#
#figure out dominant taxa for NMDS/other ms figs

#plot proportion of density that each taxon makes up (sum of all individual days)
ggplot(all_zoops, aes(x=Taxon, y=prop)) +
  theme_bw() + geom_bar(stat="identity") + 
  theme(text = element_text(size=5), 
        axis.text = element_text(size=5, color="black"), 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, hjust=1), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#count how many 0s there are
#num_zero <- all_zoops |> group_by(Taxon) |> 
#  summarise(num = sum(dens == 0)) 
#only keep taxa with <60 0 values (drop lecane, diaphanosoma, asplanchna, chydorus)

#create new df of mean prop per taxon
summary_prop <- all_zoops |> group_by(Taxon) |> 
  summarise(prop = mean(prop))

ggplot(summary_prop, aes(x=reorder(Taxon,-prop), y=prop)) +
  theme_bw() + geom_bar(stat="identity") + 
  theme(text = element_text(size=5), 
        axis.text = element_text(size=5, color="black"), 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, hjust=1), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#there are 19 (technically 16 not including clads, copes, rots) 
#taxa in common between JPD and HLW years - just using 12 for NMDS based on prevalence of 0 values

#------------------------------------------------------------------------------#
# new df for total zoop dens
zoops_total <- all_zoops |> 
  group_by(DateTime) |> 
  summarise(Total = sum(dens[Taxon %in% c("Cladocera","Copeoda","Rotifera")]),
            sd = mean(sd[Taxon %in% c("Cladocera","Copeoda","Rotifera")],na.rm=T)) |> 
  ungroup() |> 
  mutate(year = format(DateTime, "%Y"),
         month = format(DateTime, "%m")) |> 
  group_by(year, month) |>
  summarise(Total_avg = mean(Total,na.rm=T),
            Total_sd = mean(sd,na.rm=T)) |> 
  ungroup() |> 
  mutate(Dens_std =(Total_avg - min(Total_avg)) / 
           (max(Total_avg) - min(Total_avg)))

# 3 group zoop dens
zoops_3_groups <- all_zoops |> group_by(Taxon, DateTime) |> 
  filter(Taxon %in% c("Cladocera","Copepoda","Rotifera")) |> 
  mutate(year = format(DateTime, "%Y"),
         month = format(DateTime, "%m")) |> 
  ungroup() |> group_by(year, month) |> 
  summarise(Cladocera_avg = mean(dens[Taxon=="Cladocera"]),
            Cladocera_sd = mean(sd[Taxon=="Cladocera"],na.rm=T),
            Copepoda_avg = mean(dens[Taxon=="Copepoda"]),
            Copepoda_sd = mean(sd[Taxon=="Copepoda"],na.rm=T),
            Rotifera_avg = mean(dens[Taxon=="Rotifera"]),
            Rotifera_sd = mean(sd[Taxon=="Rotifera"],na.rm=T)) |> 
  pivot_longer(-c(year,month),
              names_to = c("Taxon", ".value"),
              names_sep="_" )  |> 
  ungroup() |> group_by(Taxon) |>
  mutate(min_dens = min(avg),
         max_dens = max(avg)) |> 
  mutate(standardized_dens = (avg - min_dens) / (max_dens - min_dens))
write.csv(zoops_3_groups,"Output/std_dens_3taxa.csv", row.names = F)

# 12 group zoop dens
zoops_12_groups <- all_zoops |> group_by(Taxon, DateTime) |> 
  filter(!Taxon %in% c("Cladocera","Copepoda","Rotifera")) |> 
  mutate(year = format(DateTime, "%Y"),
         month = format(DateTime, "%m")) |> 
  group_by(year, month) |> 
  summarise(Daphnia_avg = mean(dens[Taxon=="Daphnia"]),
            Daphnia_sd = mean(sd[Taxon=="Daphnia"],na.rm=T),
            Calanoida_avg = mean(dens[Taxon=="Calanoida"]),
            Calanoida_sd = mean(sd[Taxon=="Calanoida"],na.rm=T),
            Cyclopoida_avg = mean(dens[Taxon=="Cyclopoida"]),
            Cyclopoida_sd = mean(sd[Taxon=="Cyclopoida"],na.rm=T),
            Nauplii_avg = mean(dens[Taxon=="Nauplii"]),
            Nauplii_sd = mean(sd[Taxon=="Nauplii"],na.rm=T),
            Bosmina_avg = mean(dens[Taxon=="Bosmina"]),
            Bosmina_sd = mean(sd[Taxon=="Bosmina"],na.rm=T),
            Ceriodaphnia_avg = mean(dens[Taxon=="Ceriodaphnia"]),
            Ceriodaphnia_sd = mean(sd[Taxon=="Ceriodaphnia"],na.rm=T),
            Ascomorpha_avg = mean(dens[Taxon=="Ascomorpha"]),
            Ascomorpha_sd = mean(sd[Taxon=="Ascomorpha"],na.rm=T),
            Conochilus_avg = mean(dens[Taxon=="Conochilus"]),
            Conochilus_sd = mean(sd[Taxon=="Conochilus"],na.rm=T),
            Keratella_avg = mean(dens[Taxon=="Keratella"]),
            Keratella_sd = mean(sd[Taxon=="Keratella"],na.rm=T),
            Trichocerca_avg = mean(dens[Taxon=="Trichocerca"]),
            Trichocerca_sd = mean(sd[Taxon=="Trichocerca"],na.rm=T),
            Kellicottia_avg = mean(dens[Taxon=="Kellicottia"]),
            Kellicottia_sd = mean(sd[Taxon=="Kellicottia"],na.rm=T),
            Polyarthra_avg = mean(dens[Taxon=="Polyarthra"]),
            Polyarthra_sd = mean(sd[Taxon=="Polyarthra"],na.rm=T)) |> 
  pivot_longer(-c(year,month),
               names_to = c("Taxon", ".value"),
               names_sep="_" )  |> 
  ungroup() |> group_by(Taxon) |>
  mutate(min_dens = min(avg),
         max_dens = max(avg)) |> 
  mutate(standardized_dens = (avg - min_dens) / (max_dens - min_dens))

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
  annotate("text", x=as.Date("2014-07-01"), y=12000, label= "2014") +
  annotate("text", x=as.Date("2015-07-01"), y=12000, label= "2015") +
  annotate("text", x=as.Date("2016-07-01"), y=12000, label= "2016") +
  annotate("text", x=as.Date("2019-07-01"), y=12000, label= "2019") +
  annotate("text", x=as.Date("2020-07-01"), y=12000, label= "2020") +
  annotate("text", x=as.Date("2021-07-01"), y=12000, label= "2021") +
  geom_errorbar(aes(ymin = Total_avg -Total_sd, ymax = Total_avg+Total_sd), 
                color="#009999") 
ggsave("Figures/zoop_total_dens.jpg", width=6, height=3) 

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
  annotate("text", x=as.Date("2014-07-01"), y=15000, label= "2014") +
  annotate("text", x=as.Date("2015-07-01"), y=15000, label= "2015") +
  annotate("text", x=as.Date("2016-07-01"), y=15000, label= "2016") +
  annotate("text", x=as.Date("2019-07-01"), y=15000, label= "2019") +
  annotate("text", x=as.Date("2020-07-01"), y=15000, label= "2020") +
  annotate("text", x=as.Date("2021-07-01"), y=15000, label= "2021") +
  geom_errorbar(aes(ymin = avg -sd, ymax = avg+sd)) +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("Saguaro", 3))
ggsave("Figures/zoop_3taxa_dens.jpg", width=6, height=3) 

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
  annotate("text", x=as.Date("2014-07-01"), y=12000, label= "2014") +
  annotate("text", x=as.Date("2015-07-01"), y=12000, label= "2015") +
  annotate("text", x=as.Date("2016-07-01"), y=12000, label= "2016") +
  annotate("text", x=as.Date("2019-07-01"), y=12000, label= "2019") +
  annotate("text", x=as.Date("2020-07-01"), y=12000, label= "2020") +
  annotate("text", x=as.Date("2021-07-01"), y=12000, label= "2021") +
  geom_errorbar(aes(ymin = avg -sd, ymax = avg+sd)) +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("Torres", 12))
ggsave("Figures/zoop_12taxa_dens.jpg", width=6, height=3) 

#------------------------------------------------------------------------------#
#standardized density
ggplot(zoops_total, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"),Dens_std)) +
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
ggsave("Figures/zoop_total_std_dens.jpg", width=6, height=3) 

ggplot(zoops_3_groups, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"),
                           standardized_dens, color=Taxon)) +
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
ggsave("Figures/zoop_3taxa_std_dens.jpg", width=6, height=3) 

ggplot(zoops_12_groups, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"),
                            standardized_dens, color=Taxon)) +
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
ggsave("Figures/zoop_12taxa_std_dens.jpg", width=6, height=3) 

#------------------------------------------------------------------------------#
#density vs. doy
ggplot(zoops_total, aes(yday(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d")),
                        Dens_std, color=year)) +
  geom_point() + geom_line() + theme_bw() + xlab("doy") +
scale_color_manual("",values=NatParksPalettes::natparks.pals("KingsCanyon", 6))
ggsave("Figures/zoop_total_std_dens_vs_doy.jpg", width=6, height=3) 

ggplot(zoops_3_groups, aes(yday(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d")),
                        standardized_dens, color=year)) +
  facet_wrap(~Taxon, nrow=3)+
  geom_point() + geom_line() + theme_bw() + xlab("doy") +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("KingsCanyon", 6))
ggsave("Figures/zoop_3taxa_std_dens_vs_doy.jpg", width=6, height=3) 

ggplot(zoops_12_groups, aes(yday(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d")),
                           standardized_dens, color=year)) +
  facet_wrap(~Taxon, nrow=4)+
  geom_point() + geom_line() + theme_bw() + xlab("doy") +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("KingsCanyon", 6))
ggsave("Figures/zoop_12taxa_std_dens_vs_doy.jpg", width=6, height=3) 
