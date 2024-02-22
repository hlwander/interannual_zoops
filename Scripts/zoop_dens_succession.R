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
    #Calanoida_dens = sum(Density_IndPerL[
    #  Taxon %in% c("Diaptomus")]),
    #Calanoida_sd = sd(Density_IndPerL[
    #  Taxon %in% c("Diaptomus")]),
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
    #Chydorus_dens = sum(Density_IndPerL[
    #  Taxon %in% c("Chydorus")]),
    #Chydorus_sd = sd(Density_IndPerL[
    #  Taxon %in% c("Chydorus")]),
    Ceriodaphnia_dens = sum(Density_IndPerL[
      Taxon %in% c("Ceriodaphnia")]),
    Ceriodaphnia_sd = sd(Density_IndPerL[
      Taxon %in% c("Ceriodaphnia")]),
    #Diaphanosoma_dens = sum(Density_IndPerL[
    #  Taxon %in% c("Diaphanosoma")]),
    #Diaphanosoma_sd = sd(Density_IndPerL[
    #  Taxon %in% c("Diaphanosoma")]),
    Ascomorpha_dens = sum(Density_IndPerL[
      Taxon %in% c("Ascomorpha")]),
    Ascomorpha_sd = sd(Density_IndPerL[
      Taxon %in% c("Ascomorpha")]),
    #Asplanchna_dens = sum(Density_IndPerL[
    #  Taxon %in% c("Asplanchna")]),
    #Asplanchna_sd = sd(Density_IndPerL[
    #  Taxon %in% c("Asplanchna")]),
    Conochilus_dens = sum(Density_IndPerL[
      Taxon %in% c("Conochilus")]),
    Conochilus_sd = sd(Density_IndPerL[
      Taxon %in% c("Conochilus")]),
    Keratella_dens = sum(Density_IndPerL[
      Taxon %in% c("Keratella")]),
    Keratella_sd = sd(Density_IndPerL[
      Taxon %in% c("Keratella")]),
    #Trichocerca_dens = sum(Density_IndPerL[
    #  Taxon %in% c("Trichocerca")]),
    #Trichocerca_sd = sd(Density_IndPerL[
    #  Taxon %in% c("Trichocerca")]),
    Kellicottia_dens = sum(Density_IndPerL[
      Taxon %in% c("Kellicottia")]),
    Kellicottia_sd = sd(Density_IndPerL[
      Taxon %in% c("Kellicottia")]),
    #Lecane_dens = sum(Density_IndPerL[
    #  Taxon %in% c("Lecane")]),
    #Lecane_sd = sd(Density_IndPerL[
    #  Taxon %in% c("Lecane")]),
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
          "Cyclopoida", "Nauplius", 
          "Conochilus","Keratella", "Rotifera",
          "Kellicottia","Ascomorpha",
          "Polyarthra","Cladocera", "Copepoda") 
    # "Chydorus","Diaphanosoma", "Lecane","Asplanchna","Trichocerca","Calanoida"

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

#calculate proportion of total density for each taxa (dens/total dens)
all_zoop_taxa <- all_zoops |> 
  filter(!Taxon %in% c("Cladocera","Copepoda","Rotifera")) |> 
  mutate(n = sum(dens)) |> 
           group_by(Taxon) |> 
    mutate(prop = sum(dens) / n)

zoops_3taxa <- all_zoops |> 
  filter(Taxon %in% c("Cladocera","Copepoda","Rotifera")) |> 
  mutate(n = sum(dens)) |>
  group_by(Taxon) |> 
  mutate(prop = sum(dens) / n)

#order data levels
all_zoop_taxa$data <- factor(all_zoop_taxa$data, levels=c("pre", "post"))
zoops_3taxa$data <- factor(zoops_3taxa$data, levels=c("pre", "post"))

#create bar plots comparing taxa densities pre vs post
ggplot(data=subset(zoops_3taxa, Taxon %in% 
                     c("Cladocera","Copepoda","Rotifera")), 
              aes(x=data, fill=Taxon ))+
  geom_bar(width = 1, stat="identity", aes(y=prop))

#clads
ggplot(data=subset(all_zoop_taxa, Taxon %in% 
                     c("Daphnia","Ceriodaphnia","Bosmina")), 
       aes(x=data, fill=Taxon ))+theme_bw() +
  geom_bar(width = 1, stat="identity", aes(y=prop))

#copes
ggplot(data=subset(all_zoop_taxa, Taxon %in% 
                     c("Calanoida","Cyclopoida","Nauplii")), 
       aes(x=data, fill=Taxon ))+ theme_bw() +
  geom_bar(width = 1, stat="identity", aes(y=prop))

#rots
ggplot(data=subset(all_zoop_taxa, Taxon %in% 
                     c("Conochilus","Keratella","Kellicottia",
                       "Trichocerca")), 
       aes(x=data, fill=Taxon ))+theme_bw() +
  geom_bar(width = 1, stat="identity", aes(y=prop))

#------------------------------------------------------------------------------#
#figure out dominant taxa for NMDS/other ms figs

#plot proportion of density that each taxon makes up (sum of all individual days)
ggplot(all_zoop_taxa, aes(x=reorder(Taxon,-prop), y=prop)) +
  theme_bw() + geom_bar(stat="identity") + 
  theme(text = element_text(size=5), 
        axis.text = element_text(size=5, color="black"), 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, hjust=1), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#new df to average proportion for each taxon
zoop_prop_taxa <- all_zoop_taxa |> group_by(Taxon) |> 
  summarise(prop_avg = mean(prop)) 
#drop taxa that are <1% (0.01) of density 
#dropping diaphanosmoma, lecane, chydorus, asplanchna, trichocerca, calanoida

#just keeping n=10

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
zoops_3_groups <- all_zoops |> 
  filter(Taxon %in% c("Cladocera","Copepoda","Rotifera")) |> 
  group_by(Taxon, DateTime) |> 
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

# 10 group zoop dens
zoops_10_groups <- all_zoops |> 
  filter(!Taxon %in% c("Cladocera","Copepoda","Rotifera")) |> 
  group_by(Taxon, DateTime) |> 
  mutate(year = format(DateTime, "%Y"),
         month = format(DateTime, "%m")) |> 
  ungroup() |> 
  group_by(year, month) |> 
  summarise(Daphnia_avg = mean(dens[Taxon=="Daphnia"]),
            Daphnia_sd = mean(sd[Taxon=="Daphnia"],na.rm=T),
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

ggplot(zoops_10_groups, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"),
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
ggsave("Figures/zoop_10taxa_dens.jpg", width=6, height=3) 

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

ggplot(zoops_10_groups, aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"),
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
ggsave("Figures/zoop_10taxa_std_dens.jpg", width=6, height=3) 

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

ggplot(zoops_10_groups, aes(yday(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d")),
                           standardized_dens, color=year)) +
  facet_wrap(~Taxon, nrow=4)+
  geom_point() + geom_line() + theme_bw() + xlab("doy") +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("KingsCanyon", 6))
ggsave("Figures/zoop_10taxa_std_dens_vs_doy.jpg", width=6, height=3) 

#-----------------------------------------------------------------------------#
#shaded line plot for clads, copes, and rots each year

# 3 group zoop dens - standardize within a year
zoops_3_groups_years <- all_zoops |> 
  filter(Taxon %in% c("Cladocera","Copepoda","Rotifera")) |> 
  group_by(Taxon, DateTime) |> 
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
  ungroup() |> group_by(Taxon,year) |>
  mutate(min_dens = min(avg),
         max_dens = max(avg)) |> 
  mutate(standardized_dens = (avg - min_dens) / (max_dens - min_dens)) |> 
  mutate(month = as.numeric(month))

ggplot(zoops_3_groups_years, aes(yday(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d")),
                           standardized_dens, color=year)) +
  facet_wrap(~Taxon, nrow=3)+
  geom_point() + geom_line() + theme_bw() + xlab("doy") +
  scale_color_manual("",values=NatParksPalettes::natparks.pals("KingsCanyon", 6))
ggsave("Figures/zoop_3taxa_std_dens_by_year_vs_doy.jpg", width=6, height=3) 

#playing around with order/layering of taxa
zoops_3_groups_years$Taxon <- factor(zoops_3_groups_years$Taxon,
                              #levels = c("Rotifera","Cladocera","Copepoda"))
                              levels = c("Rotifera","Copepoda","Cladocera"))
                              #levels = c("Copepoda","Rotifera","Cladocera"))
                              #levels = c("Cladocera","Copepoda","Rotifera"))

#shaded line plot time
ggplot(data=subset(zoops_3_groups_years, month %in% c(5,6,7,8,9)), 
       aes(as.Date(paste0(year,"-",month,"-01"), "%Y-%m-%d"),
                                 standardized_dens, color=Taxon)) +
  geom_area(aes(color = Taxon, fill = Taxon),
            position = "stack", stat = "identity",
            alpha=0.7) +
  facet_wrap(~year, scales = "free")+
  scale_color_manual(values = NatParksPalettes::natparks.pals("KingsCanyon", 3))+
  scale_fill_manual(values = NatParksPalettes::natparks.pals("KingsCanyon", 3),
                    breaks = c("Cladocera","Copepoda","Rotifera"))+
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
ggsave("Figures/BVR_succession_3groups_subset_stacked.jpg", width=6, height=4) 

#-----------------------------------------------------------------------------#
#read in phyto csv
all_phytos_std <- read.csv("Output/phytos.csv") |> 
  select(month, year, Total_ugL) |> 
  group_by(year) |> 
  mutate(min_val = min(Total_ugL),
         max_val = max(Total_ugL)) |> 
  mutate(phyto_abund = (Total_ugL - min_val) / (max_val - min_val)) |> 
  select(month, year, phyto_abund) |> 
  arrange(year, month)
  
all_zoops_std <- all_zoops |> 
  mutate(year = format(DateTime, "%Y"),
         month = format(DateTime, "%m")) |> 
  mutate(month = as.numeric(month),
         year = as.numeric(year)) |> 
  filter(month %in% c(5,6,7,8,9)) |> 
  group_by(DateTime) |> 
  mutate(Total = sum(dens[Taxon %in% c("Cladocera","Copeoda","Rotifera")])) |> 
  group_by(year, month) |>
  summarise(Total_avg = mean(Total,na.rm=T)) |> 
  ungroup() |>   group_by(year) |> 
  mutate(min_val = min(Total_avg),
         max_val = max(Total_avg)) |> 
  mutate(zoop_dens = (Total_avg - min_val) / (max_val - min_val)) |> 
  select(month, year, zoop_dens)

#combine zoops and phytos
zoops_phytos <- left_join(all_zoops_std, all_phytos_std) |> 
  pivot_longer(-c(month,year),
               names_to = "variable")

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
  mutate(month = as.numeric(month),
         year = as.numeric(year)) |> 
  group_by(year, month) |>
  summarise(secchi_avg = mean(Secchi_m,na.rm=T)) |> 
  ungroup() |>  group_by(year) |> 
  mutate(min_val = min(secchi_avg),
         max_val = max(secchi_avg)) |> 
  mutate(secchi_std = (secchi_avg - min_val) / (max_val - min_val)) |> 
  select(month, year, secchi_std) 

#classic PEG model fig - phytos + zoops 
ggplot(zoops_phytos, aes(as.Date("2023-12-31") + yday(as.Date(
  paste0(year,"-",month,"-01"), "%Y-%m-%d")), value, color=variable)) + 
  geom_point() + geom_line() + theme_bw() + xlab("Month") +
  facet_wrap(~year) +
  geom_area(aes(color=variable, fill = variable), 
            position="identity", alpha = 0.7) +
  scale_color_manual(values = c("#809848", "#2E0014")) +
  scale_fill_manual(values = c("#809848", "#2E0014")) +
  geom_point(data=secchi, aes(as.Date("2023-12-31") + yday(as.Date(
    paste0(year,"-",month,"-01"), "%Y-%m-%d")), secchi_std,
    group=as.numeric(year)), col="red") +
  geom_line(data=secchi, aes(as.Date("2023-12-31") + yday(as.Date(
    paste0(year,"-",month,"-01"), "%Y-%m-%d")), secchi_std,
    group=as.numeric(year)), col="red")

ggsave("Figures/phyto_zoop_annual_peg_plus_secchi.jpg", width=6, height=4)


