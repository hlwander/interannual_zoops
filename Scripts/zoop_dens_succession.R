#zoop density succession figs

#read in packages
pacman::p_load(zoo, dplR, dplyr, tidyverse, ggplot2, ggpubr, sf, lubridate)

#read in zoop data from EDI
inUrl1  <-  "https://pasta.lternet.edu/package/data/eml/edi/197/3/9eb6db370194bd3b2824726d89a008a6" 
infile1 <-  tempfile()
download.file(inUrl1,infile1,method="curl")

zoops <- read.csv(infile1, header=T) |>
  filter(CollectionMethod=="Tow" & Reservoir %in% c("BVR") &
           StartDepth_m > 7.1 & EndDepth_m == 0.1) |> 
  select(-c(Site,EndDepth_m,CollectionMethod))

#split data into pre 2019 and post
zoops_2016_2018 <- zoops[as.Date(zoops$DateTime)<"2019-01-01",]
zoops_2019_2021 <- zoops[as.Date(zoops$DateTime)>="2019-01-01",]

#-----------------------------------------------------------------------#
# IMPORTANT - There are n=7 total rotifer densities that are INCORRECT  #
# In this section, I'm manually updating these values, but note that    #
# we will eventually need to change these data in the next EDI pub      #
#-----------------------------------------------------------------------#
# 14May2014, 29May2014, 4Jun2014, 2Jul2014, 23Jul2014, 13Aug2014, 4Sep2014

zoops_2016_2018$Density_IndPerL[zoops_2016_2018$DateTime=="2014-05-14 10:40:00" & 
                                zoops_2016_2018$Taxon=="Total rotifers"] <- 66.19

zoops_2016_2018$Density_IndPerL[zoops_2016_2018$DateTime=="2014-05-29 11:10:00" & 
                                zoops_2016_2018$Taxon=="Total rotifers"] <- 20.71

zoops_2016_2018$Density_IndPerL[zoops_2016_2018$DateTime=="2014-06-04 12:00:00" & 
                                zoops_2016_2018$Taxon=="Total rotifers"] <- 48.35

zoops_2016_2018$Density_IndPerL[zoops_2016_2018$DateTime=="2014-07-02 10:50:00" & 
                                zoops_2016_2018$Taxon=="Total rotifers"] <- 11.33

zoops_2016_2018$Density_IndPerL[zoops_2016_2018$DateTime=="2014-07-23 11:45:00" & 
                                  zoops_2016_2018$Taxon=="Total rotifers"] <- 63.14

zoops_2016_2018$Density_IndPerL[zoops_2016_2018$DateTime=="2014-08-13 10:05:00" & 
                                  zoops_2016_2018$Taxon=="Total rotifers"] <- 23.91

zoops_2016_2018$Density_IndPerL[zoops_2016_2018$DateTime=="2014-09-04 09:45:00" & 
                                  zoops_2016_2018$Taxon=="Total rotifers"] <- 43.55

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
#write.csv(all_zoops, paste0("Output/all_zoops_dens.csv"),row.names = FALSE)

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
            sd = mean(sd[Taxon %in% c("Cladocera","Copeoda","Rotifera")],na.rm=T))  
  #ungroup() |> 
  #mutate(year = format(DateTime, "%Y"),
  #       month = format(DateTime, "%m")) |> 
  #group_by(year, month) |>
  #summarise(Total_avg = mean(Total,na.rm=T),
  #          Total_sd = mean(sd,na.rm=T)) |> 
  #ungroup() |> 
  #mutate(Dens_std =(Total_avg - min(Total_avg)) / 
  #         (max(Total_avg) - min(Total_avg)))

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

#write.csv(zoops_3_groups,"Output/std_dens_3taxa.csv", row.names = F)

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

#order zoops by group
zoops_10_groups$Taxon <- factor(zoops_10_groups$Taxon, 
                                levels=c("Bosmina", "Ceriodaphnia", "Daphnia",
                                         "Cyclopoida","Nauplii","Ascomorpha",
                                         "Conochilus","Kellicottia","Keratella",
                                         "Polyarthra"))

#write.csv(zoops_10_groups,"Output/std_dens_10taxa.csv", row.names = F)

#create new column w/ date bc leap years are really messing up my tick marks
zoops_10_groups$date <- as.Date(paste0(zoops_10_groups$year,"-",
                                       zoops_10_groups$month,"-01"))

#shaded line plot - raw density
ggplot(data = subset(zoops_10_groups, month %in% 
                       c("05","06","07","08","09")),
       aes(x=date, y = avg, color=Taxon)) +
  geom_area(aes(color = Taxon, fill = Taxon),
            position = "fill", stat = "identity", #position = stack
            alpha=0.7) +
  facet_wrap(~year, scales = "free_x")+
  scale_color_manual(values = NatParksPalettes::natparks.pals("DeathValley", 12, direction=-1)[c(1:3,5,6,8:12)])+
  scale_fill_manual(values = NatParksPalettes::natparks.pals("DeathValley", 12, direction=-1)[c(1:3,5,6,8:12)])+
  scale_x_date(expand = c(0,0)) +#,
 #              breaks = scales::pretty_breaks(5),
 #              date_labels = "%b") +
  scale_y_continuous(expand = c(0,0))+
  xlab("") + ylab("Relative density") +
  guides(color= "none",
         fill = guide_legend(ncol=1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        text = element_text(size=10), 
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.x = element_text(face = "bold",hjust = 0),
        axis.text.x = element_text(angle=90),
        strip.background.x = element_blank(),
        axis.title.y = element_text(size = 11),
        plot.margin = unit(c(0, 1, 0, 0), "cm"),
        panel.spacing.x = unit(0.2, "in"),
        panel.background = element_rect(
          fill = "white"),
        panel.spacing = unit(0.5, "lines"))
#ggsave("Figures/BVR_succession_10groups_fill_alldens_relative.jpg", width=7, height=4) 

#----------------------------------------------------------------#
#looking at each taxon + grouping years by trajectory 

#order by trajectory year
zoops_10_groups$year <- factor(zoops_10_groups$year, 
                               levels = c("2014","2019","2021", 
                                          "2015","2016", "2020"))

#add trajectory col
zoops_10_groups$traj <- ifelse(zoops_10_groups$year %in% 
                                 c("2014","2019","2021"), 
                               "clockwise", "counterclockwise")

ggplot(data = subset(zoops_10_groups, month %in% 
                       c("05","06","07","08","09") &
                       Taxon %in% c("Cyclopoida","Nauplii",
                                    "Keratella","Kellicottia")), 
       aes(as.Date("2019-12-31") + yday(as.Date(paste0(
       year,"-",month,"-01"))), avg, color=year)) +
  geom_area(aes(color = year, fill = year),
            position = "identity", stat = "identity", #position = stack
            alpha=0.5) +
  facet_wrap(~Taxon+traj, scales = "free", ncol=2,
             labeller = labeller(.multi_line = FALSE)) +
  scale_x_date(expand = c(0,0), labels = 
                       scales::date_format("%b",tz="EST5EDT"))+
  scale_y_continuous(expand = c(0,0))+
  scale_color_manual(values = c(
    "#003366","#0099CC","#339999","#660000","#CC0000","#CC6666"))+
  scale_fill_manual(values = c(
    "#003366","#0099CC","#339999","#660000","#CC0000","#CC6666"))+
  xlab("") + ylab("Density (#/L)") +
  guides(color= "none",
         fill = guide_legend(ncol=1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        text = element_text(size=10), 
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.x = element_text(face = "bold",hjust = 0),
        axis.text.x = element_text(angle=90),
        strip.background.x = element_blank(),
        axis.title.y = element_text(size = 11),
        plot.margin = unit(c(0.4, 1, 0, 0), "cm"),
        panel.background = element_rect(
          fill = "white"),
        panel.spacing = unit(0.5, "lines"))
#ggsave("Figures/BVR_succession_4taxa.jpg", width=6, height=7) 

#order taxa depending on year 
zoops_10_groups$Taxon <- factor(zoops_10_groups$Taxon, 
                                levels = c("Conochilus","Ascomorpha","Daphnia","Keratella", "Kellicottia","Polyarthra","Nauplii","Cyclopoida","Ceriodaphnia","Bosmina")) # 2014
                                #levels = c("Kellicottia","Polyarthra","Cyclopoida","Keratella", "Daphnia","Nauplii","Ascomorpha","Conochilus","Ceriodaphnia","Bosmina")) # 2015
                                #levels = c("Ascomorpha","Keratella","Cyclopoida", "Daphnia","Conochilus", "Ceriodaphnia","Polyarthra","Nauplii","Bosmina","Kellicottia")) # 2016
                                #levels = c("Daphnia", "Polyarthra","Conochilus", "Cyclopoida","Nauplii","Kellicottia","Keratella","Ceriodaphnia","Bosmina","Ascomorpha")) # 2019
                                #levels = c("Polyarthra","Conochilus","Nauplii", "Daphnia","Cyclopoida", "Keratella","Ceriodaphnia","Bosmina","Kellicottia","Ascomorpha")) # 2020/2021

#look at "succession" for each taxon
ggplot(data = subset(zoops_10_groups, month %in% 
                       c("05","06","07","08","09") &
                       year =="2014"),
       aes(date, standardized_dens,
         color=Taxon)) +
  geom_area(aes(color = Taxon, fill = Taxon),
            stat = "identity") +
  facet_wrap(~Taxon, strip.position="left",ncol=1) +
  scale_x_date(expand = c(0,0), labels = 
                 scales::date_format("%b",tz="EST5EDT"))+
  scale_y_continuous(expand = c(0,0)) +
  scale_color_manual(values = NatParksPalettes::
                       natparks.pals("DeathValley", 12, 
                                     direction=-1)[c(1:3,5,6,8:12)])+
  scale_fill_manual(values = NatParksPalettes::
                      natparks.pals("DeathValley", 12, 
                                    direction=-1)[c(1:3,5,6,8:12)])+
  xlab("") + ylab("") +
  guides(color= "none",
         fill = "none") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        text = element_text(size=10), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background.y = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle=0),
        plot.margin = unit(c(0.4, 1, 0, 0), "cm"),
        panel.background = element_rect(
          fill = "white"),
        panel.spacing = unit(0.5, "lines"))
#ggsave("Figures/BVR_succession_alltaxa_2014.jpg", width=6, height=7) 



#order years
zoops_10_groups$year <- factor(zoops_10_groups$year, levels = c( 
                                 "2014","2019","2021", 
                                 "2015","2016", "2020"))

#line plots - colors = trajectories
ggplot(data = subset(zoops_10_groups, month %in% 
                       c("05","06","07","08","09") &
                       Taxon %in% "Bosmina"), 
       aes(as.Date("2019-12-31") + yday(as.Date(paste0(
         year,"-",month,"-01"))), avg, color=year)) +
  geom_point(cex=3) + geom_line(lwd=1) + 
  scale_x_date(expand = c(0,0), labels = 
                 scales::date_format("%b",tz="EST5EDT"))+
  scale_color_manual(values = c(rep("#01586D",3),rep("#8B0C13",3)))+
  scale_fill_manual(values = c(rep("#01586D",3),rep("#8B0C13",3)))+
  xlab("") + ylab("Density (#/L)") + 
  guides(color= "none",
         fill = guide_legend(ncol=1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        text = element_text(size=10), 
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.x = element_text(face = "bold",hjust = 0),
        axis.text.x = element_text(angle=90),
        strip.background.x = element_blank(),
        axis.title.y = element_text(size = 11),
        plot.margin = unit(c(0.4, 1, 0, 0), "cm"),
        panel.background = element_rect(
          fill = "white"),
        panel.spacing = unit(0.5, "lines"))
#ggsave("Figures/BVR_succession_lines_raw_dens_bosmina.jpg", width=7, height=4) 

#----------------------------------------------------------------#
#zooming in on clads - raw density
ggplot(data = subset(zoops_10_groups, Taxon %in% c("Bosmina","Ceriodaphnia","Daphnia")), 
       aes(as.Date("2019-12-31") + 
             yday(as.Date(paste0(year,"-",month,"-01"))), avg, color=Taxon)) +
  geom_area(aes(color = Taxon, fill = Taxon),
            position = "identity", stat = "identity",
            alpha=0.7) +
  facet_wrap(~year, scales = "free")+
  scale_color_manual(values = NatParksPalettes::natparks.pals("Glacier", 3, direction=-1))+
  scale_fill_manual(values = NatParksPalettes::natparks.pals("Glacier", 3, direction=-1))+
  scale_x_date(expand = c(0,0),
               labels = scales::date_format("%b",tz="EST5EDT")) +
  scale_y_continuous(expand = c(0,0))+
  xlab("") + ylab("Density (#/L)") +
  guides(color= "none",
         fill = guide_legend(ncol=1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        text = element_text(size=10), 
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.x = element_text(face = "bold",hjust = 0),
        axis.text.x = element_text(angle=90),
        strip.background.x = element_blank(),
        axis.title.y = element_text(size = 11),
        plot.margin = unit(c(0, 1, 0, 0), "cm"),
        panel.background = element_rect(
          fill = "white"),
        panel.spacing = unit(0.5, "lines"))
#ggsave("Figures/BVR_succession_clads_alldens_raw.jpg", width=7, height=4) 

#add functional groups 
zoops_10_groups$trophi <- ifelse(zoops_10_groups$Taxon %in% c("Ascomorpha", "Polyarthra"), "virgate",
                           ifelse(zoops_10_groups$Taxon %in% c("Keratella", "Kellicottia"), "malleate",
                                  ifelse(zoops_10_groups$Taxon %in% c("Conochilus"), "malleoramate", "NA")))
#virgate = raptorial, attached rotifers, incudate = omnivorous/predacious feeding behavior, 
#malleate = teeth for chewing and grasping prey, malleoramate = similar to malleate but more strongly toothed

zoops_10_groups$feeding_type <- ifelse(zoops_10_groups$Taxon %in% c("Cyclopoida"), "predators",
                                 ifelse(zoops_10_groups$Taxon %in% c("Ascomorpha", "Polyarthra"), "suctors",
                                        "filter feeders"))
#"Daphnia","Ceriodaphnia","Bosmina", "Nauplius", "Keratella", "Kellicottia"

zoops_10_groups$trophic_group <- ifelse(zoops_10_groups$Taxon %in% c("Cyclopoida"), "Omnicarnivore",
                                  ifelse(zoops_10_groups$Taxon %in% c(
                                    "Ceriodaphnia","Daphnia","Bosmina", "Nauplii"), "Herbivore",
                                    ifelse(zoops_10_groups$Taxon %in% c("Ascomorpha", "Polyarthra"), "Raptorial",
                                           "Microphagous"))) #conochilus, keratella, kellicottia

#order years again
zoops_10_groups$year <- factor(zoops_10_groups$year, levels = c( 
  "2014","2015","2016", 
  "2019","2020", "2021"))

#functional groups
zoops_10_groups |> group_by(year, month) |> 
  summarise(filter_feeders = sum(avg[feeding_type == "filter feeders"]),
            predators = sum(avg[feeding_type == "predators"]),
            suctors = sum(avg[feeding_type == "suctors"])) |> 
  pivot_longer(-c(year,month),
               names_to = "feeding_types")|> 
  ungroup() |> 
ggplot(aes(as.Date("2019-12-31") + 
             yday(as.Date(paste0(year,"-",month,"-01"))), 
             value, color=feeding_types)) +
  geom_area(aes(color = feeding_types, fill = feeding_types),
            position = "identity",# stat = "identity",
            alpha=0.7) +
  facet_wrap(~year, scales = "free")+
  scale_color_manual(values = NatParksPalettes::natparks.pals("Halekala", 3))+
  scale_fill_manual(values = NatParksPalettes::natparks.pals("Halekala", 3))+
  scale_x_date(expand = c(0,0),
               labels = scales::date_format("%b",tz="EST5EDT")) +
  scale_y_continuous(expand = c(0,0))+
  xlab("") + ylab("Density (#/L)") +
  guides(color= "none",
         fill = guide_legend(ncol=1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        text = element_text(size=10), 
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.x = element_text(face = "bold",hjust = 0),
        axis.text.x = element_text(angle=90),
        strip.background.x = element_blank(),
        axis.title.y = element_text(size = 11),
        plot.margin = unit(c(0, 1, 0, 0), "cm"),
        panel.background = element_rect(
          fill = "white"),
        panel.spacing = unit(0.5, "lines"))
#ggsave("Figures/BVR_succession_feeding_types_alldens_raw.jpg", width=7, height=4) 

zoops_10_groups |> group_by(year, month) |> 
  summarise(virgate = sum(avg[trophi == "virgate"]),
            malleoramate = sum(avg[trophi == "malleoramate"]),
            malleate = sum(avg[trophi == "malleate"]),
            none = sum(avg[trophi == "NA"])) |> 
  pivot_longer(-c(year,month),
               names_to = "trophi")|> 
  ungroup() |> 
  ggplot(aes(as.Date("2019-12-31") + 
               yday(as.Date(paste0(year,"-",month,"-01"))), 
             value, color=trophi)) +
  geom_area(aes(color = trophi, fill = trophi),
            position = "identity",# stat = "identity",
            alpha=0.7) +
  facet_wrap(~year, scales = "free")+
  scale_color_manual(values = NatParksPalettes::natparks.pals("Halekala", 4))+
  scale_fill_manual(values = NatParksPalettes::natparks.pals("Halekala", 4))+
  scale_x_date(expand = c(0,0),
               labels = scales::date_format("%b",tz="EST5EDT")) +
  scale_y_continuous(expand = c(0,0))+
  xlab("") + ylab("Density (#/L)") +
  guides(color= "none",
         fill = guide_legend(ncol=1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        text = element_text(size=10), 
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.x = element_text(face = "bold",hjust = 0),
        axis.text.x = element_text(angle=90),
        strip.background.x = element_blank(),
        axis.title.y = element_text(size = 11),
        plot.margin = unit(c(0, 1, 0, 0), "cm"),
        panel.background = element_rect(
          fill = "white"),
        panel.spacing = unit(0.5, "lines"))
#ggsave("Figures/BVR_succession_trophi_alldens_raw.jpg", width=7, height=4) 

zoops_10_groups |> group_by(year, month) |> 
  summarise(Herbivore = sum(avg[trophic_group == "Herbivore"]),
            Omnicarnivore = sum(avg[trophic_group == "Omnicarnivore"]),
            Microphagous = sum(avg[trophic_group == "Microphagous"]),
            Raptorial = sum(avg[trophic_group == "Raptorial"])) |> 
  pivot_longer(-c(year,month),
               names_to = "trophic_group")|> 
  ungroup() |> 
  ggplot(aes(as.Date("2019-12-31") + 
               yday(as.Date(paste0(year,"-",month,"-01"))), 
             value, color=trophic_group)) +
  geom_area(aes(color = trophic_group, fill = trophic_group),
            position = "identity",# stat = "identity",
            alpha=0.7) +
  facet_wrap(~year, scales = "free")+
  scale_color_manual(values = NatParksPalettes::natparks.pals("Halekala", 4))+
  scale_fill_manual(values = NatParksPalettes::natparks.pals("Halekala", 4))+
  scale_x_date(expand = c(0,0),
               labels = scales::date_format("%b",tz="EST5EDT")) +
  scale_y_continuous(expand = c(0,0))+
  xlab("") + ylab("Density (#/L)") +
  guides(color= "none",
         fill = guide_legend(ncol=1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        text = element_text(size=10), 
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.x = element_text(face = "bold",hjust = 0),
        axis.text.x = element_text(angle=90),
        strip.background.x = element_blank(),
        axis.title.y = element_text(size = 11),
        plot.margin = unit(c(0, 1, 0, 0), "cm"),
        panel.background = element_rect(
          fill = "white"),
        panel.spacing = unit(0.5, "lines"))
#ggsave("Figures/BVR_succession_trophic_group_alldens_raw.jpg", width=7, height=4) 

#-----------------------------------------------------------------------------#
#shaded line plot for clads, copes, and rots each year

# 3 group zoop dens - standardize within a year including all data (not just averaged by month)
zoops_3_groups_years <- all_zoops |> 
  filter(Taxon %in% c("Cladocera","Copepoda","Rotifera")) |> 
  group_by(DateTime) |>
  mutate(month = format(DateTime, "%m")) |> 
  filter(month %in% c("05","06","07","08","09")) |> 
  summarise(Cladocera_avg = mean(dens[Taxon=="Cladocera"]),
            Cladocera_sd = mean(sd[Taxon=="Cladocera"],na.rm=T),
            Copepoda_avg = mean(dens[Taxon=="Copepoda"]),
            Copepoda_sd = mean(sd[Taxon=="Copepoda"],na.rm=T),
            Rotifera_avg = mean(dens[Taxon=="Rotifera"]),
            Rotifera_sd = mean(sd[Taxon=="Rotifera"],na.rm=T)) |> 
  mutate(C_R = sum(Cladocera_avg,Copepoda_avg, na.rm=T) / Rotifera_avg) |>
  mutate(percent_clad = Cladocera_avg/ (Cladocera_avg + Copepoda_avg + Rotifera_avg) * 100,
         percent_cope = Copepoda_avg/ (Cladocera_avg + Copepoda_avg + Rotifera_avg) * 100,
         percent_rot = Rotifera_avg/ (Cladocera_avg + Copepoda_avg + Rotifera_avg) * 100) |> 
  #higher C:R means more crustaceans relative to rotifers
  mutate(year = format(DateTime, "%Y")) |> 
  pivot_longer(Cladocera_avg:Rotifera_sd,
               names_to = c("Taxon", ".value"),
               names_sep="_" )  |> 
  ungroup() |> group_by(Taxon,year) |>
  mutate(min_dens = min(avg),
         max_dens = max(avg)) |> 
  mutate(standardized_dens = (avg - min_dens) / (max_dens - min_dens)) |> 
  ungroup()

#add a month column
zoops_3_groups_years$month <- format(zoops_3_groups_years$DateTime,"%m") 

zoops_3_groups_years <- zoops_3_groups_years |> 
  group_by(month, year) |> mutate(C_R_sd = sd(C_R))

#plot crustacean:rotifer over time
ggplot(zoops_3_groups_years,aes(as.Date("2019-12-31") + 
      yday(as.Date(DateTime)), C_R, color=year)) +
  geom_point() + geom_line() + xlab("") +
  scale_color_manual(values = NatParksPalettes::
                       natparks.pals("Triglav", 6))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        text = element_text(size=10), 
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.x = element_text(face = "bold",hjust = 0),
        axis.text.x = element_text(angle=90),
        strip.background.x = element_blank(),
        axis.title.y = element_text(size = 11),
        plot.margin = unit(c(0, 1, 0, 0), "cm"),
        panel.background = element_rect(
          fill = "white"),
        panel.spacing = unit(0.5, "lines"))
#ggsave("Figures/BVR_yearly_crustacean_rotifer_ratio.jpg", width=7, height=4) 

#C:R boxplots
zoops_3_groups_years |> group_by(month,year) |> 
  summarise(C_R = mean(C_R),
            C_R_sd = sd(C_R_sd)) |> 
ggplot(aes(x=as.factor(month), y=C_R, fill=year)) + xlab("") +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = C_R - C_R_sd, 
                    ymax = C_R + C_R_sd),
    width=0.3, alpha=0.9, size=1.3,
    position = position_dodge(0.9)) +
  scale_x_discrete(labels= c("May","June","July",
                             "August","September")) +
  scale_fill_manual(values = c("#01586D", "#8B0C13", "#8B0C13",
                               "#01586D", "#8B0C13", "#01586D")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        text = element_text(size=10), 
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.x = element_text(face = "bold",hjust = 0),
        strip.background.x = element_blank(),
        axis.title.y = element_text(size = 11),
        plot.margin = unit(c(0, 1, 0, 0), "cm"),
        panel.background = element_rect(
          fill = "white"),
        panel.spacing = unit(0.5, "lines"))
#ggsave("Figures/BVR_barplot_crustacean_rotifer_ratio_months.jpg", width=7, height=4) 

#percent barplots
zoops_3_groups_years |> group_by(month,year) |> 
  summarise(cladocera_mean = mean(percent_clad),
            cladocera_sd = sd(percent_clad),
            copepoda_mean = mean(percent_cope),
            copepoda_sd = sd(percent_cope),
            rotifera_mean = mean(percent_rot),
            rotifera_sd = sd(percent_rot)) |> 
  pivot_longer(-c(month,year),
               names_to = c("Taxon", ".value"),
               names_sep="_" )  |> 
  ggplot(aes(x=as.factor(month), y=mean, fill=year)) + xlab("") +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = mean - sd, 
                    ymax = mean + sd),
                width=0.3, alpha=0.9, size=1.3,
                position = position_dodge(0.9)) +
  facet_wrap(~Taxon, nrow=3) +
  scale_x_discrete(labels= c("May","June","July",
                             "August","September")) +
  scale_fill_manual(values = c("#01586D", "#8B0C13", "#8B0C13",
                               "#01586D", "#8B0C13", "#01586D")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        text = element_text(size=10), 
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.x = element_text(face = "bold",hjust = 0),
        strip.background.x = element_blank(),
        axis.title.y = element_text(size = 11),
        plot.margin = unit(c(0, 1, 0, 0), "cm"),
        panel.background = element_rect(
          fill = "white"),
        panel.spacing = unit(0.5, "lines"))
#ggsave("Figures/BVR_barplot_percent_taxon_month.jpg", width=7, height=4) 

#percent barplots
zoops_3_groups_years |> group_by(month) |> 
  summarise(cladocera_mean = mean(percent_clad),
            cladocera_sd = sd(percent_clad),
            copepoda_mean = mean(percent_cope),
            copepoda_sd = sd(percent_cope),
            rotifera_mean = mean(percent_rot),
            rotifera_sd = sd(percent_rot)) |> 
  pivot_longer(-c(month),
               names_to = c("Taxon", ".value"),
               names_sep="_" )  |> 
  ggplot(aes(x=as.factor(month), y=mean, fill=Taxon)) + xlab("") +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = mean - sd, 
                    ymax = mean + sd),
                width=0.3, alpha=0.9, size=1.3,
                position = position_dodge(0.9)) +
  scale_x_discrete(labels= c("May","June","July",
                             "August","September")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right",
        text = element_text(size=10), 
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.x = element_text(face = "bold",hjust = 0),
        strip.background.x = element_blank(),
        axis.title.y = element_text(size = 11),
        plot.margin = unit(c(0, 1, 0, 0), "cm"),
        panel.background = element_rect(
          fill = "white"),
        panel.spacing = unit(0.5, "lines"))
#ggsave("Figures/BVR_barplot_percent_taxon_vs_month.jpg", width=7, height=4) 

#numbers
mean(zoops_3_groups_years$percent_rot[zoops_3_groups_years$month=="05"])
mean(zoops_3_groups_years$percent_rot[zoops_3_groups_years$month=="06"])
mean(zoops_3_groups_years$percent_rot[zoops_3_groups_years$month=="07"])
mean(zoops_3_groups_years$percent_rot[zoops_3_groups_years$month=="08"])
mean(zoops_3_groups_years$percent_rot[zoops_3_groups_years$month=="09"])

mean(zoops_3_groups_years$percent_clad[zoops_3_groups_years$month=="05"])
mean(zoops_3_groups_years$percent_clad[zoops_3_groups_years$month=="06"])
mean(zoops_3_groups_years$percent_clad[zoops_3_groups_years$month=="07"])
mean(zoops_3_groups_years$percent_clad[zoops_3_groups_years$month=="08"])
mean(zoops_3_groups_years$percent_clad[zoops_3_groups_years$month=="09"])

mean(zoops_3_groups_years$percent_cope[zoops_3_groups_years$month=="05"])
mean(zoops_3_groups_years$percent_cope[zoops_3_groups_years$month=="06"])
mean(zoops_3_groups_years$percent_cope[zoops_3_groups_years$month=="07"])
mean(zoops_3_groups_years$percent_cope[zoops_3_groups_years$month=="08"])
mean(zoops_3_groups_years$percent_cope[zoops_3_groups_years$month=="09"])

#percent barplots
zoops_3_groups_years |> group_by(year) |> 
  summarise(cladocera_mean = mean(percent_clad),
            cladocera_sd = sd(percent_clad),
            copepoda_mean = mean(percent_cope),
            copepoda_sd = sd(percent_cope),
            rotifera_mean = mean(percent_rot),
            rotifera_sd = sd(percent_rot)) |> 
  pivot_longer(-c(year),
               names_to = c("Taxon", ".value"),
               names_sep="_" )  |> 
  ggplot(aes(x=as.factor(year), y=mean, fill=Taxon)) + xlab("") +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = mean - sd, 
                    ymax = mean + sd),
                width=0.3, alpha=0.9, size=1.3,
                position = position_dodge(0.9)) +
  scale_x_discrete(labels= c("2014","2015","2016",
                             "2019","2020","2021")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right",
        text = element_text(size=10), 
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.x = element_text(face = "bold",hjust = 0),
        strip.background.x = element_blank(),
        axis.title.y = element_text(size = 11),
        plot.margin = unit(c(0, 1, 0, 0), "cm"),
        panel.background = element_rect(
          fill = "white"),
        panel.spacing = unit(0.5, "lines"))
#ggsave("Figures/BVR_barplot_percent_taxon_vs_year.jpg", width=7, height=4) 

#C:R barplots - years
zoops_3_groups_years |> group_by(year) |> 
  summarise(C_R = mean(C_R),
            C_R_sd = sd(C_R_sd)) |> 
  ggplot(aes(x=as.factor(year), y=C_R, fill=year)) + xlab("") +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = C_R - C_R_sd, 
                    ymax = C_R + C_R_sd),
                width=0.3, alpha=0.9, size=1.3,
                position = position_dodge(0.9)) +
  scale_x_discrete(labels= c("2014","2015","2016",
                             "2019","2020","2021")) +
  scale_fill_manual(values = c("#01586D", "#8B0C13", "#8B0C13",
                               "#01586D", "#8B0C13", "#01586D")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        text = element_text(size=10), 
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.x = element_text(face = "bold",hjust = 0),
        strip.background.x = element_blank(),
        axis.title.y = element_text(size = 11),
        plot.margin = unit(c(0, 1, 0, 0), "cm"),
        panel.background = element_rect(
          fill = "white"),
        panel.spacing = unit(0.5, "lines"))
#ggsave("Figures/BVR_barplot_crustacean_rotifer_ratio_year.jpg", width=7, height=4) 

#C:R boxplots - months
zoops_3_groups_years |> group_by(month) |> 
  summarise(C_R = mean(C_R),
            C_R_sd = sd(C_R_sd)) |> 
  ggplot(aes(x=as.factor(month), y=C_R)) + xlab("") +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = C_R - C_R_sd, 
                    ymax = C_R + C_R_sd),
                width=0.3, alpha=0.9, size=1.3,
                position = position_dodge(0.9)) +
  scale_x_discrete(labels= c("May","June","July",
                             "August","September")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        text = element_text(size=10), 
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.x = element_text(face = "bold",hjust = 0),
        strip.background.x = element_blank(),
        axis.title.y = element_text(size = 11),
        plot.margin = unit(c(0, 1, 0, 0), "cm"),
        panel.background = element_rect(
          fill = "white"),
        panel.spacing = unit(0.5, "lines"))
#ggsave("Figures/BVR_barplot_crustacean_rotifer_ratio_month.jpg", width=7, height=4) 


mean(zoops_3_groups_years$C_R[zoops_3_groups_years$year %in% c(2014,2019,2021)]) #C
mean(zoops_3_groups_years$C_R[zoops_3_groups_years$year %in% c(2015,2016,2020)]) #CC
#CC has higher mean + median C:R than C

#shaded line plot time
ggplot(zoops_3_groups_years, aes(as.Date("2019-12-31") + yday(as.Date(DateTime)), 
           standardized_dens, color=Taxon)) +
  geom_area(aes(color = Taxon, fill = Taxon),
            position = "stack", stat = "identity",
            alpha=0.7) +
  facet_wrap(~year, scales = "free")+
  scale_color_manual(values = NatParksPalettes::natparks.pals("KingsCanyon", 3, direction=-1))+
  scale_fill_manual(values = NatParksPalettes::natparks.pals("KingsCanyon", 3, direction=-1),
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
        text = element_text(size=10), 
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.x = element_text(face = "bold",hjust = 0),
        axis.text.x = element_text(angle=90),
        strip.background.x = element_blank(),
        axis.title.y = element_text(size = 11),
        plot.margin = unit(c(0, 1, 0, 0), "cm"),
        panel.background = element_rect(
          fill = "white"),
        panel.spacing = unit(0.5, "lines"))
#ggsave("Figures/BVR_succession_3groups_stacked_alldens_std.jpg", width=7, height=4) 

#shaded line plot time - raw density
ggplot(zoops_3_groups_years, aes(as.Date("2019-12-31") + yday(as.Date(DateTime)), 
                                 avg, color=Taxon)) +
  geom_area(aes(color = Taxon, fill = Taxon),
            position = "stack", stat = "identity",
            alpha=0.7) +
  facet_wrap(~year, scales = "free")+
  scale_color_manual(values = NatParksPalettes::natparks.pals("KingsCanyon", 3, direction=-1))+
  scale_fill_manual(values = NatParksPalettes::natparks.pals("KingsCanyon", 3, direction=-1),
                    breaks = c("Cladocera","Copepoda","Rotifera"))+
  scale_x_date(expand = c(0,0),
               labels = scales::date_format("%b",tz="EST5EDT")) +
  scale_y_continuous(expand = c(0,0))+
  xlab("") + ylab("Density (#/L)") +
  guides(color= "none",
         fill = guide_legend(ncol=1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        text = element_text(size=10), 
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.x = element_text(face = "bold",hjust = 0),
        axis.text.x = element_text(angle=90),
        strip.background.x = element_blank(),
        axis.title.y = element_text(size = 11),
        plot.margin = unit(c(0, 1, 0, 0), "cm"),
        panel.background = element_rect(
          fill = "white"),
        panel.spacing = unit(0.5, "lines"))
#ggsave("Figures/BVR_succession_3groups_stacked_alldens_raw.jpg", width=7, height=4) 

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

#ggsave("Figures/phyto_zoop_annual_peg_plus_secchi.jpg", width=6, height=4)

