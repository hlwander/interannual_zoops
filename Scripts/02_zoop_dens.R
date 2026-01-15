#zoop density succession figs

#read in packages
pacman::p_load(zoo, dplR, dplyr, tidyverse, ggplot2, ggpubr, lubridate, ggtext)

#cb friendly year palette
year_cols <- c("#011f51","#1f78b4","#33a02c","#fdfa66","#ff7f00","#e31a1c","#6a3d9a")

#read in zoop data from EDI (v5)
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/197/5/38fc9d1a4c8b6976c71e56bda5ff073b" 
infile1 <-  tempfile()
try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

zoops <- read.csv(infile1, header=T) |>
  filter(CollectionMethod=="Tow" & Reservoir %in% c("BVR") &
           StartDepth_m > 7.1 & EndDepth_m == 0.1) |> 
  select(-c(Site,EndDepth_m,CollectionMethod)) |>
  mutate(Taxon = ifelse(Taxon == "nauplius", "Nauplii", Taxon),
         DateTime = as.POSIXct(DateTime)) |>
  filter(!year(DateTime) %in% c("2019","2020","2021","2022","2023","2024","2025")) 
#recreated summary file for 2019-2025 and it is not yet published on EDI 

zoops_2019_2021 <- read.csv("inputs/EDI_zoop_taxa_2019-2022.csv") |>
  filter(CollectionMethod=="Tow" & Reservoir %in% c("BVR") &
           StartDepth_m > 7.1 & EndDepth_m == 0 & Site==50) |>
  select(-c(Site,EndDepth_m,CollectionMethod)) |>
  mutate(Taxon = ifelse(Taxon == "nauplius", "Nauplii", Taxon),
         DateTime = as.POSIXct(DateTime, format = "%Y-%m-%d %H:%M:%S", tz="EST")) |>
  filter(!year(DateTime) %in% c("2022")) |>
  distinct()

zoops_2023 <- read.csv("inputs/EDI_zoop_summary_2022-2025.csv") |>
  filter(CollectionMethod=="Tow" & Reservoir %in% c("BVR") &
           StartDepth_m > 7.1 & EndDepth_m == 0 & Site ==50) |>
  select(-c(Site,EndDepth_m,CollectionMethod)) |>
  mutate(Taxon = ifelse(Taxon == "nauplius", "Nauplii", Taxon),
         DateTime = as.POSIXct(DateTime, format = "%Y-%m-%d %H:%M:%S", tz="EST")) |>
  filter(year(DateTime) %in% c("2023"))
  
#add in 2022-2025 zoop data
zoops <- bind_rows(zoops, zoops_2019_2021, zoops_2023) |>
  filter(month(DateTime) %in% 4:11)

#------------------------------------------------------------------------------#
#figure out dominant taxa for NMDS/other ms figs
zoop_taxa_props <- zoops |>
  filter(!is.na(Taxon), !is.na(Density_IndPerL),
         !Taxon %in% c("Cladocera","Copepoda","Rotifera","Crustacea")) |>
  group_by(Taxon) |>
  mutate(n_days = n_distinct(DateTime[Density_IndPerL > 0])) |>
  filter(n_days > 10) |> #keep taxa with > 10 observation days
  summarise(total_dens = sum(Density_IndPerL, na.rm = TRUE)) |>
  ungroup() |>
  mutate(prop = total_dens / sum(total_dens)) |>
  arrange(desc(prop))

# Note that this fig corresponds to Table S1 of the manuscript
#plot proportion of density that each taxon makes up (sum of all individual days)
ggplot(zoop_taxa_props, aes(x=reorder(Taxon,-prop), y=prop)) +
  theme_bw() + geom_bar(stat="identity") + 
  theme(text = element_text(size=5), 
        axis.text = element_text(size=5, color="black"), 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, hjust=1), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#new df to average proportion for each taxon (Supplemental Table S1)
zoop_prop_table <- zoop_taxa_props |> group_by(Taxon) |> 
  summarise(prop_avg = mean(prop)) |>
  filter(prop_avg >= 0.01) |>
  filter(!Taxon %in% c("Conochiloides","Gastropus")) #bc NA for 2016-2018
#drop taxa that are <1% (0.01) of density 
#write.csv(zoop_prop_table,"Output/dominant_taxa_props.csv", row.names=F)

#just keeping n=8
#consider grouping as Gastropidae (Gastropus and Ascomorpha) and Conochilidae (Conochilus and Conochiloides)
#just hard to know if astropus and conochiloides were really absent from 2014-2016 samples or if they were ided differently

#---------------------------------
#split data into pre 2019 and post
zoops_2014_2016 <- zoops[as.Date(zoops$DateTime)<"2019-01-01",] |> 
  filter(!if_all(everything(), is.na))
zoops_2019_2023 <- zoops[as.Date(zoops$DateTime)>="2019-01-01",] |> 
  filter(!if_all(everything(), is.na))

#add daphnia (D. catawba, D. ambigua), calanoida (diaptomus) groups for 2014-2016 data
zoops_pre <- zoops_2014_2016 |> 
  group_by(Reservoir, DateTime, StartDepth_m) |> 
  summarise(Daphnia_dens = sum(Density_IndPerL[
    Taxon %in% c("D. catawba","D. ambigua")]),
    Daphnia_sd = sd(Density_IndPerL[
      Taxon %in% c("D. catawba","D. ambigua")]),
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
    Conochilus_dens = sum(Density_IndPerL[
      Taxon %in% c("Conochilus")]),
    Conochilus_sd = sd(Density_IndPerL[
      Taxon %in% c("Conochilus")]),
    Keratella_dens = sum(Density_IndPerL[
      Taxon %in% c("Keratella")]),
    Keratella_sd = sd(Density_IndPerL[
      Taxon %in% c("Keratella")]),
    Kellicottia_dens = sum(Density_IndPerL[
      Taxon %in% c("Kellicottia")]),
    Kellicottia_sd = sd(Density_IndPerL[
      Taxon %in% c("Kellicottia")]),
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
  filter(hour(DateTime) %in% c(8,9,10,11,12,13,14)) |> #drop nighttime samples
  mutate(DateTime = as.Date(DateTime)) |> 
  mutate(dens = dens * (1/0.031)) |>  #10m bvr neteff from 2016 (n=2) - note that 7m neteff was also 0.031
#avg from 2020 and 2021 is 0.021 for reference
  mutate(DateTime = as.Date(ifelse(DateTime %in% as.Date("2014-09-25"), 
                           as.Date("2014-10-23"), DateTime))) #choosing 10-23 bc we have ctd profiles that day
#adding 6 days to this 2014 sep point to add another month to the SS NMDS

#list dominant taxa (+3 overall zoop groups)
taxa <- c(zoop_prop_table$Taxon,"Cladocera", "Copepoda","Rotifera") 

#average reps when appropriate
zoops_post <- zoops_2019_2023 |> 
  #mutate(DateTime = as.POSIXct(DateTime, format="%Y-%m-%d %H:%M:%S", tz="UTC")) |> 
  filter(hour(DateTime) %in% c(9,10,11,12,13,14) &  #drop nighttime samples
          !year(DateTime) %in% c(2022, 2024, 2025)) |> #we lost the 17 june 2024 sample somewhere between the field and the lab so that means this whole year needs to be excluded :(
  filter(Taxon %in% c(taxa)) |> 
  mutate(DateTime = as.Date(DateTime)) |> 
  group_by(Reservoir, DateTime, StartDepth_m, Taxon) |> 
  summarise(dens = mean(Density_IndPerL),
            sd = sd(Density_IndPerL), .groups = "drop")

#combine all zoop data
all_zoops <- bind_rows(zoops_pre, zoops_post) |> 
  mutate_all(~replace(., is.nan(.), NA)) |>  #replace NAN with NA
  ungroup() |> select(-StartDepth_m)  #dropping, but note that depths range from 8-11.5m....

#write.csv(all_zoops, paste0("Output/all_zoops_dens.csv"),row.names = FALSE)

#------------------------------------------------------------------------------#
# 8 group zoop dens
zoops_8_groups <- all_zoops |> 
  filter(!Taxon %in% c("Cladocera","Copepoda","Rotifera")) |> 
  group_by(DateTime) |> 
  summarise(Daphnia_avg = mean(dens[Taxon=="Daphnia"], na.rm=T),
            Daphnia_sd = mean(sd[Taxon=="Daphnia"],na.rm=T),
            Cyclopoida_avg = mean(dens[Taxon=="Cyclopoida"], na.rm=T),
            Cyclopoida_sd = mean(sd[Taxon=="Cyclopoida"],na.rm=T),
            Nauplii_avg = mean(dens[Taxon=="Nauplii"], na.rm=T),
            Nauplii_sd = mean(sd[Taxon=="Nauplii"],na.rm=T),
            Bosmina_avg = mean(dens[Taxon=="Bosmina"], na.rm=T),
            Bosmina_sd = mean(sd[Taxon=="Bosmina"],na.rm=T),
            Conochilus_avg = mean(dens[Taxon=="Conochilus"], na.rm=T),
            Conochilus_sd = mean(sd[Taxon=="Conochilus"],na.rm=T),
            Keratella_avg = mean(dens[Taxon=="Keratella"], na.rm=T),
            Keratella_sd = mean(sd[Taxon=="Keratella"],na.rm=T),
            Kellicottia_avg = mean(dens[Taxon=="Kellicottia"], na.rm=T),
            Kellicottia_sd = mean(sd[Taxon=="Kellicottia"],na.rm=T),
            Polyarthra_avg = mean(dens[Taxon=="Polyarthra"], na.rm=T),
            Polyarthra_sd = mean(sd[Taxon=="Polyarthra"],na.rm=T)) |> 
  pivot_longer(-c(DateTime),
               names_to = c("Taxon", ".value"),
               names_sep="_" )  |> 
  mutate(total = sum(avg, na.rm=T)) |> ungroup() |> 
  mutate(p_dens = avg / total,
  year = year(DateTime),
         month = month(DateTime),
         day   = day(DateTime),
         pseudoDate = as.Date(sprintf("2000-%02d-%02d", month, day))) |>
  select(-day) 

#order zoops by group
zoops_8_groups$Taxon <- factor(zoops_8_groups$Taxon, 
                                levels=c("Bosmina","Daphnia",
                                         "Cyclopoida","Nauplii",
                                         "Conochilus","Kellicottia",
                                         "Keratella","Polyarthra")) 

taxon_labels <- c(expression(italic("Bosmina")),
  expression(italic("Daphnia")),"Cyclopoida","Nauplii",
  expression(italic("Conochilus")),
  expression(italic("Kellicottia")),
  expression(italic("Keratella")),
  expression(italic("Polyarthra")))

#shaded line plot - raw density (Manuscript Figure 1)
ggplot(data = zoops_8_groups|> filter(month %in% c(4:11)),
       aes(x=pseudoDate, y = avg, color=Taxon)) +
  geom_area(aes(color = Taxon, fill = Taxon),
            position = "stack", stat = "identity") +
  facet_wrap(~year, scales = "free")+
  scale_color_manual(values = NatParksPalettes::natparks.pals(
    "DeathValley", 14, direction=-1)[c(1,3,5,7,10,11,13,14)],
    labels = taxon_labels) +
  scale_fill_manual(values = NatParksPalettes::natparks.pals(
    "DeathValley", 14, direction=-1)[c(1,3,5,7,10,11,13,14)],
    labels = taxon_labels) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b",
               expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  xlab("") + ylab(expression("Density (individuals L"^{-1}*")")) +
  guides(color= "none", fill = guide_legend(nrow=4)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.72,0.1),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        text = element_text(size=9), 
        axis.text.y = element_text(size = 9),
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
#ggsave("Figures/BVR_8groups_fill_alldens.jpg", width=5, height=4) 

#shaded line plot - relative density (Fig. S1)
ggplot(data = subset(zoops_8_groups, month %in% c(4:11)),
       aes(x=pseudoDate, y = avg, color=Taxon)) +
  geom_area(aes(color = Taxon, fill = Taxon),
            position = "fill", stat = "identity") +
  facet_wrap(~year(DateTime), scales = "free_x")+
  scale_color_manual(values = NatParksPalettes::natparks.pals(
    "DeathValley", 14, direction=-1)[c(1,3,5,7,10,11,13,14)],
    labels = taxon_labels) +
  scale_fill_manual(values = NatParksPalettes::natparks.pals(
    "DeathValley", 14, direction=-1)[c(1,3,5,7,10,11,13,14)],
    labels = taxon_labels) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b",
               expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  xlab("") + ylab("Relative density") +
  guides(color= "none", fill = guide_legend(nrow=4)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.72,0.1),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        text = element_text(size=9), 
        axis.text.y = element_text(size = 9),
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
#ggsave("Figures/BVR_succession_8groups_fill_alldens_relative.jpg", width=5, height=4) 

#----------------------
#numbers for ms results
zoops_8_groups_max <- zoops_8_groups |>
  group_by(Taxon, year) |>
  slice_max(order_by = avg, n = 1, with_ties = FALSE) |>
  ungroup()

#bosmina mostly peak between aug-nov (except 2014; april)
#daphnia peak between may-oct
#cyclopoid peak between jun-nov
#nauplii peak between jun-sep
#conochilus peak between may-oct
#kellicottia peak between apr-aug
#keratella peak between apr-aug
#polyarthra peak between apr-oct

mean(c(zoops_8_groups_max$month[zoops_8_groups_max$year %in% "2014" &
                                  zoops_8_groups_max$Taxon %in% c("Bosmina","Daphnia",
                                                            "Cyclopoida","Nauplii")]))

mean(c(zoops_8_groups_max$month[zoops_8_groups_max$year %in% "2015" &
                                   zoops_8_groups_max$Taxon %in% c("Bosmina","Daphnia",
                                                               "Cyclopoida","Nauplii")]))

mean(c(zoops_8_groups_max$month[zoops_8_groups_max$year %in% "2016" &
                                  zoops_8_groups_max$Taxon %in% c("Bosmina","Daphnia",
                                                               "Cyclopoida","Nauplii")]))

mean(c(zoops_8_groups_max$month[zoops_8_groups_max$year %in% "2019" &
                                   zoops_8_groups_max$Taxon %in% c("Bosmina","Daphnia",
                                                               "Cyclopoida","Nauplii")]))

mean(c(zoops_8_groups_max$month[zoops_8_groups_max$year %in% "2020" &
                                  zoops_8_groups_max$Taxon %in% c("Bosmina","Daphnia",
                                                               "Cyclopoida","Nauplii")]))

mean(c(zoops_8_groups_max$month[zoops_8_groups_max$year %in% "2021" &
                                  zoops_8_groups_max$Taxon %in% c("Bosmina","Daphnia",
                                                               "Cyclopoida","Nauplii")]))

mean(c(zoops_8_groups_max$month[zoops_8_groups_max$year %in% "2023" &
                                  zoops_8_groups_max$Taxon %in% c("Bosmina","Daphnia",
                                                               "Cyclopoida","Nauplii")]))

#clads only
c(zoops_8_groups_max$month[zoops_8_groups_max$year %in% "2014" &
                              zoops_8_groups_max$Taxon %in% c("Bosmina","Daphnia")])

c(zoops_8_groups_max$month[zoops_8_groups_max$year %in% "2015" &
                             zoops_8_groups_max$Taxon %in% c("Bosmina","Daphnia")])

c(zoops_8_groups_max$month[zoops_8_groups_max$year %in% "2016" &
                             zoops_8_groups_max$Taxon %in% c("Bosmina","Daphnia")])

c(zoops_8_groups_max$month[zoops_8_groups_max$year %in% "2019" &
                              zoops_8_groups_max$Taxon %in% c("Bosmina","Daphnia")])

c(zoops_8_groups_max$month[zoops_8_groups_max$year %in% "2020" &
                             zoops_8_groups_max$Taxon %in% c("Bosmina","Daphnia")])

c(zoops_8_groups_max$month[zoops_8_groups_max$year %in% "2021" &
                             zoops_8_groups_max$Taxon %in% c("Bosmina","Daphnia")])

c(zoops_8_groups_max$month[zoops_8_groups_max$year %in% "2023" &
                              zoops_8_groups_max$Taxon %in% c("Bosmina", "Daphnia")])

