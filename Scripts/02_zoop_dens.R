#zoop density succession figs

#read in packages
pacman::p_load(zoo, dplR, dplyr, tidyverse, ggplot2, ggpubr, lubridate, ggtext)

#cb friendly year palette
year_cols <- c("#a13637","#06889b", "#facd60", "#f44034", "#011f51", "#fdfa66")

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
  filter(!year(DateTime) %in% c("2022","2023","2024","2025")) #recreated summary file for these years and it is not yet published on EDI 

zoops_2023 <- read.csv("inputs/EDI_zoop_summary_2022-2025.csv") |>
  select(-c(Site,EndDepth_m,CollectionMethod)) |>
  mutate(Taxon = ifelse(Taxon == "nauplius", "Nauplii", Taxon),
         DateTime = as.POSIXct(DateTime))
  
#add in 2022-2025 zoop data
zoops <- bind_rows(zoops, zoops_2023)

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
  filter(prop_avg >= 0.01)
#drop taxa that are <1% (0.01) of density 

#just keeping n=10

#---------------------------------
#split data into pre 2019 and post
zoops_2016_2018 <- zoops[as.Date(zoops$DateTime)<"2019-01-01",]
zoops_2019_2023 <- zoops[as.Date(zoops$DateTime)>="2019-01-01",]

#add daphnia (D. catawba, D. ambigua), calanoida (diaptomus) groups for 2014-2018 data
zoops_pre <- zoops_2016_2018 |> 
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
    Conochiloides_dens = sum(Density_IndPerL[
      Taxon %in% c("Conochiloides")]),
    Conochiloides_sd = sd(Density_IndPerL[
      Taxon %in% c("Conochiloides")]),
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
    Ploima_dens = sum(Density_IndPerL[
      Taxon %in% c("Ploima")]),
    Ploima_sd = sd(Density_IndPerL[
      Taxon %in% c("Ploima")]),
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
  ungroup() |> select(-StartDepth_m) #dropping, but note that depths range from 8-11.5m....

#write.csv(all_zoops, paste0("Output/all_zoops_dens.csv"),row.names = FALSE)

#------------------------------------------------------------------------------#
# 10 group zoop dens
zoops_10_groups <- all_zoops |> 
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
            Conochiloides_avg = mean(dens[Taxon=="Conochiloides"], na.rm=T),
            Conochiloides_sd = mean(sd[Taxon=="Conochiloides"],na.rm=T),
            Keratella_avg = mean(dens[Taxon=="Keratella"], na.rm=T),
            Keratella_sd = mean(sd[Taxon=="Keratella"],na.rm=T),
            Kellicottia_avg = mean(dens[Taxon=="Kellicottia"], na.rm=T),
            Kellicottia_sd = mean(sd[Taxon=="Kellicottia"],na.rm=T),
            Ploima_avg = mean(dens[Taxon=="Ploima"], na.rm=T),
            Ploima_sd = mean(sd[Taxon=="Ploima"],na.rm=T),
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
zoops_10_groups$Taxon <- factor(zoops_10_groups$Taxon, 
                                levels=c("Bosmina", "Daphnia","Cyclopoida",
                                         "Nauplii","Conochiloides",
                                         "Conochilus","Kellicottia","Keratella",
                                         "Ploima","Polyarthra")) 

#shaded line plot - raw density (Manuscript Figure 2)
ggplot(data = zoops_10_groups|> filter(month %in% c(4:11)),
       aes(x=pseudoDate, y = avg, color=Taxon)) +
  geom_area(aes(color = Taxon, fill = Taxon),
            position = "stack", stat = "identity") +
  facet_wrap(~year, scales = "free")+
  scale_color_manual(values = NatParksPalettes::natparks.pals(
    "DeathValley", 14, direction=-1)[c(1,3,5,7,9:14)]) +
  scale_fill_manual(values = NatParksPalettes::natparks.pals(
    "DeathValley", 14, direction=-1)[c(1,3,5,7,9:14)]) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b",
              # limits = as.Date(c("2000-04-01", "2000-11-30")),
               expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  xlab("") + ylab("Density (ind/L)") +
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
#ggsave("Figures/BVR_10groups_fill_alldens.jpg", width=5, height=4) 

#shaded line plot - relative density (Fig. S1)
ggplot(data = subset(zoops_10_groups, month(DateTime) %in% c(5:10)),
       aes(x=DateTime, y = avg, color=Taxon)) +
  geom_area(aes(color = Taxon, fill = Taxon),
            position = "fill", stat = "identity") +
  facet_wrap(~year(DateTime), scales = "free_x")+
  scale_color_manual(values = NatParksPalettes::natparks.pals(
    "DeathValley", 14, direction=-1)[c(1,3,5,7,9:14)]) +
  scale_fill_manual(values = NatParksPalettes::natparks.pals(
    "DeathValley", 14, direction=-1)[c(1,3,5,7,9:14)]) +
  scale_x_date(expand = c(0,0)) +
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
#ggsave("Figures/BVR_succession_10groups_fill_alldens_relative.jpg", width=5, height=4) 

#numbers for ms results
mean(c(sum(zoops_10_groups$p_dens[zoops_10_groups$year=="2019" & zoops_10_groups$month=="06" &
                                 zoops_10_groups$Taxon %in% c("Ploima","Conochilus", "Conochiloides",
                                                              "Keratella","Kellicottia",
                                                              "Polyarthra", "Pompholyx")]),
        sum(zoops_10_groups$p_dens[zoops_10_groups$year=="2019" & zoops_10_groups$month=="07" &
                                 zoops_10_groups$Taxon %in% c("Ploima","Conochilus", "Conochiloides",
                                                              "Keratella","Kellicottia",
                                                              "Polyarthra", "Pompholyx")]),
        sum(zoops_10_groups$p_dens[zoops_10_groups$year=="2019" & zoops_10_groups$month=="08" &
                                 zoops_10_groups$Taxon %in% c("Ploima","Conochilus", "Conochiloides",
                                                              "Keratella","Kellicottia",
                                                              "Polyarthra", "Pompholyx")])))

#kellicottia
mean(c(zoops_10_groups$p_dens[zoops_10_groups$year=="2021" & zoops_10_groups$month=="05" &
                                    zoops_10_groups$Taxon %in% c("Kellicottia")],
       zoops_10_groups$p_dens[zoops_10_groups$year=="2021" & zoops_10_groups$month=="06" &
                                    zoops_10_groups$Taxon %in% c("Kellicottia")],
       zoops_10_groups$p_dens[zoops_10_groups$year=="2021" & zoops_10_groups$month=="07" &
                                    zoops_10_groups$Taxon %in% c("Kellicottia")],
       zoops_10_groups$p_dens[zoops_10_groups$year=="2021" & zoops_10_groups$month=="08" &
                                    zoops_10_groups$Taxon %in% c("Kellicottia")],
       zoops_10_groups$p_dens[zoops_10_groups$year=="2021" & zoops_10_groups$month=="09" &
                                    zoops_10_groups$Taxon %in% c("Kellicottia")]))

#keratella
mean(c(zoops_10_groups$p_dens[zoops_10_groups$year=="2021" & zoops_10_groups$month=="05" &
                                    zoops_10_groups$Taxon %in% c("Keratella")],
       zoops_10_groups$p_dens[zoops_10_groups$year=="2021" & zoops_10_groups$month=="06" &
                                    zoops_10_groups$Taxon %in% c("Keratella")],
       zoops_10_groups$p_dens[zoops_10_groups$year=="2021" & zoops_10_groups$month=="07" &
                                    zoops_10_groups$Taxon %in% c("Keratella")],
       zoops_10_groups$p_dens[zoops_10_groups$year=="2021" & zoops_10_groups$month=="08" &
                                    zoops_10_groups$Taxon %in% c("Keratella")],
       zoops_10_groups$p_dens[zoops_10_groups$year=="2021" & zoops_10_groups$month=="09" &
                                    zoops_10_groups$Taxon %in% c("Keratella")]))

#cyclopoids
mean(c(zoops_10_groups$p_dens[zoops_10_groups$year=="2021" & zoops_10_groups$month=="05" &
                                zoops_10_groups$Taxon %in% c("Cyclopoida")],
       zoops_10_groups$p_dens[zoops_10_groups$year=="2021" & zoops_10_groups$month=="06" &
                                zoops_10_groups$Taxon %in% c("Cyclopoida")],
       zoops_10_groups$p_dens[zoops_10_groups$year=="2021" & zoops_10_groups$month=="07" &
                                zoops_10_groups$Taxon %in% c("Cyclopoida")],
       zoops_10_groups$p_dens[zoops_10_groups$year=="2021" & zoops_10_groups$month=="08" &
                                zoops_10_groups$Taxon %in% c("Cyclopoida")],
       zoops_10_groups$p_dens[zoops_10_groups$year=="2021" & zoops_10_groups$month=="09" &
                                zoops_10_groups$Taxon %in% c("Cyclopoida")]))
  
#clads
mean(c(sum(zoops_10_groups$avg[zoops_10_groups$year=="2014" & zoops_10_groups$month=="05" &
                                    zoops_10_groups$Taxon %in% c("Bosmina","Ceriodaphnia",
                                                                 "Daphnia")]),
       sum(zoops_10_groups$avg[zoops_10_groups$year=="2014" & zoops_10_groups$month=="06" &
                                    zoops_10_groups$Taxon %in% c("Bosmina","Ceriodaphnia",
                                                                 "Daphnia")]),
       sum(zoops_10_groups$avg[zoops_10_groups$year=="2014" & zoops_10_groups$month=="07" &
                                    zoops_10_groups$Taxon %in% c("Bosmina","Ceriodaphnia",
                                                                 "Daphnia")]),
       sum(zoops_10_groups$avg[zoops_10_groups$year=="2014" & zoops_10_groups$month=="07" &
                                    zoops_10_groups$Taxon %in% c("Bosmina","Ceriodaphnia",
                                                                 "Daphnia")]),
       sum(zoops_10_groups$avg[zoops_10_groups$year=="2014" & zoops_10_groups$month=="09" &
                                    zoops_10_groups$Taxon %in% c("Bosmina","Ceriodaphnia",
                                                                 "Daphnia")])))

#copes
mean(c(sum(zoops_10_groups$avg[zoops_10_groups$year=="2021" & zoops_10_groups$month=="05" &
                                 zoops_10_groups$Taxon %in% c("Cyclopoida","Nauplii")]),
       sum(zoops_10_groups$avg[zoops_10_groups$year=="2021" & zoops_10_groups$month=="06" &
                                 zoops_10_groups$Taxon %in% c("Cyclopoida","Nauplii")]),
       sum(zoops_10_groups$avg[zoops_10_groups$year=="2021" & zoops_10_groups$month=="07" &
                                 zoops_10_groups$Taxon %in% c("Cyclopoida","Nauplii")]),
       sum(zoops_10_groups$avg[zoops_10_groups$year=="2021" & zoops_10_groups$month=="07" &
                                 zoops_10_groups$Taxon %in% c("Cyclopoida","Nauplii")]),
       sum(zoops_10_groups$avg[zoops_10_groups$year=="2021" & zoops_10_groups$month=="09" &
                                 zoops_10_groups$Taxon %in% c("Cyclopoida","Nauplii")])))

#rots
mean(c(sum(zoops_10_groups$avg[zoops_10_groups$year=="2014" & zoops_10_groups$month=="05" &
                                 zoops_10_groups$Taxon %in% c("Ascomorpha","Conochilus",
                                                              "Keratella","Kellicottia",
                                                              "Polyarthra")]),
       sum(zoops_10_groups$avg[zoops_10_groups$year=="2014" & zoops_10_groups$month=="06" &
                                 zoops_10_groups$Taxon %in% c("Ascomorpha","Conochilus",
                                                              "Keratella","Kellicottia",
                                                              "Polyarthra")]),
       sum(zoops_10_groups$avg[zoops_10_groups$year=="2014" & zoops_10_groups$month=="07" &
                                 zoops_10_groups$Taxon %in% c("Ascomorpha","Conochilus",
                                                              "Keratella","Kellicottia",
                                                              "Polyarthra")]),
       sum(zoops_10_groups$avg[zoops_10_groups$year=="2014" & zoops_10_groups$month=="07" &
                                 zoops_10_groups$Taxon %in% c("Ascomorpha","Conochilus",
                                                              "Keratella","Kellicottia",
                                                              "Polyarthra")]),
       sum(zoops_10_groups$avg[zoops_10_groups$year=="2014" & zoops_10_groups$month=="09" &
                                 zoops_10_groups$Taxon %in% c("Ascomorpha","Conochilus",
                                                              "Keratella","Kellicottia",
                                                              "Polyarthra")])))

#---------------------------------------------------------------------#
#C:R vs nmds1

# different taxa ratios
zoop_ratios <- all_zoops |> 
  filter(Taxon %in% c("Cladocera","Cyclopoida","Nauplii","Keratella",
                      "Kellicottia","Conochilus","Rotifera","Copepoda")) |> 
  mutate(year = year(DateTime),
         month = format(DateTime, "%m")) |> 
  filter(month %in% c("05","06","07","08","09","10")) |> 
  group_by(DateTime) |> 
  mutate(total_dens = sum(dens),
         crust_dens = sum(dens[Taxon=="Cladocera"], 
                          dens[Taxon=="Copepoda"],na.rm=T),
         p_dens = dens / total_dens) |> 
  ungroup() |> group_by(year) |>
  summarise(C_R =  mean(crust_dens) /
               mean(dens[Taxon=="Rotifera"], na.rm=T),
            cyc_nau = mean(dens[Taxon=="Cyclopoida"], na.rm=T) /
              mean(dens[Taxon=="Nauplii"], na.rm=T),
            ker_kel = mean(dens[Taxon=="Keratella"], na.rm=T) /
              mean(dens[Taxon=="Kellicottia"], na.rm=T),
            ker_naup = mean(dens[Taxon=="Keratella"], na.rm=T) /
              mean(dens[Taxon=="Nauplii"], na.rm=T),
            kel_naup = mean(dens[Taxon=="Kellicottia"], na.rm=T) /
              mean(dens[Taxon=="Nauplii"], na.rm=T),
            clad_rot = mean(dens[Taxon=="Cladocera"], na.rm=T) /
              mean(dens[Taxon=="Rotifera"], na.rm=T),
            naup_cono = mean(dens[Taxon=="Nauplii"], na.rm=T) /
              mean(dens[Taxon=="Conochilus"], na.rm=T),
            cyc_cono = mean(dens[Taxon=="Cyclopoida"], na.rm=T) /
              mean(dens[Taxon=="Conochilus"], na.rm=T),
            ker_cono = mean(dens[Taxon=="Keratella"], na.rm=T) /
              mean(dens[Taxon=="Conochilus"], na.rm=T),
            kel_cono = mean(dens[Taxon=="Kellicottia"], na.rm=T) /
              mean(dens[Taxon=="Conochilus"], na.rm=T),
            p_cyc = mean(p_dens[Taxon=="Cyclopoida"],na.rm=T),
            p_nau = mean(p_dens[Taxon=="Nauplii"],na.rm=T),
            p_ker = mean(p_dens[Taxon=="Keratella"],na.rm=T),
            p_kel = mean(p_dens[Taxon=="Kellicottia"],na.rm=T),
            p_con = mean(p_dens[Taxon=="Conochilus"],na.rm=T),
            p_rot = mean(p_dens[Taxon=="Rotifera"],na.rm=T))
#write.csv(zoop_ratios, "Output/zoop_proportions.csv")

nmds1 <- read.csv("Output/ss_nmds1.csv")
zoop_ratios$nmds1 <- nmds1$nmds1

#convert from wide to long
zoop_ratios_long <- zoop_ratios |> 
  pivot_longer(-c(year,nmds1),
               names_to = "metric") 

metric_names <- c("C_R" = "Crustacean:Rotifer",
                  "clad_rot" = "Cladocera:Rotifer",
                  "cyc_cono" = "Cyclopoida:*Conochilus*", 
                  "cyc_nau" = "Cyclopoida:Nauplius",
                  "kel_cono" = "*Kellicottia:Conochilus*",
                  "kel_naup" = "*Kellicottia*:Nauplius",
                  "ker_cono" = "*Keratella:Conochilus*",
                  "ker_kel" = "*Keratella:Kellicottia*",
                  "ker_naup" = "*Keratella*:Nauplius",
                  "naup_cono" = "Nauplius:*Conochilus*",
                  "p_con" = "Percent *Conochilus* (%)",
                  "p_cyc" = "Percent Cyclopoida (%)",
                  "p_kel" = "Percent *Kellicottia* (%)", 
                  "p_ker" = "Percent *Keratella* (%)", 
                  "p_nau" = "Percent Nauplius (%)", 
                  "p_rot" = "Percent Rotifer (%)")

#-------------------------------------------------------#
# NMDS1 vs. zoop metrics (Supplemental Figure S6)
ggplot(zoop_ratios_long, aes(nmds1, value, color=as.factor(year))) +
  geom_point(size=4) + theme_bw() + xlab("NMDS axis 1 value") +
  scale_color_manual(values = year_cols) + ylab("Value") +
  facet_wrap(~metric, scales = "free_y", 
             labeller = as_labeller(metric_names)) +
  guides(color = guide_legend("",nrow=1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        strip.text = ggtext::element_markdown(),
        text = element_text(size=10), 
        panel.border = element_rect(colour = "black", fill = NA),
        plot.margin = unit(c(0, 1, 0, 0), "cm"),
        legend.box.margin=margin(-5,-10,-15,-10),
        panel.spacing.x = unit(0.2, "in"),
        panel.background = element_rect(
          fill = "white"),
        panel.spacing = unit(0.5, "lines"))
#ggsave("Figures/zoop_metrics_vs_nmds1.jpg", width=8, height=6)
