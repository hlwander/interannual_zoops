#zoop density succession figs

#read in packages
pacman::p_load(zoo, dplR, dplyr, tidyverse, ggplot2, ggpubr, lubridate, ggtext)

#read in zoop data from EDI (v6)
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/197/6/38fc9d1a4c8b6976c71e56bda5ff073b?key=yltMpS4UEIk12AvB9L7OL5uRiG0"
infile1 <-  tempfile()
try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

zoops <- read.csv(infile1) |>
  filter(CollectionMethod=="Tow" & Reservoir %in% c("BVR") &
           StartDepth_m > 7.1 & EndDepth_m == 0.1 & Site==50) |> 
  select(-c(Site,EndDepth_m,CollectionMethod)) |>
  mutate(Taxon = ifelse(Taxon == "nauplius", "Nauplii",
                 ifelse(Taxon == "Total Rotifers", "Total rotifers",
                 ifelse(Taxon == "Cyclopoids", "Cyclopoida", Taxon))))  |>
  filter(month(DateTime) %in% 3:12) |># filter months used in this analysis
  mutate(DateTime = ifelse(DateTime == "2014-09-25 09:50:00", "2014-10-23 12:00:00", DateTime))
  #setting the second late sep 2014 sample as oct 23 (bc we have ctd data on this day)

#------------------------------------------------------------------------------#
# data coverage plot (Figure S1)

#list all dates
all_zoop_dates <- c(unique(as.Date(zoops$DateTime)))

#create df for plotting
coverage_df <- data.frame(SampleDate = all_zoop_dates) |>
  mutate(SampleMonth = floor_date(SampleDate, "month"),
         month = month(SampleDate))

# heat map
coverage_df <- data.frame(SampleDate = all_zoop_dates) |>
  mutate(Year = year(SampleDate),
         Month = month(SampleDate)) |>
  count(Year, Month) |>
  filter(!Year %in% c("2022","2024","2025"))

# convert numeric month to factor with labels
coverage_df$Month <- factor(coverage_df$Month, levels = 3:12, 
  labels = c("Mar","Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

#plot
ggplot(coverage_df, aes(x = Month, y = factor(Year))) +
  geom_tile(aes(fill = n), color = "white") +
  scale_fill_gradient(low = "orange", high = "darkblue") +
  labs(x = "Month", y = "Year", fill = "Sample Count") +
  theme_minimal()
#ggsave("Figures/zoop_data_coverage_heatmap.jpg", width=7, height=4) 

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
      Taxon %in% c("Cyclopoida")]),
    Cyclopoida_sd = sd(Density_IndPerL[
      Taxon %in% c("Cyclopoida")]),
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
    Trichocerca_dens = sum(Density_IndPerL[
      Taxon %in% c("Trichocerca")]),
    Trichocerca_sd = sd(Density_IndPerL[
      Taxon %in% c("Trichocerca")]),
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
     Taxon %in% c("Diaptomus","Nauplii", "Cyclopoids")]),
    Chydorus_dens = sum(Density_IndPerL[
      Taxon %in% c("Chydorus")]),
    Chydorus_sd = sd(Density_IndPerL[
      Taxon %in% c("Chydorus")]),
    Ascomorpha_dens = sum(Density_IndPerL[
      Taxon %in% c("Ascomorpha")]),
    Ascomorpha_sd = sd(Density_IndPerL[
      Taxon %in% c("Ascomorpha")]),
    Asplanchna_dens = sum(Density_IndPerL[
      Taxon %in% c("Asplanchna")]),
    Asplanchna_sd = sd(Density_IndPerL[
      Taxon %in% c("Asplanchna")]),
    Lecane_dens = sum(Density_IndPerL[
      Taxon %in% c("Lecane")]),
    Lecane_sd = sd(Density_IndPerL[
      Taxon %in% c("Lecane")]),
    Ceriodaphnia_dens = sum(Density_IndPerL[
      Taxon %in% c("Ceriodaphnia")]),
    Ceriodaphnia_sd = sd(Density_IndPerL[
      Taxon %in% c("Ceriodaphnia")]),
    Diaphanosoma_dens = sum(Density_IndPerL[
      Taxon %in% c("Diaphanosoma")]),
    Diaphanosoma_sd = sd(Density_IndPerL[
      Taxon %in% c("Diaphanosoma")]),
    Calanoida_dens = sum(Density_IndPerL[
      Taxon %in% c("Diaptomus")]),
    Calanoida_sd = sd(Density_IndPerL[
      Taxon %in% c("Diaptomus")])) |> 
    pivot_longer(-c(Reservoir,DateTime,StartDepth_m),
               names_to = c("Taxon", ".value"),
               names_sep="_" )  |> 
  filter(hour(DateTime) %in% c(8,9,10,11,12,13,14)) |> #drop nighttime samples
  mutate(DateTime = as.Date(DateTime)) |> 
  mutate(dens = dens * (1/0.031))  #10m bvr neteff from 2016 (n=2) - note that 7m neteff was also 0.031
#avg from 2020 and 2021 is 0.021 for reference

#average reps when appropriate
zoops_post <- zoops_2019_2023 |> 
  #mutate(DateTime = as.POSIXct(DateTime, format="%Y-%m-%d %H:%M:%S", tz="UTC")) |> 
  filter(hour(DateTime) %in% c(9,10,11,12,13,14) &  #drop nighttime samples
          !year(DateTime) %in% c(2022, 2024, 2025)) |> #we lost the 17 june 2024 sample somewhere between the field and the lab so that means this whole year needs to be excluded :(
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
#figure out dominant taxa for NMDS/other ms figs
zoop_taxa_props <- all_zoops |>
  filter(!is.na(Taxon), !is.na(dens),
         !Taxon %in% c("Cladocera","Copepoda","Rotifera","Crustacea",
                       "Total rotifers", "Chaoborus")) |>
  group_by(Taxon) |>
  mutate(n_years = n_distinct(year(DateTime)[dens > 0])) |>
  filter(n_years == 7) |> #keep taxa present in each year
  summarise(total_dens = sum(dens, na.rm = TRUE)) |>
  ungroup() |>
  mutate(prop = total_dens / sum(total_dens)) |>
  arrange(desc(prop))

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
# n = 10 taxa

#------------------------------------------------------------------------------#
# 10 group zoop dens
zoops_10_groups <- all_zoops |> 
  filter(Taxon %in% zoop_taxa_props$Taxon) |> 
  group_by(DateTime) |> 
  summarise(Daphnia_avg = mean(dens[Taxon=="Daphnia"], na.rm=T),
            Daphnia_sd = mean(sd[Taxon=="Daphnia"],na.rm=T),
            Ceriodaphnia_avg = mean(dens[Taxon=="Ceriodaphnia"], na.rm=T),
            Ceriodaphnia_sd = mean(sd[Taxon=="Ceriodaphnia"],na.rm=T),
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
            Polyarthra_sd = mean(sd[Taxon=="Polyarthra"],na.rm=T),
            Trichocerca_avg = mean(dens[Taxon=="Trichocerca"], na.rm=T),
            Trichocerca_sd = mean(sd[Taxon=="Trichocerca"],na.rm=T)) |>
  pivot_longer(-c(DateTime),
               names_to = c("Taxon", ".value"),
               names_sep="_" ) |>
  rename(dens = avg)
#write.csv(zoops_10_groups, paste0("Output/zoop_dens_10groups_nmds.csv"),row.names = FALSE)

#Table S2 - percent contribution to total density and frequency
n_total_dates <- n_distinct(all_zoops$DateTime)

zoop_summary_table <- all_zoops |>
  filter(!Taxon %in% c("Rotifera","Crustacea","Copepoda","Cladocera")) |>
  group_by(DateTime, Taxon) |>
  summarise(dens_avg = mean(dens, na.rm = TRUE), .groups = "drop") |>
  group_by(Taxon) |>
  summarise(total_dens = sum(dens_avg, na.rm = TRUE),
            n_dates_present = sum(dens_avg > 0, na.rm = TRUE)) |>
  mutate(percent_contribution = round(100 * total_dens / sum(total_dens), 2),
         frequency_occurrence = round(100 * n_dates_present / n_total_dates, 1)) |>
  arrange(desc(percent_contribution)) |>
  select(Taxon, frequency_occurrence, percent_contribution)
#write.csv(zoop_summary_table, "Output/zoop_percent_contribution.csv", row.names = FALSE)

#for plotting
zoops_10_groups <- zoops_10_groups |>
  rename(avg = dens) |>
  mutate(total = sum(avg, na.rm=T)) |> ungroup() |> 
  mutate(p_dens = avg / total,
  year = year(DateTime),
         month = month(DateTime),
         day   = day(DateTime),
         pseudoDate = as.Date(sprintf("2000-%02d-%02d", month, day))) |>
  select(-day) 

#order zoops by group
zoops_10_groups$Taxon <- factor(zoops_10_groups$Taxon, 
                                levels=c("Bosmina","Daphnia", "Ceriodaphnia",
                                         "Cyclopoida","Nauplii",
                                         "Conochilus","Kellicottia",
                                         "Keratella","Polyarthra", "Trichocerca")) 

taxon_labels <- c(expression(italic("Bosmina")),
  expression(italic("Daphnia")),
  expression(italic("Ceriodaphnia")),
  "Cyclopoida","Nauplii",
  expression(italic("Conochilus")),
  expression(italic("Kellicottia")),
  expression(italic("Keratella")),
  expression(italic("Polyarthra")),
  expression(italic("Trichocerca")))

#shaded line plot - raw density (Manuscript Figure 1)
ggplot(data = zoops_10_groups|> filter(month %in% c(4:11)),
       aes(x=pseudoDate, y = avg, color=Taxon)) +
  geom_area(aes(color = Taxon, fill = Taxon),
            position = "stack", stat = "identity") +
  facet_wrap(~year, scales = "free")+
  scale_color_manual(values = NatParksPalettes::natparks.pals(
    "DeathValley", 14, direction=-1)[c(1,2,3,5,7,10,11,12,13,14)],
    labels = taxon_labels) +
  scale_fill_manual(values = NatParksPalettes::natparks.pals(
    "DeathValley", 14, direction=-1)[c(1,2,3,5,7,10,11,12,13,14)],
    labels = taxon_labels) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b",
               expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  xlab("") + ylab("Density (individuals/L)") +
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

#shaded line plot - relative density (Figure. S2)
ggplot(data = subset(zoops_10_groups, month %in% c(4:11)),
       aes(x=pseudoDate, y = avg, color=Taxon)) +
  geom_area(aes(color = Taxon, fill = Taxon),
            position = "fill", stat = "identity") +
  facet_wrap(~year(DateTime), scales = "free_x")+
  scale_color_manual(values = NatParksPalettes::natparks.pals(
    "DeathValley", 14, direction=-1)[c(1,2,3,5,7,10,11,12,13,14)],
    labels = taxon_labels) +
  scale_fill_manual(values = NatParksPalettes::natparks.pals(
    "DeathValley", 14, direction=-1)[c(1,2,3,5,7,10,11,12,13,14)],
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
#ggsave("Figures/BVR_succession_10groups_fill_alldens_relative.jpg", width=5, height=4) 

#----------------------
#numbers for ms results
zoops_10_groups_max <- zoops_10_groups |>
  group_by(Taxon, year) |>
  slice_max(order_by = avg, n = 1, with_ties = FALSE) |>
  ungroup()

#bosmina mostly peak between aug-nov (except 2014; end of april)
#daphnia peak between may-oct
#ceriodaphnia peak between jun-sep
#cyclopoid peak between jun-nov
#nauplii peak between jun-sep
#conochilus peak between may-oct
#kellicottia peak between apr-aug
#keratella peak between jun-aug (except 2023; april)
#polyarthra peak between may-dec (super variable)
#trichocerca peak between jul-oct (except 2023; april)

#rotifers
mean(c(zoops_10_groups_max$month[zoops_10_groups_max$Taxon %in% c("Conochilus")]))   # 6.9
mean(c(zoops_10_groups_max$month[zoops_10_groups_max$Taxon %in% c("Keratella")]))    # 6.4
mean(c(zoops_10_groups_max$month[zoops_10_groups_max$Taxon %in% c("Kellicottia")]))  # 6.1
mean(c(zoops_10_groups_max$month[zoops_10_groups_max$Taxon %in% c("Polyarthra")]))   # 8
mean(c(zoops_10_groups_max$month[zoops_10_groups_max$Taxon %in% c("Trichocerca")]))  # 8.3

#crustaceans
mean(c(zoops_10_groups_max$month[zoops_10_groups_max$Taxon %in% c("Bosmina")]))      # 9.3
mean(c(zoops_10_groups_max$month[zoops_10_groups_max$Taxon %in% c("Daphnia")]))      # 6.9
mean(c(zoops_10_groups_max$month[zoops_10_groups_max$Taxon %in% c("Ceriodaphnia")])) # 7.6

mean(c(zoops_10_groups_max$month[zoops_10_groups_max$Taxon %in% c("Cyclopoida")]))   # 8
mean(c(zoops_10_groups_max$month[zoops_10_groups_max$Taxon %in% c("Nauplii")]))      # 7.3
