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
  summarise(Daphnia = sum(Biomass_ugL[
    Taxon %in% c("D. catawba","D. ambigua")]),
    Calanoida = sum(Biomass_ugL[
      Taxon %in% c("Diaptomus")]),
    Cyclopoida = sum(Biomass_ugL[
      Taxon %in% c("Cyclopoids")]),
    Nauplii = sum(Biomass_ugL[
      Taxon %in% c("Nauplii")]),
    Bosmina = sum(Biomass_ugL[
      Taxon %in% c("Bosmina")]),
    Ceriodaphnia = sum(Biomass_ugL[
      Taxon %in% c("Ceriodaphnia")]),
    Conochilus = sum(Biomass_ugL[
      Taxon %in% c("Conochilus")]),
    Keratella = sum(Biomass_ugL[
      Taxon %in% c("Keratella")]),
    Trichocerca = sum(Biomass_ugL[
      Taxon %in% c("Trichocerca")]),
    Kellicottia = sum(Biomass_ugL[
      Taxon %in% c("Kellicottia")]),
    Lecane = sum(Biomass_ugL[
      Taxon %in% c("Lecane")]),
    Rotifera = sum(Biomass_ugL[
      Taxon %in% c("Total Rotifers")]),
    Cladocera = sum(Biomass_ugL[
      Taxon %in% c("Bosmina","D. catawba", "Chydorus","D. ambigua",
                   "Diaphanosoma","Ceriodaphnia")]),
    Copepoda = sum(Biomass_ugL[
      Taxon %in% c("Diaptomus","Nauplii", "Cyclopoids")]))

#convert back to long
zoops_final_pre <- zoops_pre |> 
  group_by(Reservoir, DateTime, StartDepth_m) |> 
  pivot_longer(cols=Daphnia:Copepoda,
               names_to = c("Taxon"),
               values_to = "Biomass_ugL") |> 
  filter(hour(DateTime) %in% c(9,10,11,12,13,14)) |> #drop nighttime samples
  mutate(DateTime = as.Date(DateTime))

#list common taxa between pre and post
taxa <- c("Bosmina", "Daphnia", "Ceriodaphnia",
          "Cyclopoida","Calanoida", "nauplius", 
          "Conochilus","Keratella", "Rotifera",
          "Trichocerca","Kellicottia", "Lecane",
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
  summarise(Biomass_ugL = mean(Biomass_ugL))

#combine all zoop data
all_zoops <- bind_rows(zoops_final_pre, zoops_final_post) |> 
  mutate_all(~replace(., is.nan(.), NA)) |>  #replace NAN with NA
  ungroup() |> select(-StartDepth_m) #dropping, but note that depths range from 7.5-11.5m....

#add column for pre vs post
all_zoops$data <- ifelse(all_zoops$DateTime<="2019-01-01","pre","post")

#order data levels
all_zoops$data <- factor(all_zoops$data, levels=c("pre", "post"))

#get data into format for wavelet transformation
all_zoops_wavelet <- all_zoops |> 
  select(DateTime, Taxon, Biomass_ugL) |> 
  filter(Taxon %in% c("Cladocera","Copepoda","Rotifera")) |> 
  group_by(DateTime) |> 
  mutate(Total_biom = mean(Biomass_ugL, na.rm=T)) |> #should I sum???
  ungroup() |> 
  pivot_wider(names_from = Taxon, values_from = Biomass_ugL) |> 
  mutate(year = format(DateTime, "%Y"),
         month = format(DateTime, "%m")) |> 
  group_by(year, month) |> 
  summarise(Total_biom = mean(Total_biom),
            Cladocera = mean(Cladocera),
            Copepoda = mean(Copepoda),
            Rotifera = mean(Rotifera)) |> 
  ungroup() |> 
  select(year, month, Total_biom, Cladocera, Copepoda, Rotifera) |> 
  mutate(date = as.Date(paste0(year,"-",month,"-01")))

#now drop 2022 and only keep May-Sep
all_zoops_wavelet <- all_zoops_wavelet |> 
  filter(!year %in% "2022",
         month %in% c("05","06","07","08","09")) |> 
  mutate(Number = row_number())

#visualize data
ggplot(data=subset(all_zoops_wavelet, year< 2019), 
       aes(month, Total_biom, group=year, color=year))+
  geom_point() + geom_path() + theme_bw()

ggplot(data=subset(all_zoops_wavelet, year>= 2019), 
       aes(month, Total_biom, group=year, color=year))+
  geom_point() + geom_path() + theme_bw()
