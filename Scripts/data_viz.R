#code to visualize seasonal and multiannual zoop patterns

#read in packages
pacman::p_load(tidyverse)

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
  group_by(Reservoir, StartDepth_m, DateTime, Taxon) |> 
  summarise(Density_IndPerL = mean(Density_IndPerL))

#combine all zoop data
all_zoops <- bind_rows(zoops_final_pre, zoops_final_post) |> 
  mutate_all(~replace(., is.nan(.), NA)) |>  #replace NAN with NA
  ungroup() |> select(-StartDepth_m) #dropping, but note that depths range from 7.5-11.5m....

#add doy and year column
all_zoops$doy <- yday(all_zoops$DateTime)
all_zoops$year <- year(all_zoops$DateTime)

#look at doy on x and year by color
ggplot(all_zoops, aes(doy, Density_IndPerL, color=as.factor(year))) + 
  geom_point() + theme_bw() + geom_line() +
  facet_wrap(~Taxon+Reservoir, scales="free_y", nrow=3) 

#not sure if I should look at all 11 taxa or just cladocerans, copepods, and rotifers

ggplot(data=subset(all_zoops, Taxon %in% c("Cladocera", "Copepoda", "Rotifera")),
                   aes(doy, Density_IndPerL, color=as.factor(year))) + 
  geom_point() + theme_bw() + geom_line() +
  facet_wrap(~Taxon+Reservoir, scales="free_y") 

#------------------------------------------------------------------------------#
#Pull in environmental data to make sure we have data associated with each zoop sample





