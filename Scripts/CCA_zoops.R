#trying out another multivariate approach for zoop data

#library(mvabund)
#read in packages
pacman::p_load(zoo,dplR,dplyr,tidyverse,ggplot2,ggpubr,sf)

#read in zoop data from EDI
inUrl1  <- "https://pasta-s.lternet.edu/package/data/eml/edi/1090/14/c7a04035b0a99adc489f5b6daec1cd52" 
infile1 <-  tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")

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
      Taxon %in% c("Lecane")]))

#convert back to long
zoops_final_pre <- zoops_pre |> 
  group_by(Reservoir, DateTime, StartDepth_m) |> 
  pivot_longer(cols=Daphnia:Lecane,
               names_to = c("Taxon"),
               values_to = "Density_IndPerL") |> 
  filter(hour(DateTime) %in% c(9,10,11,12,13,14)) |> #drop nighttime samples
  mutate(DateTime = as.Date(DateTime)) |> 
  mutate(Density_IndPerL = Density_IndPerL * (1/0.031))  #10m bvr neteff from 2016 (n=2) - note that 7m neteff was also 0.31
#avg from 2020 and 2021 is 0.021...

#list common taxa between pre and post
taxa <- c("Bosmina", "Daphnia", "Ceriodaphnia",
          "Cyclopoida","Calanoida", "nauplius", 
          "Conochilidae","Keratella", "Lecane",
          "Trichocercidae","Kellicottia")

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
  summarise(Density_IndPerL = mean(Density_IndPerL))

#combine all zoop data
all_zoops <- bind_rows(zoops_final_pre, zoops_final_post) |> 
  mutate_all(~replace(., is.nan(.), NA)) |>  #replace NAN with NA
  ungroup() |> select(-StartDepth_m) #dropping, but note that depths range from 7.5-11.5m....

#add column for pre vs post
all_zoops$data <- ifelse(all_zoops$DateTime<="2019-01-01","pre","post")

#convert to wide for CCA
all_zoops_wide <- all_zoops |> 
pivot_wider(names_from = Taxon, values_from = Density_IndPerL)


#-----------------------------------------------------------------------------#
#read in env data
dates <- unique(all_zoops_wide$DateTime)

#read in ctd
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/13/27ceda6bc7fdec2e7d79a6e4fe16ffdf" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

ctd <-read.csv(infile1,header=T) |> 
  mutate(DateTime = as.Date(DateTime)) |>
  filter(DateTime %in% dates & 
           Depth_m > 0 & Site ==50 & Reservoir %in% c("BVR")) |> 
  select(Reservoir, Site, DateTime, Depth_m, 
         Temp_C, DO_mgL, DOsat_percent, Chla_ugL)

#round ctd to nearest m
ctd_final <- ctd |>
  dplyr::mutate(rdepth = plyr::round_any(Depth_m, 1)) |>
  dplyr::group_by(Reservoir, DateTime, rdepth) |>
  dplyr::summarise(Temp_C = mean(Temp_C),
                   DO_mgL = mean(DO_mgL),
                   DOsat_percent = mean(DOsat_percent),
                   Chla_ugL = mean(Chla_ugL)) |>
  dplyr::rename(Depth_m = rdepth) |> 
  dplyr::filter(Depth_m %in% ifelse(Reservoir=="BVR", c(1, last(Depth_m[Reservoir=="BVR"])),
                                    c(1,9)))

#read in chem from edi
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/199/11/509f39850b6f95628d10889d66885b76" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

chem <-read.csv(infile1,header=T) |> 
  dplyr::mutate(DateTime = as.Date(DateTime)) |>
  dplyr::filter(DateTime %in% dates & 
                  Site ==50 & Reservoir %in% c("BVR")) |> 
  dplyr::select(Reservoir, DateTime, Depth_m, Rep,
                TN_ugL, TP_ugL,NH4_ugL ,NO3NO2_ugL, SRP_ugL) |> 
  dplyr::group_by(DateTime, Reservoir) |> 
  dplyr::filter(Depth_m %in% ifelse(Reservoir=="BVR", c(0.1,last(Depth_m[Reservoir=="BVR"])),
                                    c(0.1,last(Depth_m[Reservoir=="FCR"])))) |> 
  dplyr::filter(Depth_m <=0.1 | Depth_m > 7) |>  #drop data when only one random depth was analyzed
  dplyr::ungroup() |> dplyr::group_by(DateTime, Depth_m) |> 
  dplyr::summarise(TN_ugL = mean(TN_ugL,na.rm=T),
                   TP_ugL = mean(TP_ugL,na.rm=T),
                   NH4_ugL = mean(NH4_ugL,na.rm=T),
                   NO3NO2_ugL = mean(NO3NO2_ugL,na.rm=T),
                   SRP_ugL = mean(SRP_ugL,na.rm=T))

#list of matching dates (within 1 day)
matching_dates <- c("2014-05-14","2014-05-29","2014-06-04", "2014-06-18", 
                    "2014-06-25","2014-07-02", "2014-07-09", "2014-07-23",
                    "2014-08-13", "2014-09-04", "2014-09-25", "2015-05-14",
                    "2015-05-28", "2015-06-11", "2015-06-25", "2015-07-09",
                    "2015-07-23", "2015-07-30","2015-09-04", "2015-09-18",
                    "2015-10-02", "2015-10-23","2016-06-02", "2016-06-23",
                    "2016-07-07", "2016-07-21", "2016-08-11", "2016-08-23",
                    "2016-09-06", "2016-09-27", "2019-06-13", "2019-06-27",
                    "2019-07-11","2019-07-24", "2019-08-14", "2019-09-04",
                    "2019-10-04", "2019-11-15", "2020-06-25", "2020-07-02",
                    "2020-08-06","2020-08-12", "2020-09-03", "2020-10-01",
                    "2020-11-10","2021-03-08", "2021-04-05","2021-05-14",
                    "2021-06-28","2021-07-26","2022-03-22","2022-04-12", "2022-05-02")

#convert to wide and only select matching dates
chem_final_wide <- chem |> filter(DateTime %in% matching_dates) |> 
  filter(Depth_m %in% c(0.1)) |> 
  pivot_longer(cols=TN_ugL:SRP_ugL,
               names_to = c("var"),
               values_to = "value") |> 
  pivot_wider(names_from = c(var,Depth_m), values_from = value,
              names_glue = "{var}_{Depth_m}")

ctd_final_wide <- ctd_final |> 
  filter(DateTime %in% matching_dates) |>
  select(!DOsat_percent) |> 
  filter(Depth_m %in% c(1)) |> 
  pivot_longer(cols=Temp_C:Chla_ugL,
               names_to = c("var"),
               values_to = "value") |> 
  pivot_wider(names_from = c(var,Depth_m), values_from = value,
              names_glue = "{var}_{Depth_m}") |> ungroup() |> 
  select(!Reservoir)
  

#combine env data
env <- dplyr::bind_cols(ctd_final_wide,chem_final_wide) |> 
  select(!DateTime...5) |> 
  rename(DateTime = DateTime...1)

#just select data
env_cca <- env |> select(Temp_C_1:SRP_ugL_0.1) 

#replace NAN with 0
env_cca$NH4_ugL_0.1[is.nan(env_cca$NH4_ugL_0.1)]<-0
env_cca$NO3NO2_ugL_0.1[is.nan(env_cca$NO3NO2_ugL_0.1)]<-0
env_cca$SRP_ugL_0.1[is.nan(env_cca$SRP_ugL_0.1)]<-0
env_cca$DO_mgL_1[is.na(env_cca$DO_mgL_1)]<-0
  #NOTE that I can't do this, but just a placeholder for now



#transform env data
env_cca_trans <- log(env_cca + 1)

#only select matching dates in zoop data
all_zoops_wide <- all_zoops_wide |> filter(DateTime %in% matching_dates)
  
#only select data cols
all_zoops_cca <- all_zoops_wide |> select(Daphnia:Lecane)

#hellinger transform
all_zoops_cca_trans <- labdsv::hellinger(all_zoops_cca)


#correlation for species
cor(all_zoops_cca_trans)
cor(env_cca_trans)

pairs(all_zoops_cca_trans,env_cca_trans)   
#pairs(all_zoops_cca_trans)

#  Running the CCA
CCA <- cca(all_zoops_cca_trans, env_cca_trans)   
#  RDA <- rda(all_zoops_cca_trans, env_cca_trans)

####  Basic Graphing Platform ####
#  Produces the results from Ter Braak (well, mostly)
plot(CCA, choices=c(1,2), display=c('lc','sp','bp'), scaling=-1)   
plot(CCA, choices=c(1,2), display=c('wa','sp','bp'), scaling=-1)  

#### Utilities #####

#  Permutation (Monte Carlo) test for CCA   
permutest(CCA, permutations=500)

#  Variance Inflation factor for CCA
vif.cca(CCA)

# Species-Environment correlations
inertcomp(CCA)

# Interset correlations
intersetcor(CCA)


#####  Advanced plotting    #####

# making a scaling variable so that we don't have to change it every time
scale.par <- -1

#  Drawing the plot axes with no data
plot(CCA, choices=c(1,2), display=c('wa','sp','bp'), scaling=scale.par, 
     type='none', xlim=c(-5, 3), ylim=c(-3,4))
#  Adding the site scores
points(CCA, display='lc', pch=24, col='black', bg='wheat3', scaling=scale.par)
#  Adding the species scores    
text(CCA, display='sp', col='red', scaling=scale.par)
#  Adding the arrows with labels (Environmental vars)
#  You can use 'points' instead of 'text' to give arrows alone
text(CCA, display='bp', col='black', scaling=scale.par)    



