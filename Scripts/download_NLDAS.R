###########################################################
### Downloading NLDAS2 data for meteorological hourly forcing
### http://ldas.gsfc.nasa.gov/nldas/NLDAS2forcing.php
### Author: Hilary Dugan hilarydugan@gmail.com and Mary Lofton melofton@vt.edu
### Date last modified: 2024-04-01 (adapted for bvr + extended dates - HLW)
###########################################################

#### step 1: download nc files for each day ####
# v 04Mar2021: BGS updated to remove hardcoding: only first 43 lines need to be edited

library(RCurl)
library(lubridate)
library(raster)
library(ncdf4)
library(sf)
library(httr)
library(curl)
library(stringr)
library(dplyr)
library(readr)

###########################################################
### Point to dump directory where data will be saved
###########################################################
dumpdir_nc = "./NLDAS"

###########################################################
### Enter password information
###########################################################
#https://urs.earthdata.nasa.gov/profile <-- GET A EARTHDATA LOGIN
username = "hwander"
password = "Zooplankton2024!"
#in addition, make sure you have authorized your account access to the GEODISC archives:
# https://disc.gsfc.nasa.gov/earthdata-login

###########################################################
### Use shapefile of lake to set bounding box
###########################################################
# read in lake file (as a .shp file) to get bounding box
# lakeShape = st_read('shapefile.shp') 
extent = c(-79.82,37.31,-79.81,37.32) 
#if the extent is loaded from the shapefile (above), make sure they are in decimal degrees, otherwise this code will not work

###########################################################
### Set timeframe
###########################################################
startdatetime = '2021-12-15 02:00:00'
enddatetime = '2021-12-15 03:00:00'
loc_tz = 'UTC'

# sequence the datetime over your desired time period
out.ts = seq.POSIXt(as.POSIXct(startdatetime, tz = loc_tz),as.POSIXct(enddatetime,tz=loc_tz), by = 'hour')

# Create output list of tables
output = list()

###########################################################
### Run hourly loop
###########################################################
# Start the clock!
ptm <- proc.time()

for (i in 1:length(out.ts)) {
  print(out.ts[i])
  yearOut = year(out.ts[i])
  monthOut = format(out.ts[i], "%m")
  dayOut = format(out.ts[i], "%d")
  hourOut = format(out.ts[i], "%H%M")
  doyOut = format(out.ts[i],'%j')
  
  filename = format(out.ts[i], "%Y%m%d%H%M")
  
  # URL3 = paste('http://',username,':',password,'@hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?',
  #              'FILENAME=%2Fdata%2FNLDAS%2FNLDAS_FORA0125_H.002%2F',yearOut,'%2F',doyOut,'%2FNLDAS_FORA0125_H.A',yearOut,monthOut,dayOut,'.',hourOut,'.002.grb&',
  #              'FORMAT=bmV0Q0RGLw&BBOX=',extent[2],'%2C',extent[1],'%2C',extent[4],'%2C',extent[3],'&',
  #              'LABEL=NLDAS_FORA0125_H.A',yearOut,monthOut,dayOut,'.',hourOut,'.002.2017013163409.pss.nc&',
  #              'SHORTNAME=NLDAS_FORA0125_H&SERVICE=SUBSET_GRIB&VERSION=1.02&DATASET_VERSION=002',sep='')
  #
  # URL <- paste('http://hydro1.sci.gsfc.nasa.gov/daac-bin/OTF/HTTP_services.cgi?',
  #              'FILENAME=%2Fdata%2Fs4pa%2FNLDAS%2FNLDAS_FORA0125_H.002%2F',yearOut,
  #              '%2F',doyOut,
  #              '%2FNLDAS_FORA0125_H.A',yearOut,monthOut,dayOut,'.',
  #              hourOut,'.002.grb&',
  #              'FORMAT=bmV0Q0RGLw&BBOX=',
  #              extent[2],'%2C',
  #              extent[1],'%2C',
  #              extent[4],'%2C',
  #              extent[3],'&',
  #              'LABEL=NLDAS_FORA0125_H.A',yearOut,monthOut,dayOut,'.',
  #              hourOut,'.002.2016116144611.pss.nc&',
  #              'SHORTNAME=NLDAS_FORA0125_H&SERVICE=SUBSET_GRIB&VERSION=1.02&DATASET_VERSION=002',sep='')
  
  ## Patricia Tran note (2020-12-03 : I updated the link the webpage)
  
  URL <- paste('https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FNLDAS%2FNLDAS_FORA0125_H.002%2F',
               yearOut, '%2F',
               str_pad(as.numeric(yday(as.Date(paste0(yearOut,"-", monthOut,'-' ,dayOut)))), 3, pad = "0"), ## The URL changes for every chunk of 24 hours
               '%2FNLDAS_FORA0125_H.A',
               yearOut, monthOut, dayOut, '.',
               hourOut,  '.002.grb&FORMAT=bmM0Lw&BBOX=', 
               round(extent[2], 2),'%2C', # In the new version of the URL, the coordinates are only up to 2 digits
               round(extent[1], 2),'%2C',
               round(extent[4], 2),'%2C',
               round(extent[3], 2),
               '&LABEL=NLDAS_FORA0125_H.A',
               yearOut,monthOut,dayOut,'.',
               hourOut,
               '.002.grb.SUB.nc4&SHORTNAME=NLDAS_FORA0125_H&SERVICE=L34RS_LDAS&VERSION=1.02&DATASET_VERSION=002',
               sep='')
  
  
  # IMPORTANT MESSAGE Dec 05, 2016    The GES DISC will be migrating from http to https throughout December
  # As part of our ongoing migration to HTTPS, the GES DISC will begin redirecting all HTTP traffic to HTTPS.
  # We expect to have all GES DISC sites redirecting traffic by January 4th. For most access methods, the redirect will be transparent to the user.
  # However, users with locally developed scripts or utilities that do not support an HTTP code 301 redirect may find that the scripts will fail.
  # If you access our servers non-interactively (i.e. via a mechanism other than a modern web browser), you will want to modify your scripts to
  # point to the HTTPS addresses to avoid the enforced redirect.
  
  # x = download.file(URL3,destfile = paste(filename,'.nc',sep=''),mode = 'wb',quiet = T)
  # x = download.file(URL,destfile = paste(filename,'.nc',sep=''),mode = 'wb',quiet = T)
  
  
  lk <- URL
  
  #wget:
  #r <- GET(lk,
  #          authenticate("ptran5@wisc.edu", "Earthdata1"),
  #          path = "~/Documents/MendotaRawData/")
  
  # or this with curl
  h <- curl::new_handle()
  
  curl::handle_setopt(
    handle = h,
    httpauth = 1,
    userpwd = paste0(username, ':', password)
  )
  
  # resp <- curl::curl_fetch_memory(lk, handle = h)
  resp <- curl::curl_fetch_disk(url = lk, 
                                path = paste(dumpdir_nc, '/NLDAS_Data_2020_2021/', filename, '_', loc_tz, '.nc',sep=''), 
                                handle = h)
  
  #Sys.sleep(2)
  
}
# Stop the clock
proc.time() - ptm

#------------------------------------------------------------------#
#### step 2: save each variable as a new csv ####

###########################################################
### set dump directory for .csv files and lake name
###########################################################
#where your .nc files are stored -- make sure all .nc files have a size >0, otherwise your loop will get hung up!
dumpdir_nc = './NLDAS/NLDAS_Data_2020_2021'


lake_name = 'BVR'

###########################################################
### Need to know how many cells your lake falls within
### Can download one instance of data from the earthdata site and see how many columns there are
### use 'nc_open(filename)' to see how many cells there are
###########################################################
cellNum = 1 
#How many output cells will there be? Need to check this beforehand by downloading a single netcdf file for your location

loc_tz = 'UTC'  #this needs to match your previous entry

###########################################################
### Set up the output data frame
###########################################################

#save the variable names in nc files
vars_nc = c('TMP','SPFH', 'PRES', 'UGRD', 'VGRD', 
            'DLWRF', 'CONVfrac', 'CAPE', 'PEVAP', 'APCP', 'DSWRF')

#set up output dataframe for the number of cells above and the number of columns of data
output <- NULL
for (l in 1:length(vars_nc)){
  colClasses = c("POSIXct", rep("numeric",cellNum))
  col.names = c('dateTime',rep(vars_nc[l],cellNum))
  output[[l]] = read.table(text = "",colClasses = colClasses,col.names = col.names)
  attributes(output[[l]]$dateTime)$tzone = 'GMT'
}

###########################################################
### Run file list loop
###########################################################
# Start the clock!
ptm <- proc.time()

nc_files <- list.files(dumpdir_nc)

for (i in 1:length(nc_files)) {
  print(nc_files[i])
  
  for (v in 1:length(vars_nc)) {
    nldasvar <- vars_nc[v]
    br = nc_open(paste0(dumpdir_nc, '/', nc_files[i]))
    output[[v]][i,1] =  as.POSIXct(paste0(substr(nc_files[i], 1, 4),'-', substr(nc_files[i], 5,6), '-', substr(nc_files[i], 7,8), ' ', substr(nc_files[i], 9,10), ':', substr(nc_files[i], 11,12)), tz=loc_tz)
    output[[v]][i,-1] = ncvar_get(br, nldasvar)
    nc_close(br)
  }
  rm(br)
}

# Stop the clock
proc.time() - ptm

###########################################################
### save each variable in a .csv
###########################################################
for (f in 1:length(vars_nc)){
  write.csv(output[[f]],paste0(dumpdir_nc, '/', lake_name, vars_nc[f],'.csv'),row.names=F)
}

#------------------------------------------------------------------#
#### step 3: merge NLDAS csv files together ####

#where you stored your .csv's
dumpdir_csv = './NLDAS/NLDAS_Data_2020_2021/'
# define the dump directory for your final .csv files
dumpdir_final = './NLDAS/'

#define the lakename, and the box number
LakeName = 'BVR'

cellNum = 1 #number of cells in your area of interest
box = 1 # Chosen cell of 'cellNum' from combineNLDAS.R, you'll have to look at these to figure out which one is best. use nc_get(filename.nc, 'lat'), etc to find bounding boxes and choose which one you want

# define the time range
startdatetime = '2020-01-01 00:00:00'
enddatetime = '2021-12-31 23:00:00'
loc_tz = 'UTC' #this should be the same as previous entries

#enter the timezone you would like to have the final data in. see OlsonNames() for options
local_tz = 'EST'

#save the variable names in nc files
vars_nc = c('TMP','SPFH', 'PRES', 'UGRD', 
            'VGRD', 'DLWRF', 'CONVfrac', 'CAPE', 
            'PEVAP', 'APCP', 'DSWRF')

###########################################################
### run loop to collate all data
###########################################################

# make a list of the files previously collated
files = list.files(dumpdir_csv, pattern = '.csv')

#make a null dataframme with the sequence of datetimes from above
final.box = data.frame(dateTime = seq.POSIXt(
  as.POSIXct(startdatetime, tz= loc_tz),
  as.POSIXct(enddatetime,tz=loc_tz),by = 'hour'))

#index each box csv to break out each of the cells
for (i in 1:11){
  fileIndx = grep(vars_nc[i],files)
  
  df = read_csv(paste0(dumpdir_csv,files[fileIndx[1]]),
                col_types = c('cn')) |> 
    dplyr::mutate(dateTime = ifelse(nchar(dateTime) == 10,
           paste(dateTime, "00:00:00"),
           dateTime)) |> 
    dplyr::mutate(dateTime = as.POSIXct(dateTime, tz=loc_tz,
                                        "%Y-%m-%d %H:%M:%S")) |> 
    arrange(dateTime) # chronological order   
  
  if(length(fileIndx) >1) {
    for (f in 2:length(fileIndx)){
      df2 = read.csv(paste0(dumpdir_csv, files[fileIndx[f]]))
      df = rbind(df,df2)
    }
  }
  
  # Total time series
  out = data.frame(dateTime = seq.POSIXt(
    as.POSIXct(startdatetime, tz= loc_tz),
    as.POSIXct(enddatetime,tz=loc_tz),by = 'hour')) 
  
  missingDates = out |> 
    anti_join(df)
  print(nrow(missingDates)) # Check for missing dates. 
  
  out = out %>% 
    left_join(df)
  print(nrow(out))
  # out <- distinct(out) #check for duplicate time stamps
  
  out %>% 
    mutate(dateTime = as.character(dateTime)) %>% 
    write_csv(.,paste0(dumpdir_final,LakeName,
                       format(as.POSIXct(startdatetime), '%Y-%m-%d'),
                       '_', format(as.POSIXct(enddatetime), '%Y-%m-%d'),
                       '_', vars_nc[i],'.csv'))
  
  final.box <- final.box %>% 
    left_join(out)
}

####### Create a Single Dataframe and adjust time zone###########
head(final.box)
tail(final.box)
which(duplicated(final.box)) #check for duplicate time stamps - if this list is long, something is wrong!! There should be ZERO duplicated timestamps.
# final.box <- distinct(final.box)
which(is.na(final.box$TMP)) # check for NA values

# adjust to local timezone #
final.box <- final.box %>% 
  mutate(local_dateTime = with_tz(dateTime, tzone = local_tz))
head(final.box)

# Air saturation as a function of temperature and pressure
# Used to calculate relative humidity 
qsat = function(Ta, Pa){
  ew = 6.1121*(1.0007+3.46e-6*Pa)*exp((17.502*Ta)/(240.97+Ta)) # in mb
  q  = 0.62197*(ew/(Pa-0.378*ew))                              # mb -> kg/kg
  return(q)
}


# Variable names for NLDAS2 forcing file:
# PDS_IDs:Short_Name:Full_Name [Unit]
# 63:ACPCPsfc:Convective precipitation hourly total [kg/m^2]
# 61:APCPsfc:Precipitation hourly total [kg/m^2]
# 118:BRTMPsfc:Surface brightness temperature from GOES-UMD Pinker [K]
# 157:CAPEsfc:Convective Available Potential Energy [J/kg]
# 205:DLWRFsfc:LW radiation flux downwards (surface) [W/m^2]
# 204:DSWRFsfc:SW radiation flux downwards (surface) [W/m^2]
# 101:PARsfc:PAR Photosynthetically Active Radiation from GOES-UMD Pinker [W/m^2]
# 201:PEDASsfc:Precipitation hourly total from EDAS [kg/m^2]
# 202:PRDARsfc:Precipitation hourly total from StageII [kg/m^2]
# 1:PRESsfc:Surface pressure [Pa]
# 206:RGOESsfc:SW radiation flux downwards (surface) from GOES-UMD Pinker [W/m^2]
# 51:SPFH2m:2-m above ground Specific humidity [kg/kg]
# 11:TMP2m:2-m above ground Temperature [K]
# 33:UGRD10m:10-m above ground Zonal wind speed [m/s]
# 34:VGRD10m:10-m above ground Meridional wind speed [m/s]

# Following code used to reformat dataframe to format used with GLM-AED

drivers <- final.box %>% 
  dplyr::rename(PotentialEvap = PEVAP,
                LongWave.W_m2=DLWRF,
                ShortWave.W_m2=DSWRF,
                ConvectivePrecip = CONVfrac,
                ConvectivePotentialEnergy = CAPE,
                Precipitation = APCP,
                SpecHumidity.kg_kg=SPFH,
                WindSpeed_Zonal = VGRD, 
                WindSpeed_Meridional = UGRD,
                AirTemp2m = TMP,
                SurfPressure.Pa = PRES) %>% 
  dplyr::mutate(RelHum = 100*SpecHumidity.kg_kg/qsat(AirTemp2m-273.15, SurfPressure.Pa*0.01),
                WindSpeed.m_s=sqrt(WindSpeed_Zonal^2+WindSpeed_Meridional^2),
                AirTemp.C = AirTemp2m - 273.15, 
                Rain.m_day = Precipitation*24/1000) %>% 
  dplyr::select(local_dateTime,AirTemp.C,ShortWave.W_m2,LongWave.W_m2,
                SpecHumidity.kg_kg,RelHum,WindSpeed.m_s,Rain.m_day,SurfPressure.Pa)
drivers %>% 
  mutate(local_dateTime = as.character(local_dateTime)) %>% 
  write_csv(.,paste0(dumpdir_final,LakeName, '_', 
                     format(as.POSIXct(startdatetime), '%Y_%m_%d'),
                     '_', format(as.POSIXct(enddatetime), '%Y_%m_%d'),'_hourly.csv'))

drivers |> 
  group_by(local_dateTime) %>% 
  filter(n()>1) 

drivers <- drivers |> 
  rename(time = local_dateTime,
         ShortWave = ShortWave.W_m2,
         LongWave = LongWave.W_m2,
         AirTemp = AirTemp.C,
         WindSpeed = WindSpeed.m_s,
         Rain = Rain.m_day) |> 
  select(-c(SurfPressure.Pa, SpecHumidity.kg_kg))

plot(drivers$time,drivers$Rain,type = 'l')
plot(drivers$time,drivers$ShortWave,type = 'l')

#combine with 2013-2019 NLDAS csv
nldas_13_19 <- read.csv(paste0(getwd(),"/NLDAS/BVR_GLM_NLDAS_010113_123119_GMTadjusted.csv")) |> 
  mutate(time = as.POSIXct(time, tz=loc_tz,
                           "%Y-%m-%d %H:%M"))

all_NLDAS <- bind_rows(drivers, nldas_13_19)

#save combined nldas file
write.csv(all_NLDAS, "./inputs/BVR_GLM_NLDAS_010113_123121_GMTadjusted.csv", rownames=F)

#quick plots to make sure things look okay
plot(all_NLDAS$time,all_NLDAS$Rain,type = 'l')
plot(all_NLDAS$time,all_NLDAS$ShortWave,type = 'l')
plot(all_NLDAS$time,all_NLDAS$LongWave,type = 'l')
plot(all_NLDAS$time,all_NLDAS$RelHum,type = 'l')
plot(all_NLDAS$time,all_NLDAS$WindSpeed,type = 'l')
