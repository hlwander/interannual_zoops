#####################################################################
#Script for downloading NLDAS2 data for meteorological hourly forcing
#Note that as of 1Aug2024, GRIB-1 files are no longer available
#this workflow downloads, subsets, and merges met data from .nc files
#Created 29 Oct 2025 by HLW
#####################################################################

library(ncdf4)
library(lubridate)
library(dplyr)
library(data.table)

# -----------------------------
# Setup
# -----------------------------
dumpdir_nc <- "./NLDAS/NLDAS_Data_2020_2023"
dir.create(dumpdir_nc, recursive = TRUE, showWarnings = FALSE)

username <- "hwander"
password <- "Zooplankton2025!"

startdatetime <- '2020-01-01 00:00:00'
enddatetime   <- '2024-01-01 23:00:00'
loc_tz <- 'UTC'

lat_range <- c(37.31, 37.32)
lon_range <- c(-79.82, -79.81)

var_names <- c("Tair", "Qair","Wind_E", "PSurf",
               "Wind_N", "LWdown", "Rainf", "SWdown")  

#list of all the vars
#Tair = 2-m above ground temperature; K
#Qair = 2-m above ground specific humidity; kg/kg
#PSurf = surface pressure; Pa
#Wind_E = 10-m above ground zonal wind speed; m/s
#Wind_N = 10-m above ground meridional wind speed; m/s
#LWdown = longwave radiation flux downwards (surface); W/m2
#CRainf_frac = fraction of total precipitation that is convective
#CAPE = convective available potential energy; J/kg
#PotEvap = potential evaporation; kg/m2
#Rainf = total precipitation; kg/m2
#SWdown = shortwave radiation flux downwards (surface); W/m2

out.ts <- seq.POSIXt(as.POSIXct(startdatetime, tz = loc_tz),
                     as.POSIXct(enddatetime, tz = loc_tz), by = 'hour')

all_data <- data.frame()

# -----------------------------
# Loop over each hour
# -----------------------------
for (i in seq_along(out.ts)) {
  
  yearOut <- year(out.ts[i])
  monthOut <- format(out.ts[i], "%m")
  dayOut <- format(out.ts[i], "%d")
  hourOut <- format(out.ts[i], "%H%M")
  doyOut <- format(out.ts[i], '%j')
  
  filename <- paste0("NLDAS_FORA0125_H.A", yearOut, monthOut, dayOut, ".", hourOut, ".020.nc")
  remote_url <- paste0(
    "https://hydro1.gesdisc.eosdis.nasa.gov/data/NLDAS/NLDAS_FORA0125_H.2.0/",
    yearOut, "/", doyOut, "/", filename)
  
  destfile <- file.path(dumpdir_nc, filename)
  
  # Skip if file already exists and is readable
  if (file.exists(destfile)) {
    tryCatch({
      nc_test <- nc_open(destfile)
      nc_close(nc_test)
      message("✅ Already downloaded and valid: ", filename)
      next
    }, error = function(e) {
      message("⚠️ File exists but is corrupt, re-downloading: ", filename)
      file.remove(destfile)
    })}
  
  # Download URL
  remote_url <- paste0(
    "https://hydro1.gesdisc.eosdis.nasa.gov/data/NLDAS/NLDAS_FORA0125_H.2.0/",
    yearOut, "/", doyOut, "/", filename)
  
  # Retry download up to 3 times
  success <- FALSE
  for (attempt in 1:3) {
    message("⬇️ Attempt ", attempt, " to download ", filename)
    cmd <- paste(
      "wget --auth-no-challenge",
      paste0("--user=", username),
      paste0("--password=", password),
      "-O", shQuote(destfile),
      shQuote(remote_url))
    
    res <- tryCatch(system(cmd, intern = TRUE), warning = function(w) NULL, error = function(e) NULL)
    
    if (file.exists(destfile) && file.info(destfile)$size > 10000) {  # sanity check: file >10KB
      message("✅ Downloaded successfully: ", filename)
      success <- TRUE
      break
    } else {
      message("⚠️ Download failed, retrying...")
      Sys.sleep(3)
    }
  }
  
  if (!success) {
    message("❌ Skipping after 3 failed attempts: ", filename)
    next}
  
  # Try opening to ensure it's valid NetCDF
  tryCatch({
    nc <- nc_open(destfile)
    message("📦 Valid NetCDF: ", filename)
    nc_close(nc)
  }, error = function(e) {
    message("🚫 Invalid NetCDF file, removing: ", filename)
    file.remove(destfile)
  })

  # Open NetCDF
  nc <- nc_open(destfile)
  lon_idx <- which(nc$dim$lon$vals >= lon_range[1] & nc$dim$lon$vals <= lon_range[2])
  lat_idx <- which(nc$dim$lat$vals >= lat_range[1] & nc$dim$lat$vals <= lat_range[2])
  
  lon_sub <- nc$dim$lon$vals[lon_idx]
  lat_sub <- nc$dim$lat$vals[lat_idx]
  
  # Loop over variables
  for (var_name in var_names) {
    
    full_data <- ncvar_get(nc, var_name)
    dims <- dim(full_data)
    
    # subset based on dimensions
    if (length(dims) == 2) {
      subset_data <- full_data[lon_idx, lat_idx]
    } else if (length(dims) == 3) {
      subset_data <- full_data[lon_idx, lat_idx, 1]  # time = 1 since each file = 1 hour
    } else {
      stop("Unexpected variable dimensions")
    }
    
    # expand grid and reshape
    df <- expand.grid(lon = lon_sub, lat = lat_sub)
    df$time <- out.ts[i]
    df$variable <- var_name
    df$value <- as.vector(subset_data)
    
    # append
    all_data <- bind_rows(all_data, df)
  }
  
  nc_close(nc)
}

# Generate expected filenames
expected_files <- sapply(out.ts, function(t) {
  yearOut <- year(t)
  monthOut <- format(t, "%m")
  dayOut <- format(t, "%d")
  hourOut <- format(t, "%H%M")
  paste0("NLDAS_FORA0125_H.A", yearOut, monthOut, dayOut, ".", hourOut, ".020.nc")
})

# List files actually present
existing_files <- list.files(dumpdir_nc, pattern = "\\.nc$", full.names = FALSE)

# Find missing ones
missing_files <- setdiff(expected_files, existing_files)

# -------------------------------
# rebuild all_data from .nc files
#-------------------------------
# Ensure local_tz used for output 
local_tz <- "EST"

# Create helper to parse time from filename
parse_time_from_filename <- function(fn) {
  base <- basename(fn)
  # Accept several possible filename patterns using regex:
  # pattern A: ...A20211215.0200.020.nc  (your usual)
  mA <- regexec("A(\\d{8})\\.(\\d{4})\\.020", base)
  gA <- regmatches(base, mA)[[1]]
  if (length(gA) == 3) {
    ymd <- gA[2]
    hhmm <- gA[3]
    out <- tryCatch({
      as.POSIXct(paste0(substr(ymd,1,4),"-",substr(ymd,5,6),"-",substr(ymd,7,8)," ",
                        substr(hhmm,1,2), ":", substr(hhmm,3,4), ":00"),
                 tz = "UTC")
    }, error = function(e) NA)
    if (!is.na(out)) return(out)
  }
  # pattern B: fallback for filenames like ..._20211215_0200.nc
  mB <- regexec("(\\d{4})(\\d{2})(\\d{2})[._-](\\d{2})(\\d{2})", base)
  gB <- regmatches(base, mB)[[1]]
  if (length(gB) >= 6) {
    out <- tryCatch({
      as.POSIXct(sprintf("%s-%s-%s %s:%s:00", gB[1], gB[2], gB[3], gB[4], gB[5]), tz = "UTC")
    }, error = function(e) NA)
    if (!is.na(out)) return(out)
  }
  return(NA)
}

nc_time_to_posix <- function(nc) {
  if (!("time" %in% names(nc$dim))) return(NA)
  tvals <- tryCatch(ncvar_get(nc, "time"), error = function(e) NULL)
  if (is.null(tvals) || length(tvals) < 1) return(NA)
  tu <- tryCatch(nc$dim$time$units, error = function(e) NULL)
  if (is.null(tu)) tu <- tryCatch(ncatt_get(nc, "time", "units")$value, error = function(e) NULL)
  if (is.null(tu) || !is.character(tu)) return(NA)
  # parse units: "<units> since YYYY-MM-DD HH:MM:SS"
  m <- regexec("^(\\w+) since (\\d{4}-\\d{2}-\\d{2}(?:[ T]\\d{2}:\\d{2}:?\\d{0,2})?)", tu)
  mm <- regmatches(tu, m)[[1]]
  if (length(mm) < 3) return(NA)
  gran <- mm[2]  
  origin <- mm[3]
  # normalize origin string
  origin <- sub("T", " ", origin)
  origin <- ifelse(nchar(origin) == 10, paste0(origin, " 00:00:00"), origin)
  origin_posix <- tryCatch(as.POSIXct(origin, tz = "UTC"), error = function(e) NA)
  if (is.na(origin_posix)) return(NA)
  # convert based on granularity
  if (grepl("^hour", gran, ignore.case = TRUE)) {
    out_time <- origin_posix + tvals * 3600
  } else if (grepl("^day", gran, ignore.case = TRUE)) {
    out_time <- origin_posix + tvals * 86400
  } else if (grepl("^minute", gran, ignore.case = TRUE)) {
    out_time <- origin_posix + tvals * 60
  } else {
    # fallback: try as seconds
    out_time <- origin_posix + tvals
  }
  return(out_time)
}

get_file_time <- function(fname, nc = NULL) {
  # try filename
  ft <- parse_time_from_filename(fname)
  if (!is.na(ft)) return(ft)
  # try netcdf time variable if nc is provided
  if (!is.null(nc)) {
    nc_times <- nc_time_to_posix(nc)
    if (!identical(nc_times, NA) && length(nc_times) >= 1) {
      # if the file contains multiple times, typically we want the first (files are hourly)
      return(as.POSIXct(nc_times[1], tz = "UTC"))
    }
  }
  warning("Could not determine time for file: ", fname)
  return(NA)
}

# Build expected filename list and get existing ones
expected_files <- sapply(out.ts, function(t) {
  yearOut <- year(t)
  monthOut <- format(t, "%m")
  dayOut <- format(t, "%d")
  hourOut <- format(t, "%H%M")
  paste0("NLDAS_FORA0125_H.A", yearOut, monthOut, dayOut, ".", hourOut, ".020.nc")})

existing_files <- list.files(dumpdir_nc, pattern = "\\.nc$", full.names = FALSE)

# Work only on the intersection of expected_files & existing_files, in time order
to_process <- intersect(expected_files, existing_files)
to_process <- sort(to_process)

if(length(to_process) == 0){
  stop("No expected downloaded files found in ", dumpdir_nc)
}

# container to accumulate per-file/variable data (use list for memory efficiency)
acc <- list()
acc_i <- 1

for(fname in to_process){
  fpath <- file.path(dumpdir_nc, fname)
  # quick sanity checks
  if(!file.exists(fpath) || file.info(fpath)$size < 2000){
    message("SKIP (missing or tiny): ", fname); next
    }
  nc <- tryCatch(nc_open(fpath), error = function(e){ message("BAD NC (skip): ", fname, " -> ", e$message); NULL })
  if(is.null(nc)) next
  
  # find indices within bounding box (if none, skip file)
  lon_vals <- nc$dim$lon$vals
  lat_vals <- nc$dim$lat$vals
  lon_idx <- which(lon_vals >= lon_range[1] & lon_vals <= lon_range[2])
  lat_idx <- which(lat_vals >= lat_range[1] & lat_vals <= lat_range[2])
  if(length(lon_idx) == 0 || length(lat_idx) == 0){
    message("NO GRID CELLS IN BBOX for file: ", fname)
    nc_close(nc); next
  }
  
  file_time <- get_file_time(fname, nc)
  
  # loop over variables requested
  for(var in var_names){
    if(!(var %in% names(nc$var))){
      message("Variable not present, skipping: ", var, " in ", fname)
      next
    }
    full_data <- tryCatch(ncvar_get(nc, var), error = function(e){ message("Can't read var ", var, " in ", fname); NULL })
    if(is.null(full_data)) next
    
    # handle typical shapes: [lon,lat] or [lon,lat,time]
    if(length(dim(full_data)) == 2){
      sub_mat <- full_data[lon_idx, lat_idx]
    } else if(length(dim(full_data)) == 3){
      # use first time slice (each file is one hour)
      sub_mat <- full_data[lon_idx, lat_idx, 1]
    } else {
      message("Unexpected dims for var ", var, " in ", fname, ": ", paste(dim(full_data), collapse="x"))
      next
    }
    
    # build small df for this var/file
    df <- expand.grid(lon = lon_vals[lon_idx], lat = lat_vals[lat_idx])
    df$time <- file_time            # UTC POSIXct
    df$variable <- var
    df$value <- as.vector(sub_mat)
    
    acc[[acc_i]] <- df
    acc_i <- acc_i + 1
  }
  
  nc_close(nc)
  message("Processed: ", fname)
}

# bind all pieces
final_all_data <- dplyr::bind_rows(acc)

# Optional: save intermediate combined all_data (if you want)
write.csv(final_all_data, file.path(dumpdir_nc, "NLDAS_subset_combined_multivar.csv"), row.names = FALSE)

# -----------------------------
# reformat and save csv
# -----------------------------

qsat <- function(Ta, Pa){
  ew <- 6.1121*(1.0007+3.46e-6*Pa)*exp((17.502*Ta)/(240.97+Ta)) # in mb
  q  <- 0.62197*(ew/(Pa-0.378*ew))                              # mb -> kg/kg
  return(q)}

# pivot to wide so variables are columns
drivers_wide <- final_all_data |>
  tidyr::pivot_wider(names_from = variable, values_from = value) |>
  # ensure time is POSIXct and in UTC
  mutate(time = as.POSIXct(time, tz = "UTC"))

# create local_dateTime column (convert UTC to desired local tz)
drivers_wide <- drivers_wide |>
  mutate(local_dateTime = lubridate::with_tz(time, tzone = local_tz))

#reformat dataframe to match GLM-AED driver file
drivers <- drivers_wide |>
  dplyr::select(-c(lat,lon,time)) |>
  dplyr::rename(time = local_dateTime,
                LongWave.W_m2=LWdown,
                ShortWave.W_m2=SWdown,
                Precipitation = Rainf,
                SpecHumidity.kg_kg=Qair,
                WindSpeed_Zonal = Wind_E, 
                WindSpeed_Meridional = Wind_N,
                AirTemp2m = Tair,
                SurfPressure.Pa = PSurf) |>
  dplyr::mutate(RelHum = 100*SpecHumidity.kg_kg/qsat(AirTemp2m-273.15, SurfPressure.Pa*0.01),
                WindSpeed.m_s=sqrt(WindSpeed_Zonal^2+WindSpeed_Meridional^2),
                AirTemp.C = AirTemp2m - 273.15, 
                Rain.m_day = Precipitation*24/1000) |> #m/d
  dplyr::select(time,AirTemp.C,ShortWave.W_m2,LongWave.W_m2,
                SpecHumidity.kg_kg,RelHum,WindSpeed.m_s,Rain.m_day) |>
  filter(time >= as.POSIXct("2020-01-01") &
           time < as.POSIXct("2024-01-01"))

#rename drivers to match glm aed col names
drivers <- drivers |>
  dplyr::rename(ShortWave = ShortWave.W_m2,
                LongWave = LongWave.W_m2,
                AirTemp = AirTemp.C,
                WindSpeed = WindSpeed.m_s,
                Rain = Rain.m_day)
  
#combine with 2013-2019 NLDAS csv
nldas_13_19 <- read.csv(paste0(getwd(),"/NLDAS/BVR_GLM_NLDAS_010113_123119_GMTadjusted.csv")) |> 
  mutate(time = as.POSIXct(time, tz=loc_tz,
                           "%Y-%m-%d %H:%M"))

all_NLDAS <- bind_rows(drivers, nldas_13_19) |>
  arrange(time)

#write new houry driver file
write.csv(all_NLDAS, "inputs/BVR_GLM_NLDAS_2013-2023_hourly.csv", row.names = FALSE)

