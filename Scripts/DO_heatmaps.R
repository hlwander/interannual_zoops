# Makes heatmaps of the CTD data in Beaverdam reservoir for all 6 years of zoop data

# load libraries
pacman::p_load(akima, dplyr, ggplot2, tidyverse, reshape2, lubridate, ggpubr,
               gridExtra, grid, colorRamps, RColorBrewer, rLakeAnalyzer)

#list of dates between May-Sep for each year
date_list <- c(seq(as.Date("2014-05-01"), as.Date("2014-09-30"), by="days"),
                seq(as.Date("2015-05-01"), as.Date("2015-09-30"), by="days"),
                seq(as.Date("2016-05-01"), as.Date("2016-09-30"), by="days"),
                seq(as.Date("2019-05-01"), as.Date("2019-09-30"), by="days"),
                seq(as.Date("2020-05-01"), as.Date("2020-09-30"), by="days"),
                seq(as.Date("2021-05-01"), as.Date("2021-09-30"), by="days"))


#read in ctd
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/13/27ceda6bc7fdec2e7d79a6e4fe16ffdf" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1, timeout = max(300, getOption("timeout"))))

ctd <-read.csv(infile1,header=T) |> 
  mutate(DateTime = as.Date(DateTime)) |>
  filter(DateTime %in% date_list & 
           Depth_m > 0 & Site ==50 & Reservoir %in% c("BVR") &
           !is.na(DO_mgL)) |> 
  select(Reservoir, Site, DateTime, Depth_m, DO_mgL) 

#Trim dataset
depths <- seq(0,12, by = .1)
newDepths <- depths
ctd_final<- ctd |>  group_by(DateTime, Reservoir) |> 
  slice(which.min(abs(as.numeric(Depth_m) - depths[1]))) #Create a new dataframe
ctd_final$Depth_m <- newDepths[1]
for (i in 2:length(depths)){ #loop through all depths and add the closest values to the final dataframe
  ctd_atThisDepth <- ctd |>  group_by(DateTime, Reservoir) |> 
    slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
  ctd_atThisDepth$Depth_m <- newDepths[i]
  ctd_final <- rbind(ctd_final,ctd_atThisDepth)
}


inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/198/11/6e5a0344231de7fcebbe6dc2bed0a1c3" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))

ysi <- read.csv(infile1,header=T) |> 
  mutate(DateTime = as.Date(DateTime)) 

all_do <- ctd_final |> 
  full_join(ysi |> select("Reservoir",     
                         "Site",     
                         "DateTime",     
                         "Depth_m",  
                         "DO_mgL") |> 
              filter(!is.na(DO_mgL) &
                     Site == 50 &
                     Reservoir == "BVR" &
                     DateTime %in% date_list))

ysi%>%
  filter(Site==50)%>%
  ggplot(aes(x = as.Date(DateTime), y = DO_mgL))+
  geom_point()

interp_do = data.frame(x = NA, y = NA, z = NA, Reservoir = NA)
years = c(2014:2016,2019:2021) 

for(i in 1:length(years)){ 
  do_bvr<- all_do  |> 
    mutate(Year = year(DateTime))  |> 
    filter(Year == years[i],
           Reservoir == "BVR",
           Depth_m <= 12)

  interp_do <- interp2xyz(interp(do_bvr$DateTime, do_bvr$Depth_m, do_bvr$DO_mgL,
                                     xo = seq(min(do_bvr$DateTime), max(do_bvr$DateTime), 1),
                                     yo = seq(min(do_bvr$Depth_m), max(do_bvr$Depth_m), 
                                              by = .01), duplicate = "mean"),data.frame = T) |> 
    filter(!is.na(z))

#add manual groups for heatmap 
interp_do$groups <- cut(interp_do$z,    
                        breaks = c(0, 2, 4, 6, 8, 10, 12, 14, 16))

#DO heatmap
do_temp <- interp_do  |> 
  mutate(x = as.Date(x, origin = "1970-01-01"),
         Year = year(x))  |> 
  ggplot(aes(x=x, y=y,fill=groups)) +
  geom_tile()+ 
  scale_y_reverse(expand = c(0,0))+
  labs(x = "", y = "Depth (m)", title = "Dissolved Oxygen ",fill=expression('(mg/L)'))+
  scale_fill_manual(breaks=levels(interp_do$groups),
                    values=c("#CC3300","#FF9900","#CCFF00","#66FF66","#66FF99","#33FFFF","#3399FF","#0066CC"))+
  scale_x_date(expand = c(0,0))+
  theme(panel.border = element_rect(fill = NA))

  ggsave(sprintf("Figures/BVR_DO_heatmap_%s.jpg", years[i]), width=6, height=6)

}
