#script to calculate quantitative values for each trajectory based on De Caceres et al. 2019

#load packages
pacman::p_load(ecotraj, tidyverse)

#read in zoop dens csv
all_zoops_dens <- read.csv("Output/all_zoops_dens.csv",header = TRUE)

taxa <- unique(all_zoops_dens$Taxon)

#taxa as cols, dates as rows, average by month
all_zoops_nmds <- all_zoops_dens |> 
  select(DateTime, Taxon, dens) |> 
  filter(Taxon %in% taxa) |> 
  mutate(DateTime = as.Date(DateTime)) |> 
  pivot_wider(names_from = Taxon, values_from = dens) |> 
  mutate(year = format(DateTime, "%Y"),
         month = format(DateTime, "%m")) |> 
  mutate_all(~replace(., is.na(.), 0)) |>  #replace NA with 0
  ungroup() |> group_by(year, month) |> 
  summarise(Bosmina = mean(Bosmina),
            Ceriodaphnia = mean(Ceriodaphnia),
            Daphnia = mean(Daphnia),
            Cyclopoida = mean(Cyclopoida),
            Nauplii = mean(Nauplii),
            Ascomorpha = mean(Ascomorpha),
            Conochilus = mean(Conochilus),
            Keratella = mean(Keratella),
            Kellicottia = mean(Kellicottia),
            Polyarthra = mean(Polyarthra)) |> 
  ungroup()

#only keep may-sep samples and drop 2022
all_zoops_nmds <- all_zoops_nmds |> 
  filter(month %in% c("05","06","07","08","09"), 
         !year %in% c("2022"))
#only 5 months per year bc needs to be equal for pairwise correlation matrix

#select only data cols
zoops_dens <- all_zoops_nmds |> select(Bosmina:Polyarthra)

#hellinger transform data
zoop_dens_trans <- labdsv::hellinger(zoops_dens)

#turn transformed community data into b-c distance matrix 
zoop_bray <- as.matrix(vegan::vegdist(zoop_dens_trans, method='euclidean'))

#calculate length of each trajectory
lengths <-trajectoryLengths(zoop_bray, all_zoops_nmds$year) |> 
  mutate(years = c(2014,2015,2016,2019,2020,2021)) |> 
  pivot_longer(cols = S1:S4, names_to = "segment") 

ggplot(lengths, aes(segment, value, color = as.factor(years), cex = value)) +
  #geom_line(aes(group=years)) + 
  geom_point() +
  ylab("Length") + theme_bw() + xlab("") +
  guides(cex = "none") +
  scale_x_discrete(labels=c(S1 = "May-Jun", 
                            S2 = "Jun-Jul",
                            S3 = "Jul-Aug",
                            S4 = "Aug-Sep")) +
  theme(text = element_text(size=14), 
        axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "vertical",
        legend.title = element_blank(),
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "right", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,1,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) 

#calculate angle of each trajectory
trajectoryAngles(zoop_bray, all_zoops_nmds$year)

#calculate speed of each trajectory (I think average speed is just avg length / (150-30 days)
lengths$speed <- lengths$value / 30

ggplot(lengths, aes(segment, speed, color = as.factor(years))) +
  geom_line(aes(group=years)) + geom_point() +
  ylab("Speed") + theme_bw() + xlab("") +
  scale_x_discrete(labels=c(S1 = "May-Jun", 
                            S2 = "Jun-Jul",
                            S3 = "Jul-Aug",
                            S4 = "Aug-Sep")) +
  theme(text = element_text(size=14), 
        axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "vertical",
        legend.title = element_blank(),
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = c(0.93, 0.86), 
        legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,1,0,0), 'lines'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) 

#summarize across years
yearly_stats <- data.frame("year" = c(2014,2015,2016,2019,2020,2021))
yearly_stats$length <- c(sum(lengths$Trajectory[lengths$years==2014]),
                         sum(lengths$Trajectory[lengths$years==2015]),
                         sum(lengths$Trajectory[lengths$years==2016]),
                         sum(lengths$Trajectory[lengths$years==2019]),
                         sum(lengths$Trajectory[lengths$years==2020]),
                         sum(lengths$Trajectory[lengths$years==2021]))
yearly_stats$speed <- c(mean(lengths$speed[lengths$years==2014]),
                         mean(lengths$speed[lengths$years==2015]),
                         mean(lengths$speed[lengths$years==2016]),
                         mean(lengths$speed[lengths$years==2019]),
                         mean(lengths$speed[lengths$years==2020]),
                         mean(lengths$speed[lengths$years==2021]))

barplot(yearly_stats$length, ylab = "Length",
        names.arg = c("2014","2015","2016",
                      "2019","2020","2021"))

barplot(yearly_stats$speed, ylab = "Speed",
        names.arg = c("2014","2015","2016",
                      "2019","2020","2021"))

#calculate direction of each trajectory
dir <- trajectoryDirectionality(zoop_bray, all_zoops_nmds$year)

barplot(dir, ylab = "Direction")


