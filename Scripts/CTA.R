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

#turn transformed community data into b-c distance matrix 
zoop_hell <- as.matrix(vegan::vegdist(zoops_dens, method='hellinger'))

#check to make sure the dissimilarity matrix satisfies triangle inequality
is.metric(zoop_hell)
#note - if I hellinger transform then use bc, does not satisfy triangle inequality

#now do NMDS w/ 3 dimensions 
NMDS_bray <- vegan::metaMDS(zoop_hell, distance='bray', k=2, trymax=20, 
                                   autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_bray$stress # 0.08 for k=3
#stress=0.14 for k=2

ord <- vegan::ordiplot(NMDS_bray,display = c('sites','species'),
                       choices = c(1,2),type = "n")

coords <- data_frame(x = ord$sites[,1], y = ord$sites[,2])
surveys=rep(c(1,2,3,4,5),6)

#calculate length of each trajectory
lengths <-trajectoryLengths2D(coords, all_zoops_nmds$year, surveys) |> 
  mutate(year = c(2014,2015,2016,2019,2020,2021)) |> 
  pivot_longer(cols = S1:S4, names_to = "segment") 

ggplot(data = subset(lengths, year %in% c(2020)), 
       aes(segment, value, color = as.factor(year))) +
  geom_line(aes(group=year)) + 
  geom_point(size = 4) +
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

ggplot(lengths, aes(as.factor(year), value, group=as.factor(year))) + 
  geom_boxplot() + theme_bw() + ylab("mean length") + xlab("year")

ggplot(lengths, aes(as.factor(year), Trajectory, group=as.factor(year))) + 
  geom_point(size=5) + theme_bw() + ylab("total length") + xlab("year")

#anova to see if any of these are sig different (4 lengths/yr)
res_aov <- aov(value ~ segment,
               data = lengths
)
summary(res_aov)
#no sig differences... (unless I look across segments...)

stats::TukeyHSD(res_aov)
#significant difference between may-jun vs. jun-jul lengths

#calculate angle of each trajectory
angles <- trajectoryAngles2D(coords, all_zoops_nmds$year, surveys) |> 
  mutate(year = c("2014","2015","2016","2019","2020","2021")) |> 
  pivot_longer(cols = c("t1-t2":"t3-t4"),names_to = "segment")

ggplot(angles, aes(year, as.numeric(value), group=year)) + 
  geom_boxplot() + theme_bw() + ylab("angle")

#calculate speed of each trajectory (I think average speed is just avg length / (150-30 days)
lengths$speed <- lengths$value / 30

ggplot(lengths, aes(segment, speed, color = as.factor(year))) +
  geom_line(aes(group=year)) + geom_point() +
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

ggplot(lengths, aes(as.factor(year), speed, group=as.factor(year))) + 
  geom_boxplot() + theme_bw() + ylab("speed") + xlab("year")


#summarize across years
yearly_stats <- data.frame("year" = c(2014,2015,2016,2019,2020,2021))
yearly_stats$length <- c(sum(lengths$value[lengths$year==2014]),
                         sum(lengths$value[lengths$year==2015]),
                         sum(lengths$value[lengths$year==2016]),
                         sum(lengths$value[lengths$year==2019]),
                         sum(lengths$value[lengths$year==2020]),
                         sum(lengths$value[lengths$year==2021]))
yearly_stats$speed <- c(mean(lengths$speed[lengths$year==2014]),
                         mean(lengths$speed[lengths$year==2015]),
                         mean(lengths$speed[lengths$year==2016]),
                         mean(lengths$speed[lengths$year==2019]),
                         mean(lengths$speed[lengths$year==2020]),
                         mean(lengths$speed[lengths$year==2021]))

barplot(yearly_stats$length, ylab = "Length",
        names.arg = c("2014","2015","2016",
                      "2019","2020","2021"))

barplot(yearly_stats$speed, ylab = "Speed",
        names.arg = c("2014","2015","2016",
                      "2019","2020","2021"))

#calculate direction of each trajectory
dir <- trajectoryDirectionality(zoop_hell, all_zoops_nmds$year)

barplot(dir, ylab = "Direction")

#convergence vs. divergence
trajectoryConvergence(zoop_hell, all_zoops_nmds$year)
#p-values are all > 0.05 so there is divergence among all trajectories

#now check across months
trajectoryConvergence(zoop_hell, all_zoops_nmds$month)
#may and jul are converging? p=0.02 but not sure what that really means

ord <- vegan::ordiplot(NMDS_bray,display = c('sites','species'),
                       choices = c(1,2),type = "n")

year <- ggordiplots::gg_ordiplot(ord, all_zoops_nmds$year,
                                  kind = "sd", ellipse=FALSE, hull = FALSE, 
                                  plot = FALSE, pt.size=0.9) 
year$plot + geom_point() + geom_line() + theme_bw() + 
  scale_fill_manual("",values=viridis::viridis(6, option="D"))+
  scale_color_manual("",values=viridis::viridis(6, option="D"),
                     label=c('2014','2015',"2016","2019","2020","2021")) +
  geom_path(data=year$df_ord,
               aes(x=x, y=y, group = Group),
               arrow=arrow(length=unit(0.3,"cm")),
               color=c(rep("#440154FF",5),rep("#414487FF",5),
                       rep("#2A788EFF",5),rep("#22A884FF",5),
                       rep("#7AD151FF",5),rep("#FDE725FF",5)))

#----------------------------------------------------------------------#
#jpeg("Figures/NMDS_k=2_yearly_dens.jpg")
par(mfrow = c(2, 3))
par(cex = 0.6)
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

for (i in 1:6) {
  ord <- vegan::ordiplot(NMDS_bray,display = c('sites','species'),
                         choices = c(1,2),type = "none", 
                         xlim= c(-0.7, 0.7), ylim= c(-0.3, 0.6),
                         xaxt = 'n', yaxt = 'n')
  
  lines(NMDS_bray$points[all_zoops_nmds$year == 
                                 unique(all_zoops_nmds$year)[i] ,1], 
        NMDS_bray$points[(all_zoops_nmds$year ==
                                  unique(all_zoops_nmds$year)[i]),2], 
        col=viridis::viridis(6, option="D")[i])
  
  points(NMDS_bray$points[all_zoops_nmds$year == 
                                  unique(all_zoops_nmds$year)[i] ,1], 
         NMDS_bray$points[(all_zoops_nmds$year == 
                                   unique(all_zoops_nmds$year)[i]),2], 
         col=viridis::viridis(6, option="D")[i],
         pch = as.character(1:5), font=2, cex=3)
  
  mtext(c("2014","2015","2016","2019","2020","2021")[i], side = 3, line = -2, 
        adj = 0.05, cex = 1.5, col = "black")
  
  mtext("NMDS1", side = 1, outer = TRUE, cex = 1.2, line = 2.2,
        col = "black")
  mtext("NMDS2", side = 2, outer = TRUE, cex = 1.2, line = 2.2,
        col = "black")
  
}

#dev.off()    




#----------------------------------------------------------------------#
#try length vs. env variables?

env <- read.csv("Output/env.csv",header=T) |> 
  arrange(year, month)

env_years <- env |> 
  group_by(year) |> 
  summarise_all(list(mean = mean), na.rm=T) |> 
  select(-month_mean)

env_traj <- inner_join(yearly_stats, env, by="year")
env_traj_final <- inner_join(env_traj, env_years, by="year")

env_traj_final$dir <- ifelse(env_traj_final$year==2014, dir[1], 
                      ifelse(env_traj_final$year==2015, dir[2],
                      ifelse(env_traj_final$year==2016, dir[3],
                      ifelse(env_traj_final$year==2019, dir[4],
                      ifelse(env_traj_final$year==2020, dir[5],
                      ifelse(env_traj_final$year==2021, dir[6], NA))))))

ggplot(env_traj_final, aes(Temp_C_epi, length, col=as.factor(year))) +
  geom_point() + theme_bw() + geom_line() +
  geom_point(aes(Temp_C_epi_mean, length,col = as.factor(year), size = 3))+
  guides(size = "none") 

ggplot(env_traj_final, aes(Temp_C_hypo, length, col=as.factor(year))) +
  geom_point() + theme_bw() + geom_line() +
  geom_point(aes(Temp_C_hypo_mean, length,col = as.factor(year), size = 3))+
  guides(size = "none") 

ggplot(env_traj_final, aes(DO_mgL_epi, length, col=as.factor(year))) +
  geom_point() + theme_bw() + geom_line() +
  geom_point(aes(DO_mgL_epi_mean, length,col = as.factor(year), size = 3))+
  guides(size = "none") 

ggplot(env_traj_final, aes(TN_ugL_epi, length, col=as.factor(year))) +
  geom_point() + theme_bw() + geom_line() +
  geom_point(aes(TN_ugL_epi_mean, length,col = as.factor(year), size = 3))+
  guides(size = "none") 

ggplot(env_traj_final, aes(TN_ugL_hypo, length, col=as.factor(year))) +
  geom_point() + theme_bw() + geom_line() +
  geom_point(aes(TN_ugL_hypo_mean, length,col = as.factor(year), size = 3))+
  guides(size = "none") 

ggplot(env_traj_final, aes(TP_ugL_epi, length, col=as.factor(year))) +
  geom_point() + theme_bw() + geom_line() +
  geom_point(aes(TP_ugL_epi_mean, length,col = as.factor(year), size = 3))+
  guides(size = "none") 

ggplot(env_traj_final, aes(TP_ugL_hypo, length, col=as.factor(year))) +
  geom_point() + theme_bw() + geom_line() +
  geom_point(aes(TP_ugL_hypo_mean, length,col = as.factor(year), size = 3))+
  guides(size = "none") 

ggplot(env_traj_final, aes(waterlevel, length, col=as.factor(year))) +
  geom_point() + theme_bw() + geom_line() +
  geom_point(aes(waterlevel_mean, length,col = as.factor(year), size = 3))+
  guides(size = "none") 

ggplot(env_traj_final, aes(therm_depth, length, col=as.factor(year))) +
  geom_point() + theme_bw() + geom_line() +
  geom_point(aes(therm_depth_mean, length,col = as.factor(year), size = 3))+
  guides(size = "none") 

ggplot(env_traj_final, aes(oxy_depth, length, col=as.factor(year))) +
  geom_point() + theme_bw() + geom_line() +
  geom_point(aes(oxy_depth_mean, length,col = as.factor(year), size = 3))+
  guides(size = "none") 

ggplot(env_traj_final, aes(SS, length, col=as.factor(year))) +
  geom_point() + theme_bw() + geom_line() +
  geom_point(aes(SS_mean, length,col = as.factor(year), size = 3))+
  guides(size = "none") 

ggplot(env_traj_final, aes(BF, length, col=as.factor(year))) +
  geom_point() + theme_bw() + geom_line() +
  geom_point(aes(BF_mean, length,col = as.factor(year), size = 3))+
  guides(size = "none") 

ggplot(env_traj_final, aes(AirTemp, length, col=as.factor(year))) +
  geom_point() + theme_bw() + geom_line() +
  geom_point(aes(AirTemp_mean, length,col = as.factor(year), size = 3))+
  guides(size = "none") 

ggplot(env_traj_final, aes(WindSpeed, length, col=as.factor(year))) +
  geom_point() + theme_bw() + geom_line() +
  geom_point(aes(WindSpeed_mean, length,col = as.factor(year), size = 3))+
  guides(size = "none") 
