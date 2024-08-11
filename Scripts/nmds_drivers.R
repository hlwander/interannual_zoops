# Env drivers of nmds cluster years
# 23 Jul 2024

pacman::p_load(tidyverse, Hmisc)

#cb friendly year palette
year_cols <- c("#a13637","#06889b", "#facd60", "#f44034", "#011f51", "#fdfa66")

# Which drivers explain patterns in zooplankton community patterns among years?

#read in driver df
zoop_drivers <- read.csv("Output/all_drivers.csv") |> 
  rename("wind speed" = WindSpeed,
         "air temp" = AirTemp,
         "residence time" = res_time_d,
         "water level" = waterlevel,
         "thermocline depth" = therm_depth,
         "oxycline depth" = oxy_depth,
         "buoyancy frequency" = BF,
         "Schmidt stability" = SS,
         "epilimnetic DO" = DO_mgL_epi,
         "epilimnetic TN" = TN_ugL_epi,
         "epilimnetic TP" = TP_ugL_epi,
         "Secchi depth" = secchi,
         "hypolimnetic temp" = Temp_C_hypo,
         "total phytoplankton biomass" = Total_ugL)

#group years based on second-stage nmds groupings
zoop_drivers$traj <- ifelse(zoop_drivers$year %in% c("2014","2019"), 
                           "g1", ifelse(zoop_drivers$year %in% c(
                             "2016", "2021"), "g2","g3")) 

#and also add nmds axis 1 values from second stage nmds
nmds1 <- read.csv("Output/ss_nmds1.csv")

zoop_drivers$nmds1 <- ifelse(zoop_drivers$year=="2014",nmds1$nmds1[1],
                      ifelse(zoop_drivers$year=="2015",nmds1$nmds1[2],
                      ifelse(zoop_drivers$year=="2016",nmds1$nmds1[3],
                      ifelse(zoop_drivers$year=="2019",nmds1$nmds1[4],
                      ifelse(zoop_drivers$year=="2020",nmds1$nmds1[5],
                                                       nmds1$nmds1[6])))))

#convert from wide to long
zoop_drivers_long <- zoop_drivers |> 
  select(-Shortwave) |> 
  pivot_longer(-c(month,year,nmds1,traj),
               names_to = "variable") |> 
  mutate(year = as.factor(year)) 

#change order of years
zoop_drivers_long$year <- factor(zoop_drivers_long$year , 
                                levels=c("2014", "2015" ,"2016", 
                                         "2019", "2020", "2021"))


#median vs nmds for all vars
zoop_drivers_long |> group_by(year, nmds1, variable, traj) |> 
  summarize(median = median(value, na.rm=T)) |> 
  ungroup() |> 
  ggplot(aes(x=nmds1, y=median, group = traj)) +
  geom_point(aes(color=year, alpha = 0.95), cex=3) + 
  facet_wrap(~variable, nrow=5, scales = "free_y") +
  theme_bw() + xlab("NMDS1") + 
  guides(alpha = "none", fill = "none",
         color = guide_legend(ncol=2)) +
  scale_color_manual("",values=year_cols)+
  theme(text = element_text(size=8), 
        axis.text = element_text(size=7, color="black"), 
        axis.text.x = element_text(angle=45, vjust=0.7, hjust=0.6),
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        plot.margin = unit(c(0.5,0.3,-0,0), 'lines'),
        legend.position = c(0.85, 0.05),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-15,-0,-0,-0),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#ggsave("Figures/all_drivers_vs_nmds1.jpg", width=7, height=5)

#drop rows with NA
zoop_drivers_long <- zoop_drivers_long |> filter(!is.na(value))

#correlation for all n=27 vars (dropping sw obvi)
var <- unique(zoop_drivers_long$variable)

#initialize correlation matrix
cor_df <- data.frame("cor" = rep(NA,27),
                     "p-val" = rep(NA,27))

for(i in 1:length(var)){
  cor_df[i,1] <- rcorr(zoop_drivers_long$nmds1[zoop_drivers_long$month==5 & 
                            zoop_drivers_long$variable==var[i]],
                zoop_drivers_long$value[zoop_drivers_long$month==5 & 
                            zoop_drivers_long$variable==var[i]],
                type = "spearman")$r[2]
  
  cor_df[i,2] <- rcorr(zoop_drivers_long$nmds1[zoop_drivers_long$month==5 & 
                                                 zoop_drivers_long$variable==var[i]],
                       zoop_drivers_long$value[zoop_drivers_long$month==5 & 
                                                 zoop_drivers_long$variable==var[i]],
                       type = "spearman")$P[2]
}
 
#add row names
rownames(cor_df) <- var
#p-vals < 0.002 are significant (0.05/27); which means that NO correlations are sig
 
#hand-picking the vars that are most different
vars <- c("air temp", "Schmidt stability",
          "Secchi depth",
          "total phytoplankton biomass",
          "epilimnetic DO", "epilimnetic TN",
          "wind speed", "residence time")

#order variables
zoop_drivers_long$variable <- factor(zoop_drivers_long$variable, 
                                      levels=c(vars, unique(zoop_drivers_long$variable[
                                        !zoop_drivers_long$variable %in% vars]))) 

#median vs trajectory
zoop_drivers_long |> group_by(year, nmds1, variable, traj) |> 
  summarize(median = median(value, na.rm=T)) |> 
  ungroup() |> 
  filter(variable %in% vars) |> 
  ggplot(aes(x=nmds1, y=median, group = traj)) +
  geom_point(aes(color=year, alpha = 0.95), cex=3) + 
  facet_wrap(~variable, nrow=5, scales = "free_y") +
  theme_bw() + xlab("NMDS1") + 
  guides(alpha = "none", fill = "none", color = guide_legend(nrow=1)) +
  scale_color_manual("",values=year_cols)+
  theme(text = element_text(size=12), 
        axis.text = element_text(size=11, color="black"), 
        axis.text.x = element_text(angle=45, vjust=0.7, hjust=0.6),
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        plot.margin = unit(c(0.5,0.5,0,0), 'lines'),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(0,-0,-0,-0),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#ggsave("Figures/8drivers_vs_clusters_vertical_points.jpg", width=7, height=5)

#month on x and colored points for years
zoop_drivers_long |> group_by(year, month, variable) |> 
  summarize(median = median(value)) |> 
  ungroup() |> 
  filter(variable %in% vars) |> 
  ggplot(aes(x=as.factor(month), y=median, group = year)) +
  geom_point(aes(color=year, alpha = 0.95), cex=3) + 
  geom_line(aes(color=year, alpha = 0.95)) +
  facet_wrap(~variable, nrow=5, scales = "free_y") +
  theme_bw() + xlab("") + guides(alpha = "none", fill = "none",
                                 colour = guide_legend(nrow = 1)) +
  scale_color_manual("",values= year_cols)+
  scale_x_discrete(labels= c("May","June","July",
                             "August","September")) +
  theme(text = element_text(size=15), 
        axis.text = element_text(size=12, color="black"), 
        axis.text.x = element_text(angle=45, vjust=0.7, hjust=0.6),
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        plot.margin = unit(c(0.5,0.3,1,0), 'lines'),
        legend.position = "bottom",
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-15,-0,-0,-0),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#ggsave("Figures/8drivers_vs_months.jpg", width=7, height=5)

mean(zoop_drivers_long$value[zoop_drivers_long$traj=="g1" &
                                 zoop_drivers_long$variable=="epilimnetic TP"])

mean(zoop_drivers_long$value[zoop_drivers_long$traj=="g2" &
                               zoop_drivers_long$variable=="epilimnetic TP"])


#month on x and colored points for years - all drivers
zoop_drivers_long |> group_by(month, year, variable) |> 
  summarize(median = median(value, na.rm=T)) |> 
  ungroup() |> 
  filter(!variable %in% "Shortwave") |> 
  ggplot(aes(x=as.factor(month), y=median, group = year)) +
  geom_point(aes(color=year, alpha = 0.95), cex=3) + 
  geom_line(aes(color=year, alpha = 0.95)) +
  facet_wrap(~variable, ncol=5, scales = "free_y") +
  theme_bw() + xlab("") + guides(alpha = "none", fill = "none",
                                 colour = guide_legend(nrow = 1)) +
  #scale_color_manual("",values= c("#003366","#0099CC","#339999",
  #                                "#660000","#CC0000","#CC6666"))+
  scale_color_manual("",values= year_cols)+
  scale_x_discrete(labels= c("May","June","July",
                             "August","September")) +
  theme(text = element_text(size=10), 
        axis.text = element_text(size=9, color="black"), 
        axis.text.x = element_text(angle=45, vjust=0.7, hjust=0.6),
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        plot.margin = unit(c(0.5,0.3,1,0), 'lines'),
        legend.position = "bottom",
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-15,-0,-0,-0),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#ggsave("Figures/all_drivers_vs_months.jpg", width=7, height=5)

#---------------------------------------------------------------------#
#plot median zoop dens + env drivers for copepod/nauplii and keratella/kellicottia

#read in zoop csv
zoops <- read.csv("Output/std_dens_10taxa.csv")

#plot represnetative drivers associated with zoop density changes
do <- zoop_drivers_long |> 
  group_by(year, month, variable) |> 
  summarize(median = median(value)) |> 
  ungroup() |> 
  filter(variable %in% c("epilimnetic DO")) |> 
  mutate(g = ifelse(year %in% c("2014","2015","2021"),"1","2")) |> 
  group_by(g, month, variable) |> 
  summarize(median = median(median)) |> 
  ggplot(aes(x=as.factor(month), y=median, group = g)) + 
  geom_line(aes(y=median),linetype=1) +
  theme_bw() + xlab("") +
  scale_x_discrete(expand = c(0,0), breaks = c("5","6","7","8","9"),
                   labels = c("May","June","Jul","Aug","Sep")) +
  scale_y_continuous(expand = c(0,0), position = "right")+ 
  facet_wrap(~g,
             labeller = label_wrap_gen(multi_line=FALSE))+ 
  scale_linetype_manual("",values= c(1))+
  coord_cartesian(ylim = c(4, 9), clip = "off") +
  ylab("") +
  #ylab("DO (mg/L)") +
  theme(text = element_text(size=10), 
        strip.text.x = element_blank(),
        axis.text = element_text(size=9, color="black"), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.right = element_text(vjust = +3),
        strip.background = element_rect(fill = "transparent"), 
        plot.margin = unit(c(1.1,1.9,-0.5,1), 'lines'),
        legend.position = "none",
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(0,-0,-0,-0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.3, "in"))

phy <- zoop_drivers_long |> 
  group_by(year, month, variable) |> 
  summarize(median = median(value,na.rm=T)) |> 
  ungroup() |> 
  filter(variable %in% c("total phytoplankton biomass")) |> 
  mutate(g = ifelse(year %in% c("2014","2015","2021"),"1","2")) |> 
  group_by(g, month, variable) |> 
  summarize(median = median(median, na.rm=T)) |> 
  ggplot(aes(x=as.factor(month), y=median, group = g)) + 
  geom_line(aes(y=median),linetype=2) +
  theme_bw() + xlab("") +
  scale_x_discrete(expand = c(0,0), breaks = c("5","6","7","8","9"),
                   labels = c("May","June","Jul","Aug","Sep")) +
  scale_y_continuous(expand = c(0,0), position = "right")+ 
  facet_wrap(~g,
             labeller = label_wrap_gen(multi_line=FALSE))+ 
  scale_linetype_manual("",values= c(1))+
  coord_cartesian(ylim = c(2, 18), clip = "off") +
  ylab("median value") +
  #ylab(bquote(phytoplankton~biomass~(mu*g/L))) +
  theme(text = element_text(size=10), 
        strip.text.x = element_blank(),
        axis.text = element_text(size=9, color="black"), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.right = element_text(vjust = +4),
        strip.background = element_rect(fill = "transparent"), 
        plot.margin = unit(c(0,1.6,0.7,1), 'lines'),
        legend.position = "none",
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(0,-0,-0,-0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.3, "in"))

bf <- zoop_drivers_long |> 
  group_by(year, month, variable) |> 
  summarize(median = median(value)) |> 
  ungroup() |> 
  filter(variable %in% c("buoyancy frequency")) |> 
  mutate(g = ifelse(year %in% c("2014","2015","2021"),"1","2")) |> 
  group_by(g, month, variable) |> 
  summarize(median = median(median)) |> 
  ggplot(aes(x=as.factor(month), y=median, group = g)) + 
  geom_line(aes(y=median), linetype=3) +
  theme_bw() + xlab("") +
  scale_x_discrete(expand = c(0,0), breaks = c("5","6","7","8","9"),
                   labels = c("May","June","Jul","Aug","Sep")) +
  scale_y_continuous(expand = c(0,0), position = "right")+ 
  facet_wrap(~g,
             labeller = label_wrap_gen(multi_line=FALSE)) + 
  geom_line(aes(x=as.factor(month),y=50, 
                linetype= "epiimnetic DO")) +
  geom_line(aes(x=as.factor(month),y=50, 
                linetype= "total phytoplankton biomass")) +
  geom_line(aes(x=as.factor(month),y=50, 
                linetype= "buoyancy frequency")) +
  scale_linetype_manual("", values = c("epiimnetic DO" = 1,
                                "total phytoplankton biomass" = 2,
                                "buoyancy frequency" = 3)) +
  coord_cartesian(ylim = c(0, 0.015), clip = "off") +
  ylab("") +
  #ylab(bquote(buoyancy~frequency~(1/sec^2))) +
  theme(text = element_text(size=10), 
        strip.text.x = element_blank(),
        axis.text = element_text(size=9, color="black"), 
        axis.ticks.x = element_line(colour = c(rep("black",5), "transparent")), 
        axis.title.y = element_text(vjust = +3),
        strip.background = element_rect(fill = "transparent"), 
        plot.margin = unit(c(-1.2,0.8,1,1), 'lines'),
        legend.box.spacing = unit(0.1, "cm"),
        legend.position = "bottom",
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(0,-0,-0,-0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.3, "in"))

ggpubr::ggarrange(do, phy, bf, nrow=3)
#ggsave("Figures/env_vars_conceptual_monthly.jpg", width=7, height=5)
  
#----------------------------------------------------------------#
#plot representative zoop densities
daph <- zoops |> select(year, month, Taxon, avg) |> 
  filter(Taxon %in% c("Daphnia"), #two representative taxa
         month %in% c(5:9),
         year %in% c("2014","2019")) |> #two representative years
  mutate(g = ifelse(year %in% c("2014"), "1","2")) |> 
  group_by(Taxon, g, month) |> 
  summarize(median = median(avg)) |> 
  ggplot(aes(x=as.factor(month), y=median, group = g)) + 
  geom_area(aes(color = Taxon, fill = Taxon),
            position = "identity", stat = "identity",
            alpha=0.7) + theme_bw() + xlab("") +
  scale_x_discrete(expand = c(0,0), breaks = c("5","6","7","8","9"),
                   labels = c("May","June","Jul","Aug","Sep")) +
  scale_y_continuous(expand = c(0,0))+ 
  facet_wrap(~g,
             labeller = label_wrap_gen(multi_line=FALSE))+ 
  scale_color_manual("",values= c("#644554","#E99848","#AB4C1E"))+
  scale_fill_manual("",values= c("#644554","#E99848","#AB4C1E"),
                    labels = c("Daphnia",
                               "Conochilus",
                               "Keratella"))+
  geom_text(x = 3, y = 250, 
            aes(label = c("2014",NA,NA,NA,NA,
                          "2019",NA,NA,NA,NA))) +
  coord_cartesian(ylim = c(0, 230), clip = "off") +
  guides(color = "none") +
  ylab("") +
  theme(text = element_text(size=10), 
        strip.text.x = element_blank(),
        axis.text = element_text(size=9, color="black"), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(vjust = +3),
        strip.background = element_rect(fill = "transparent"), 
        plot.margin = unit(c(1,1,-0.5,0.6), 'lines'),
        legend.position = "none",
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(0,-0,-0,-0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.3, "in"))

con <- zoops |> select(year, month, Taxon, avg) |> 
  filter(Taxon %in% c("Conochilus"), #two representative taxa
         month %in% c(5:9),
         year %in% c("2014","2019")) |> #two representative years
  mutate(g = ifelse(year %in% c("2014"), "1","2")) |> 
  group_by(Taxon, g, month) |> 
  summarize(median = median(avg)) |> 
  ggplot(aes(x=as.factor(month), y=median, group = g)) + 
  geom_area(aes(color = Taxon, fill = Taxon),
            position = "identity", stat = "identity",
            alpha=0.7) + theme_bw() + xlab("") +
  scale_x_discrete(expand = c(0,0), breaks = c("5","6","7","8","9"),
                   labels = c("May","June","Jul","Aug","Sep")) +
  scale_y_continuous(expand = c(0,0))+ 
  facet_wrap(~g,
             labeller = label_wrap_gen(multi_line=FALSE))+ 
  scale_color_manual("",values= c("#E99848","#AB4C1E", "#644554"))+
  scale_fill_manual("",values= c("#E99848","#AB4C1E", "#644554"),
                    labels = c("Conochilus",
                               "Keratella",
                               "Daphnia"))+
  coord_cartesian(ylim = c(0, 1200), clip = "off") +
  guides(color = "none") +
  ylab("zooplankton density (#/L)") +
  theme(text = element_text(size=10), 
        strip.text.x = element_blank(),
        axis.text = element_text(size=9, color="black"), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(vjust = +3),
        strip.background = element_rect(fill = "transparent"), 
        plot.margin = unit(c(0,1,0.5,0.3), 'lines'),
        legend.position = "none",
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(0,-0,-0,-0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.3, "in"))


ker <- zoops |> select(year, month, Taxon, avg) |> 
  filter(Taxon %in% c("Keratella"), #representative taxa
         month %in% c(5:9),
         year %in% c("2014","2019")) |> #two representative years
  mutate(g = ifelse(year %in% c("2014"), "1","2")) |> 
  group_by(Taxon, g, month) |> 
  summarize(median = median(avg)) |> 
  ggplot(aes(x=as.factor(month), y=median, group = g)) + 
  geom_area(aes(color = Taxon, fill = Taxon),
            position = "identity", stat = "identity",
            alpha=0.7) + theme_bw() + xlab("") +
  scale_x_discrete(expand = c(0,0), breaks = c("5","6","7","8","9"),
                   labels = c("May","June","Jul","Aug","Sep")) +
  scale_y_continuous(expand = c(0,0))+ 
  coord_cartesian(ylim = c(0, 1200), clip = "off") +
  facet_wrap(~g,
             labeller = label_wrap_gen(multi_line=FALSE))+ 
  scale_color_manual("",values= c("#AB4C1E","#AB4C1E", "#AB4C1E"),
                     labels = "Keratella")+
  scale_fill_manual("",values= c("#644554","#E99848", "#AB4C1E"),
                    labels = c("Daphnia",
                               "Conochilus",
                               "Keratella"))+
  geom_point(aes(x=as.factor(month),y=100, fill = "#644554"), alpha = 0) +
  geom_point(aes(x=as.factor(month),y=100, fill = "#E99848"), alpha = 0) +
  guides(color = "none") +
  ylab("") +
  theme(text = element_text(size=10), 
        strip.text.x = element_blank(),
        axis.text = element_text(size=9, color="black"), 
        axis.ticks.x = element_line(colour = c(rep("black",5), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        plot.margin = unit(c(-1,1,0.7,0.3), 'lines'),
        legend.position = "bottom",
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.box.spacing = unit(0.1, "cm"),
        legend.margin=margin(-1,-0,-0,-0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.3, "in"))

ggpubr::ggarrange(daph, con, ker, nrow=3)
#ggsave("Figures/zoop_dens_conceptual_monthly.jpg", width=7, height=5)


