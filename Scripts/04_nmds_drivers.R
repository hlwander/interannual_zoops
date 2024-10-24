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
         "residence time" = res_time,
         "water level" = waterlevel,
         "thermocline depth" = therm_depth,
         "oxycline depth" = oxy_depth,
         "Schmidt stability" = SS,
         "epilimnetic DO" = DO_mgL_epi,
         "epilimnetic TN" = TN_ugL_epi,
         "epilimnetic TP" = TP_ugL_epi,
         "hypolimnetic TN" = TN_ugL_hypo,
         "hypolimnetic TP" = TP_ugL_hypo,
         "Secchi depth" = secchi,
         "hypolimnetic temp" = Temp_C_hypo,
         "epilimnetic temp" = Temp_C_epi,
         "phytoplankton biomass" = Total_ugL)

#and also add nmds axis 1 values from second stage nmds
nmds1 <- read.csv("Output/ss_nmds1.csv")

zoop_drivers$nmds1 <- ifelse(zoop_drivers$year=="2014",nmds1$nmds1[1],
                      ifelse(zoop_drivers$year=="2015",nmds1$nmds1[2],
                      ifelse(zoop_drivers$year=="2016",nmds1$nmds1[3],
                      ifelse(zoop_drivers$year=="2019",nmds1$nmds1[4],
                      ifelse(zoop_drivers$year=="2020",nmds1$nmds1[5],
                                                       nmds1$nmds1[6])))))

#read in proportion metrics
prop <- read.csv("Output/zoop_proportions.csv") 

#convert from wide to long
zoop_drivers_long <- zoop_drivers |> 
  select(-Shortwave) |> 
  mutate(prop_cyclopoid = ifelse(year%in% "2014",prop$p_cyc[1],
                                ifelse(year %in% "2015", prop$p_cyc[2],
                                       ifelse(year %in% "2016", prop$p_cyc[3],
                                              ifelse(year %in% "2019", prop$p_cyc[4],
                                                     ifelse(year %in% "2020", prop$p_cyc[5], prop$p_cyc[6])))))) |> 
  pivot_longer(-c(month,year,nmds1,prop_cyclopoid),
               names_to = "variable") |> 
  mutate(year = as.factor(year)) 

#change order of years
zoop_drivers_long$year <- factor(zoop_drivers_long$year , 
                                levels=c("2014", "2015" ,"2016", 
                                         "2019", "2020", "2021"))

#some weird glitch where you need to unload then reload dplyr before getting next chunk of code to run
detach("package:dplyr", unload = TRUE)
library(dplyr)

#median vs nmds for all vars (Supplemental Figure S5)
zoop_drivers_long |> group_by(year, prop_cyclopoid, variable) |> 
  summarize(median = median(value, na.rm=T)) |> 
  ungroup() |> 
  mutate(variable = factor(variable, labels = c(
    bquote('air temp ('*degree*'C)'), bquote('bluegreen ('*mu*g~L^-1~')'),
    bquote('brown ('*mu*g~L^-1~')'), "epilimnetic~DO~(mg~L^{-1})",
    bquote('epilimnetic temp ('*degree*'C)'), bquote('epilimnetic TN ('*mu*g~L^-1~')'),
    bquote('epilimnetic TP ('*mu*g~L^-1~')'), bquote('green ('*mu*g~L^-1~')'),
    bquote('hypolimnetic temp ('*degree*'C)'), bquote('hypolimnetic TN ('*mu*g~L^-1~')'),
    bquote('hypolimnetic TP ('*mu*g~L^-1~')'), "longwave~(W~m^{-2})",
    bquote('mixed ('*mu*g~L^-1~')'), "oxycline~depth~(m)",
    bquote('phytoplankton biomass ('*mu*g~L^-1~')'), "precipitation~(m~d^{-1})",
    "relative~humidity~('%')", "residence~time (d)",
    "Schmidt~stability~(J~m^{-2})", "Secchi~depth~(m)",
    "thermocline~depth~(m)", "water~level~(m)",
    "wind~speed~(m~s^{-1})", "water~level~cv"))) |> 
  ggplot(aes(x=median, y=prop_cyclopoid*100, group = year)) +
  geom_point(aes(color=year, alpha = 0.95), cex=3) + 
  facet_wrap(~variable, nrow=5, scales = "free_x",
             labeller = label_parsed) +
  theme_bw() + 
  ylab("Percent Cyclopoid (%)") + 
  guides(alpha = "none", fill = "none",
         color = guide_legend(ncol=2)) +
  scale_color_manual("",values=year_cols)+
  theme(text = element_text(size=7), 
        axis.text = element_text(size=7, color="black"), 
        axis.text.x = element_text(angle=45, vjust=0.7, hjust=0.6),
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        plot.margin = unit(c(0.5,0.3,-0,0), 'lines'),
        legend.position = c(0.9, 0.06),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-15,-0,-0,-0),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#ggsave("Figures/p_cyc_vs_all_drivers.jpg", width=7, height=5)

#hand-picking the vars that are most different
vars <- c("air temp", 
           "Rain",
          "Schmidt stability",
          "Secchi depth", 
          "epilimnetic TN",
          "epilimnetic TP")

#order variables
zoop_drivers_long$variable <- factor(zoop_drivers_long$variable, 
                                      levels=c(vars, unique(zoop_drivers_long$variable[
                                        !zoop_drivers_long$variable %in% vars]))) 

#median vs trajectory (Manuscript Figure 5)
zoop_drivers_long |> group_by(year, prop_cyclopoid, variable) |> 
  summarize(median = median(value, na.rm=T)) |> 
  ungroup() |> 
  filter(variable %in% vars) |> 
  mutate(variable = factor(variable,
                          labels = c(bquote('air temp ('*degree*'C)'),
                                     "precipitation~(m~d^{-1})",
                                     "Schmidt~stability~(J~m^{-2})",
                                     "Secchi~depth~(m)",
                                     bquote('epilimnetic TN ('*mu*g~L^-1~')'),
                                     bquote('epilimnetic TP ('*mu*g~L^-1~')')))) |> 
  ggplot(aes(x=median, y=prop_cyclopoid*100, group = year)) +
  geom_point(aes(color=year, alpha = 0.95), cex=3) + 
  facet_wrap(~variable, nrow=3, scales = "free_x", 
             labeller = label_parsed) +
  theme_bw() + ylab("Percent Cyclopoid (%)") + 
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
#ggsave("Figures/p_cyc_vs_6drivers.jpg", width=7, height=5)

#ms values for results text
res <- zoop_drivers_long |> group_by(year, variable) |> 
          summarize(median = median(value, na.rm=T)) |> 
          filter(variable %in% vars)

(mean(res$median[res$variable %in% 'air temp' & 
             res$year %in% c("2014","2019")]) -
mean(res$median[res$variable %in% 'air temp' & 
                  res$year %in% c("2015","2016","2020","2021")])) /
  mean(res$median[res$variable %in% 'air temp' & 
                    res$year %in% c("2015","2016","2020","2021")]) * 100

(mean(res$median[res$variable %in% 'wind speed' & 
                   res$year %in% c("2014","2019")]) -
    mean(res$median[res$variable %in% 'wind speed' & 
                      res$year %in% c("2015","2016","2020","2021")])) /
  mean(res$median[res$variable %in% 'wind speed' & 
                    res$year %in% c("2015","2016","2020","2021")]) * 100

(mean(res$median[res$variable %in% 'Schmidt stability' & 
                   res$year %in% c("2014","2019")]) -
    mean(res$median[res$variable %in% 'Schmidt stability' & 
                      res$year %in% c("2015","2016","2020","2021")])) /
  mean(res$median[res$variable %in% 'Schmidt stability' & 
                    res$year %in% c("2015","2016","2020","2021")]) * 100

(mean(res$median[res$variable %in% 'Secchi depth' & 
                   res$year %in% c("2014","2019")]) -
    mean(res$median[res$variable %in% 'Secchi depth' & 
                      res$year %in% c("2015","2016","2020","2021")])) /
  mean(res$median[res$variable %in% 'Secchi depth' & 
                    res$year %in% c("2015","2016","2020","2021")]) * 100

(mean(res$median[res$variable %in% 'epilimnetic TN' & 
                   res$year %in% c("2014","2019")]) -
    mean(res$median[res$variable %in% 'epilimnetic TN' & 
                      res$year %in% c("2015","2016","2020","2021")])) /
  mean(res$median[res$variable %in% 'epilimnetic TN' & 
                    res$year %in% c("2015","2016","2020","2021")]) * 100

(mean(res$median[res$variable %in% 'epilimnetic TP' & 
                   res$year %in% c("2014","2019")]) -
    mean(res$median[res$variable %in% 'epilimnetic TP' & 
                      res$year %in% c("2015","2016","2020","2021")])) /
  mean(res$median[res$variable %in% 'epilimnetic TP' & 
                    res$year %in% c("2015","2016","2020","2021")]) * 100

#month on x and colored points for years - all drivers (Supplemental Figure S7)
zoop_drivers_long |> group_by(month, year, variable) |> 
  summarize(median = median(value, na.rm=T)) |> 
  ungroup() |> 
  mutate(variable = factor(variable, labels = c(
    bquote('air temp ('*degree*'C)'), "precipitation~(m~d^{-1})", 
    "Schmidt~stability~(J~m^{-2})", "Secchi~depth~(m)",
    bquote('epilimnetic TN ('*mu*g~L^-1~')'), bquote('epilimnetic TP ('*mu*g~L^-1~')'),
    bquote('hypolimnetic TN ('*mu*g~L^-1~')'), bquote('hypolimnetic TP ('*mu*g~L^-1~')'),
    bquote('epilimnetic temp ('*degree*'C)'), bquote('hypolimnetic temp ('*degree*'C)'),
    "epilimnetic~DO~(mg~L^{-1})", "water~level~(m)", "water~level~cv",
    "thermocline~depth~(m)", "oxycline~depth~(m)",
    bquote('green ('*mu*g~L^-1~')'), bquote('bluegreen ('*mu*g~L^-1~')'),
    bquote('brown ('*mu*g~L^-1~')'), bquote('mixed ('*mu*g~L^-1~')'),
    bquote('phytoplankton biomass ('*mu*g~L^-1~')'),
    "longwave~(W~m^{-2})", "relative~humidity~('%')", 
    "wind~speed~(m~s^{-1})","residence~time (d)"))) |> 
  ggplot(aes(x=as.factor(month), y=median, group = year)) +
  geom_point(aes(color=year, alpha = 0.95), cex=3) + 
  geom_line(aes(color=year, alpha = 0.95)) +
  facet_wrap(~variable, ncol=5, scales = "free_y",
             labeller = label_parsed) +
  theme_bw() + xlab("") + guides(alpha = "none", fill = "none",
                                 colour = guide_legend(nrow = 1)) +
  scale_color_manual("",values= year_cols)+
  scale_x_discrete(labels= c("May","June","July",
                             "August","September")) +
  theme(text = element_text(size=6), 
        axis.text = element_text(size=7, color="black"), 
        axis.text.x = element_text(angle=45, vjust=0.7, hjust=0.6),
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        plot.margin = unit(c(0.5,0.3,1,0), 'lines'),
        legend.position = "bottom",
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-15,-0,-0,-0),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#ggsave("Figures/all_drivers_vs_months.jpg", width=7, height=5)