#script for second-stage NMDS

pacman::p_load(zoo,dplR,dplyr,tidyverse,ggplot2,ggpubr,
               sf,vegan,FSA,rcompanion,NatParksPalettes,ggrepel)

#read in all_zoops df
all_zoops_dens <- read.csv("Output/all_zoops_dens.csv",header = TRUE)

#cb friendly year palette
#year_cols <- c("#f93414","#fe6c31","#06889b","#006b8f","#011f51","#feae6c")

year_cols <- c("#a13637","#06889b", "#facd60", "#f44034", "#011f51", "#fdfa66")
# 2014, 2019, 2016, 2021, 2015, 2020 
#------------------------------------------------------------------------------#
#NMDS 
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

#initialize diversity df
div <- data.frame(year = c(rep(2014,5),rep(2015,5),rep(2016,5),
                           rep(2019,5),rep(2020,5),rep(2021,5)),
                  month = rep(c(5,6,7,8,9),6))

#calculate diversity
simp <- diversity(zoops_dens, index="simpson")
shan <- diversity(zoops_dens, index="shannon")

div$simpsonD <- simp
div$shannonD <- shan

ggplot(div, aes(as.Date("2019-12-31") + 
                 yday(as.Date(paste0(year,"-0",month,"-01"))), 
                                     shannonD, color=as.factor(year))) +
  geom_point() + geom_line() + xlab("") +
  scale_color_manual(values = year_cols)+
  scale_x_date(labels = scales::date_format("%b",tz="EST5EDT"),
              breaks = "1 month") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        text = element_text(size=10), 
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.x = element_text(face = "bold",hjust = 0),
        axis.text.x = element_text(angle=90),
        strip.background.x = element_blank(),
        axis.title.y = element_text(size = 11),
        plot.margin = unit(c(0, 1, 0, 0), "cm"),
        panel.background = element_rect(
          fill = "white"),
        panel.spacing = unit(0.5, "lines"))
#ggsave("Figures/BVR_yearly_shannon.jpg", width=7, height=4) 

#hellinger transform data
zoop_dens_trans <- labdsv::hellinger(zoops_dens)

#turn transformed community data into b-c distance matrix 
zoop_bray <- as.matrix(vegan::vegdist(zoop_dens_trans, method='bray'))

#------------------------------------------------------------------------------#
#first-stage NMDS
#scree plot to choose dimension 
#jpeg("figures/scree.jpg") 
goeveg::dimcheckMDS(zoop_bray, distance = "bray", 
                    k = 6, trymax = 20, autotransform = TRUE)
#dev.off()

set.seed(11)
#now do NMDS w/ 3 dimensions 
NMDS_bray_first <- vegan::metaMDS(zoop_bray, distance='bray', k=3, trymax=20, 
                                  autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_bray_first$stress
# 0.10

#plot
ord <- vegan::ordiplot(NMDS_bray_first,display = c('sites','species'),
                       choices = c(1,2),type = "n")
year <- ggordiplots::gg_ordiplot(ord, all_zoops_nmds$year,
                                 kind = "ehull", ellipse=FALSE, hull = TRUE, 
                                 plot = FALSE, pt.size=0.9) 
year <- year$plot + geom_point() + theme_bw() + 
  geom_polygon(data = year$df_hull, aes(x = x, y = y, fill = Group), 
               alpha=0.2) +
  geom_point(data=year$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, 
             fill=year_cols) +
  theme(text = element_text(size=10), 
        axis.text = element_text(size=6, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "vertical",
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = c(0.88,0.18), 
        legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,0.1,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) +
  scale_fill_manual("",values=year_cols)+
  scale_color_manual("",values=year_cols,
                     label=c('2014','2015',"2016","2019","2020","2021")) 

month <- ggordiplots::gg_ordiplot(ord, all_zoops_nmds$month,
                                  kind = "ehull", ellipse=FALSE, hull = TRUE, 
                                  plot = FALSE, pt.size=0.9) 
month <- month$plot + geom_point() + theme_bw() + 
  geom_polygon(data = month$df_hull, aes(x = x, y = y, fill = Group), 
               alpha=0.2) +
  geom_point(data=month$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, 
             fill=viridis::viridis(5, option="F")) +
  theme(text = element_text(size=10), 
        axis.text = element_text(size=6, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "vertical",
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = c(0.83,0.16),
        legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,0.1,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + guides(fill='none') +
  scale_fill_manual("",values=viridis::viridis(5, option="F"))+
  scale_color_manual("",values=viridis::viridis(5, option="F"),
                     label=c('May',"June","July","August","September")) 

ggpubr::ggarrange(year,month,ncol=2, common.legend = F)
#ggsave("Figures/first_stage_NMDS_2v1_dens.jpg", width=5, height=3)

#------------------------------------------------------------------------------#
#track months in each year and connect
#jpeg("Figures/first_stage_NMDS_2v1_dens.jpg")
par(mfrow = c(2, 3))
par(cex = 0.6)
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))

par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
for (i in 1:6) {
  ord <- vegan::ordiplot(NMDS_bray_first,display = c('sites','species'),
                         choices = c(1,2),type = "n", xaxt='n', yaxt='n')
  
  lines(NMDS_bray_first$points[all_zoops_nmds$year == 
                                 unique(all_zoops_nmds$year)[i] ,1], 
        NMDS_bray_first$points[(all_zoops_nmds$year ==
                                  unique(all_zoops_nmds$year)[i]),2], 
        col=year_cols[i])
  
  points(NMDS_bray_first$points[all_zoops_nmds$year == 
                                  unique(all_zoops_nmds$year)[i] ,1], 
         NMDS_bray_first$points[(all_zoops_nmds$year == 
                                   unique(all_zoops_nmds$year)[i]),2], 
         col=year_cols[i],
         pch = as.character(5:9), font=2, cex=2)
  
  mtext(c("2014","2015","2016","2019","2020","2021")[i], side = 3, line = -2, 
        adj = 0.04, cex = 1, col = "black")
  #if (i %in% c(4, 5, 6))
  #  axis(1, col = "black", col.axis = "black")
  #if (i %in% c(1, 4))
  #  axis(2, col = "black", col.axis = "black")
  #box(col = "black")
  #if (i %in% c(4))
  #  legend("bottomleft", legend=c('May','June','July','August','September'),
  #         pch=c(as.character(1:5)) ,bty = "n", cex=1.5) 
}
mtext("NMDS1", side = 1, outer = TRUE, cex = 0.8, line = 1,
      col = "black")
mtext("NMDS2", side = 2, outer = TRUE, cex = 0.8, line = 1,
      col = "black")

#legend("topright", legend=c('2014','2015','2016','2019','2020','2021'),
#       pt.bg=c(viridis::viridis(5, option="D")) ,bty = "n", cex=1.2, pch=21) 

#dev.off()
#------------------------------------------------------------------------------#
#create distance matrices for all years (first stage pairwise dissimilarities)
dist2014 <- as.matrix(vegan::vegdist(zoop_dens_trans[all_zoops_nmds$year=="2014",],
                                     method='bray'))
dist2015 <- as.matrix(vegan::vegdist(zoop_dens_trans[all_zoops_nmds$year=="2015",],
                                     method='bray'))
dist2016 <- as.matrix(vegan::vegdist(zoop_dens_trans[all_zoops_nmds$year=="2016",],
                                     method='bray'))
dist2019 <- as.matrix(vegan::vegdist(zoop_dens_trans[all_zoops_nmds$year=="2019",],
                                     method='bray'))
dist2020 <- as.matrix(vegan::vegdist(zoop_dens_trans[all_zoops_nmds$year=="2020",],
                                     method='bray'))
dist2021 <- as.matrix(vegan::vegdist(zoop_dens_trans[all_zoops_nmds$year=="2021",],
                                     method='bray'))
#combine all distance matrices - each col is a year/month combo (5X30)
#alldist <- as.data.frame(cbind(dist2014, dist2015, dist2016,
#                              dist2019, dist2020, dist2021))
#perform pairwise correlations among dissimilarity matrices
#and convert from similarity to distance 
#stage2 <- as.matrix(vegan::vegdist(1-cor(alldist)), method="bray")

# Function to calculate pairwise correlation between two matrices
pairwise_correlation <- function(mat1, mat2) {
  vec1 <- as.vector(mat1)
  vec2 <- as.vector(mat2)
  return(cor(vec1, vec2)) 
}

# List of matrices
matrices <- list(dist2014, dist2015, dist2016, dist2019, dist2020, dist2021)

# Initialize an empty matrix to store correlations
n <- length(matrices)
correlation_matrix <- matrix(NA, n, n)

# Calculate pairwise correlations and fill the matrix
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    correlation_matrix[i, j] <- pairwise_correlation(matrices[[i]], matrices[[j]])
    correlation_matrix[j, i] <- correlation_matrix[i, j] # Since correlation is symmetric
  }
}

# Set diagonal to 1 (correlation of a matrix with itself)
diag(correlation_matrix) <- 1

#run NMDS again - clustering will indicate years where changes through time are similar
set.seed(11)
#now do NMDS w/ 3 dimensions 
NMDS_bray_second <- vegan::metaMDS(correlation_matrix, distance='bray', k=3, trymax=20, 
                                   autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_bray_second$stress
# 0.05

#--------------------------------------------------------------------------#
#NMDS plot - second-stage
#Note that warnings are okay here because there is only one point per year (can't determine different df_ellipse values for years)
ord <- vegan::ordiplot(NMDS_bray_second,display = c('sites','species'),
                       choices = c(1,2),type = "n")
year1 <- ggordiplots::gg_ordiplot(ord, unique(all_zoops_nmds$year), kind = "ehull", 
                                 spiders = FALSE, ellipse = FALSE,
                                 label = FALSE, hull = FALSE, 
                                 plot = FALSE, pt.size=0.9) 
plot1 <- year1$plot + geom_point() + theme_bw() + 
  geom_point(data=year1$df_mean.ord, aes(x, y), 
             pch=21, size=3, 
             fill=year_cols) +
#  ylim(-0.4,0.4) + xlim(-0.6,0.6) +
  theme(text = element_text(size=14), 
        axis.text = element_text(size=7, color="black"), 
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "none", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,0,0,0), 'lines'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(
    color = year_cols))) +
  scale_fill_manual("",values=year_cols,
                    label=c('2014','2015',"2016","2019","2020","2021")) +
  scale_color_manual("",values=rep("gray",6) )

ord <- vegan::ordiplot(NMDS_bray_second,display = c('sites','species'),
                       choices = c(1,3),type = "n")
year2 <- ggordiplots::gg_ordiplot(ord, unique(all_zoops_nmds$year), kind = "ehull", 
                                  spiders = FALSE, ellipse = FALSE,
                                  label = FALSE, hull = FALSE, 
                                  plot = FALSE, pt.size=0.9) 
plot2 <- year2$plot + geom_point() + theme_bw() + 
  geom_point(data=year2$df_mean.ord, aes(x, y), 
             pch=21, size=3, 
             fill=year_cols) +
  #ylim(-0.4,0.4) + xlim(-0.6,0.6) +
  theme(text = element_text(size=14), 
        axis.text = element_text(size=7, color="black"), 
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "right",
        plot.margin = unit(c(0,0,0,0), 'lines'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(
    color = year_cols))) +
  scale_fill_manual("",values = year_cols,
                    label=c('2014','2015',"2016","2019","2020","2021")) +
  scale_color_manual("",values=rep("gray",6) )


ggpubr::ggarrange(plot1,plot2,ncol=2, common.legend = T)
#ggsave("Figures/second_stage_NMDS_dens.jpg", width=6, height=3) 

#create new df for export
nmds1 <- data.frame("year" = c(2014,2015,2016,2019,2020,2021),
                    "nmds1" = ord$sites[,1])

#export ss nmds1 for driver analysis
write.csv(nmds1,"Output/ss_nmds1.csv", row.names = F)

#-----------------------------------------------------------#
#read in env csv
env_drivers <- read.csv("./Output/env.csv") |> 
  rename("epi TN" = "TN_ugL_epi",
         "hypo TN" = "TN_ugL_hypo",
         "epi TP" = "TP_ugL_epi",
         "hypo TP" = "TP_ugL_hypo",
         "epi temp" = "Temp_C_epi",
         "hypo temp" = "Temp_C_hypo",
         "epi DO" = "DO_mgL_epi",
         "water level" = "waterlevel",
         "wl cv" = "wl_cv",
         "thermocline depth" = "therm_depth",
         "oxycline depth" = "oxy_depth",
         "Schmidt stability" = "SS",
         "buoyancy frequency" = "BF",
         "air temp" = "AirTemp",
         "shortwave" = "Shortwave",
         "longwave" = "Longwave",
         "relative humidity" = "RelHum",
         "wind speed" = "WindSpeed",
         "rain" = "Rain")

#get zoops into same order as drivers
all_zoops_nmds <- all_zoops_nmds |> arrange(month, year)

#join driver and env dfs
zoops_plus_drivers <- bind_cols(all_zoops_nmds, env_drivers[
  !colnames(env_drivers) %in% c("month", "year")])

ord <- vegan::ordiplot(NMDS_bray_first,display = c('sites','species'),
                       choices = c(1,2),type = "n")
#fit environmental drivers onto ordination
fit_env <- envfit(ord$sites, zoops_plus_drivers[,c(13:23, 25, 29:31)])
# 1 vs. 2 - dropping vars with p > 0.9; SS, airtemp, sw, lw
#fit_env <- envfit(ord$sites, zoops_plus_drivers[,c(13:15, 17:23,25,26, 28:31)])
# 1 vs. 3 - dropping sw, hypo tp, ss

#pull out vectors - need to multiply by the sqrt of r2 to get magnitude!
scores <- data.frame((fit_env$vectors)$arrows * sqrt(fit_env$vectors$r), 
                     pvals=(fit_env$vectors)$pvals)
scores <- cbind(scores, env = rownames(scores))

#supplmental table w/ r2 and p values for ms
driver_NMDS_correlation <- data.frame("variable" = scores$env,
                                      "R2" = fit_env$vectors$r,
                                      "p-value" = fit_env$vectors$pvals)
#write.csv(driver_NMDS_correlation, "Output/driver_NMDS_correlation.csv", row.names=F)

#look for correlations between all pairs of env variables
#driver_correlation <- data.frame(cor(zoops_plus_drivers[,c(13:26, 28:31)], method = "spearman"))
#write.csv(driver_correlation, "Output/driver_correlation.csv", row.names=F)

#plot drivers w/ NMDS
ord <- vegan::ordiplot(NMDS_bray_first,display = c('sites'),
                       choices = c(1,2),type = "n")
month <- ggordiplots::gg_ordiplot(ord, all_zoops_nmds$month,
                                  kind = "ehull", ellipse=FALSE, hull = TRUE, 
                                  plot = FALSE, pt.size=0.9) 
month$plot + geom_point() + theme_bw() + 
  geom_polygon(data = month$df_hull, aes(x = x, y = y, fill = Group), 
               alpha=0.2) +
  geom_point(data=month$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, 
             fill=viridis::viridis(5, option="F")) +
  theme(text = element_text(size=10), 
        axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "vertical",
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = c(0.86,0.15), 
        legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,-0.1,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + guides(fill='none') +
  scale_fill_manual("",values=viridis::viridis(5, option="F"))+
  scale_color_manual("",values=viridis::viridis(5, option="F"),
                     label=c('May',"June","July","August","September")) +
  #xlim(-0.4,0.4) +
  geom_segment(data = scores,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), linewidth= 0.3,
               arrow = arrow(length = unit(0.1, "cm")), colour = "black") +
  geom_text_repel(data = scores, aes(x = NMDS1, y = NMDS2, label = env), 
                  size = 1.5, box.padding = 0.1)
#ggsave("Figures/first_stage_NMDS_2v1_months_envfit.jpg", width=6, height=4)


ord <- vegan::ordiplot(NMDS_bray_first,display = c('sites'),
                       choices = c(1,2),type = "n")
year <- ggordiplots::gg_ordiplot(ord, all_zoops_nmds$year,
                                 kind = "ehull", ellipse=FALSE, hull = TRUE, 
                                 plot = FALSE, pt.size=0.9) 
year$plot + geom_point() + theme_bw() + 
  geom_polygon(data = year$df_hull, aes(x = x, y = y, fill = Group), 
               alpha=0.2) +
  geom_point(data=year$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, 
             fill=year_cols) +
  theme(text = element_text(size=10), 
        axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.text = element_text(size=6),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "vertical",
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = c(0.11,0.14), legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,-0.1,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) +
  guides(fill='none') +
  scale_fill_manual("",values=year_cols)+
  xlim(-0.41,0.41) +
  scale_color_manual("",values=year_cols,
                     label=c('2014','2015',"2016","2019","2020","2021")) +
  geom_segment(data = scores,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), linewidth= 0.3,
               arrow = arrow(length = unit(0.1, "cm")), colour = "black") +
  geom_text_repel(data = scores, aes(x = NMDS1, y = NMDS2, label = env), 
                  size = 1.5, box.padding = 0.1)
#ggsave("Figures/first_stage_NMDS_2v1_years_envfit.jpg", width=4, height=3)


#-----------------------------------------------------------------#
#calculate dispersion between years vs. months - monte carlo approach

#order zoop epi tows by year and month
all_zoops_nmds <- all_zoops_nmds |> dplyr::arrange(year, month)

#convert matrix back into distance structure for next steps
zoop_bray <- vegdist(zoop_dens_trans, method='bray', 
                     upper = TRUE)

#Now calculate the centroids of each polygon AND the avg distance of each point to its polygon centroid
centroids_year <- betadisper(zoop_bray, group = as.factor(all_zoops_nmds$year), 
                             type="centroid")
centroids_month <- betadisper(zoop_bray, group = as.factor(all_zoops_nmds$month), 
                              type="centroid")

#calculate dispersion as average distance of each point to polygon centroid
#randomly select 5 points per group and 5 groups

#create a df to store variability values
var_results <- data.frame("year_disp"=rep(NA,500))

set.seed(1)

for (i in 1:500){ 
  #randomly select 5 points in each group
  year_sub <- sample(unique(all_zoops_nmds$year), 5)
  month_sub <-  sample(unique(all_zoops_nmds$month), 5)
  
  zoop_sub <-  all_zoops_nmds |> group_by(year) |> 
    filter(year %in% c(year_sub), month %in% c(month_sub)) |>
    slice_sample(n=5)
  
  #only select data cols
  zoop_dens_sub <- zoop_sub |> ungroup() |> select(-c(year,month))
  
  #hellinger transformation
  zoop_dens_sub_trans <- labdsv::hellinger(zoop_dens_sub)
  
  #convert matrix back into distance structure for next steps
  zoop_bray_sub <- vegdist(zoop_dens_sub_trans, method='bray', 
                           upper = TRUE)
  
  #Now calculate the centroids of each polygon AND the avg distance of each point to its polygon centroid
  centroids_year_sub <- betadisper(zoop_bray_sub, group = as.factor(zoop_sub$year), 
                                   type="centroid")
  centroids_month_sub <-  betadisper(zoop_bray_sub, group = as.factor(zoop_sub$month), 
                                     type="centroid")
  
  #distance calcs
  var_results$year_disp[i] <- mean(centroids_year_sub$group.distances)
  var_results$month_disp[i] <- mean(centroids_month_sub$group.distances)
  
}

#Kruskal-wallis test to determine if group means are significant

#first convert wide to long
disp_df <- var_results[,grepl("disp",colnames(var_results))] |>  
  pivot_longer(everything(), names_to="group")

#now kw test
kw_disp <- kruskal.test(value ~ group, data = disp_df) #sig
#dispersion among years is greater than months

#plot
ggboxplot(disp_df, x = "group", y = "value", 
          fill = "group", palette = c("#A4C6B8", "#5E435D"),
          order = c("year_disp", "month_disp"),
          ylab = "Dispersion") +
  theme(text = element_text(size=9),
        plot.margin = unit(c(0.2,0,-0.5,0), 'lines')) +
  annotate("text",label=c("a","b"), x=c(1.1,2.1),
           y=c(mean(disp_df$value[disp_df$group=="year_disp"]) + 
                 sd(disp_df$value[disp_df$group=="year_disp"])/1.3,
               mean(disp_df$value[disp_df$group=="month_disp"]) + 
                 sd(disp_df$value[disp_df$group=="month_disp"])),
           size = 6) +
  guides(fill = "none") +
  scale_x_discrete(name ="", 
                   labels=c("year_disp"="year",
                            "month_disp"="month"))
#ggsave("Figures/among_variability_boxplots_dens.jpg", width=4, height=4)

#create table for kw test results
kw_results <- data.frame("Scale" = c("Year", "Month"),
                         "n" = c(rep(500,2)),
                         "mean" = c(round(mean(var_results$year_disp),2), 
                                    round(mean(var_results$month_disp),2)),
                         "sd" = c(round(sd(var_results$year_disp),2), 
                                  round(sd(var_results$month_disp),2)),
                         "df" = c(kw_disp$parameter, " "),
                         "χ2" = c(round(kw_disp$statistic,3), " "),
                         "p-value" = c(kw_disp$p.value, " "))
#write.csv(kw_results, "Output/Bray_distances_bootstrapped_kw_results_dens.csv", row.names = FALSE)

#dfs to calculate significance within years, and months
within_year_dist <- data.frame("group" = c(rep("2014",5),rep("2015",5),rep("2016",5),
                                           rep("2019",5),rep("2020",5),rep("2021",5)),
                               "dist" = c(centroids_year$distances[all_zoops_nmds$year=="2014"],
                                          centroids_year$distances[all_zoops_nmds$year=="2015"],
                                          centroids_year$distances[all_zoops_nmds$year=="2016"],
                                          centroids_year$distances[all_zoops_nmds$year=="2019"],
                                          centroids_year$distances[all_zoops_nmds$year=="2020"],
                                          centroids_year$distances[all_zoops_nmds$year=="2021"]))

within_month_dist <- data.frame("group" = c(rep("May",6),rep("June",6),rep("July",6),
                                            rep("August",6),rep("September",6)),
                                "dist" = c(centroids_month$distances[all_zoops_nmds$month=="05"],
                                           centroids_month$distances[all_zoops_nmds$month=="06"],
                                           centroids_month$distances[all_zoops_nmds$month=="07"],
                                           centroids_month$distances[all_zoops_nmds$month=="08"],
                                           centroids_month$distances[all_zoops_nmds$month=="09"]))

#now kw tests for significance 
kw_year <- kruskal.test(dist ~ group, data = within_year_dist) #not significant
kw_month <- kruskal.test(dist ~ group, data = within_month_dist) #sig

#dunn post-hoc test
dunn_within_month <- dunnTest(dist ~ as.factor(group),
                              data=within_month_dist,
                              method="bonferroni")

cldList(P.adj ~ Comparison, data=dunn_within_month$res, threshold = 0.05)

#within boxplot for years and months
year_box <- ggboxplot(within_year_dist, x = "group", y = "dist", 
                      fill = "group", 
                      palette = year_cols,
                      order = c("2014", "2015","2016","2019","2020","2021"),
                      ylab = "Dispersion", xlab = "") +
  ylim(c(0,0.5)) +
  theme(text = element_text(size=10),
        plot.margin = unit(c(0,-0.5,0,0), 'lines'),
        axis.text.x = element_text(angle=45, vjust=0.8, hjust=0.8)) +
  annotate("text", x=1.5, y=0.5, label= "a: years",
           fontface = "italic", size=5) +
  guides (fill = "none")

month_box <- ggboxplot(within_month_dist, x = "group", y = "dist", 
                       fill = "group", 
                       palette = viridis::viridis(5, option="F"),
                       order = c("May", "June", "July", "August", "September"),
                       ylab = "", xlab = "") + ylim(c(0,0.5)) +
  theme(text = element_text(size=10),
        plot.margin = unit(c(0,0.2,0,-0.4), 'lines'),
        axis.text.x = element_text(angle=45, vjust=0.8, hjust=0.8),
        axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  annotate("text",label=c("a","b","b","b","ab"), x=c(1.2, 2.2, 3.2, 4.2, 5.2), size=4,
           y=c(0.33, 0.24, 0.23, 0.25, 0.25)) +
  annotate("text", x=1.5, y=0.5, label= "b: months",
           fontface = "italic", size=5) +
  guides (fill = "none")

within_scales <- egg::ggarrange(year_box, month_box, nrow=1, widths = c(2, 1.9))
#ggsave("Figures/within_variability_boxplots_dens.jpg", within_scales, width=5, height=4)

#create table for within scale kw test results
kw_results_disp <- data.frame("Group" = c("2014", "2015", "2016", "2019", 
                                          "2020", "2021", "May", "June",
                                          "July", "August", "September"),
                              "n" = c(6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5),
                              "mean" = c(round(mean(within_year_dist$dist[within_year_dist$group=="2014"]),2), 
                                         round(mean(within_year_dist$dist[within_year_dist$group=="2015"]),2),
                                         round(mean(within_year_dist$dist[within_year_dist$group=="2016"]),2),
                                         round(mean(within_year_dist$dist[within_year_dist$group=="2019"]),2),
                                         round(mean(within_year_dist$dist[within_year_dist$group=="2020"]),2),
                                         round(mean(within_year_dist$dist[within_year_dist$group=="2021"]),2),
                                         round(mean(within_month_dist$dist[within_month_dist$group=="May"]),2),
                                         round(mean(within_month_dist$dist[within_month_dist$group=="June"]),2),
                                         round(mean(within_month_dist$dist[within_month_dist$group=="July"]),2),
                                         round(mean(within_month_dist$dist[within_month_dist$group=="August"]),2),
                                         round(mean(within_month_dist$dist[within_month_dist$group=="September"]),2)),
                              "sd" = c(round(sd(within_year_dist$dist[within_year_dist$group=="2014"]),2), 
                                       round(sd(within_year_dist$dist[within_year_dist$group=="2015"]),2),
                                       round(sd(within_year_dist$dist[within_year_dist$group=="2016"]),2),
                                       round(sd(within_year_dist$dist[within_year_dist$group=="2019"]),2),
                                       round(sd(within_year_dist$dist[within_year_dist$group=="2020"]),2),
                                       round(sd(within_year_dist$dist[within_year_dist$group=="2021"]),2),
                                       round(sd(within_month_dist$dist[within_month_dist$group=="May"]),2),
                                       round(sd(within_month_dist$dist[within_month_dist$group=="June"]),2),
                                       round(sd(within_month_dist$dist[within_month_dist$group=="July"]),2),
                                       round(sd(within_month_dist$dist[within_month_dist$group=="August"]),2),
                                       round(sd(within_month_dist$dist[within_month_dist$group=="September"]),2)),
                              "df" = c(rep(" ",2), kw_year$parameter, rep(" ",5), kw_month$parameter, rep(" ",2)),
                              "χ2" = c(rep(" ",2), kw_year$statistic, rep(" ",5), kw_month$statistic, rep(" ",2)),
                              "p-value" = c(rep(" ",2), kw_year$p.value, rep(" ",5), kw_month$p.value, rep(" ",2)))
#write.csv(kw_results_disp, "Output/within_group_dispersion_kw_results_dens.csv",row.names = FALSE)
