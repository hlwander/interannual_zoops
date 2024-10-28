#script for second-stage NMDS

pacman::p_load(zoo,dplR,dplyr,tidyverse,ggplot2,ggpubr,
               vegan,FSA,rcompanion,NatParksPalettes,ggrepel,labdsv,
               goeveg, ggordiplots, cowplot)

#read in all_zoops df
all_zoops_dens <- read.csv("Output/all_zoops_dens.csv",header = TRUE)

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

#hellinger transform data
zoop_dens_trans <- labdsv::hellinger(zoops_dens)

#turn transformed community data into b-c distance matrix 
zoop_bray <- as.matrix(vegan::vegdist(zoop_dens_trans, method='bray'))

#------------------------------------------------------------------------------#
#first-stage NMDS
#scree plot to choose dimension 
#jpeg("figures/scree.jpg") (Supplemental Figure S1)
goeveg::dimcheckMDS(zoop_bray, distance = "bray", 
                    k = 6, trymax = 20, autotransform = TRUE)
#dev.off()

set.seed(11)
#now do NMDS w/ 3 dimensions 
NMDS_bray_first <- vegan::metaMDS(zoop_bray, distance='bray', k=3, trymax=20, 
                                  autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_bray_first$stress
# 0.10

#plot axis 2 vs. 1 (Supplemental Figure S3)
ord <- vegan::ordiplot(NMDS_bray_first,display = c('sites'),
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

#plot axis 3 vs. 1 (Supplemental Figure S4)
ord <- vegan::ordiplot(NMDS_bray_first,display = c('sites'),
                       choices = c(1,3),type = "n")
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
#ggsave("Figures/first_stage_NMDS_3v1_dens.jpg", width=5, height=3)

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
#NMDS plot - second-stage (Manuscript Figure 3)
#Note that warnings are okay here because there is only one point per year (can't determine different df_ellipse values for years)
ord <- vegan::ordiplot(NMDS_bray_second,display = c('sites'),
                       choices = c(1,2),type = "n")
year1 <- ggordiplots::gg_ordiplot(ord, unique(all_zoops_nmds$year), kind = "ehull", 
                                 spiders = FALSE, ellipse = FALSE,
                                 label = FALSE, hull = FALSE, 
                                 plot = FALSE, pt.size=NA) 
plot1 <- year1$plot + geom_point() + theme_bw() + xlim(c(-0.45,0.45)) +
  geom_text(data = year1$df_mean.ord, 
            aes(x, y, label = as.factor(Group),
                color = as.factor(Group)), 
            size = 4) + ylim(c(-0.5,0.5)) +
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
  scale_color_manual("",values=year_cols,
                    label=c('2014','2015',"2016","2019","2020","2021")) 

ord <- vegan::ordiplot(NMDS_bray_second,display = c('sites','species'),
                       choices = c(1,3),type = "n")
year2 <- ggordiplots::gg_ordiplot(ord, unique(all_zoops_nmds$year), kind = "ehull", 
                                  spiders = FALSE, ellipse = FALSE,
                                  label = FALSE, hull = FALSE, 
                                  plot = FALSE, pt.size=NA) 
plot2 <- year2$plot + theme_bw() + xlim(c(-0.45,0.45)) +
  geom_text(data = year2$df_mean.ord, 
            aes(x, y, label = as.factor(Group),
                     color = as.factor(Group)), 
            size = 4) + ylim(c(-0.5,0.5)) +
  theme(text = element_text(size=14), 
        axis.text = element_text(size=7, color="black"), 
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0.1), 'lines'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual("",values = year_cols,
                    label=c('2014','2015',"2016","2019","2020","2021")) 


ggpubr::ggarrange(plot1,plot2,ncol=2, common.legend = F)
#ggsave("Figures/second_stage_NMDS_dens.jpg", width=6, height=3) 

#create new df for export
nmds1 <- data.frame("year" = c(2014,2015,2016,2019,2020,2021),
                    "nmds1" = ord$sites[,1])

#export ss nmds1 for driver analysis
#write.csv(nmds1,"Output/ss_nmds1.csv", row.names = F)

#-----------------------------------------------------------#
#read in env csv
env_drivers <- read.csv("./Output/all_drivers.csv") |> 
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
         "Secchi depth" = "secchi",
         "air temp" = "AirTemp",
         "shortwave" = "Shortwave",
         "longwave" = "Longwave",
         "relative humidity" = "RelHum",
         "wind speed" = "WindSpeed",
         "rain" = "Rain",
         "residence time" = "res_time",
         "phytoplankton biomass" = "Total_ugL") |> 
  select(-c(Mixed_ugL, Bluegreen_ugL, Brown_ugL, Green_ugL, shortwave))

#get zoops into same order as drivers
all_zoops_nmds <- all_zoops_nmds |> arrange(month, year)

#join driver and env dfs
zoops_plus_drivers <- bind_cols(all_zoops_nmds, env_drivers[
  !colnames(env_drivers) %in% c("month", "year")])

#---------------------------------------------------------------#
# summarize by year
zoops_plus_drivers_yearly <- zoops_plus_drivers |> 
  group_by(year) |> 
  summarise_at(vars(-month), funs(mean(., na.rm=TRUE)))

ord <- vegan::ordiplot(NMDS_bray_second,display = c('sites'),
                       choices = c(1,2),type = "n")
#fit environmental drivers onto ordination
fit_env <- envfit(ord['sites'], zoops_plus_drivers_yearly[,c(12:31)])

#pull out vectors - need to multiply by the sqrt of r2 to get magnitude!
scores <- data.frame((fit_env$vectors)$arrows * sqrt(fit_env$vectors$r), 
                     pvals=(fit_env$vectors)$pvals)
scores <- cbind(scores, env = rownames(scores))
#---------------------------------------------------------------#
#plot drivers w/ second stage NMDS (Manuscript Figure 4)
ss_year <- ggordiplots::gg_ordiplot(ord, unique(all_zoops_nmds$year),
                                  kind = "ehull", ellipse=FALSE, hull = TRUE, 
                                  plot = FALSE, pt.size=0.9) 
env_plot1 <- ss_year$plot + geom_point() + theme_bw() + 
  geom_polygon(data = ss_year$df_hull, aes(x = x, y = y, fill = Group), 
               alpha=0.2) +
  geom_point(data=ss_year$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, 
             fill=year_cols) +
  theme(text = element_text(size=10), 
        axis.text = element_text(size=7, color="black"), 
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "none",
        plot.margin = unit(c(0,-1,0,-2), 'lines'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(nrow=1, override.aes = list(
    color = year_cols)), fill = "none") +
  scale_fill_manual("",values=year_cols)+
  scale_color_manual("",values=year_cols,
                     label=c("2014","2015","2016","2019","2020","2021")) +
  xlim(-0.7,0.9) + ylim(-1,0.9) +
  geom_segment(data = scores,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), linewidth= 0.3,
               arrow = arrow(length = unit(0.1, "cm")), colour = "black") +
  geom_segment(data = scores[scores$env %in% c("epi TP"),],
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), linewidth= 0.3,
               arrow = arrow(length = unit(0.1, "cm")), colour = "lightgray") +
  geom_text_repel(data = scores, aes(x = NMDS1, y = NMDS2, label = env), 
                  size = 1.5, box.padding = 0.2, max.overlaps=Inf)


#axis 1 vs. 3
ord <- vegan::ordiplot(NMDS_bray_second,display = c('sites'),
                       choices = c(1,3),type = "n")
#fit environmental drivers onto ordination
fit_env <- envfit(ord['sites'], zoops_plus_drivers_yearly[,c(12:31)])

#pull out vectors - need to multiply by the sqrt of r2 to get magnitude!
scores <- data.frame((fit_env$vectors)$arrows * sqrt(fit_env$vectors$r), 
                     pvals=(fit_env$vectors)$pvals)
scores <- cbind(scores, env = rownames(scores))


ss_year <- ggordiplots::gg_ordiplot(ord, unique(all_zoops_nmds$year),
                                    kind = "ehull", ellipse=FALSE, hull = TRUE, 
                                    plot = FALSE, pt.size=0.9) 
env_plot2 <- ss_year$plot + geom_point() + theme_bw() + 
  geom_polygon(data = ss_year$df_hull, aes(x = x, y = y, fill = Group), 
               alpha=0.2) +
  geom_point(data=ss_year$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, 
             fill=year_cols) +
  theme(text = element_text(size=10), 
        axis.text = element_text(size=7, color="black"), 
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "right",
        plot.margin = unit(c(0,0,0,-2), 'lines'),
        legend.margin = margin(-10,-10,-10,-10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(
    color = year_cols)), fill = "none") +
  scale_fill_manual("",values=year_cols)+
  scale_color_manual("",values=year_cols,
                     label=c("2014","2015","2016","2019","2020","2021")) +
  xlim(-0.7,0.9) + ylim(-1,0.9) +
  geom_segment(data = scores,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS3), linewidth= 0.3,
               arrow = arrow(length = unit(0.1, "cm")), colour = "black") +
  geom_segment(data = scores[scores$env %in% c("Schmidt stability", 
                                                "hypo temp", "Secchi depth"),],
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS3), linewidth= 0.3,
               arrow = arrow(length = unit(0.1, "cm")), colour = "lightgray") +
  geom_text_repel(data = scores, aes(x = NMDS1, y = NMDS3, label = env), 
                  size = 1.5, box.padding = 0.2, max.overlaps=Inf)

ggpubr::ggarrange(env_plot1,env_plot2,ncol=2, common.legend = F)
#ggsave("Figures/second_stage_NMDS_envfit.jpg", width=6, height=3)

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

#plot (Manuscript Figure 1 A)
YvsM <-ggboxplot(disp_df, x = "group", y = "value", 
          fill = "group", palette = c("#A4C6B8", "#5E435D"),
          order = c("year_disp", "month_disp"),
          ylab = "Dispersion") +
  theme(text = element_text(size=10),
        plot.margin = unit(c(0.2,0,-0.2,0.4), 'lines')) +
  annotate("text",label=c("a","b"), x=c(1.1,2.1),
           y=c(mean(disp_df$value[disp_df$group=="year_disp"]) + 
                 sd(disp_df$value[disp_df$group=="year_disp"]) - 0.0005 ,
               mean(disp_df$value[disp_df$group=="month_disp"]) + 
                 sd(disp_df$value[disp_df$group=="month_disp"]) + 0.001),
           size = 6) +
  guides(fill = "none") +
  scale_x_discrete(name ="", 
                   labels=c("year_disp"="Year",
                            "month_disp"="Month"))

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

#within boxplot for years and months (Manuscript Figure 1 B+C)
year_box <- ggboxplot(within_year_dist, x = "group", y = "dist", 
                      fill = "group", 
                      palette = year_cols,
                      order = c("2014", "2015","2016","2019","2020","2021"),
                      ylab = "Dispersion", xlab = "") +
  ylim(c(0,0.5)) +
  theme(text = element_text(size=10),
        plot.margin = unit(c(-0.3,-0.2,0.4,0.8), 'lines'),
        axis.text.x = element_text(angle=45, vjust=0.8, hjust=0.8)) +
  annotate("text",label=c("ns"), x=c(3.5), size=4, 
           y=c(0.04), fontface = "italic") +
  guides (fill = "none")

month_box <- ggboxplot(within_month_dist, x = "group", y = "dist", 
                       fill = "group", 
                       palette = viridis::viridis(5, option="F"),
                       order = c("May", "June", "July", "August", "September"),
                       ylab = "", xlab = "") + ylim(c(0,0.5)) +
  theme(text = element_text(size=10),
        plot.margin = unit(c(0,0.4,-0.9,0.2), 'lines'),
        axis.text.x = element_text(angle=45, vjust=0.8, hjust=0.8),
        axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  annotate("text",label=c("a","b","b","b","ab"), x=c(1.2, 2.2, 3.2, 4.2, 5.2), size=4,
           y=c(0.35, 0.26, 0.25, 0.27, 0.27)) +
  guides (fill = "none")

# Combine year_box and month_box in a single row
second_row <- plot_grid(year_box, month_box, labels = c("B", "C"), ncol = 2)

# Combine YvsM and the second row vertically
final_plot <- plot_grid(YvsM, second_row, labels = c("A", ""), ncol = 1)

# Print the final plot
print(final_plot)
#ggsave("Figures/variability_boxplots_dens.jpg", width=5, height=4)

#create table for within scale kw test results (Supplemental Table S2)
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
