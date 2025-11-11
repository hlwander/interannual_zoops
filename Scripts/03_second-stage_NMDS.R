#script for second-stage NMDS

pacman::p_load(zoo,dplR,dplyr,tidyverse,ggplot2,ggpubr,
               vegan,FSA,rcompanion,NatParksPalettes,ggrepel,labdsv,
               goeveg, ggordiplots, cowplot)

#read in all_zoops df
all_zoops_dens <- read.csv("Output/all_zoops_dens.csv",header = TRUE)

year_cols <- c("#011f51","#06889b","#2E8B57","#fdfa66","#facd60","#f44034","#a13637")
# 2014, 2015, 2016, 2019, 2020, 2021, 2023

#---------------------------------------------------------------------------------#
#NMDS: 89 total samples; 16 removed that fall outside of may-oct range for SS NMDS
taxa <- unique(all_zoops_dens$Taxon)

#taxa as cols, dates as rows, average by month
all_zoops_nmds <- all_zoops_dens |> 
  dplyr::select(DateTime, Taxon, dens) |> 
  filter(Taxon %in% taxa) |> 
  mutate(DateTime = as.Date(DateTime)) |> 
  pivot_wider(names_from = Taxon, values_from = dens) |> 
  mutate_all(~replace(., is.na(.), 0)) |>  #replace NA with 0
  ungroup() |>
  group_by(DateTime) |> 
  summarise(Bosmina = mean(Bosmina),
            Daphnia = mean(Daphnia),
            Cyclopoida = mean(Cyclopoida),
            Nauplii = mean(Nauplii),
            Conochilus = mean(Conochilus),
            Keratella = mean(Keratella),
            Kellicottia = mean(Kellicottia),
            Polyarthra = mean(Polyarthra)) |> 
  ungroup() |>
  mutate(year = year(DateTime),
         month = month(DateTime)) |>
  filter(!month %in% c(3,12)) #drop the months that could not be imputed
#write.csv(all_zoops_nmds, "./Output/zoop_raw_dens.csv", row.names=FALSE)

#select only data cols
zoops_dens_all <- all_zoops_nmds |> select(Bosmina:Polyarthra)

#hellinger transform data
zoop_dens_trans_all <- labdsv::hellinger(zoops_dens_all)
#write.csv(zoop_dens_trans_all, "./Output/zoop_dens_trans.csv", row.names=FALSE)

#turn transformed community data into b-c distance matrix 
zoop_bray_first <- as.matrix(vegan::vegdist(zoop_dens_trans_all, method='bray'))

#------------------------------------------------------------------------------#
#first-stage NMDS
#scree plot to choose dimension 
#jpeg("figures/scree.jpg") (Supplemental Figure S1)
goeveg::dimcheckMDS(zoop_bray_first, distance = "bray", 
                    k = 6, trymax = 20, autotransform = TRUE)
#dev.off()

set.seed(11)
#now do NMDS w/ 3 dimensions 
NMDS_bray_first <- vegan::metaMDS(zoop_bray_first, distance='bray', k=3, trymax=20, 
                                  autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_bray_first$stress
# 0.13

#plot axis 2 vs. 1 (Supplemental Figure S3)
ord <- vegan::ordiplot(NMDS_bray_first,display = c('sites'),
                       choices = c(1,2),type = "n")
year12 <- ggordiplots::gg_ordiplot(ord, all_zoops_nmds$year,
                                 kind = "ehull", ellipse=FALSE, hull = TRUE, 
                                 plot = FALSE, pt.size=0.9) 
year12 <- year12$plot + geom_point() + theme_bw() + 
  geom_polygon(data = year12$df_hull, aes(x = x, y = y, fill = Group), 
               alpha=0.2) +
  geom_point(data=year12$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, 
             fill=year_cols) +
  theme(text = element_text(size=9), 
        axis.text = element_text(size=6, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "horizontal",
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "top", 
        legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,0.1,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) +
  guides(fill = guide_legend(nrow = 2)) +
  scale_fill_manual("",values=year_cols)+
  scale_color_manual("",values=year_cols,
                     label=c('2014','2015',"2016","2019","2020","2021","2023")) 

month12 <- ggordiplots::gg_ordiplot(ord, all_zoops_nmds$month,
                                  kind = "ehull", ellipse=FALSE, hull = TRUE, 
                                  plot = FALSE, pt.size=0.9) 
month12 <- month12$plot + geom_point() + theme_bw() + 
  geom_polygon(data = month12$df_hull, aes(x = x, y = y, fill = Group), 
               alpha=0.2) +
  geom_point(data=month12$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, 
             fill=viridis::viridis(8, option="F")) +
  theme(text = element_text(size=9), 
        axis.text = element_text(size=6, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "horizontal",
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "top",
        legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,0.1,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + guides(fill='none') +
  guides(color = guide_legend(nrow = 2)) +
  scale_fill_manual("",values=viridis::viridis(8, option="F"))+
  scale_color_manual("",values=viridis::viridis(8, option="F"),
                     label=c("April",'May',"June","July","August",
                             "September","October","November")) 

ggpubr::ggarrange(year12,month12,ncol=2, common.legend = F)
#ggsave("Figures/first_stage_NMDS_2v1_dens_all.jpg", width=5, height=3)

#plot axis 3 vs. 1 (Supplemental Figure S4)
ord <- vegan::ordiplot(NMDS_bray_first,display = c('sites'),
                       choices = c(1,3),type = "n")
year13 <- ggordiplots::gg_ordiplot(ord, all_zoops_nmds$year,
                                 kind = "ehull", ellipse=FALSE, hull = TRUE, 
                                 plot = FALSE, pt.size=0.9) 
year13 <- year13$plot + geom_point() + theme_bw() + 
  geom_polygon(data = year13$df_hull, aes(x = x, y = y, fill = Group), 
               alpha=0.2) +
  geom_point(data=year13$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, 
             fill=year_cols) +
  theme(text = element_text(size=9), 
        axis.text = element_text(size=6, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "horizontal",
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "top", 
        legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,0.1,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) +
  guides(color = guide_legend(ncol = 1)) +
  scale_fill_manual("",values=year_cols)+
  scale_color_manual("",values=year_cols,
                     label=c('2014','2015',"2016","2019","2020","2021","2023")) 

month13 <- ggordiplots::gg_ordiplot(ord, all_zoops_nmds$month,
                                  kind = "ehull", ellipse=FALSE, hull = TRUE, 
                                  plot = FALSE, pt.size=0.9) 
month13 <- month13$plot + geom_point() + theme_bw() + 
  geom_polygon(data = month13$df_hull, aes(x = x, y = y, fill = Group), 
               alpha=0.2) +
  geom_point(data=month13$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, 
             fill=viridis::viridis(8, option="F")) +
  theme(text = element_text(size=9), 
        axis.text = element_text(size=6, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "horizontal",
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "top",
        legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,0.1,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + guides(fill='none') +
  scale_fill_manual("",values=viridis::viridis(8, option="F"))+
  scale_color_manual("",values=viridis::viridis(8, option="F"),
                     label=c("April",'May',"June","July","August",
                             "September","October","November")) 

ggpubr::ggarrange(year13,month13,ncol=2, common.legend = F)
#ggsave("Figures/first_stage_NMDS_3v1_dens_all.jpg", width=5, height=3)

#------------------------------------------------------------------------------#
# prepare dataset for SS NMDS

#only keep may-sep samples and drop 2022
monthly_zoops_nmds <- all_zoops_nmds |> 
  mutate(year = format(DateTime, "%Y"),
         month = format(DateTime, "%m")) |>
  group_by(month, year) |>
  summarise(Bosmina = mean(Bosmina),
            Daphnia = mean(Daphnia),
            Cyclopoida = mean(Cyclopoida),
            Nauplii = mean(Nauplii),
            Conochilus = mean(Conochilus),
            Keratella = mean(Keratella),
            Kellicottia = mean(Kellicottia),
            Polyarthra = mean(Polyarthra)) |> 
  ungroup() |>
  filter(month %in% c("05","06","07","08","09","10"), 
         !year %in% c("2022"))

#select only data cols
zoops_dens <- monthly_zoops_nmds |> select(Bosmina:Polyarthra)

#hellinger transform data
zoop_dens_trans <- labdsv::hellinger(zoops_dens)

#turn transformed community data into b-c distance matrix 
zoop_bray_second <- as.matrix(vegan::vegdist(zoop_dens_trans, method='bray'))

#------------------------------------------------------------------------------#
#create distance matrices for all years (first stage pairwise dissimilarities)
dist2014 <- as.matrix(vegan::vegdist(zoop_dens_trans[monthly_zoops_nmds$year=="2014",],
                                     method='bray'))
dist2015 <- as.matrix(vegan::vegdist(zoop_dens_trans[monthly_zoops_nmds$year=="2015",],
                                     method='bray'))
dist2016 <- as.matrix(vegan::vegdist(zoop_dens_trans[monthly_zoops_nmds$year=="2016",],
                                     method='bray'))
dist2019 <- as.matrix(vegan::vegdist(zoop_dens_trans[monthly_zoops_nmds$year=="2019",],
                                     method='bray'))
dist2020 <- as.matrix(vegan::vegdist(zoop_dens_trans[monthly_zoops_nmds$year=="2020",],
                                     method='bray'))
dist2021 <- as.matrix(vegan::vegdist(zoop_dens_trans[monthly_zoops_nmds$year=="2021",],
                                     method='bray'))
dist2023 <- as.matrix(vegan::vegdist(zoop_dens_trans[monthly_zoops_nmds$year=="2023",],
                                     method='bray'))

# Function to calculate pairwise correlation between two matrices
pairwise_correlation <- function(mat1, mat2) {
  vec1 <- as.vector(mat1)
  vec2 <- as.vector(mat2)
  return(cor(vec1, vec2)) 
}

# List of matrices
matrices <- list(dist2014, dist2015, dist2016, dist2019, 
                 dist2020, dist2021, dist2023)

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
# 0.08

#--------------------------------------------------------------------------#
#NMDS plot - second-stage (Manuscript Figure 3)
#Note that warnings are okay here because there is only one point per year (can't determine different df_ellipse values for years)
ord <- vegan::ordiplot(NMDS_bray_second,display = c('sites'),
                       choices = c(1,2),type = "n")
year1_ss <- ggordiplots::gg_ordiplot(ord, unique(monthly_zoops_nmds$year), kind = "sd", 
                                 spiders = FALSE, ellipse = FALSE,
                                 label = FALSE, hull = FALSE, 
                                 plot = FALSE, pt.size=NA) 
plot1 <- year1_ss$plot + geom_point() + theme_bw() +
  geom_text(data = year1_ss$df_mean.ord, 
            aes(x, y, label = as.factor(Group),
                color = as.factor(Group)), 
            size = 4) + ylim(c(-0.5,0.5)) +
  theme(text = element_text(size=14), 
        axis.text = element_text(size=7, color="black"), 
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",5), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "none", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,0,0,0), 'lines'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(
    color = year_cols))) +
  scale_x_continuous(limits = c(-0.55, 0.55),
                     breaks = c(-0.5, -0.25, 0, 0.25, 0.5),
                     labels = c("-0.5", "-0.25", "0", "0.25", "0.5")) +
  scale_color_manual("",values=year_cols,
                    label=c('2014','2015',"2016","2019","2020","2021","2023")) 

ord <- vegan::ordiplot(NMDS_bray_second,display = c('sites'),
                       choices = c(1,3),type = "n")
year2_ss <- ggordiplots::gg_ordiplot(ord, unique(monthly_zoops_nmds$year), kind = "sd", 
                                  spiders = FALSE, ellipse = FALSE,
                                  label = FALSE, hull = FALSE, 
                                  plot = FALSE, pt.size=NA) 
plot2 <- year2_ss$plot + theme_bw() + 
  geom_text(data = year2_ss$df_mean.ord, 
            aes(x, y, label = as.factor(Group),
                     color = as.factor(Group)), size = 4) +
  ylim(c(-0.5,0.5)) +
  theme(text = element_text(size=14), 
        axis.text = element_text(size=7, color="black"), 
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",5), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0.1), 'lines'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(-0.55, 0.55),
                     breaks = c(-0.5, -0.25, 0, 0.25, 0.5),
                     labels = c("-0.5", "-0.25", "0", "0.25", "0.5")) +
  scale_color_manual("",values = year_cols,
                    label=c('2014','2015',"2016","2019","2020","2021","2023")) 

ggpubr::ggarrange(plot1,plot2,ncol=2, common.legend = F)
#ggsave("Figures/second_stage_NMDS_dens.jpg", width=6, height=3) 

#create new df for export
#nmds1 <- data.frame("year" = c(2014,2015,2016,2019,2020,2021,2023),
#                    "nmds1" = ord$sites[,1])

#export ss nmds1 for driver analysis
#write.csv(nmds1,"Output/ss_nmds1.csv", row.names = F)

# extract 3D site scores
site_scores <- as.data.frame(vegan::scores(NMDS_bray_second, display = "sites"))

#assign years
site_scores$year <- c(unique(monthly_zoops_nmds$year))

# compute Euclidean distance matrix across the first 3 NMDS dims
dist_mat <- as.matrix(dist(site_scores[, c(1,2,3)], method = "euclidean"))

# find the most similar pair (smallest non-zero distance)
# set diagonal to NA so it won't be picked
diag(dist_mat) <- NA
min_val <- min(dist_mat, na.rm = TRUE)
which_min <- which(dist_mat == min_val, arr.ind = TRUE)
most_similar_pair <- data.frame(
  year1 = site_scores$year[which_min[1, "row"]],
  year2 = site_scores$year[which_min[1, "col"]],
  distance = min_val)

# after computing dist_mat
rownames(dist_mat) <- site_scores$year
colnames(dist_mat) <- site_scores$year

#wide to long
dist_long <- as.data.frame(dist_mat) |>
  mutate(year1 = site_scores$year) |>
  pivot_longer(-year1, names_to = "year2", values_to = "distance") |>
  mutate(year2 = as.integer(year2)) |>
  filter(!is.na(distance)) |>
  arrange(distance) |>
  filter(year1 < year2)

ggplot(dist_long, aes(x = factor(year1), y = factor(year2), fill = distance)) +
  geom_tile() +
  geom_text(aes(label = round(distance, 2)), size = 3) +
  scale_fill_gradient(low = "steelblue", high = "white") +
  labs(x = "", y = "", fill = "Euclidean\ndistance",
       title = "Pairwise distances between years (k = 3)") +
  theme_minimal() +
  coord_fixed()
#ggsave("Figures/second_stage_NMDS_distance_heatmap.jpg", width=6, height=3) 

# Use the distance matrix (as.dist)
hc <- hclust(as.dist(dist(site_scores[,1:3])), method = "average")

#jpeg("Figures/hierarchical_yar_clusters_ss.jpg", width = 6, height = 4, units = "in", res = 300)
plot(hc, labels = site_scores$year, xlab = NA, ylab = "Mean euclidean distance", 
     main = "Hierarchical clustering of years (k=3 NMDS distances)")
#dev.off()

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
         "thermocline depth" = "therm_depth",
         "oxycline depth" = "oxy_depth",
         "Schmidt stability" = "SS",
         "Secchi depth" = "secchi",
         "air temp" = "AirTemp",
         "shortwave" = "Shortwave",
         "longwave" = "Longwave",
         "relative humidity" = "RelHum",
         "wind speed" = "WindSpeed",
         "rain" = "Rain") |> 
  select(-c(Total_ugL, shortwave)) |>
  mutate(month = month(DateTime),
         year = year(DateTime)) 

#get zoops into same order as drivers
all_zoops_nmds <- all_zoops_nmds |> arrange(month, year)

#join zoop and env dfs
zoops_plus_drivers <- bind_cols(all_zoops_nmds, env_drivers[
  !colnames(env_drivers) %in% c("month", "year","DateTime")]) 

#---------------------------------------------------------------#
# plot env vectors on first stage nmds
env_drivers_only <- env_drivers |> dplyr::select(-c(month,year,DateTime))

# run envfit (use permutations if you want p-values)
set.seed(3)
ord <- vegan::ordiplot(NMDS_bray_first,display = c('sites'),
                       choices = c(1,2),type = "n")
ef <- envfit(ord, env_drivers_only, permutations = 999, na.rm = TRUE)

#pull out vectors; multiply by sqrt(r2) to get correct magnitudes
scores <- data.frame(ef$vectors$arrows * sqrt(ef$vectors$r),
                     pvals = ef$vectors$pvals,
                     env = rownames(ef$vectors$arrows))

#scale arrows by significance
min_mult <- 0.3  # shorter arrows for non-significant variables
scores$sig_mult <- ifelse(scores$pvals <= 0.05, 1, min_mult)
scores[, c("NMDS1", "NMDS2")] <- scores[, c("NMDS1", "NMDS2")] * scores$sig_mult * 0.5 
#halving the vectors to make them fit on the plot better

#extract site scores
sites_scores <- as.data.frame(vegan::scores(NMDS_bray_first, 
                                            display = "sites"))[, c(1,2)]
colnames(sites_scores) <- c("NMDS1", "NMDS2")

#plot
year_with_env_12 <- year12 +
  geom_segment(data = filter(scores, pvals <= 0.05),
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.02, "npc")), 
               color = "black") +
  geom_text_repel(data = filter(scores, pvals <= 0.05),
                  aes(x = NMDS1, y = NMDS2, label = env), 
                  size = 2) 

month_with_env_12 <- month12 +
  geom_segment(data = filter(scores, pvals <= 0.05),
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.02, "npc")), 
               color = "black") +
  geom_text_repel(data = filter(scores, pvals <= 0.05),
                  aes(x = NMDS1, y = NMDS2, label = env), 
                  size = 2) 

ggpubr::ggarrange(year_with_env_12, month_with_env_12, ncol=2, common.legend = F)
#ggsave("Figures/NMDS_2v1_all_dens_env.jpg", width=6, height=3) 

#----------------------------------------------------------------------------#
# same for axis 3 vs 1

# run envfit 
set.seed(3)
ord <- vegan::ordiplot(NMDS_bray_first,display = c('sites'),
                       choices = c(1,3),type = "n")
ef <- envfit(ord, env_drivers_only, permutations = 999, na.rm = TRUE)

#pull out vectors; multiply by sqrt(r2) to get correct magnitudes
scores <- data.frame(ef$vectors$arrows * sqrt(ef$vectors$r),
                     pvals = ef$vectors$pvals,
                     env = rownames(ef$vectors$arrows))

#scale arrows by significance
min_mult <- 0.3  # shorter arrows for non-significant variables
scores$sig_mult <- ifelse(scores$pvals <= 0.05, 1, min_mult)
scores[, c("NMDS1", "NMDS3")] <- scores[, c("NMDS1", "NMDS3")] * scores$sig_mult * 0.5
#halving the vectors to make them fit on the plot better

#extract site scores
sites_scores <- as.data.frame(vegan::scores(NMDS_bray_first, 
                                            display = "sites"))[, c(1,3)]
colnames(sites_scores) <- c("NMDS1", "NMDS3")

#plot
year_with_env_13 <- year13 +
  geom_segment(data = filter(scores, pvals <= 0.05),
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS3),
               arrow = arrow(length = unit(0.02, "npc")), 
               color = "black") +
  geom_text_repel(data = filter(scores, pvals <= 0.05),
                  aes(x = NMDS1, y = NMDS3, label = env), 
                  size = 2) 

month_with_env_13 <- month13 +
  geom_segment(data = filter(scores, pvals <= 0.05),
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS3),
               arrow = arrow(length = unit(0.02, "npc")), 
               color = "black") +
  geom_text_repel(data = filter(scores, pvals <= 0.05),
                  aes(x = NMDS1, y = NMDS3, label = env), 
                  size = 2) 

ggpubr::ggarrange(year_with_env_13, month_with_env_13, ncol=2, common.legend = F)
#ggsave("Figures/NMDS_3v1_all_dens_env.jpg", width=6, height=4) 

#---------------------------------------------------------------#
# summarize by year
ss_env <- env_drivers |>
  dplyr::select(-DateTime) |>
  group_by(month, year) |>
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), 
            .groups = "drop") |>
  filter(month %in% c(5:10))
  
zoops_plus_drivers_yearly <- bind_cols(monthly_zoops_nmds, ss_env[
  !colnames(ss_env) %in% c("month","year")]) |> 
  group_by(year) |> 
  summarise_at(vars(-month), funs(mean(., na.rm=TRUE)))

ord <- vegan::ordiplot(NMDS_bray_second,display = c('sites'),
                       choices = c(1,2),type = "n")
#fit environmental drivers onto ordination
fit_env <- envfit(ord, zoops_plus_drivers_yearly[,c(10:30)])

#pull out vectors - need to multiply by the sqrt of r2 to get magnitude!
scores <- data.frame((fit_env$vectors)$arrows * sqrt(fit_env$vectors$r), 
                     pvals=(fit_env$vectors)$pvals)
scores <- cbind(scores, env = rownames(scores))

#plot drivers w/ second stage NMDS (Manuscript Figure 4)
ss_year <- ggordiplots::gg_ordiplot(ord, unique(monthly_zoops_nmds$year),
                                  kind = "sd", ellipse=FALSE, hull = TRUE, 
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
        plot.margin = unit(c(0,1.3,0,0), 'lines'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(nrow=1, override.aes = list(
    color = year_cols)), fill = "none") +
  scale_fill_manual("",values=year_cols)+
  scale_color_manual("",values=year_cols,
                     label=c("2014","2015","2016",
                             "2019","2020","2021","2023")) +
  xlim(-0.7,0.9) + ylim(-1,0.9) +
  #geom_segment(data = scores,
  #             aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), linewidth= 0.3,
  #             arrow = arrow(length = unit(0.1, "cm")), colour = "lightgray") +
  geom_segment(data = filter(scores, pvals <= 0.05),
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), linewidth= 0.3,
               arrow = arrow(length = unit(0.1, "cm")), colour = "black") +
  geom_text_repel(data = filter(scores, pvals <= 0.05),
                  aes(x = NMDS1, y = NMDS2, label = env), 
                  size = 1.5, box.padding = 0.2, max.overlaps=Inf)

#axis 1 vs. 3
ord <- vegan::ordiplot(NMDS_bray_second,display = c('sites'),
                       choices = c(1,3),type = "n")
#fit environmental drivers onto ordination
fit_env <- envfit(ord, zoops_plus_drivers_yearly[,c(10:30)])

#pull out vectors - need to multiply by the sqrt of r2 to get magnitude!
scores <- data.frame((fit_env$vectors)$arrows * sqrt(fit_env$vectors$r), 
                     pvals=(fit_env$vectors)$pvals)
scores <- cbind(scores, env = rownames(scores))

ss_year <- ggordiplots::gg_ordiplot(ord, unique(all_zoops_nmds$year),
                                    kind = "sd", ellipse=FALSE, hull = TRUE, 
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
        legend.key.width =unit(0.1,"line"),
        plot.margin = unit(c(0,0.1,0,-1), 'lines'),
        legend.margin = margin(0,0,0,0),
        legend.spacing = unit(-0.5, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(
    color = year_cols)), fill = "none") +
  scale_fill_manual("",values=year_cols)+
  scale_color_manual("",values=year_cols,
                     label=c("2014","2015","2016",
                             "2019","2020","2021","2023")) +
  xlim(-0.7,0.9) + ylim(-1,0.9) +
  geom_segment(data = filter(scores, pvals <= 0.05),
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS3), linewidth= 0.3,
               arrow = arrow(length = unit(0.1, "cm")), colour = "black") +
  geom_text_repel(data = filter(scores, pvals <= 0.05),
                  aes(x = NMDS1, y = NMDS3, label = env), 
                  size = 1.5, box.padding = 0.2, max.overlaps=Inf)

ggpubr::ggarrange(env_plot1, env_plot2, ncol=2, common.legend = F)
#ggsave("Figures/second_stage_NMDS_envfit.jpg", width=6, height=3)

#---------------------------------------------------------------------------
# SIGNIFICANCE TESTS --> Is zoop community structure different among years?

#check assumptions
dist_bray <- vegdist(zoop_dens_trans_all, method = "bray")
bd_year <- betadisper(dist_bray, group = all_zoops_nmds$year)
plot(bd_year)           
permutest(bd_year, permutations = 999) #dispersion is similar

set.seed(3)
adonis2(zoop_bray_first ~ year, data = all_zoops_nmds, 
        permutations = 999, by = "margin") #not sig, no year differences

#same for SS NMDS
year_df <- data.frame(year = as.integer(rownames(zoop_bray_second)))

#look at SS pairwise year correlations
mean(correlation_matrix[lower.tri(correlation_matrix)])

# Assign names to the list of matrices
names(matrices) <- c("2014","2015","2016","2019","2020","2021","2023")

year_pairs <- combn(names(matrices), 2, simplify = FALSE)
cor_summary <- sapply(year_pairs, function(p) {
  i <- which(names(matrices) == p[1])
  j <- which(names(matrices) == p[2])
  correlation_matrix[i,j]
})
names(cor_summary) <- sapply(year_pairs, paste, collapse = "_vs_")
cor_summary # 2014 and 2021 have most similar summer mean dispersion

#test within-year variability
# Make a factor of year per row in your community matrix
year_factor <- monthly_zoops_nmds$year

# Bray–Curtis distances on Hellinger-transformed data
dist_matrix <- vegdist(zoop_dens_trans, method = "bray")

# test homogeneity of dispersion
disp <- betadisper(dist_matrix, year_factor)
anova(disp) #not significant so dispersion is homogenous
boxplot(disp)

