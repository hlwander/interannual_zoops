# first stage NMDS without aggregating data

pacman::p_load(zoo,dplR,dplyr,tidyverse,ggplot2,ggpubr,
               vegan,FSA,rcompanion,NatParksPalettes,ggrepel,labdsv,
               goeveg, ggordiplots, cowplot)

#year colors
year_cols <- viridis(6, option = "D")  # 6 years

#read in all_zoops df
all_zoops_dens <- read.csv("Output/all_zoops_dens.csv",header = TRUE)

#NMDS 
taxa <- unique(all_zoops_dens$Taxon)

#taxa as cols, dates as rows 
all_zoops_nmds <- all_zoops_dens |> 
  select(DateTime, Taxon, dens) |> 
  filter(Taxon %in% taxa) |> 
  mutate(DateTime = as.Date(DateTime)) |> 
  pivot_wider(names_from = Taxon, values_from = dens) |> 
  mutate_all(~replace(., is.na(.), 0)) |>  
  ungroup() |> group_by(DateTime) |>
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
  ungroup() |>
  mutate(year = format(DateTime, "%Y"),
         month = format(DateTime, "%m"))

#only keep may-sep samples and drop 2022
all_zoops_nmds <- all_zoops_nmds |> 
  filter(year %in% c("2014","2015","2016",
                     "2019","2020", "2021")) |>
  group_by(month) |>
  filter(n() >= 3) |>
  ungroup()
# 11-15 samples collected per year
#write.csv(all_zoops_nmds, "./Output/zoop_raw_dens.csv", row.names=FALSE)

#check to make sure unequal sampling isn't problematic
table(all_zoops_nmds$year)
table(all_zoops_nmds$month)

all_zoops_nmds %>%
  group_by(year, month) %>% 
  summarise(n_samples = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = month, values_from = n_samples, values_fill = 0)

ggplot(all_zoops_nmds, aes(x = factor(year))) +
  geom_bar() +
  labs(x = "Year", y = "Number of samples",
       title = "Sampling effort by year") +
  theme_minimal()

#select only data cols
zoops_dens <- all_zoops_nmds |> select(Bosmina:Polyarthra)

#hellinger transform data
zoop_dens_trans <- labdsv::hellinger(zoops_dens)
#write.csv(zoop_dens_trans, "./Output/zoop_dens_trans.csv", row.names=FALSE)

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
# 0.15

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
             fill= viridis::viridis(6, option="D")) +
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
  scale_fill_manual("",values=viridis::viridis(6, option="D")) +#NatParksPalettes::natparks.pals('DeathValley', 6))+
  scale_color_manual("",values=viridis::viridis(6, option="D"),
                     label=c('2014','2015','2016','2019','2020','2021')) 

month <- ggordiplots::gg_ordiplot(ord, all_zoops_nmds$month,
                                  kind = "ehull", ellipse=FALSE, hull = TRUE, 
                                  plot = FALSE, pt.size=0.9) 
month <- month$plot + geom_point() + theme_bw() + 
  geom_polygon(data = month$df_hull, aes(x = x, y = y, fill = Group), 
               alpha=0.2) +
  geom_point(data=month$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, 
             fill=viridis::viridis(8, option="F")) +
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
  scale_fill_manual("",values=viridis::viridis(8, option="F"))+
  scale_color_manual("",values=viridis::viridis(8, option="F"),
                     label=c("April",'May',"June","July","August",
                             "September","October","November")) 

ggpubr::ggarrange(year,month,ncol=2, common.legend = F)
#ggsave("Figures/first_stage_NMDS_2v1_dens_alldata.jpg", width=5, height=3)

#plot axis 3 vs. 1 
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
             fill= viridis::viridis(6, option="D")) +
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
  scale_fill_manual("",values=viridis::viridis(6, option="D")) +#NatParksPalettes::natparks.pals('DeathValley', 6))+
  scale_color_manual("",values=viridis::viridis(6, option="D"),
                     label=c('2014','2015','2016','2019','2020','2021')) 

month <- ggordiplots::gg_ordiplot(ord, all_zoops_nmds$month,
                                  kind = "ehull", ellipse=FALSE, hull = TRUE, 
                                  plot = FALSE, pt.size=0.9) 
month <- month$plot + geom_point() + theme_bw() + 
  geom_polygon(data = month$df_hull, aes(x = x, y = y, fill = Group), 
               alpha=0.2) +
  geom_point(data=month$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, 
             fill=viridis::viridis(8, option="F")) +
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
  scale_fill_manual("",values=viridis::viridis(8, option="F"))+
  scale_color_manual("",values=viridis::viridis(8, option="F"),
                     label=c("April",'May',"June","July","August",
                             "September","October","November")) 

ggpubr::ggarrange(year,month,ncol=2, common.legend = F)
#ggsave("Figures/first_stage_NMDS_3v1_dens_alldata.jpg", width=5, height=3)


# check to see if there are year effects after controlling for month
adonis2(zoop_bray ~ year + month, data = all_zoops_nmds, method = "bray", 
        permutations = 999, by = "margin")
# year explains 17% of the variation and month explains 21%
# SO year effects are robust after controlling for months 

#make sure zoop_bray is a dist object
zoop_bray_dist <- vegan::vegdist(zoop_dens_trans, method='bray')

# variability differences among years
bd_year <- betadisper(zoop_bray_dist, all_zoops_nmds$year)
#highest dispersion --> 2020; lowest dispersion --> 2015

permutest(bd_year, permutations = 999)
# no significant differences in dispersion among years

# mean distance to centroid (variability per year)
disp_df <- data.frame(
  year = bd_year$group,
  dist_to_centroid = bd_year$distances
) %>%
  group_by(year) %>%
  summarise(mean_disp = mean(dist_to_centroid),
            sd_disp = sd(dist_to_centroid),
            n = n(),
            .groups = "drop")

# Plot within-year dispersion
ggplot(disp_df, aes(x = year, y = mean_disp)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_disp - sd_disp,
                    ymax = mean_disp + sd_disp),
                width = 0.2) +
  theme_bw() +
  labs(x = "",
       y = "Mean distance to centroid ± SD") +
  theme(text = element_text(size=12), 
        axis.text = element_text(size=10, color="black"), 
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
#ggsave("Figures/zoop_dispersion.jpg", width=6, height=3) 

#-----------------------------------------------------------------------------#
# SAME CODE EXCEPT NOW LOOKING AT MONTHLY AGGREGATIONS
#taxa as cols, dates as rows, 
all_zoops_nmds_monthly <- all_zoops_dens |> 
  select(DateTime, Taxon, dens) |> 
  filter(Taxon %in% taxa) |> 
  mutate(DateTime = as.Date(DateTime),
         year = format(DateTime, "%Y"),
         month = format(DateTime, "%m")) |> 
  pivot_wider(names_from = Taxon, values_from = dens) |> 
  mutate_all(~replace(., is.na(.), 0)) |>
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

#only keep months with more than 3 points
all_zoops_nmds_monthly <- all_zoops_nmds_monthly |> 
  filter(month %in% c("05","06","07","08","09"), #need to have the same number of months per year for ss nmds
         year %in% c("2014","2015","2016",
                     "2019","2020", "2021")) |>
  group_by(month) |>
  filter(n() >= 3) |>
  ungroup()
# 1-17 samples for each month grouping

#select only data cols
zoops_dens_monthly <- all_zoops_nmds_monthly |> select(Bosmina:Polyarthra)

#hellinger transform data
zoop_dens_trans_monthly <- labdsv::hellinger(zoops_dens_monthly)

#turn transformed community data into b-c distance matrix 
zoop_bray_monthly <- as.matrix(vegan::vegdist(zoop_dens_trans_monthly, method='bray'))

#------------------------------------------------------------------------------#
#first-stage NMDS
#scree plot to choose dimension 
#jpeg("figures/scree.jpg") (Supplemental Figure S1)
goeveg::dimcheckMDS(zoop_bray_monthly, distance = "bray", 
                    k = 6, trymax = 20, autotransform = TRUE)
#dev.off()

set.seed(11)
#now do NMDS w/ 3 dimensions 
NMDS_bray_first_monthly <- vegan::metaMDS(zoop_bray_monthly, distance='bray', k=3, trymax=20, 
                                  autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_bray_first_monthly$stress
# 0.11

#plot axis 2 vs. 1 monthly aggregated NMDS
ord <- vegan::ordiplot(NMDS_bray_first_monthly,display = c('sites'),
                       choices = c(1,2),type = "n")
year <- ggordiplots::gg_ordiplot(ord, all_zoops_nmds_monthly$year,
                                 kind = "ehull", ellipse=FALSE, hull = TRUE, 
                                 plot = FALSE, pt.size=0.9) 
year <- year$plot + geom_point() + theme_bw() + 
  geom_polygon(data = year$df_hull, aes(x = x, y = y, fill = Group), 
               alpha=0.2) +
  geom_point(data=year$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, 
             fill= NatParksPalettes::natparks.pals('DeathValley', 6)) +
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
  scale_fill_manual("",values=NatParksPalettes::natparks.pals('DeathValley', 6))+
  scale_color_manual("",values=NatParksPalettes::natparks.pals('DeathValley', 6),
                     label=c('2014','2015','2016','2019','2020','2021')) 

month <- ggordiplots::gg_ordiplot(ord, all_zoops_nmds_monthly$month,
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
                     label=c('May',"June","July","August",
                             "September")) 

ggpubr::ggarrange(year,month,ncol=2, common.legend = F)
#ggsave("Figures/first_stage_NMDS_2v1_dens_monthly.jpg", width=5, height=3)

#plot axis 3 vs. 1 monthly aggregated NMDS
ord <- vegan::ordiplot(NMDS_bray_first_monthly,display = c('sites'),
                       choices = c(1,3),type = "n")
year <- ggordiplots::gg_ordiplot(ord, all_zoops_nmds_monthly$year,
                                 kind = "ehull", ellipse=FALSE, hull = TRUE, 
                                 plot = FALSE, pt.size=0.9) 
year <- year$plot + geom_point() + theme_bw() + 
  geom_polygon(data = year$df_hull, aes(x = x, y = y, fill = Group), 
               alpha=0.2) +
  geom_point(data=year$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, 
             fill= NatParksPalettes::natparks.pals('DeathValley', 6)) +
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
  scale_fill_manual("",values=NatParksPalettes::natparks.pals('DeathValley', 6))+
  scale_color_manual("",values=NatParksPalettes::natparks.pals('DeathValley', 6),
                     label=c('2014','2015','2016','2019','2020','2021')) 

month <- ggordiplots::gg_ordiplot(ord, all_zoops_nmds_monthly$month,
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
                     label=c('May',"June","July","August",
                             "September")) 

ggpubr::ggarrange(year,month,ncol=2, common.legend = F)
#ggsave("Figures/first_stage_NMDS_3v1_dens_monthly.jpg", width=5, height=3)


#------------------------------------------------------------------------------#
# Second-stage NMDS 

#create distance matrices for all years (first stage pairwise dissimilarities)
dist2014 <- as.matrix(vegan::vegdist(zoop_dens_trans_monthly[all_zoops_nmds_monthly$year=="2014",],
                                     method='bray'))
dist2015 <- as.matrix(vegan::vegdist(zoop_dens_trans_monthly[all_zoops_nmds_monthly$year=="2015",],
                                     method='bray'))
dist2016 <- as.matrix(vegan::vegdist(zoop_dens_trans_monthly[all_zoops_nmds_monthly$year=="2016",],
                                     method='bray'))
dist2019 <- as.matrix(vegan::vegdist(zoop_dens_trans_monthly[all_zoops_nmds_monthly$year=="2019",],
                                     method='bray'))
dist2020 <- as.matrix(vegan::vegdist(zoop_dens_trans_monthly[all_zoops_nmds_monthly$year=="2020",],
                                     method='bray'))
dist2021 <- as.matrix(vegan::vegdist(zoop_dens_trans_monthly[all_zoops_nmds_monthly$year=="2021",],
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
year1 <- ggordiplots::gg_ordiplot(ord, unique(all_zoops_nmds_monthly$year), kind = "sd", 
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
year2 <- ggordiplots::gg_ordiplot(ord, unique(all_zoops_nmds_monthly$year), kind = "sd", 
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

