# code for second-stage NMDS
# 20 December 2023

#read in packages
pacman::p_load(zoo,dplR,dplyr,tidyverse,ggplot2,ggpubr,
               sf,vegan,FSA,rcompanion,NatParksPalettes,ggrepel)

#read in all_zoops df
all_zoops_dens <- read.csv("Output/all_zoops_dens.csv",header = TRUE)
#all_zoops_biom <- read.csv("Output/all_zoops_biom.csv",header = TRUE)

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

#center the dataset means that multiplying with a vector has the same effect as subtracting the mean of the components of the vector from every component of the vector
#zoop_dens_trans_centered <- scale(zoop_dens_trans, scale=FALSE)
#not doing this bc the metaMDS automatically centers data

#turn transformed community data into b-c distance matrix 
zoop_bray <- as.matrix(vegan::vegdist(zoop_dens_trans, method='bray'))
  
#------------------------------------------------------------------------------#
#first-stage NMDS

#scree plot to choose dimension 
#jpeg("figures/scree_first_stage.jpg") 
goeveg::dimcheckMDS(zoop_bray, distance = "bray", 
                    k = 6, trymax = 20, autotransform = TRUE)
#dev.off()

set.seed(11)

#now do NMDS w/ 4 dimensions 
NMDS_bray_first <- vegan::metaMDS(zoop_bray, distance='bray', k=4, trymax=20, 
                            autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_bray_first$stress
# 0.076

#plot
ord <- vegan::ordiplot(NMDS_bray_first,display = c('sites','species'),
                       choices = c(1,2),type = "n")
year <- ggordiplots::gg_ordiplot(ord, all_zoops_nmds$year,
                                  kind = "ehull", ellipse=FALSE, hull = TRUE, 
                                  plot = FALSE, pt.size=0.9) 
year$plot + geom_point() + theme_bw() + 
  geom_polygon(data = year$df_hull, aes(x = x, y = y, fill = Group), 
               alpha=0.2) +
  geom_point(data=year$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, 
             fill=viridis::viridis(6, option="D")) +
  theme(text = element_text(size=14), 
        axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "vertical",
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "right", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,-0.1,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) +
  scale_fill_manual("",values=viridis::viridis(6, option="D"))+
  scale_color_manual("",values=viridis::viridis(6, option="D"),
                     label=c('2014','2015',"2016","2019","2020","2021")) 
#ggsave("Figures/first_stage_NMDS_2v1_years_dens.jpg", width=5, height=3)

month <- ggordiplots::gg_ordiplot(ord, all_zoops_nmds$month,
                                 kind = "ehull", ellipse=FALSE, hull = TRUE, 
                                 plot = FALSE, pt.size=0.9) 
month$plot + geom_point() + theme_bw() + 
  geom_polygon(data = month$df_hull, aes(x = x, y = y, fill = Group), 
               alpha=0.2) +
  geom_point(data=month$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, 
             fill=viridis::viridis(5, option="F")) +
  theme(text = element_text(size=14), 
        axis.text = element_text(size=7, color="black"), 
        legend.background = element_blank(), 
        legend.key.height=unit(0.3,"line"),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "vertical",
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = "right", legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,-0.1,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) + guides(fill='none') +
  scale_fill_manual("",values=viridis::viridis(5, option="F"))+
  scale_color_manual("",values=viridis::viridis(5, option="F"),
                     label=c('May',"June","July","August","September")) 
#ggsave("Figures/first_stage_NMDS_2v1_months_dens.jpg", width=5, height=3)

#track months in each year and connect
#jpeg("Figures/first_stage_NMDS_2v1_dens.jpg")
par(mfrow = c(2, 3))
par(cex = 0.6)
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

for (i in 1:6) {
  ord <- vegan::ordiplot(NMDS_bray_first,display = c('sites','species'),
                         choices = c(1,2),type = "n", xaxt = 'n', yaxt = 'n')
  
  lines(NMDS_bray_first$points[all_zoops_nmds$year == 
                                 unique(all_zoops_nmds$year)[i] ,1], 
        NMDS_bray_first$points[(all_zoops_nmds$year ==
                                  unique(all_zoops_nmds$year)[i]),2], 
        col=viridis::viridis(6, option="D")[i])
  
  points(NMDS_bray_first$points[all_zoops_nmds$year == 
                                  unique(all_zoops_nmds$year)[i] ,1], 
         NMDS_bray_first$points[(all_zoops_nmds$year == 
                                   unique(all_zoops_nmds$year)[i]),2], 
         col=viridis::viridis(6, option="D")[i],
         pch = as.character(1:5), font=2, cex=3)
  
  mtext(c("2014","2015","2016","2019","2020","2021")[i], side = 3, line = -2, 
         adj = 0.05, cex = 1.5, col = "black")
 # mtext(c("clockwise","??","counterclockwise","clockwise",
 #         "counterclockwise","clockwise")[i], side = 3, line = -4, 
 #       adj = 0.05, cex = 1.3, col = c("#01586D","black", "#8B0C13", 
 #                                     "#01586D", "#8B0C13", "#01586D")[i])
   #if (i %in% c(4, 5, 6))
    # axis(1, col = "black", col.axis = "black")
   #if (i %in% c(1, 4))
   #  axis(2, col = "black", col.axis = "black")
   #  box(col = "black")
  #if (i %in% c(4))
    #legend("bottomleft", legend=c('May','June','July','August','September'),
     #      pch=c(as.character(1:5)) ,bty = "n", cex=1.5) 
}

mtext("NMDS1", side = 1, outer = TRUE, cex = 1.5, line = 0.9,
       col = "black")
mtext("NMDS2", side = 2, outer = TRUE, cex = 1.5, line = 0.5,
         col = "black")
  
#legend("topright", legend=c('2014','2015','2016','2019','2020','2021'),
#       pt.bg=c(viridis::viridis(5, option="D")) ,bty = "n", cex=1.2, pch=21) 

#dev.off()

#------------------------------------------------------------------------------#
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
kw_disp <- kruskal.test(value ~ group, data = disp_df) #ns
#month and year are NOT different from each other

#plot
ggboxplot(disp_df, x = "group", y = "value", 
                      fill = "group", palette = c("#A4C6B8", "#5E435D"),
                      order = c("year_disp", "month_disp"),
                      ylab = "Dispersion") +
  theme(text = element_text(size=9),
        plot.margin = unit(c(0.2,0,-0.5,0), 'lines')) +
  #annotate("text",label=c("a","b"), x=c(1.1,2.1),
  #         y=c(mean(disp_df$value[disp_df$group=="year_disp"]) + 
  #               sd(disp_df$value[disp_df$group=="year_disp"])/2,
  #             mean(disp_df$value[disp_df$group=="month_disp"]) + 
  #               sd(disp_df$value[disp_df$group=="month_disp"]))) +
  guides(fill = "none") +
  scale_x_discrete(name ="", 
                   labels=c("year_disp"="year",
                            "month_disp"="month"))
#ggsave("Figures/among_variability_boxplots_dens.jpg", width=3, height=4)

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
                      fill = "group", palette = viridis::viridis(6, option="D"),
                      #palette = natparks.pals("CapitolReef",6),
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
                     fill = "group", palette = viridis::viridis(6, option="F"),
                     #palette = natparks.pals("DeathValley",5),
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
                              "n" = c(5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6),
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
alldist <- as.data.frame(cbind(dist2014, dist2015, dist2016,
                              dist2019, dist2020, dist2021))

#perform pairwise correlations among dissimilarity matrices
#and convert from similarity to distance 
stage2 <- as.matrix(vegan::vegdist(1-cor(alldist)), method="bray")

#run NMDS again - clustering will indicate years where changes through time are similar
#scree plot to choose dimension 
#jpeg("figures/scree.jpg") 
goeveg::dimcheckMDS(stage2, distance = "bray", 
                    k = 4, trymax = 20, autotransform = TRUE)
#dev.off()

set.seed(11)

#now do NMDS w/ 4 dimensions for consistency
NMDS_bray_second <- vegan::metaMDS(stage2, distance='bray', k=4, trymax=20, 
                            autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_bray_second$stress
# 0.014

#--------------------------------------------------------------------------#
#NMDS plot - second-stage

ord <- vegan::ordiplot(NMDS_bray_second,display = c('sites','species'),
                       choices = c(1,2),type = "n")
year <- ggordiplots::gg_ordiplot(ord, all_zoops_nmds$year, kind = "ehull", 
                                 spiders = TRUE, ellipse = FALSE,
                                 label = FALSE, hull = FALSE, 
                                 plot = FALSE, pt.size=0.9) 

year$plot + geom_point() + theme_bw() + 
      geom_point(data=year$df_mean.ord, aes(x, y), 
                pch=21, size=2, 
                 fill=viridis::viridis(6, option="D")) +
      theme(text = element_text(size=8), 
            axis.text = element_text(size=6, color="black"), 
            axis.title = element_text(size=6),
            legend.background = element_blank(), 
            legend.key.height=unit(0.3,"line"),
            legend.box.margin=margin(-10,-10,-10,-10),
            legend.margin=margin(-0,-0,-0,-0),
            legend.direction = "vertical",
            axis.text.x = element_text(vjust = 0.5), 
            axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
            strip.background = element_rect(fill = "transparent"), 
            legend.position = "right", legend.spacing = unit(-0.5, 'cm'),
            plot.margin = unit(c(0,-0.1,0,0), 'lines'),
            panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
            legend.key.width =unit(0.1,"line")) +
  guides(color = guide_legend(override.aes = list(
    color = viridis::viridis(6, option="D")))) +
      scale_fill_manual("",values=viridis::viridis(6, option="D"),
                        label=c('2014','2015',"2016","2019","2020","2021")) +
      scale_color_manual("",values=rep("gray",6) )
#ggsave("Figures/second_stage_NMDS_2v1_dens.jpg", width=3, height=1) 


#ANOSIM test on second-stage similarities
ano_year = anosim(stage2, all_zoops_nmds$year, distance = "bray", permutations = 9999)
#p > 0.05 so no sig difference among years

ano_month = anosim(stage2, all_zoops_nmds$month, distance = "bray", permutations = 9999)
#p < 0.01 so there is a difference among months!

#------------------------------------------------------------------------------#
#read in env csv
env_drivers <- read.csv("./Output/env.csv")

#get zoops into same order as drivers
all_zoops_nmds <- all_zoops_nmds |> arrange(month, year)

#join driver and env dfs
zoops_plus_drivers <- bind_cols(all_zoops_nmds, env_drivers[
  !colnames(env_drivers) %in% c("month", "year")])

#fit environmental drivers onto ordination
fit_env <- envfit(ord$sites, zoops_plus_drivers[,c(13:23, 25:26)]) #dropping oxycline depth bc 2 NAs

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
driver_correlation <- data.frame(cor(zoops_plus_drivers[,c(13:23, 25:26)], method = "spearman"))
#write.csv(driver_correlation, "Output/driver_correlation.csv", row.names=F)

#plot drivers w/ NMDS
ord <- vegan::ordiplot(NMDS_bray_first,display = c('sites','species'),
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
  geom_segment(data = scores,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), linewidth= 0.3,
               arrow = arrow(length = unit(0.1, "cm")), colour = "black") +
  geom_text_repel(data = scores, aes(x = NMDS1, y = NMDS2, label = env), size = 1.5)
#ggsave("Figures/first_stage_NMDS_2v1_months_envfit.jpg", width=4, height=4)
  

ord <- vegan::ordiplot(NMDS_bray_first,display = c('sites','species'),
                       choices = c(1,2),type = "n")
year <- ggordiplots::gg_ordiplot(ord, all_zoops_nmds$year,
                                 kind = "ehull", ellipse=FALSE, hull = TRUE, 
                                 plot = FALSE, pt.size=0.9) 
year$plot + geom_point() + theme_bw() + 
  geom_polygon(data = year$df_hull, aes(x = x, y = y, fill = Group), 
               alpha=0.2) +
  geom_point(data=year$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, 
             fill=viridis::viridis(6, option="D")) +
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
        legend.position = c(0.86,0.2), legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,-0.1,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.1,"line")) +
  scale_fill_manual("",values=viridis::viridis(6, option="D"))+
  scale_color_manual("",values=viridis::viridis(6, option="D"),
                     label=c('2014','2015',"2016","2019","2020","2021")) +
  geom_segment(data = scores,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), linewidth= 0.3,
               arrow = arrow(length = unit(0.1, "cm")), colour = "black") +
  geom_text_repel(data = scores, aes(x = NMDS1, y = NMDS2, label = env), size = 1.5)
#ggsave("Figures/first_stage_NMDS_2v1_years_envfit.jpg", width=4, height=3)

#------------------------------------------------------------------------------#
# Which drivers explain patterns in zooplankton seasoanl succession among years?

#read in driver df
zoop_drivers <- read.csv("Output/all_drivers.csv") |> select(-diff) |> 
  rename(buoyancy_frequency = BF,
         schmidt_stability = SS,
         water_level = waterlevel,
         dissolved_oxygen_epi = DO_mgL_epi,
         thermocline_depth = therm_depth,
         total_nitrogen_epi = TN_ugL_epi,
         total_phosphorus_epi = TP_ugL_epi)
#delta wl is the same for every month so dropping this var

#group years based on hysteresis direction
zoop_drivers$box <- ifelse(zoop_drivers$year %in% c("2014", "2019", "2021"),
                           "clockwise", ifelse(zoop_drivers$year %in% c(
                             "2016", "2020"), "counterclockwise","mixed"))

#convert from wide to long
zoop_drivers_long <- zoop_drivers |> 
  pivot_longer(-c(month,year,box),
               names_to = "variable") |> 
  mutate(year = as.factor(year))

#change order of boxplots
zoop_drivers_long$year <- factor(zoop_drivers_long$year , 
                                levels=c("2014", "2019", "2021", 
                                         "2016", "2020", "2015"))

#now look at boxplots for each driver
ggplot(data=subset(zoop_drivers_long, !box %in% c("mixed")),
       aes(x=year, y=value, group = year)) +
  geom_boxplot(aes(fill=box, alpha = 0.95)) + 
  facet_wrap(~variable, nrow=5, scales = "free_y") +
  theme_bw() + xlab("") + guides(alpha = "none") +
  scale_fill_manual("",values=c("#01586D", "#8B0C13"))+
  theme(text = element_text(size=10), 
        axis.text = element_text(size=7, color="black"), 
        legend.key.height=unit(0.3,"line"),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(-0,-0,-0,-0),
        legend.direction = "vertical",
        legend.title = element_blank(),
        axis.text.x = element_text(angle=45, vjust=0.8, hjust=0.8),
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        legend.position = c(0.86,0.03), legend.spacing = unit(-0.5, 'cm'),
        plot.margin = unit(c(0,0,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.key.width =unit(0.3,"line"))
#ggsave("Figures/drivers_vs_hysteresis_boxplot_years.jpg", width=6, height=5)

#hand-picking the vars that are most different
vars <- c("buoyancy_frequency","dissolved_oxygen_epi", "thermocline_depth",
           "total_nitrogen_epi","total_phosphorus_epi","water_level")

#just group by direction
ggplot(data=subset(zoop_drivers_long, !box %in% c("mixed") &
                     variable %in% vars),
       aes(x=box, y=value, group = box)) +
  geom_boxplot(aes(fill=box, alpha = 0.95)) + 
  facet_wrap(~variable, nrow=5, scales = "free_y") +
  theme_bw() + xlab("") + guides(alpha = "none", fill = "none") +
  scale_fill_manual("",values=c("#01586D", "#8B0C13"))+
  theme(text = element_text(size=15), 
        axis.text = element_text(size=12, color="black"), 
        axis.text.x = element_text(angle=45, vjust=0.7, hjust=0.6),
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        plot.margin = unit(c(0,0,0,0), 'lines'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#ggsave("Figures/drivers_vs_hysteresis_boxplot_6.jpg", width=6, height=5)

#kw tests to see if the vars are significantly different 
kruskal_test_hyst_bf <- 
  kruskal.test(zoop_drivers$buoyancy_frequency[
    !zoop_drivers$box %in% c("mixed")] ~ zoop_drivers$box[
      !zoop_drivers$box %in% c("mixed")]) #ns

kruskal_test_hyst_do <- 
  kruskal.test(zoop_drivers$dissolved_oxygen_epi[
    !zoop_drivers$box %in% c("mixed")] ~ zoop_drivers$box[
      !zoop_drivers$box %in% c("mixed")]) #sig (but maybe not after bf correction..)

kruskal_test_hyst_td <- 
  kruskal.test(zoop_drivers$thermocline_depth[
    !zoop_drivers$box %in% c("mixed")] ~ zoop_drivers$box[
      !zoop_drivers$box %in% c("mixed")]) #ns

kruskal_test_hyst_tn <- 
  kruskal.test(zoop_drivers$total_nitrogen_epi[
    !zoop_drivers$box %in% c("mixed")] ~ zoop_drivers$box[
      !zoop_drivers$box %in% c("mixed")]) #sig (but maybe not after bf correction..)

kruskal_test_hyst_tp <- 
  kruskal.test(zoop_drivers$total_phosphorus_epi[
    !zoop_drivers$box %in% c("mixed")] ~ zoop_drivers$box[
      !zoop_drivers$box %in% c("mixed")]) #ns

kruskal_test_hyst_wl <- 
  kruskal.test(zoop_drivers$water_level[
    !zoop_drivers$box %in% c("mixed")] ~ zoop_drivers$box[
      !zoop_drivers$box %in% c("mixed")]) #ns
