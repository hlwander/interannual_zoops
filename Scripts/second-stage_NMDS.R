# code for second-stage NMDS
# 20 December 2023

#read in packages
pacman::p_load(zoo,dplR,dplyr,tidyverse,ggplot2,ggpubr,sf,vegan)

#read in zoop data from EDI
inUrl1  <-"https://pasta-s.lternet.edu/package/data/eml/edi/1090/16/c7a04035b0a99adc489f5b6daec1cd52"
infile1 <-  tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")

zoops <- read.csv(infile1, header=T) |>
  filter(CollectionMethod=="Tow" & Reservoir %in% c("BVR") &
           StartDepth_m > 7.1) |> 
  select(-c(Site,EndDepth_m,CollectionMethod))

#split data into pre 2019 and post
zoops_2016_2018 <- zoops[as.Date(zoops$DateTime)<"2019-01-01",]
zoops_2019_2021 <- zoops[as.Date(zoops$DateTime)>="2019-01-01",]

#add daphnia (D. catawba, D. ambigua), calanoida (diaptomus) for 2014-2018 data
zoops_pre <- zoops_2016_2018 |> 
  group_by(Reservoir, DateTime, StartDepth_m) |> 
  summarise(Daphnia = sum(Density_IndPerL[
    Taxon %in% c("D. catawba","D. ambigua")]),
    Calanoida = sum(Density_IndPerL[
      Taxon %in% c("Diaptomus")]),
    Cyclopoida = sum(Density_IndPerL[
      Taxon %in% c("Cyclopoids")]),
    Nauplii = sum(Density_IndPerL[
      Taxon %in% c("Nauplii")]),
    Bosmina = sum(Density_IndPerL[
      Taxon %in% c("Bosmina")]),
    Ceriodaphnia = sum(Density_IndPerL[
      Taxon %in% c("Ceriodaphnia")]),
    Conochilus = sum(Density_IndPerL[
      Taxon %in% c("Conochilus")]),
    Keratella = sum(Density_IndPerL[
      Taxon %in% c("Keratella")]),
    Trichocerca = sum(Density_IndPerL[
      Taxon %in% c("Trichocerca")]),
    Kellicottia = sum(Density_IndPerL[
      Taxon %in% c("Kellicottia")]),
    Lecane = sum(Density_IndPerL[
      Taxon %in% c("Lecane")]))

#convert back to long
zoops_final_pre <- zoops_pre |> 
  group_by(Reservoir, DateTime, StartDepth_m) |> 
  pivot_longer(cols=Daphnia:Lecane,
               names_to = c("Taxon"),
               values_to = "Density_IndPerL") |> 
  filter(hour(DateTime) %in% c(9,10,11,12,13,14)) |> #drop nighttime samples
  mutate(DateTime = as.Date(DateTime)) |> 
  mutate(Density_IndPerL = Density_IndPerL * (1/0.031))  #10m bvr neteff from 2016 (n=2) - note that 7m neteff was also 0.31
#avg from 2020 and 2021 is 0.021...

#list common taxa between pre and post
taxa <- c("Daphnia","Calanoida","Cyclopoida","nauplius",
          "Bosmina","Ceriodaphnia","Conochilus","Keratella",
          "Trichocerca","Kellicottia")

#average reps when appropriate
zoops_final_post <- zoops_2019_2021 |> 
  mutate(DateTime = as.POSIXct(DateTime, format="%Y-%m-%d %H:%M:%S", tz="UTC")) |> 
  filter(hour(DateTime) %in% c(9,10,11,12,13,14)) |> #drop nighttime samples
  filter(Taxon %in% c(taxa)) |> 
  mutate(DateTime = as.Date(DateTime)) |> 
  mutate(Taxon = ifelse(Taxon=="nauplius", "Nauplii",Taxon)) |> 
  group_by(Reservoir, DateTime, StartDepth_m, Taxon) |> 
  summarise(Density_IndPerL = mean(Density_IndPerL))

#combine all zoop data
all_zoops <- bind_rows(zoops_final_pre, zoops_final_post) |> 
  mutate_all(~replace(., is.nan(.), NA)) |>  #replace NAN with NA
  ungroup() |> select(-StartDepth_m) #dropping, but note that depths range from 7.5-11.5m....

#add column for pre vs post
all_zoops$data <- ifelse(all_zoops$DateTime<="2019-01-01","pre","post")

#------------------------------------------------------------------------------#
#NMDS 

taxa <- unique(all_zoops$Taxon)

#taxa as cols, dates as rows, average by month
all_zoops_nmds <- all_zoops |> 
  select(DateTime, data, Taxon, Density_IndPerL) |> 
  filter(Taxon %in% taxa) |> 
  pivot_wider(names_from = Taxon, values_from = Density_IndPerL) |> 
  mutate(year = format(DateTime, "%Y"),
         month = format(DateTime, "%m")) |> 
  ungroup() |> group_by(year, month, data) |> 
  summarise(Bosmina = mean(Bosmina),
            Ceriodaphnia = mean(Ceriodaphnia),
            Daphnia = mean(Daphnia),
            Calanoida = mean(Calanoida),
            Cyclopoida = mean(Cyclopoida),
            Nauplii = mean(Nauplii),
            Conochilus = mean(Conochilus),
            Keratella = mean(Keratella),
            Kellicottia = mean(Kellicottia),
            Trichocerca = mean(Trichocerca)) |> 
  ungroup()

#only keep may-sep samples and drop 2022
all_zoops_nmds <- all_zoops_nmds |> 
  filter(month %in% c("05","06","07","08","09"), 
         !year %in% c("2022"))
#only 5 months per year bc needs to be equal for pairwise correlation matrix

#select only data cols
zoops_dens <- all_zoops_nmds |> select(Bosmina:Trichocerca)

#hellinger transform data
zoop_dens_trans <- labdsv::hellinger(zoops_dens)

#center the dataset - is this necessary for all multivariate approaches??
#zoop_dens_trans_centered <- scale(zoop_dens_trans, scale=FALSE)
#not doing this bc the metaMDS automatically centers data

#turn transformed community data into Euclidean distance matrix (bc negative values after centering)
#centered means that multiplying with a vector has the same effect as subtracting the mean of the components of the vector from every component of the vector
zoop_bray <- as.matrix(vegan::vegdist(zoop_dens_trans, method='bray'))
  
#------------------------------------------------------------------------------#
#first-stage NMDS

#scree plot to choose dimension 
#jpeg("figures/scree.jpg") 
goeveg::dimcheckMDS(zoop_bray, distance = "bray", 
                    k = 6, trymax = 20, autotransform = TRUE)
#dev.off()

set.seed(1)

#now do NMDS w/ 4 dimensions 
NMDS_bray_first <- vegan::metaMDS(zoop_bray, distance='bray', k=4, trymax=20, 
                            autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_bray_first$stress

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
#ggsave("Figures/first_stage_NMDS_2v1_years.jpg", width=5, height=3)

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
#ggsave("Figures/first_stage_NMDS_2v1_months.jpg", width=5, height=3)

#track months in each year and connect
#jpeg("Figures/first_stage_NMDS_2v1_2021.jpg")
ord <- vegan::ordiplot(NMDS_bray_first,display = c('sites','species'),
                       choices = c(1,2),type = "n")

lines(NMDS_bray_first$points[all_zoops_nmds$year== "2014" ,1], 
      NMDS_bray_first$points[(all_zoops_nmds$year== "2014" ),2], 
      col=viridis::viridis(6, option="D")[1])

points(NMDS_bray_first$points[all_zoops_nmds$year== "2014" ,1], 
       NMDS_bray_first$points[(all_zoops_nmds$year== "2014" ),2], 
       col=viridis::viridis(6, option="D")[1],
       pch = as.character(1:5), font=2, cex=1.5)

lines(NMDS_bray_first$points[all_zoops_nmds$year== "2015" ,1], 
      NMDS_bray_first$points[(all_zoops_nmds$year== "2015" ),2], 
      col=viridis::viridis(6, option="D")[2])

points(NMDS_bray_first$points[all_zoops_nmds$year== "2015" ,1], 
       NMDS_bray_first$points[(all_zoops_nmds$year== "2015" ),2], 
       col=viridis::viridis(6, option="D")[2],
       pch = as.character(1:5), font=2, cex=1.5)

lines(NMDS_bray_first$points[all_zoops_nmds$year== "2016" ,1], 
      NMDS_bray_first$points[(all_zoops_nmds$year== "2016" ),2], 
      col=viridis::viridis(6, option="D")[3])

points(NMDS_bray_first$points[all_zoops_nmds$year== "2016" ,1], 
       NMDS_bray_first$points[(all_zoops_nmds$year== "2016" ),2], 
       col=viridis::viridis(6, option="D")[3],
       pch = as.character(1:5), font=2, cex=1.5)

lines(NMDS_bray_first$points[all_zoops_nmds$year== "2019" ,1], 
      NMDS_bray_first$points[(all_zoops_nmds$year== "2019" ),2], 
      col=viridis::viridis(6, option="D")[4])

points(NMDS_bray_first$points[all_zoops_nmds$year== "2019" ,1], 
       NMDS_bray_first$points[(all_zoops_nmds$year== "2019" ),2], 
       col=viridis::viridis(6, option="D")[4],
       pch = as.character(1:5), font=2, cex=1.5)

lines(NMDS_bray_first$points[all_zoops_nmds$year== "2020" ,1], 
      NMDS_bray_first$points[(all_zoops_nmds$year== "2020" ),2], 
      col=viridis::viridis(6, option="D")[5])

points(NMDS_bray_first$points[all_zoops_nmds$year== "2020" ,1], 
       NMDS_bray_first$points[(all_zoops_nmds$year== "2020" ),2], 
       col=viridis::viridis(6, option="D")[5],
       pch = as.character(1:5), font=2, cex=1.5)

lines(NMDS_bray_first$points[all_zoops_nmds$year== "2021" ,1], 
      NMDS_bray_first$points[(all_zoops_nmds$year== "2021" ),2], 
      col=viridis::viridis(6, option="D")[6])

points(NMDS_bray_first$points[all_zoops_nmds$year== "2021" ,1], 
       NMDS_bray_first$points[(all_zoops_nmds$year== "2021" ),2], 
       col=viridis::viridis(6, option="D")[6],
       pch = as.character(1:5), font=2, cex=1.5)

#legend("topright", legend=c('2014','2015','2016','2019','2020','2021'),
#       pt.bg=c(viridis::viridis(5, option="D")) ,bty = "n", cex=1.2, pch=21) 
legend("topleft", legend=c('May','June','July','August','September'),
       pch=c(as.character(1:5)) ,bty = "n", cex=1.2) 
#dev.off()

#------------------------------------------------------------------------------#
#calculate dispersion between years vs. months - monte carlo approach

#DISP CODE HERE

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

set.seed(1)

#now do NMDS w/ 4 dimensions for consistency
NMDS_bray_second <- vegan::metaMDS(stage2, distance='bray', k=4, trymax=20, 
                            autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_bray_second$stress
# 0.02

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
  guides(color = guide_legend(override.aes = list(
    color = viridis::viridis(6, option="D")))) +
      scale_fill_manual("",values=viridis::viridis(6, option="D"),
                        label=c('2014','2015',"2016","2019","2020","2021")) +
      scale_color_manual("",values=rep("gray",6) )
#ggsave("Figures/second_stage_NMDS_2v1.jpg", width=5, height=3) 


#ANOSIM test on second-stage similarities
ano_year = anosim(stage2, all_zoops_nmds$year, distance = "bray", permutations = 9999)
#p > 0.05 so no sig difference among years

ano_month = anosim(stage2, all_zoops_nmds$month, distance = "bray", permutations = 9999)
#p < 0.01 so there is a difference between months!

#indicator species analysis to see which species are more abundant in one group vs another?
