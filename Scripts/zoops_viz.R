#nmds for density data grouped by year

#read in packages
pacman::p_load(zoo,dplR,dplyr,tidyverse,ggplot2,ggpubr,sf)

#read in zoop data from EDI
inUrl1  <- "https://pasta-s.lternet.edu/package/data/eml/edi/1090/14/c7a04035b0a99adc489f5b6daec1cd52" 
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
      Taxon %in% c("Lecane")]),
    Rotifera = sum(Density_IndPerL[
      Taxon %in% c("Total Rotifers")]),
    Cladocera = sum(Density_IndPerL[
      Taxon %in% c("Bosmina","D. catawba", "Chydorus","D. ambigua",
                   "Diaphanosoma","Ceriodaphnia")]),
    Copepoda = sum(Density_IndPerL[
      Taxon %in% c("Diaptomus","Nauplii", "Cyclopoids")]))

#convert back to long
zoops_final_pre <- zoops_pre |> 
  group_by(Reservoir, DateTime, StartDepth_m) |> 
  pivot_longer(cols=Daphnia:Copepoda,
               names_to = c("Taxon"),
               values_to = "Density_IndPerL") |> 
  filter(hour(DateTime) %in% c(9,10,11,12,13,14)) |> #drop nighttime samples
  mutate(DateTime = as.Date(DateTime)) |> 
  mutate(Density_IndPerL = Density_IndPerL * (1/0.031))  #10m bvr neteff from 2016 (n=2) - note that 7m neteff was also 0.31
#avg from 2020 and 2021 is 0.021...

#list common taxa between pre and post
taxa <- c("Bosmina", "Daphnia", "Ceriodaphnia",
          "Cyclopoida","Calanoida", "nauplius", 
          "Conochilidae","Keratella", "Rotifera",
          "Trichocercidae","Kellicottia", "Lecane",
          "Cladocera", "Copepoda")

#average reps when appropriate
zoops_final_post <- zoops_2019_2021 |> 
  mutate(DateTime = as.POSIXct(DateTime, format="%Y-%m-%d %H:%M:%S", tz="UTC")) |> 
  filter(hour(DateTime) %in% c(9,10,11,12,13,14)) |> #drop nighttime samples
  filter(Taxon %in% c(taxa)) |> 
  mutate(DateTime = as.Date(DateTime)) |> 
  mutate(Taxon = ifelse(Taxon=="nauplius", "Nauplii",
                        ifelse(Taxon=="Trichocercidae", "Trichocerca", 
                               ifelse(Taxon=="Conochilidae", "Conochilus", Taxon)))) |> 
  group_by(Reservoir, DateTime, StartDepth_m, Taxon) |> 
  summarise(Density_IndPerL = mean(Density_IndPerL))

#combine all zoop data
all_zoops <- bind_rows(zoops_final_pre, zoops_final_post) |> 
  mutate_all(~replace(., is.nan(.), NA)) |>  #replace NAN with NA
  ungroup() |> select(-StartDepth_m) #dropping, but note that depths range from 7.5-11.5m....

#add column for pre vs post
all_zoops$data <- ifelse(all_zoops$DateTime<="2019-01-01","pre","post")

#calculate proportion of total density for each taxa by day
all_zoops <- all_zoops |> group_by(DateTime) |> 
  mutate(n = sum(Density_IndPerL)) |> ungroup() |> 
           group_by(DateTime,Taxon) |> 
    mutate(prop = Density_IndPerL / sum(n))

#order data levels
all_zoops$data <- factor(all_zoops$data, levels=c("pre", "post"))

#create bar plots comparing taxa densities pre vs post
ggplot(data=subset(all_zoops, Taxon %in% 
                     c("Cladocera","Copepoda","Rotifera")), 
              aes(x=data, fill=Taxon ))+
  geom_bar(width = 1, stat="identity", aes(y=Density_IndPerL))
  #coord_polar("y")

#clads
ggplot(data=subset(all_zoops, Taxon %in% 
                     c("Daphnia","Ceriodaphnia","Bosmina")), 
       aes(x=data, fill=Taxon ))+theme_bw() +
  geom_bar(width = 1, stat="identity", aes(y=prop))

#copes
ggplot(data=subset(all_zoops, Taxon %in% 
                     c("Calanoida","Cyclopoida","Nauplii")), 
       aes(x=data, fill=Taxon ))+ theme_bw() +
  geom_bar(width = 1, stat="identity", aes(y=prop))

#rots
ggplot(data=subset(all_zoops, Taxon %in% 
                     c("Conochilus","Keratella","Kellicottia",
                       "Lecan","Trichocerca")), 
       aes(x=data, fill=Taxon ))+theme_bw() +
  geom_bar(width = 1, stat="identity", aes(y=prop))

#----------------------------------------------------------------------#
#NMDS 

taxa <- c("Bosmina", "Daphnia", "Ceriodaphnia",
          "Cyclopoida","Calanoida", "nauplius", 
          "Conochilidae","Keratella","Lecane",
          "Trichocercidae","Kellicottia")

#taxa as cols, dates as rows, average by month
all_zoops_nmds <- all_zoops |> 
  select(DateTime, data, Taxon, Density_IndPerL) |> 
  filter(Taxon %in% taxa) |> 
  pivot_wider(names_from = Taxon, values_from = Density_IndPerL) |> 
  mutate(year = format(DateTime, "%Y"),
          month = format(DateTime, "%m")) |> 
  ungroup() |> group_by(year, month, data) |> 
  summarise(Daphnia = mean(Daphnia),
            Calanoida = mean(Calanoida),
            Cyclopoida = mean(Cyclopoida),
            Bosmina = mean(Bosmina),
            Ceriodaphnia = mean(Ceriodaphnia),
            Keratella = mean(Keratella),
            Kellicottia = mean(Kellicottia),
            Lecane = mean(Lecane)) |> 
  ungroup()

#only keep apr-sep samples and drop 2022
all_zoops_nmds <- all_zoops_nmds |> 
  filter(month %in% c("04","05","06","07","08","09","10"),
         !year %in% c("2022"))

#select only data cols
zoops_dens <- all_zoops_nmds |> select(Daphnia:Lecane)

#hellinger transform data
zoop_dens_trans <- labdsv::hellinger(zoops_dens)

#turn transformed community data into Euclidean distance matrix
zoop_euc <- as.matrix(vegan::vegdist(zoop_dens_trans, method='euclidean'))

#scree plot to choose dimension 
#jpeg("figures/scree.jpg") 
goeveg::dimcheckMDS(zoop_euc, distance = "bray", 
                    k = 6, trymax = 20, autotransform = TRUE)
#dev.off()

set.seed(15)

#now do NMDS using averages w/ 4 dimensions for consistency
NMDS_bray <- vegan::metaMDS(zoop_euc, distance='bray', k=4, trymax=20, 
                                  autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_bray$stress

#--------------------------------------------------------------------------#
#NMDS plot

ord <- vegan::ordiplot(NMDS_bray,display = c('sites','species'),
                choices = c(1,2),type = "n")
year <- ggordiplots::gg_ordiplot(ord, all_zoops_nmds$year, kind = "ehull", 
                     ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

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

#month
month <- ggordiplots::gg_ordiplot(ord, all_zoops_nmds$month, kind = "ehull", 
                                 ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

month$plot + geom_point() + theme_bw() + 
  geom_polygon(data = month$df_hull, aes(x = x, y = y, fill = Group), 
               alpha=0.2) +
  geom_point(data=month$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, 
             fill=viridis::viridis(7, option="F")) +
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
  scale_fill_manual("",values=viridis::viridis(7, option="F"))+
  scale_color_manual("",values=viridis::viridis(7, option="F"),
                     label=c('April','May',"June","July","August","September","October")) 

#order post then pre
all_zoops_nmds$data <- factor(all_zoops_nmds$data, levels=c("post", "pre"))

#pre vs post
data <- ggordiplots::gg_ordiplot(ord, all_zoops_nmds$data, kind = "ehull", 
                                  ellipse=FALSE, hull = TRUE, plot = FALSE, pt.size=0.9) 

data$plot + geom_point() + theme_bw() + 
  geom_polygon(data = data$df_hull, aes(x = x, y = y, fill = Group), 
               alpha=0.2) +
  geom_point(data=data$df_mean.ord, aes(x, y), 
             color="black", pch=21, size=2, 
             fill=viridis::viridis(2, option="C")) +
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
  scale_fill_manual("",values=viridis::viridis(2, option="C"))+
  scale_color_manual("",values=viridis::viridis(2, option="C")) 

#track months in each year and connect
ord <- vegan::ordiplot(NMDS_bray,display = c('sites','species'),
                       choices = c(1,2),type = "n")

vegan::ordihull(ord, all_zoops_nmds$year, display = "sites", draw = c("polygon"),
         col = viridis::viridis(6, option="D"), alpha = 75,cex = 2)
lines(NMDS_bray$points[all_zoops_nmds$year== "2014" ,1], 
      NMDS_bray$points[(all_zoops_nmds$year== "2014" ),2], 
      col=viridis::viridis(6, option="D")[1])

points(NMDS_bray$points[all_zoops_nmds$year== "2014" ,1], 
      NMDS_bray$points[(all_zoops_nmds$year== "2014" ),2], 
      col=viridis::viridis(6, option="D")[1],
      pch = as.character(1:7), font=2, cex=1.5)

lines(NMDS_bray$points[all_zoops_nmds$year== "2015" ,1], 
      NMDS_bray$points[(all_zoops_nmds$year== "2015" ),2], 
      col=viridis::viridis(6, option="D")[2])

points(NMDS_bray$points[all_zoops_nmds$year== "2015" ,1], 
       NMDS_bray$points[(all_zoops_nmds$year== "2015" ),2], 
       col=viridis::viridis(6, option="D")[2],
       pch = as.character(2:7), font=2, cex=1.5)

lines(NMDS_bray$points[all_zoops_nmds$year== "2016" ,1], 
      NMDS_bray$points[(all_zoops_nmds$year== "2016" ),2], 
      col=viridis::viridis(6, option="D")[3])

points(NMDS_bray$points[all_zoops_nmds$year== "2016" ,1], 
       NMDS_bray$points[(all_zoops_nmds$year== "2016" ),2], 
       col=viridis::viridis(6, option="D")[3],
       pch = as.character(1:7), font=2, cex=1.5)

lines(NMDS_bray$points[all_zoops_nmds$year== "2019" ,1], 
      NMDS_bray$points[(all_zoops_nmds$year== "2019" ),2], 
      col=viridis::viridis(6, option="D")[4])

points(NMDS_bray$points[all_zoops_nmds$year== "2019" ,1], 
       NMDS_bray$points[(all_zoops_nmds$year== "2019" ),2], 
       col=viridis::viridis(6, option="D")[4],
       pch = as.character(1:7), font=2, cex=1.5)

lines(NMDS_bray$points[all_zoops_nmds$year== "2020" ,1], 
      NMDS_bray$points[(all_zoops_nmds$year== "2020" ),2], 
      col=viridis::viridis(6, option="D")[5])

points(NMDS_bray$points[all_zoops_nmds$year== "2020" ,1], 
       NMDS_bray$points[(all_zoops_nmds$year== "2020" ),2], 
       col=viridis::viridis(6, option="D")[5],
       pch = as.character(2:7), font=2, cex=1.5)

lines(NMDS_bray$points[all_zoops_nmds$year== "2021" ,1], 
      NMDS_bray$points[(all_zoops_nmds$year== "2021" ),2], 
      col=viridis::viridis(6, option="D")[6])

points(NMDS_bray$points[all_zoops_nmds$year== "2021" ,1], 
       NMDS_bray$points[(all_zoops_nmds$year== "2021" ),2], 
       col=viridis::viridis(6, option="D")[6],
       pch = as.character(1:7), font=2, cex=1.5)

legend("topright", legend=c('2014','2015','2016','2019','2020','2021'),
       pt.bg=c(viridis::viridis(6, option="D")) ,bty = "n", cex=1.2, pch=21) 
legend("topleft", legend=c('April','May','June','July','August','September','October'),
       pch=c(as.character(1:7)) ,bty = "n", cex=1.2) 

