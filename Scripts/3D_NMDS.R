#trying a 3D plot to visualize all NMDS dimensions

pacman::p_load(scatterplot3d, plotly)

#NMDS 
set.seed(11)

#now do NMDS w/ 4 dimensions 
NMDS_bray_first <- vegan::metaMDS(zoop_bray, distance='bray', k=3, trymax=20, 
                                  autotransform=FALSE, pc=FALSE, plot=FALSE)
NMDS_bray_first$stress
# 0.1

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
                         choices = c(1,2),type = "none", xaxt = 'n', yaxt = 'n')
  
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
  #if (i %in% c(4)){
   # legend(-0.3, 0.5, legend=c('May','June','July','August','September'),
    #       pch=c(as.character(1:5)) ,bty = "n", cex=1.8) 
  }
#}

mtext("NMDS1", side = 1, outer = TRUE, cex = 1.5, line = 0.9,
      col = "black")
mtext("NMDS2", side = 2, outer = TRUE, cex = 1.5, line = 0.5,
      col = "black")

#legend("topright", legend=c('2014','2015','2016','2019','2020','2021'),
#       pt.bg=c(viridis::viridis(5, option="D")) ,bty = "n", cex=1.2, pch=21) 

dev.off()

#3d plot

nmds_3d <- as.data.frame(NMDS_bray_first$points)
nmds_3d$year <- c(rep("2014",5),rep("2015",5),rep("2016",5),
                  rep("2019",5),rep("2020",5),rep("2021",5))
nmds_3d$month <- c(rep(c("5","6","7","8","9"),6))

#scatterplot3d(nmds_3d[,1:3], color=nmds_3d$year, type="b")

#plotly plot
fig <- plot_ly(nmds_3d, x = ~MDS1[year=="2014"], 
               y = ~MDS2[year=="2014"],
               z = ~MDS3[year=="2014"],
               text = ~month[year== "2014"],
               type = 'scatter3d', mode ='text+lines',
               textfont = list(color = '#003366', width = 1),
               line = list(color = '#003366', width = 1),
               name = "2014") 
fig <- fig |> layout(scene = list(xaxis = list(title = "x"),
                                   yaxis = list(title = "y"), 
                                   zaxis = list(title = "z"))) |> 
     add_trace(x = ~MDS1[year=="2015"], 
               y = ~MDS2[year=="2015"], 
               z = ~MDS3[year=="2015"],
               text = ~month[year== "2015"],
               type = 'scatter3d', mode ='text+lines',
               textfont = list(color = '#660000', width = 1),
               line = list(color = '#660000', width = 1),
               name = "2015")
fig <- fig |> add_trace(x = ~MDS1[year=="2016"], 
               y = ~MDS2[year=="2016"], 
               z = ~MDS3[year=="2016"],
               text = ~month[year== "2016"],
               type = 'scatter3d', mode ='text+lines',
               textfont = list(color = '#CC0000', width = 1),
               line = list(color = '#CC0000', width = 1),
               name = "2016")
fig <- fig |> add_trace(x = ~MDS1[year=="2019"], 
                        y = ~MDS2[year=="2019"], 
                        z = ~MDS3[year=="2019"],
                        text = ~month[year== "2019"],
                        type = 'scatter3d', mode ='text+lines',
                        textfont = list(color = '#0099CC', width = 1),
                        line = list(color = '#0099CC', width = 1),
                        name = "2019")
fig <- fig |> add_trace(x = ~MDS1[year=="2020"], 
                        y = ~MDS2[year=="2020"], 
                        z = ~MDS3[year=="2020"],
                        text = ~month[year== "2020"],
                        type = 'scatter3d', mode ='text+lines',
                        textfont = list(color = '#CC6666', width = 1),
                        line = list(color = '#CC6666', width = 1),
                        name = "2020")
fig <- fig |> add_trace(x = ~MDS1[year=="2021"], 
                        y = ~MDS2[year=="2021"], 
                        z = ~MDS3[year=="2021"],
                        text = ~month[year== "2021"],
                        type = 'scatter3d', mode ='text+lines',
                        textfont = list(color = '#339999', width = 1),
                        line = list(color = '#339999', width = 1),
                        name = "2021")

