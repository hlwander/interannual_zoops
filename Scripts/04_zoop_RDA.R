# Redundancy analysis for all BVR zoop data 2014-2021

pacman::p_load(vegan, tidyr, data.table, lubridate,
               rLakeAnalyzer, car, ggnewscale, dplyr)

year_cols <- c("#011f51","#06889b","#2E8B57","#fdfa66","#facd60","#f44034","#a13637")
# 2014, 2015, 2016, 2019, 2020, 2021, 2023

#read in zoop data
all_zoops_dens <- read.csv("Output/all_zoops_dens.csv",header = TRUE)

#list of all taxa
taxa <- unique(all_zoops_dens$Taxon)

#taxa as cols, dates as rows (n = 85)
all_zoops <- all_zoops_dens |> 
  select(DateTime, Taxon, dens) |> 
  filter(Taxon %in% taxa) |> 
  mutate(DateTime = as.Date(DateTime)) |> 
  pivot_wider(names_from = Taxon, values_from = dens) |> 
  mutate_all(~replace(., is.na(.), 0)) |> 
  ungroup() |> group_by(DateTime) |>
  summarise(Bosmina = mean(Bosmina),
            Daphnia = mean(Daphnia),
            Cyclopoida = mean(Cyclopoida),
            Nauplii = mean(Nauplii),
            Conochilus = mean(Conochilus),
            Conochiloides = mean(Conochiloides),
            Keratella = mean(Keratella),
            Kellicottia = mean(Kellicottia),
            Ploima = mean(Ploima),
            Polyarthra = mean(Polyarthra)) |> 
  ungroup() |>
  filter(!month(DateTime) %in% c(3,12)) #removing edge months with low sample size bc could not be imputed
        
#select only data cols
zoops_dens <- all_zoops |> select(Bosmina:Polyarthra)

#hellinger transform data
zoop_dens_trans <- labdsv::hellinger(zoops_dens)

#list of dates to match up with env data (n=85)
dates <- unique(all_zoops$DateTime)

#read in env drivers
env_drivers <- read.csv("Output/all_drivers.csv") |>
  dplyr::select(-c(Total_ugL, AirTemp, Longwave, Temp_C_epi)) 
#too many NAs that could not be imputed for secchi; total VIF too high

# Make sure the rows match between predictors and responses
stopifnot(nrow(env_drivers) == nrow(all_zoops))

#check gradient lengths to determine whether to do RDA vs. CCA
dca <- decorana(zoop_dens_trans)   
# axis lengths are < 3, so RDA should be good

#drop datetime cols
all_drivers_num <- env_drivers[ , !names(env_drivers) %in% c("DateTime")]

# Run the RDA
rda_model <- rda(zoop_dens_trans ~ ., data = all_drivers_num)

#next see whether env vars are colinear (VIF>5)
vif.cca(rda_model) #dropping epi temp, air temp, longwave, and total bc VIF > 7

# 42% of total zooplankton variation is explained by env variables (0.1236/0.2945)
# 58% is unexplained (0.1709/0.2945)
summary(rda_model)

# ANOVA to test significance
anova(rda_model, by = "axis", permutations = 999) 
#RDA1 and RDA2 explain a significant fraction of the variance!

#term
rda_anova <- anova(rda_model, by = "term", permutations = 999)
#hypo tn, epi tp, hypo temp, epi DO, wl, ss, secchi, and brown are significant drivers of zoop community structure

# Site (sample) scores
site_scores <- vegan::scores(rda_model, display = "sites", choices = c(1,2)) |>
  as.data.frame() |>
  mutate(year = format(all_zoops$DateTime, "%Y"))

# Species scores
species_scores <- vegan::scores(rda_model, display = "species", choices = c(1,2)) |>
  as.data.frame() |>
  rownames_to_column(var = "Species")

# List of significant drivers from anova(rda_model, by="term")
sig_drivers <- rownames(rda_anova)[!is.na(rda_anova$`Pr(>F)`) & 
                                     rda_anova$`Pr(>F)` < 0.05]

# Environmental (biplot) scores
env_scores <- vegan::scores(rda_model, display = "bp", choices = c(1,2)) |>
  as.data.frame() |>
  rownames_to_column(var = "variable") |>
  mutate(Significant = ifelse(variable %in% sig_drivers, "yes", "no"),
         RDA1_end = RDA1 * 2, #multiply arrows by a constant for visibility
         RDA2_end = RDA2 * 2)

ggplot() +
  geom_point(data = site_scores, aes(x = RDA1, y = RDA2, color=year), size = 2) +
  scale_color_manual(values = year_cols, name = NULL) + 
  ggnewscale::new_scale_color() +
  geom_segment(data = env_scores |> filter(Significant == "yes"), 
               aes(x = 0, y = 0, xend = RDA1_end, yend = RDA2_end), color = "black",
               arrow = arrow(length = unit(0.3,"cm")), show.legend = FALSE) +
  geom_text_repel(data = env_scores |> filter(Significant == "yes"), 
            aes(x = RDA1_end, y = RDA2_end, label = variable),
            color =  "black", size = 2, vjust = -0.5, show.legend = FALSE) +
  xlab(paste0("RDA1 (", round(summary(rda_model)$cont$importance[2,1]*100,1), "%)")) +
  ylab(paste0("RDA2 (", round(summary(rda_model)$cont$importance[2,2]*100,1), "%)")) +
  theme_minimal() + 
  theme(text = element_text(size=10, color = "black"), 
        legend.position = "right",
        axis.text = element_text(size=6, color="black"), 
        axis.text.x = element_text(vjust = 0.5), 
        axis.ticks.x = element_line(colour = c(rep("black",4), "transparent")), 
        strip.background = element_rect(fill = "transparent"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)) 
#ggsave("Figures/zoop_RDA.jpg", width=5, height=4) 

#now see how much variation the 4 sig env vars make up
rda_model_sig <- rda(zoop_dens_trans ~ TN_ugL_hypo + TP_ugL_epi +
                       Temp_C_hypo + DO_mgL_epi + waterlevel + SS +
                       secchi + Brown_ugL,
                 data = all_drivers_num)

# Overall RDA results
summary(rda_model_sig) #these vars explain 31% of the total variation (0.09/0.295)
