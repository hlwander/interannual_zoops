# Redundancy analysis for all BVR zoop data 2014-2021

pacman::p_load(vegan, tidyr, data.table, lubridate, tibble, ggrepel, cols4all,
               rLakeAnalyzer, car, ggnewscale, dplyr, ggplot2, patchwork)

#cb friendly year palette (2014, 2015, 2016, 2019, 2020, 2021, 2023)
year_cols <- c4a("cols4all.friendly7", n = 7)

#read in zoop data
all_zoops_dens <- read.csv("Output/all_zoops_dens.csv",header = TRUE) |>
  select(-c(Reservoir, sd)) |>
  filter(!Taxon %in% c("Cladocera","Copepoda","Rotifera","Crustacea"))

#list of all taxa
taxa <- unique(all_zoops_dens$Taxon)

#taxa as cols, dates as rows (n = 85)
all_zoops <- all_zoops_dens |> 
  mutate(DateTime = as.Date(DateTime)) |> 
  pivot_wider(names_from = Taxon, values_from = dens) |> 
  mutate_all(~replace(., is.na(.), 0)) |> 
  ungroup() |> 
  filter(!month(DateTime) %in% c(3,12)) #removing edge months with low sample size bc could not be imputed
#write.csv(all_zoops, "./Output/zoop_raw_dens.csv", row.names=FALSE)

#select only data cols
zoops_dens <- all_zoops |> select(Daphnia:Synchaeta)

#hellinger transform data
zoop_dens_trans <- labdsv::hellinger(zoops_dens)
#write.csv(zoop_dens_trans, "./Output/zoop_dens_trans.csv", row.names=FALSE)

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
rda_mod <- rda(zoop_dens_trans ~ ., data = all_drivers_num)

#next see whether env vars are colinear (VIF>5)
vif.cca(rda_mod) #dropping epi temp, air temp, longwave, and total bc VIF > 7

# 41% of total zooplankton variation is explained by env variables (0.1380/0.3405)
# 59% is unexplained (0.2025/0.3405)
summary(rda_mod)

# ANOVA to test significance
anova_rda_axis <- anova(rda_mod, by = "axis", permutations = 999) 
#RDA1 and RDA2 explain a significant fraction of the variance
#this is raw canonical variance though not total variance

axis_rda_df <- as.data.frame(anova_rda_axis) |>
  rownames_to_column("axis") |>
  rename(F_value = F, P_value = `Pr(>F)`) |>
  mutate(axis = str_replace(axis, "^Axis ", "Axis_")) |>
  filter(!axis %in% "Residual")  |>
  dplyr::select(-c(Df,Variance)) |>
  mutate(F_value = round(F_value, 2))

#add %var explained by each axis (Table S6)
eig <- eigenvals(rda_mod)
axis_rda_df <- axis_rda_df |>
  mutate(Variance_explained = round(100 * eig[1:nrow(axis_rda_df)] / sum(eig), 1))
#write.csv(axis_rda_df, "Output/RDA_axis_ANOVA.csv", row.names=FALSE)

#term (Table 1)
anova_rda_term <- anova(rda_mod, by = "term", permutations = 999)
#hypo tn, epi tp, hypo temp, epi DO, ss, secchi, and brown are significant drivers of zoop community structure

term_rda_df <- as.data.frame(anova_rda_term) |>
  rownames_to_column("term") |>
  rename(F_value = F, P_value = `Pr(>F)`) |>
  filter(!term %in% "Residual")  |>
  dplyr::select(-c(Df)) |>
  mutate(F_value = round(F_value, 2))

#add %var explained by each term
total_SS_rda <- sum(term_rda_df$Variance)  # sum of all constrained variance
term_rda_df <- term_rda_df |>
  mutate(Variance_pct = round(100 * Variance / total_SS_rda, 1)) |>
  select(-Variance) |>
  arrange(-Variance_pct)
#write.csv(term_rda_df, "Output/RDA_term_ANOVA.csv", row.names=FALSE)

# envfit to get arrow directions and p-values (uses same predictor table)
envfit_rda <- envfit(rda_mod, all_drivers_num, permutations = 999)

#------------------------------------------------------------------------------#
#now db rda
cap_mod <- capscale(zoop_dens_trans ~ ., data = all_drivers_num, 
                    distance = "bray")

summary(cap_mod)
# 38.7% of variability explained by env variables (4.382/11.312)
# 61.3% is unexplained (6.930/11.312) 

#permutation anovas to test sig
anova_dbrda_axis <- anova(cap_mod, by = "axis", permutations = 999)
#axes 1-3 are sig

#Table S7
axis_dbrda_df <- as.data.frame(anova_dbrda_axis) |>
  rownames_to_column("axis") |>
  rename(F_value = F, P_value = `Pr(>F)`) |>
  mutate(axis = str_replace(axis, "^(CAP|CAP\\s*|Axis)\\s*-?", "Axis_")) |>
  filter(!str_detect(axis, regex("^Residual$", ignore_case = TRUE))) |>
  mutate(Variance_pct = round(100 * SumOfSqs / sum(SumOfSqs, na.rm = TRUE), 1))|>
  dplyr::select(axis, F_value, P_value, Variance_pct)|>
  mutate(F_value = round(F_value, 2))
#write.csv(axis_dbrda_df, "Output/dbRDA_axis_ANOVA.csv", row.names=FALSE)

#term (Table 2)
anova_dbrda_term <- anova(cap_mod, by = "margin", permutations = 999) 
#sig drivers: TN epi and hypo, TP epi, hypo temp, epi DO, wl, thermo_depth, ss, bluegreen, brown, secchi

term_dbrda_df <- as.data.frame(anova_dbrda_term) |>
  rownames_to_column("term") |>
  rename(F_value = F, P_value = `Pr(>F)`) |>
  filter(!term %in% "Residual")  |>
  dplyr::select(-c(Df)) |>
  mutate(F_value = round(F_value, 2))

#add %var explained by each term
total_SS_dbrda <- sum(term_dbrda_df$SumOfSqs)  # sum of all constrained variance
term_dbrda_df <- term_dbrda_df |>
  mutate(Variance_pct = round(100 * SumOfSqs / total_SS_dbrda, 1)) |>
  select(-SumOfSqs) |>
  arrange(-Variance_pct)
#write.csv(term_dbrda_df, "Output/dbRDA_term_ANOVA.csv", row.names=FALSE)

cap_R2 <- RsquareAdj(cap_mod)$r.squared
cap_R2adj <- RsquareAdj(cap_mod)$adj.r.squared

envfit_cap <- envfit(cap_mod, all_drivers_num, permutations = 999)

# extract site scores
rda_sites <- vegan::scores(rda_mod, display = "sites", 
                           choices = 1:2, scaling = 2) |> as.data.frame()
cap_sites <- vegan::scores(cap_mod, display = "sites", 
                           choices = 1:2, scaling = 2) |> as.data.frame()

#extract species
rda_species <- vegan::scores(rda_mod, display = "species", 
                             choices = 1:2, scaling = 2) |> 
  as.data.frame() |> rownames_to_column("Taxon")
cap_species <- vegan::scores(cap_mod, display = "species", 
                             choices = 1:2, scaling = 2) |> 
  as.data.frame() |> rownames_to_column("Taxon") |>
  rename(RDA1 = CAP1, RDA2 = CAP2)

# rename columns consistently for plotting
colnames(rda_sites)[1:2] <- c("RDA1","RDA2")
colnames(cap_sites)[1:2] <- c("RDA1","RDA2")

# add year for plotting
rda_sites <- rda_sites |> rownames_to_column("Sample") |> 
  mutate(Year = format(all_zoops$DateTime, "%Y"))
cap_sites <- cap_sites |> rownames_to_column("Sample") |> 
  mutate(Year = format(all_zoops$DateTime, "%Y"))

# ensure same row order
stopifnot(nrow(rda_sites) == nrow(cap_sites))

#are db and normal rda ordinations similar? yes
pro <- procrustes(rda_sites[,c("RDA1","RDA2")], cap_sites[,c("RDA1","RDA2")])
summary(pro)
# RMSE of 0.11 is small relative to the ordination x-y ranges, so similar ordinations

protest_res <- protest(rda_sites[,c("RDA1","RDA2")], cap_sites[,c("RDA1","RDA2")], 
                       permutations = 999)
#significant procrustus correlation (0.9797) so ordinations are similar

# envfit
env_vec_rda <- as.data.frame(vegan::scores(envfit_rda, display = "vectors")) |> 
  rownames_to_column("variable") |> 
  mutate(Significant = ifelse(variable %in% rownames(anova_rda_term)[
    which(anova_rda_term$`Pr(>F)` < 0.05)], "yes", "no"),
         RDA1_end = RDA1 * 2, RDA2_end = RDA2 * 2)

env_vec_cap <- as.data.frame(vegan::scores(envfit_cap, display = "vectors")) |> 
  rownames_to_column("variable") |> 
  rename(RDA1 = CAP1, RDA2 = CAP2) |>
  mutate(Significant = ifelse(variable %in% rownames(anova_dbrda_term)[
    which(anova_dbrda_term$`Pr(>F)` < 0.05)], "yes", "no"),
    RDA1_end = RDA1 * 2, RDA2_end = RDA2 * 2)

rda_pct1 <- axis_rda_df$Variance_explained[1] 
rda_pct2 <- axis_rda_df$Variance_explained[2]  
cap_pct1 <- axis_dbrda_df$Variance_pct[1] 
cap_pct2 <- axis_dbrda_df$Variance_pct[2]  

all_x <- c(rda_sites$RDA1, cap_sites$RDA1)
all_y <- c(rda_sites$RDA2, cap_sites$RDA2)
xlim <- range(all_x, na.rm = TRUE) * 1.1
ylim <- range(all_y, na.rm = TRUE) * 1.1

# plotting function for manuscript Figure 4
make_plot <- function(sites_df, species_df, env_df, xlab_text, ylab_text) {
  ggplot() +
    geom_point(data = sites_df, aes(x = RDA1, y = RDA2, color = Year), size = 2) +
    scale_color_manual(values = year_cols, name = NULL,
                       guide = guide_legend(nrow = 1, byrow = TRUE, 
                                            override.aes = list(size = 3))) +
    geom_segment(data = env_df |> filter(Significant == "yes"),
                 aes(x = 0, y = 0, xend = RDA1_end, yend = RDA2_end),
                 arrow = arrow(length = unit(0.25, "cm")), color = "black") +
    geom_text_repel(data = env_df |> filter(Significant == "yes"),
                    aes(x = RDA1_end, y = RDA2_end, label = variable),
                    size = 2, color = "black", box.padding = 0.25, 
                    point.padding = 0.25, max.overlaps = 30, force = 1) +
    geom_segment(data = species_df,
                 aes(x = 0, y = 0, xend = RDA1 * 2, yend = RDA2 * 2),
                 arrow = arrow(length = unit(0.2, "cm")),
                 color = "darkred", alpha = 0.8) +
    geom_text_repel(data = species_df,
                    aes(x = RDA1 * 2, y = RDA2 * 2, label = Taxon),
                    size = 2, color = "darkred",
                    box.padding = 0.25, point.padding = 0.2, 
                    max.overlaps = 30, force = 1) +
    xlab(xlab_text) + ylab(ylab_text) + theme_minimal() +
    theme(text = element_text(size = 10),
          axis.text = element_text(size = 6),
          legend.position = "top", panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5))}

p_rda <- make_plot(rda_sites, rda_species, env_vec_rda,
                   xlab_text = paste0("RDA1 (", rda_pct1, "%)"),
                   ylab_text = paste0("RDA2 (", rda_pct2, "%)"))

p_cap <- make_plot(cap_sites, cap_species, env_vec_cap,
                   xlab_text = paste0("dbRDA1 (", cap_pct1, "%)"),
                   ylab_text = paste0("dbRDA2 (", cap_pct2, "%)"))

# combine side-by-side
(p_rda + p_cap) + plot_layout(guides = "collect") & 
  theme(legend.position = "top") 
#ggsave("Figures/zoop_RDA_vs_dbRDA.jpg", width=6, height=4) 

#-----------------------------------------------------------------------------
#variance partitioning to see whether physical or chemcial drivers explain more variance
# Table S8
# Define variable groups
physical_vars <- c("Temp_C_hypo", "waterlevel", "therm_depth", "oxy_depth", 
                   "SS", "secchi")  

chemical_vars <- c("TN_ugL_epi", "TN_ugL_hypo", "TP_ugL_epi", "TP_ugL_hypo", 
                   "DO_mgL_epi", "Green_ugL", "Bluegreen_ugL", "Brown_ugL",
                   "Mixed_ugL")  

met_vars <- c("Shortwave", "RelHum", "WindSpeed", "Rain")

# Function to calculate % variance explained by group (varpart-style)
group_varpart <- function(comm, env_data, group1_vars, group2_vars) {
  # Filter to only vars present in data
  g1 <- intersect(group1_vars, colnames(env_data))
  g2 <- intersect(group2_vars, colnames(env_data))
  
  X1 <- env_data[, g1, drop = FALSE]
  X2 <- env_data[, g2, drop = FALSE]
  
  vp <- varpart(comm, X1, X2)
  return(vp)}

# Subset each group to only significant drivers first
phys_rda <- all_drivers_num[, intersect(physical_vars, colnames(all_drivers_num)), drop = FALSE]
chem_rda <- all_drivers_num[, intersect(chemical_vars, colnames(all_drivers_num)), drop = FALSE]
met_rda  <- all_drivers_num[, intersect(met_vars,      colnames(all_drivers_num)), drop = FALSE]

# Check what's in each group before running
cat("Physical vars in RDA:", names(phys_rda), "\n")
cat("Chemical vars in RDA:", names(chem_rda), "\n")
cat("Met vars in RDA:",      names(met_rda),  "\n")

# Run varpart with pre-subsetted objects
rda_vp <- varpart(zoop_dens_trans, phys_rda, chem_rda, met_rda)
plot(rda_vp, digits = 2, Xnames = c("Physical", "Chemical", "Meteorological"))
summary(rda_vp)

varpart_table <- data.frame(
  Component = c("Physical (unique)", "Chemical (unique)", "Shared (physical + chemical)", "Residual"),
  Adj_R2 = c(
    round(rda_vp$part$indfract$Adj.R.square[1], 3),  # [a] unique physical
    round(rda_vp$part$indfract$Adj.R.square[2], 3),  # [b] unique chemical
    round(rda_vp$part$indfract$Adj.R.square[4], 3),  # [d] shared
    round(rda_vp$part$indfract$Adj.R.square[8], 3)   # [h] residual
  ))
#write.csv(varpart_table, "Output/RDA_varpart_grouped.csv", row.names = FALSE)
