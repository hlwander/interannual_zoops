# Redundancy analysis for all BVR zoop data 2014-2021

pacman::p_load(vegan, tidyr, data.table, lubridate, tibble, ggrepel,
               rLakeAnalyzer, car, ggnewscale, dplyr, ggplot2, patchwork)

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
            Keratella = mean(Keratella),
            Kellicottia = mean(Kellicottia),
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
rda_mod <- rda(zoop_dens_trans ~ ., data = all_drivers_num)

#next see whether env vars are colinear (VIF>5)
vif.cca(rda_mod) #dropping epi temp, air temp, longwave, and total bc VIF > 7

# 40.9% of total zooplankton variation is explained by env variables (0.1115/0.2730)
# 59.1% is unexplained (0.1616/0.2730)
summary(rda_mod)

# ANOVA to test significance
anova_rda_axis <- anova(rda_mod, by = "axis", permutations = 999) 
#RDA1 and RDA2 explain a significant fraction of the variance!
#this is raw canonical variance though not total variance

axis_rda_df <- as.data.frame(anova_rda_axis) |>
  rownames_to_column("axis") |>
  rename(F_value = F, P_value = `Pr(>F)`) |>
  mutate(axis = str_replace(axis, "^Axis ", "Axis_")) |>
  filter(!axis %in% "Residual")  |>
  dplyr::select(-c(Df,Variance)) |>
  mutate(F_value = round(F_value, 2))

#add %var explained by each axis
eig <- eigenvals(rda_mod)
axis_rda_df <- axis_rda_df |>
  mutate(Variance_explained = round(100 * eig[1:nrow(axis_rda_df)] / sum(eig), 1))
#write.csv(axis_rda_df, "Output/RDA_axis_ANOVA.csv", row.names=FALSE)

#term
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
  select(-Variance)
#write.csv(term_rda_df, "Output/RDA_term_ANOVA.csv", row.names=FALSE)

# envfit to get arrow directions and p-values (uses same predictor table)
envfit_rda <- envfit(rda_mod, all_drivers_num, permutations = 999)

#------------------------------------------------------------------------------#
#now db rda
cap_mod <- capscale(zoop_dens_trans ~ ., data = all_drivers_num, 
                    distance = "bray")

summary(cap_mod)
# 37.8% of variability explained by env variables (3.150/8.329)
# 62.2% is unexplained (5.179/8.329) 

#permutation anovas to test sig
anova_dbrda_axis <- anova(cap_mod, by = "axis", permutations = 999)
#axes 1 and 2 are sig

axis_dbrda_df <- as.data.frame(anova_dbrda_axis) |>
  rownames_to_column("axis") |>
  rename(F_value = F, P_value = `Pr(>F)`) |>
  mutate(axis = str_replace(axis, "^(CAP|CAP\\s*|Axis)\\s*-?", "Axis_")) |>
  filter(!str_detect(axis, regex("^Residual$", ignore_case = TRUE))) |>
  mutate(Variance_pct = round(100 * SumOfSqs / sum(SumOfSqs, na.rm = TRUE), 1))|>
  dplyr::select(axis, F_value, P_value, Variance_pct)|>
  mutate(F_value = round(F_value, 2))
#write.csv(axis_dbrda_df, "Output/dbRDA_axis_ANOVA.csv", row.names=FALSE)

#term
anova_dbrda_term <- anova(cap_mod, by = "term", permutations = 999)
#sig drivers: TN epi and hypo, TP epi, hypo temp, epi DO, wl, ss, bluegreen, brown, secchi

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
  select(-SumOfSqs)
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
# RMSE of 0.093 is small relative to the ordination x-y ranges, so similar ordinations

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

# plotting function
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
#ggsave("Figures/zoop_RDA_vs_dbRDA.jpg", width=5, height=4) 

#variance partitioning over significant drivers
sig_drivers_rda <- env_vec_rda$variable[env_vec_rda$Significant=="yes"]
sig_drivers_dbrda <- env_vec_cap$variable[env_vec_cap$Significant=="yes"]

# Function to calculate % variance explained by each variable
partial_rda_var <- function(drivers, comm, env_data) {
  results <- data.frame(variable = drivers, R2adj = NA, p = NA)
  
  for (var in drivers) {
    cond_vars <- setdiff(drivers, var)
    X <- env_data[, var, drop = FALSE]
    Z <- if(length(cond_vars) > 0) env_data[, cond_vars, drop = FALSE] else NULL
    
    partial <- rda(comm, X, Z)
    an <- anova(partial)
    
    # adjusted R^2
    results[results$variable == var, "R2adj"] <- RsquareAdj(partial)$adj.r.squared
    results[results$variable == var, "p"] <- an$`Pr(>F)`[1]
  }
  return(results)
}

# RDA
rda_varpart <- partial_rda_var(sig_drivers_rda, zoop_dens_trans, all_drivers_num) |>
  mutate(R2adj = round(R2adj,3))
#write.csv(rda_varpart, "Output/RDA_varpart.csv", row.names=FALSE)

# dbRDA
dbrda_varpart <- partial_rda_var(sig_drivers_dbrda, zoop_dens_trans, all_drivers_num) |>
  mutate(R2adj = round(R2adj,3))
#write.csv(dbrda_varpart, "Output/dbRDA_varpart.csv", row.names=FALSE)
