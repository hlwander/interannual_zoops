#density vs. drivers

pacman::p_load(ggplot2, tidyverse)

#read in env drivers
env <- read.csv("Output/env.csv",header=T) |> 
  arrange(year, month)

#read in density df
std_dens <- read.csv("Output/std_dens_3taxa.csv", header=T) |> 
  filter(month %in% c(5:9)) |> 
  select(year, month, Taxon, standardized_dens) |> 
  pivot_wider(names_from = Taxon, values_from = standardized_dens,
              names_glue = "{Taxon}_density_std")

#combine dens and driver dfs
density_plus_drivers <- left_join(std_dens, env)

#quick plots to look at drivers of density for each taxon
plot(density_plus_drivers$Cladocera_density_std~density_plus_drivers$DO_mgL_epi)

plot(density_plus_drivers$Copepoda_density_std~density_plus_drivers$NH4_ugL_hypo)

plot(density_plus_drivers$Rotifera_density_std~density_plus_drivers$Temp_C_epi)
