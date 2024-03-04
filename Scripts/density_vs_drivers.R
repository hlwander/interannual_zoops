#density vs. drivers

pacman::p_load(ggplot2, tidyverse, Hmisc, corrplot, boot)

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

#nonparametric analysis where months are numbers and drivers are zoop density and other env variables

#Q: Is zooplankton seasonal succession related to environmental drivers?

names <- c("cladocera_std_dens", "copepoda_std_dens", "rotifera_std_dens",
           "epi_tn", "hypo_tn", "epi_tp", "hypo_tp", "anoxic_depth",
           "epi_temp", "hypo_temp", "water_level", 
           "schmidt_stability", "buoyancy_frequency")

#spearman rank correlations for nonparametric data

#split up correlations by year
corr14 <- density_plus_drivers |> 
  filter(year %in% 2014) |> 
  select(-c(year, diff))
corr14_coef <- cor(corr14)

corr15 <- density_plus_drivers |> 
  filter(year %in% 2015) |> 
  select(-c(year, diff))
corr15_coef <- cor(corr15)

corr16 <- density_plus_drivers |> 
  filter(year %in% 2016) |> 
  select(-c(year, diff))
corr16_coef <- cor(corr16)

corr19 <- density_plus_drivers |> 
  filter(year %in% 2019) |> 
  select(-c(year, diff))
corr19_coef <- cor(corr19)

corr20 <- density_plus_drivers |> 
  filter(year %in% 2020) |> 
  select(-c(year, diff))
corr20_coef <- cor(corr20)

corr21 <- density_plus_drivers |> 
  filter(year %in% 2021) |> 
  select(-c(year, diff))
corr21_coef <- cor(corr21)


#rename rows + cols
rownames(corr14_coef) <- c("month",names)
rownames(corr15_coef) <- c("month",names)
rownames(corr16_coef) <- c("month",names)
rownames(corr19_coef) <- c("month",names)
rownames(corr20_coef) <- c("month",names)
rownames(corr21_coef) <- c("month",names)

colnames(corr14_coef) <- c("month",names)
colnames(corr15_coef) <- c("month",names)
colnames(corr16_coef) <- c("month",names)
colnames(corr19_coef) <- c("month",names)
colnames(corr20_coef) <- c("month",names)
colnames(corr21_coef) <- c("month",names)

#visualize corelations
corrplot(corr14_coef, tl.col = "black")
corrplot(corr15_coef, tl.col = "black")
corrplot(corr16_coef, tl.col = "black")
corrplot(corr19_coef, tl.col = "black")
corrplot(corr20_coef, tl.col = "black")
corrplot(corr21_coef, tl.col = "black")

#create a dataframe for each year and sum the positive and negative values
corr.df <- data.frame("year" = c(2014,2015,2016,2019,2020,2021),
                      "pos" = NA,
                      "neg" = NA)

#playing around with the cumulative "strength" of the correlation coefficients
corr.df$pos <- c(sum(corr14_coef[corr14_coef > 0]),
                 sum(corr15_coef[corr15_coef > 0]),
                 sum(corr16_coef[corr16_coef > 0]),
                 sum(corr19_coef[corr19_coef > 0]),
                 sum(corr20_coef[corr20_coef > 0]),
                 sum(corr21_coef[corr21_coef > 0]))

corr.df$neg <- c(sum(corr14_coef[corr14_coef < 0]),
                 sum(corr15_coef[corr15_coef < 0]),
                 sum(corr16_coef[corr16_coef < 0]),
                 sum(corr19_coef[corr19_coef < 0]),
                 sum(corr20_coef[corr20_coef < 0]),
                 sum(corr21_coef[corr21_coef < 0]))

# NEED TO WORK ON CODE BELOW

# Number of bootstrap samples
n_boot <- 1000

# Function to compute p-value for correlation between two variables
compute_p_value <- function(x, y) {
  observed_cor <- cor(x, y)
  if (length(unique(x)) <= 1 | length(unique(y)) <= 1) {
    p_value <- NA
  } else {
    boot_cor <- replicate(n_boot, {
      sample_indices <- sample(length(x), replace = TRUE)
      cor(x[sample_indices], y[sample_indices])
    })
    p_value <- mean(abs(boot_cor) >= abs(observed_cor))
  }
  return(p_value)
}

# Compute p-values for each pairwise correlation
p_values <- matrix(NA, ncol(corr14), ncol(corr14))
for (i in 1:(ncol(corr14) - 1)) {
  for (j in (i + 1):ncol(corr14)) {
    p_values[i, j] <- compute_p_value(corr14[[i]], corr14[[j]])
  }
}
