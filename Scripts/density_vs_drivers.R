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
corr14_coef <- cor(corr14, method = "spearman")

corr15 <- density_plus_drivers |> 
  filter(year %in% 2015) |> 
  select(-c(year, diff))
corr15_coef <- cor(corr15, method = "spearman")

corr16 <- density_plus_drivers |> 
  filter(year %in% 2016) |> 
  select(-c(year, diff))
corr16_coef <- cor(corr16, method = "spearman")

corr19 <- density_plus_drivers |> 
  filter(year %in% 2019) |> 
  select(-c(year, diff))
corr19_coef <- cor(corr19, method = "spearman")

corr20 <- density_plus_drivers |> 
  filter(year %in% 2020) |> 
  select(-c(year, diff))
corr20_coef <- cor(corr20, method = "spearman")

corr21 <- density_plus_drivers |> 
  filter(year %in% 2021) |> 
  select(-c(year, diff))
corr21_coef <- cor(corr21, method = "spearman")


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


# Initialize an empty matrix to store p-values
p_values <- matrix(nrow = ncol(density_plus_drivers), 
                   ncol = ncol(density_plus_drivers))

# Convert dataframe to numeric
density_plus_drivers <- as.data.frame(sapply(density_plus_drivers, as.numeric))

# Number of bootstrap samples
num_bootstraps <- 1000

# Loop through each pair of variables and calculate bootstraped correlation and adjusted p-value
for (i in 1:(ncol(density_plus_drivers)-1)) {
  for (j in (i+1):ncol(density_plus_drivers)) {
    # Bootstrap resampling
    boot_cor <- function(data, indices) {
      cor(data[indices, i], data[indices, j], method = "spearman")
    }
    
    boot_result <- boot(density_plus_drivers, boot_cor, R = num_bootstraps)
    
    # Calculate p-value from bootstrap distribution
    p_value <- mean(boot_result$t >= 0)  # Two-tailed test
    
    # Adjust p-value for multiple testing (Bonferroni correction)
    p_adj <- p_value * (ncol(density_plus_drivers) * 
                          (ncol(density_plus_drivers) - 1)) / 2
    
    # Store adjusted p-value
    p_values[i, j] <- p_adj
    p_values[j, i] <- p_adj
  }
}


names <- c("year","month","cladocera_std_dens", "copepoda_std_dens", 
           "rotifera_std_dens", "epi_tn", "hypo_tn", "epi_tp", 
           "hypo_tp", "anoxic_depth", "epi_temp", "hypo_temp",
           "water_level", "diff", "schmidt_stability", "buoyancy_frequency")

#add col and row names
rownames(p_values) <- c(names)
colnames(p_values) <- c(names)

#visualize correlation matrix
#jpeg("Figures/spearman_matrix.jpg")
corrplot(cor(density_plus_drivers, method = "spearman"), tl.col = "black")
#dev.off()

#initialize df for kw p vals
kw_pvals <- matrix(nrow = 16)

#Kruskal-wallis test to look for differences among months
for (i in 3:(ncol(density_plus_drivers))) {
    kruskal_test_result <- kruskal.test(density_plus_drivers[, i] ~ 
                                      density_plus_drivers$month)
    
    kw_pvals[i] <- kruskal_test_result$p.value
}

#convert to dataframe
kw_pvals <- as.data.frame(kw_pvals)

#adjusted p-value and variable name col for 14 tests
kw_pvals <- kw_pvals |> rename(p_val = V1) |> 
  mutate(p_adj = p_val  * 14,
         var = names)

#add var col
kw_pvals$var <- names

#change col names
colnames(kw_pvals) <- c("p_val","p_adj", "var")


