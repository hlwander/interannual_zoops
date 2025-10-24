# playing around with additional analyses

pacman::p_load(dplyr, vegan, ggplot2, akima, lubridate, purr, viridis, mgcv)

# create wide zoop df: one row per year, columns for each taxon
zoop_wide <- all_zoop_taxa |>
  #filter(month %in% c(5,6,7,8,9)) |>
  group_by(DateTime, Taxon) |>
  summarise(Density = sum(dens), .groups = "drop") |>
  tidyr::pivot_wider(names_from = Taxon, values_from = Density, values_fill = 0) |>
  mutate(month = month(DateTime),
         year = year(DateTime)) 

#now add in simpson col
zoop_wide$simpson_div <- 1 - vegan::diversity(zoop_wide[,-c(1,12,13)], index = "simpson")
zoop_wide$shannon_div <- diversity(zoop_wide[,-c(1,12,13)], index = "shannon")

ggplot(zoop_wide, aes(x = DateTime, y = simpson_div)) +
  geom_line() + geom_point(aes(color = as.factor(month)), size = 3) + xlab("") +
  ylab("Simpson's Diversity Index (1 - D)") +
  theme_minimal() # highest in 2019, lowest in 2014 when summarizing only by year (n=6 popints)
#some really interesting 

#look across months and color by year
ggplot(zoop_wide, aes(x = as.factor(month), y = simpson_div, group=year, 
                      color=as.factor(year))) + 
  geom_line() + geom_point(aes(color = as.factor(year)), size = 3) + 
  ylab("Simpson's Diversity Index (1 - D)") + xlab("") +
  theme_minimal()

#try barplots with SE
zoop_summary <- zoop_wide |>
  filter(month %in% c(5,6,7,8,9)) |>
  group_by(year = as.factor(year(DateTime)),
           month = factor(month(DateTime, label = TRUE, abbr = TRUE),
                          levels = c("May", "Jun", "Jul", "Aug", "Sep"))) |>
  summarise(
    mean_simpson = mean(simpson_div, na.rm = TRUE),
    se_simpson = sd(simpson_div, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

ggplot(zoop_summary, aes(x = month, y = mean_simpson, fill = year)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_simpson - se_simpson,
                    ymax = mean_simpson + se_simpson),
                position = position_dodge(width = 0.8),
                width = 0.2) +
  ylab("Simpson's Diversity Index (1 - D)") +
  xlab("") + theme_minimal() +
  scale_fill_brewer(palette = "Set1", name = "Year")

#and violin plots
ggplot(zoop_wide, aes(x = DateTime, y = simpson_div, fill = as.factor(year))) +
  geom_violin(position = position_dodge(width = 0.8), width = 0.7, trim = FALSE) +
  stat_summary(fun = mean, geom = "point", 
               position = position_dodge(width = 0.8), 
               shape = 21, size = 2, color = "black", fill = "white") +
  ylab("Simpson's Diversity Index (1 - D)") +
  xlab("") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1", name = "Year")

#----------------------------------------------------------------------------#
# temperature heat maps for each year w/ interpolation between profiles and depths
yr <- unique(year(temp_final$DateTime))
interp_list <- vector("list", length(yr))  # create empty list to store results

for(i in seq_along(yr)) {
  # Filter for the year
  df <- temp_final %>%
    filter(year(DateTime) == yr[i]) %>%
    mutate(time_num = as.numeric(DateTime)) %>%
    group_by(time_num, Depth_m) %>%
    summarise(Temp_C = mean(Temp_C, na.rm = TRUE), .groups = "drop")
  
  # Interpolation grids for May 1 to Sept 30
  start_date <- as.Date(paste0(yr[i], "-05-01"))
  end_date <- as.Date(paste0(yr[i], "-09-30"))
  
  time_seq <- seq(
    from = as.numeric(start_date),
    to = as.numeric(end_date),
    length.out = as.numeric(end_date - start_date) + 1
  )
  
  depth_seq <- seq(
    min(df$Depth_m),
    max(df$Depth_m),
    length.out = 100
  )
  
  # Perform interpolation
  interp_result <- akima::interp(
    x = df$time_num,
    y = df$Depth_m,
    z = df$Temp_C,
    xo = time_seq,
    yo = depth_seq,
    duplicate = "mean"
  )
  
  # Create dataframe from interpolation result
  interp_df <- expand.grid(
    time_num = interp_result$x,
    Depth_m = interp_result$y
  )
  interp_df$Temp_C <- as.vector(interp_result$z)
  interp_df$DateTime_orig <- as.Date(interp_df$time_num, origin = "1970-01-01")
  
  # Convert to standardized year (2000) for aligned x-axis
  interp_df$DateTime <- as.Date(
    paste0(
      "2000-",
      sprintf("%02d", month(interp_df$DateTime_orig)),
      "-",
      sprintf("%02d", day(interp_df$DateTime_orig))
    )
  )
  
  interp_df$Year <- yr[i]
  
  # Save to list
  interp_list[[i]] <- interp_df
}

# Combine all years into one dataframe
interp_all <- bind_rows(interp_list)

interp_all$Year <- as.factor(interp_all$Year)

ggplot(interp_all, aes(x = DateTime, y = Depth_m, fill = Temp_C)) +
  geom_raster() +
  scale_y_reverse(limits = c(max(interp_all$Depth_m), 0)) +
  facet_wrap(~Year, ncol = 3, scales = "fixed") +
  #scale_fill_gradient(low = "#66CCFF", high = "#990000", na.value = "grey50") +
  scale_fill_gradientn(
    colors = c("#006699", "#CCCCCC", "#990000"),
    values = scales::rescale(c(10, 20, 30)),  
    na.value = "grey50",
    name = "Temp (°C)"
  ) +
  labs(
    x = "",
    y = "Depth (m)",
    fill = "Temperature (°C)",
  ) +
  theme_minimal() +
  scale_x_date(date_labels = "%b", date_breaks = "1 month")

#-----------------------------------------------------------------------------#
# General Additive Model using the 80 tows

#make sure rows are ordered by date
zoop_wide <- zoop_wide[order(zoop_wide$DateTime), ]

# Create a numeric index
zoop_wide$time_index <- seq_len(nrow(zoop_wide))

gam_mod <- gam(simpson_div ~ s(time_index), data = zoop_wide)
summary(gam_mod) 
# this model explains <1% of variation in shannon and simpson diversity (deviance)

#-----------------------------------------------------------------------------#
# ANOVA - not significant
summary(aov(shannon_div ~ factor(year(DateTime)), data = zoop_wide))

ggplot(zoop_wide, aes(x = DateTime, y = simpson_div)) +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE)


# Redundancy analysis
rda_model <- rda(zoop_dens_trans ~ Year + Temp + TP + DO, data = env_drivers)
#need to create a new env df with datetimme not just month and year

#-----------------------------------------------------------------------------#
# PERMANOVA

# Extract community data (columns for zooplankton taxa only)
zoop_comm <- zoop_wide |>
  select(Ascomorpha:Polyarthra) 

# Extract year as grouping variable
zoop_meta <- zoop_wide |>
  select(year)

# Run PERMANOVA
permanova_result <- adonis2(zoop_comm ~ year,
                            data = zoop_meta,
                            method = "bray",   # Bray-Curtis dissimilarity is standard
                            permutations = 999)

print(permanova_result)

#test assumptions of homogeneity
bray_dist <- vegdist(zoop_comm, method = "bray")
group_dispersion <- betadisper(bray_dist, zoop_meta$year)
anova(group_dispersion)  # dispersion does differ among years so permanova may not be appropriate

#visualize dispersion
plot(group_dispersion)
boxplot(group_dispersion) # 2015 and 2016 seem like the outliers

# PERMDISP using hellinger transformed data
permanova_result <- adonis2(zoop_dens_trans ~ year, data = all_zoops_nmds, method = "euclidean")
bray_dist <- vegdist(zoop_dens_trans, method = "euclidean")
group_dispersion <- betadisper(bray_dist, all_zoops_nmds$year)
anova(group_dispersion)



# Test dispersion assumption again
disp_hel <- betadisper(vegdist(zoop_dens_trans, method = "euclidean"), year(all_zoops_nmds$DateTime))
anova(disp_hel) #dispersion is homogenous so can run PERMANOVA p>0.05

# Run PERMANOVA with Euclidean distance on Hellinger-transformed data
adonis2(zoop_dens_trans ~ year(DateTime), data = all_zoops_nmds, method = "euclidean", permutations = 999)
# valid even with unequal sample size bc of homogenous dispersion
# Zooplankton community composition significantly differs across year (p=0.001) 
# year accounts for 5.8% of the variation in zooplankton community structure
# n = 81 samples from the full year (not just may-sep which would be 61)

# Custom pairwise PERMANOVA using adonis2
library(rstatix)

# Extract grouping variable
grouping_var <- as.factor(all_zoops_nmds$year)
levels <- levels(grouping_var)

# Prepare results container
results <- data.frame()

# Loop over all pairwise combinations
for(i in 1:(length(levels)-1)){
  for(j in (i+1):length(levels)){
    
    # Subset data for the two groups
    group1 <- levels[i]
    group2 <- levels[j]
    
    idx <- grouping_var %in% c(group1, group2)
    data_sub <- zoop_dens_trans[idx, ]
    group_sub <- droplevels(grouping_var[idx])
    
    # Run adonis2
    a <- adonis2(data_sub ~ group_sub, method = "euclidean", permutations = 999)
    
    # Store results
    results <- rbind(results, data.frame(
      group1 = group1,
      group2 = group2,
      F = a$F[1],
      R2 = a$R2[1],
      p = a$`Pr(>F)`[1]
    ))
  }
}

# Adjust p-values for multiple comparisons
results <- results %>%
  mutate(p_adjusted = p.adjust(p, method = "BH"))

print(results)

pairwise_results <- results %>%
  mutate(
    p_category = case_when(
      p_adjusted < 0.05 ~ "< 0.05",
      p_adjusted == 0.05 ~ "= 0.05",
      p_adjusted > 0.05 ~ "> 0.05"
    ),
    pair = paste(group1, group2, sep = "-")
  )

ggplot(pairwise_results, aes(x = reorder(pair, p_adjusted), y = p_adjusted, fill = p_category)) +
  geom_col() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  scale_fill_manual(
    values = c(
      "< 0.05" = "red",
      "= 0.05" = "orange",
      "> 0.05" = "gray"
    )
  ) +
  labs(
    title = "Pairwise PERMANOVA Adjusted p-values",
    x = "Year Pair",
    y = "Adjusted p-value",
    fill = "p-value category"
  ) +
  coord_flip() +
  theme_minimal()
# the only 4 year pairs that are NOT significant are 2014-2015, 2014-2016, 2015-2016, and 2020-2021

# SIMPER analysis to see which taxa drive differences across years
comm <- zoop_dens_trans  # taxa abundance matrix
group <- all_zoops_nmds$year  # this should be a vector/factor, not a data frame column

# Make sure it's a factor
group <- as.factor(group)

# Run SIMPER
simper_results <- simper(comm, group, permutations = 999)
summary(simper_results)

#new df with simper results
simper_full_df <- imap_dfr(simper_results, ~{
  pair_name <- .y
  simper_pair <- .x
  
  # Get taxa names
  taxa <- simper_pair$species
  
  # Create tibble with stats extracted by taxa
  tibble(
    year_pair = pair_name,
    taxon = taxa,
    average = simper_pair$average[taxa],
    sd = simper_pair$sd[taxa],
    #ratio = simper_pair$ratio[taxa],
    #ava = simper_pair$ava[taxa],
    #avb = simper_pair$avb[taxa],
    #cumsum = simper_pair$cusum[taxa],
    p = simper_pair$p[taxa]
  )
})

# View results
print(simper_full_df)

simper_full_df %>%
  group_by(year_pair) %>%
  slice_max(order_by = average, n = 5) %>% # top 5 taxa per pair
  ungroup() %>%
  mutate(signif = ifelse(p < 0.05, "Significant", "Not significant")) %>%
  ggplot(aes(x = reorder(taxon, average), y = average, fill = signif)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ year_pair, scales = "free_y") +
  labs(title = "Top Taxa Contributions by Year Pair (SIMPER)",
       x = "Taxon",
       y = "Average Contribution") +
  theme_minimal() +
  scale_fill_manual(values = c("Significant" = "firebrick", "Not significant" = "gray70"))
#ascomorpha sig in 5/15 comps, kellicottia in 4/15, keratella 2/15, bosmina 3/15, 
#cyclopoid 2/15, nauplii and ceriodaphnia 1 each / 15

# heatmap to show contribution of each taxon
simper_full_df %>%
  mutate(signif = ifelse(p < 0.05, "*", "")) %>%
  ggplot(aes(x = year_pair, y = taxon, fill = average)) +
  geom_tile() +
  geom_text(aes(label = signif), color = "white", size = 5) +
  scale_fill_viridis_c(option = "plasma") +
  labs(title = "Heatmap of Taxa Contributions Across Year Pairs",
       x = "Year Pair",
       y = "Taxon",
       fill = "Avg Contribution") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# peak timing extraction

peak_months <- zoops_10_groups %>%
  group_by(year, Taxon) %>%
  filter(avg == max(avg)) %>%
  summarise(peak_month = first(month),
            trough_month = last(month))


# general additive model
library(mgcv)

gam_mod <- gam(as.numeric(peak_month) ~ year + Taxon, data = peak_months)
summary(gam_mod)
# no year or taxon have significant effect on peak month timing

#plot
ggplot(peak_months, aes(x = year, y = as.numeric(peak_month), color = Taxon)) +
  geom_point() + geom_line(aes(group=Taxon)) +
  theme_minimal()

# calculate spread between max and min each year
#peak_range <- zoop_data %>%
#  group_by(year, taxon) %>%
#  summarise(
#    peak_month = month[which.max(density)],
#    trough_month = month[which.min(density)],
#    range_months = abs(peak_month - trough_month)
#  )

# coverage figure of zoop data from 2014-2024

new_samples <- read.csv("/Users/heatherwander/Documents/VirginiaTech/research/zooplankton/Summer2022-2025-DataAnalysis/RawData/BVR2022-2025_Zooplanktoncounting_Density_DataEntry.csv")
new_samples <- unique(new_samples$collect_date)
new_samples_date <- as.Date(new_samples, format = "%d-%b-%y")

#list all dates
all_zoop_dates <- c(all_zoops_nmds$DateTime)

#create df for plotting
coverage_df <- data.frame(SampleDate = all_zoop_dates) |>
  mutate(SampleMonth = floor_date(SampleDate, "month"),
         month = month(SampleDate))

ggplot(coverage_df, aes(x = SampleMonth, fill = as.factor(month))) +
  geom_bar(color = "black") +  
  scale_fill_viridis_d(option = "C") +
  labs(
    title = "Monthly Zooplankton Sampling Coverage",
    x = "Date",
    y = "Number of Samples"
  ) +
  theme_minimal() 
#ggsave("Figures/zoop_data_coverage_barplot.jpg", width=7, height=4) 

ggplot(coverage_df, aes(x = SampleDate, y = 1)) +
  geom_point(size = 2, alpha = 0.8, color = "darkblue") +
  scale_y_continuous(breaks = NULL) +
  labs(
    title = "Zooplankton Sampling Date Coverage",
    x = "Date",
    y = NULL
  ) +
  theme_minimal()
#ggsave("Figures/zoop_data_coverage_point.jpg", width=7, height=4) 

# heat map
coverage_df <- data.frame(SampleDate = all_zoop_dates) %>%
  mutate(Year = year(SampleDate),
         Month = month(SampleDate, label = TRUE)) %>%
  count(Year, Month) |>
  filter(!Year %in% "2025")

#Fig SX
ggplot(coverage_df, aes(x = Month, y = factor(Year))) +
  geom_tile(aes(fill = n), color = "white") +
  scale_fill_gradient(low = "orange", high = "darkblue") +
  labs(x = "Month", y = "Year", fill = "Sample Count") +
  theme_minimal()
#ggsave("Figures/zoop_data_coverage_heatmap.jpg", width=7, height=4) 
