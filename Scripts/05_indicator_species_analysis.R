# indicator species analysis

#read in packages
pacman::p_load(indicspecies, pheatmap, ggplot2, dplyr, tidyr, 
               igraph, ggraph, viridis)

year_cols <- c("#011f51","#1f78b4","#33a02c","#fdfa66","#ff7f00","#e31a1c","#6a3d9a")
# 2014, 2015, 2016, 2019, 2020, 2021, 2023

#read in zoop data
zoop_dens <- read.csv("Output/zoop_raw_dens.csv", header=TRUE) |> 
  mutate(DateTime = as.Date(DateTime),
         year = format(DateTime, "%Y")) 
year <- zoop_dens$year

zoop_dens_trans <- read.csv("Output/zoop_dens_trans.csv", header=TRUE)

indval <- multipatt(zoop_dens_trans, year, control = how(nperm=999))
summary(indval)
# only bosmina is significantly associated with certain combinations of years
# 2014+2020+2021+2023

#create a table summarizing this
indicators <- list(
  "Bosmina" = c(2014, 2020, 2021, 2023))

zoop_dens_trans <- zoop_dens_trans |> mutate(rowid = row_number()) |> 
                   left_join(zoop_dens |> select(DateTime, year, month) |>
                               mutate(rowid = row_number()), by = "rowid") |>
  mutate(DateTime = as.Date(DateTime))

#list taxa
taxa_cols <- setdiff(names(zoop_dens_trans), c("DateTime","year","month","rowid"))

#convert to long
zoop_dens_trans_long <- zoop_dens_trans |>
  pivot_longer(cols = all_of(taxa_cols),
               names_to = "Taxon",
               values_to = "TransDensity") |>
  mutate(Taxon = as.character(Taxon),
         year = as.integer(year),
         month = month(DateTime),
         day   = day(DateTime),
         pseudoDate = as.Date(sprintf("2000-%02d-%02d", month, day))) |>
  dplyr::select(-rowid) |>
  rowwise() |>
  mutate(isa_alpha = if_else(year %in% indicators[[Taxon]], 1, 0.3)) |>
  filter(Taxon %in% rownames(indval$sign)[!is.na(indval$sign$p.value)]) #only select ind taxa

# create a named vector of facet labels including significant years
facet_labels <- sapply(names(indicators), function(taxon) {
  sig_years <- paste(indicators[[taxon]], collapse = ", ")
  paste0(taxon, " (", sig_years, ")")})
names(facet_labels) <- names(indicators)

#plot!
ggplot(zoop_dens_trans_long, aes(x = pseudoDate, y = TransDensity, 
                                 color = factor(year), group = factor(year))) +
  geom_line(aes(alpha = isa_alpha), size = 1) +
  geom_point(aes(alpha = isa_alpha), size = 1.5) +
  scale_alpha_identity() +
  scale_color_manual(values = year_cols) +
  scale_x_date(limits = as.Date(c("2000-04-01", "2000-11-30")),
                breaks = seq(as.Date("2000-04-01"), as.Date("2000-11-01"), by = "1 month"),
                date_labels = "%b") +
  facet_wrap(~Taxon, scales = "free_y", ncol = 1,
             labeller = labeller(Taxon = facet_labels)) +
  labs(x = "", y = "Hellinger transformed density", color = "") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(nrow = 1))
#ggsave("Figures/ind_sp_dens.jpg", width = 6, height = 5)

