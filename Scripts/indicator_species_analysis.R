# indicator species analysis

#read in packages
pacman::p_load(indicspecies, pheatmap, ggplot2, dplyr, tidyr, 
               igraph, ggraph, viridis)

year_cols <- c("#a13637","#06889b", "#facd60", "#f44034", "#011f51", "#fdfa66")
# 2014, 2019, 2016, 2021, 2015, 2020

#read in zoop data
zoop_dens <- read.csv("Output/zoop_raw_dens.csv", header=TRUE) |> 
  mutate(DateTime = as.Date(DateTime),
         year = format(DateTime, "%Y")) 
year <- zoop_dens$year

zoop_dens_trans <- read.csv("Output/zoop_dens_trans.csv", header=TRUE)

indval <- multipatt(zoop_dens_trans, year, control = how(nperm=999))
summary(indval)
# 3 taxa were significantly associated with certain combinations of years
# Ascomorpha is a strong indicator of early years (2014-2016)    p = 0.001
# Bosmina characterizes 2014 and 2019-2021                       p = 0.022
# Ceriodaphnia is more typical of later years (2016, 2019-2021); p = 0.015

# Manually create a data frame from results
indval_df <- tibble(
  species = c("Ascomorpha", "Bosmina", "Ceriodaphnia"),
  group = c("2014+2015+2016", "2014+2019+2020+2021", "2016+2019+2020+2021"),
  stat = c(0.796, 0.709, 0.699),
  p_value = c(0.001, 0.022, 0.015)
)

# Plot
ggplot(indval_df, aes(x = reorder(species, stat), y = stat, fill = group)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = paste0("p=", p_value)), vjust = -0.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Indicator Species Analysis",
    x = "",
    y = "Indicator Value (IndVal)"
  ) +
  theme_minimal()


# Species and the years they are associated with
species <- c("Ascomorpha", "Bosmina", "Ceriodaphnia")
edges_list <- list(
  Ascomorpha = c("2014","2015","2016"),
  Bosmina = c("2014","2019","2020","2021"),
  Ceriodaphnia = c("2016","2019","2020","2021")
)

# Flatten into edges data frame
edges <- do.call(rbind, lapply(names(edges_list), function(sp) {
  data.frame(from = sp, to = edges_list[[sp]], species = sp)
}))

# Nodes: combine all unique species and years
nodes <- data.frame(
  name = c(species, unique(unlist(edges_list))),
  type = c(rep("Species", length(species)), rep("Year", length(unique(unlist(edges_list)))))
)

g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

species_colors <- viridis(length(species), option = "plasma")
names(species_colors) <- species

year_cols <- viridis(6, option = "D")  # 6 years
names(year_cols) <- unique(nodes$name[nodes$type == "Year"])

ggraph(g, layout = "fr") +
  geom_edge_link(aes(color = species), width = 1.5, alpha = 0.8) +
  geom_node_point(aes(filter = type == "Year", fill = name), 
                  shape = 21, size = 8, color = "black", show.legend=FALSE) +
  geom_node_text(aes(label = name, filter = type == "Year"), repel = TRUE, size = 4, fontface = "bold") +
  scale_edge_color_manual(values = species_colors) +
  scale_fill_manual(values = year_cols, name = "Year") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_blank())
#ggsave("Figures/zoop_ISA.jpg", width=5, height=4)
