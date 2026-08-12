# taxa-specific analyses (ISA and SIMPER)

#read in packages
pacman::p_load(indicspecies, pheatmap, ggplot2, dplyr, tidyr, tidytext, 
               lubridate, NatParksPalettes, cols4all, patchwork,
               igraph, ggraph, viridis, vegan, purrr, tibble, stringr)

#ISA analysis

#read in zoop data
zoop_dens <- read.csv("Output/zoop_raw_dens.csv", header=TRUE) |> 
  mutate(DateTime = as.Date(DateTime),
         year = format(DateTime, "%Y"),
         month = format(DateTime, "%b")) 
year <- zoop_dens$year

zoop_dens_trans <- read.csv("Output/zoop_dens_trans.csv", header=TRUE)

zoop_dens_trans <- zoop_dens_trans |> mutate(rowid = row_number()) |> 
  left_join(zoop_dens |> select(DateTime, year, month) |>
              mutate(rowid = row_number()), by = "rowid") |>
  mutate(DateTime = as.Date(DateTime)) 

# community matrix only — no metadata columns
comm_matrix <- zoop_dens_trans |>
  dplyr::select(-rowid, -DateTime, -month, -year)

indval <- multipatt(comm_matrix, year, control = how(nperm=999))
summary(indval)

#create a table summarizing ISA results
indicators <- list(
  "Chydorus" = c(2014),
  "Synchaeta" = c(2019, 2020),
  "Ascomorpha" = c(2014, 2015, 2016),
  "Conochiloides" = c(2020, 2021, 2023),
  "Bosmina" = c(2014, 2020, 2021, 2023),
  "Gastropus" = c(2019, 2020, 2021, 2023),
  "Trichocerca" = c(2014, 2019, 2020, 2021, 2023),
  "Ceriodaphnia" = c(2016, 2019, 2020, 2021, 2023))

#-----------------------------------------------------------------------------#
# SIMPER analysis

comm_matrix_raw <- zoop_dens |>
  dplyr::select(-DateTime, -year, -month) 

# grouping factor
grp <- factor(year)

# run SIMPER
set.seed(123)
sim <- simper(comm_matrix_raw, group = grp, permutations = 999)

# inspect all pairwise comparisons
summary(sim)

# Tidy SIMPER output into one table
simper_tbl <- bind_rows(
  lapply(names(sim), function(comp){
    x <- as.data.frame(sim[[comp]])
    x$Taxon <- rownames(x)
    x$Comparison <- comp
    x})) |>
  as_tibble() |>
  relocate(Comparison, Taxon)

# Top taxa for each year-pair comparison
top_simper <- simper_tbl |>
  group_by(Comparison) |>
  arrange(desc(average), .by_group = TRUE) |>
  dplyr::slice_head(n = 11) |>
  ungroup()

#------------------------------------------------------------------------------#
# Joint ISA + SIMPER figure

# ── colour map ──────────────────────────
taxa_cols <- c(
  "Ascomorpha"    = "#212E52",
  "Bosmina"       = "#444E7E",
  "Ceriodaphnia"  = "#B7ABBC",
  "Chydorus"      = "#F9ECE8",
  "Conochiloides" = "#FCC893",
  "Diaphanosoma"  = "grey95", 
  "Gastropus"     = "#FEB424",
  "Synchaeta"     = "grey95",
  "Trichocerca"   = "#D8511D",
  "Polyarthra"    = "grey95",
  "Nauplii"       = "grey95",
  "Keratella"     = "grey95",
  "Kellicottia"  = "grey95",
  "Daphnia"       = "grey95",
  "Cyclopoida"    = "grey95",
  "Calanoida"     = "grey95",
  "Conochilus"    = "grey95",
  "Asplanchna"    = "grey95") 

taxa_order <- names(taxa_cols) 

# ── PANEL A: ISA dot plot ─────────────────────────────────────────────────────
# Build a long table: one row per taxon × year, flagged as associated or not

all_years <- c(2014, 2015, 2016, 2019, 2020, 2021, 2023)

overlap_taxa <- intersect(names(indicators), unique(top_simper$Taxon))

isa_long <- tibble(
  Taxon = rep(names(indicators), lengths(indicators)),
  Year  = as.integer(unlist(indicators))) |>
  mutate(Associated = TRUE) |>
  right_join(
    expand.grid(Taxon = names(indicators), Year = all_years,
                stringsAsFactors = FALSE),
    by = c("Taxon", "Year")) |>
  replace_na(list(Associated = FALSE)) |>
  mutate(Year  = factor(Year),
         colour_group = ifelse(Taxon %in% overlap_taxa, Taxon, "Other"))

# pull p-values from indval for label annotation
pvals <- indval$sign |>
  as.data.frame() |>
  rownames_to_column("Taxon") |>
  filter(!is.na(p.value)) |>
  select(Taxon, stat, p.value)

isa_long <- isa_long |>
  left_join(pvals, by = "Taxon") 

pA <- ggplot(isa_long, aes(x = Year, y = Taxon)) +
  geom_tile(fill = "white", colour = "grey95", linewidth = 0.4) +
  geom_point(data = filter(isa_long, Associated), aes(fill = colour_group),
               shape = 21, size = 5, colour = "black", stroke = 0.3) +
  scale_fill_manual(values = c(taxa_cols, Other = "grey95"), guide = "none") +
  geom_point(data = filter(isa_long, !Associated), shape = 21, 
                size = 5, fill = "white", color = "grey", stroke = 0.5) +
  scale_colour_manual(values = c(taxa_cols), guide = "none") +
  scale_x_discrete(expand = expansion(add = c(0.5, 1.2))) +
  labs(x = "", y = "", title = "A  Indicator Species Analysis",
       subtitle = "Filled = significant indicator year(s)") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "bold"),
    plot.title    = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 8, colour = "grey50"))

# ── PANEL B: SIMPER — mean contribution across all pairwise comparisons ───────
# Aggregate to avoid showing all 21 pairwise bars

simper_mean <- top_simper |>
  group_by(Taxon) |>
  summarise(
    mean_contrib = mean(average, na.rm = TRUE),
    se_contrib   = sd(average, na.rm = TRUE) / sqrt(n()),
    .groups = "drop") 

pB <- ggplot(simper_mean, aes(x = mean_contrib, y = Taxon, fill = Taxon)) +
  geom_col(width = 0.65) +
  geom_errorbar(aes(xmin = mean_contrib - se_contrib,
                    xmax = mean_contrib + se_contrib),
                width = 0.25, colour = "grey40", linewidth = 0.4) +
  scale_fill_manual(values = c(taxa_cols, Other = "grey95"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15)),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Mean contribution to Bray–Curtis dissimilarity",
       y = "", title = "B  SIMPER",
       subtitle = "Mean ± SE across all pairwise year comparisons  ★ = ISA indicator") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y = element_text(face = "bold"),
    plot.title         = element_text(size = 10, face = "bold"),
    plot.subtitle      = element_text(size = 8, colour = "grey50"))

# Combine with patchwork (Figure 5)
joint_fig <- pA + pB +
  plot_layout(widths = c(1.4, 1)) &
  theme(plot.margin = margin(6, 10, 6, 6))
#ggsave("Figures/isa_simper_joint.jpg", joint_fig, width = 8, height = 4.5, dpi = 300)

#------------------------------------------------------------------------------
#ISA Table S9
isa_table <- indval$sign |>
  as.data.frame() |>
  rownames_to_column("Taxon") |>
  filter(!is.na(p.value)) |>
  mutate(
    `Associated years` = sapply(Taxon, function(t) paste(indicators[[t]], collapse = ", ")),
    IndVal = round(stat, 3),
    A      = round(mapply(function(tax, idx) indval$A[tax, idx], Taxon, index), 3),
    B      = round(mapply(function(tax, idx) indval$B[tax, idx], Taxon, index), 3),
    p      = round(p.value, 3)) |>
  filter(A > 0.75 | B > 0.75,
         p < 0.05) |> 
  select(Taxon, `Associated years`, IndVal, A, B, p) |>
  arrange(match(Taxon, taxa_order))
#write.csv(isa_table, "Output/isa_table.csv", row.names = FALSE)
