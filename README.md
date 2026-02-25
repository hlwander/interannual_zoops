# interannual_zoops

Code for the analysis of zooplankton community data in Beaverdam Reservoir over seven years using an integrative multivariate framework.

# Repo structure

This repository contains three main folders: 1) Output contains zooplankton and environmental data files used for analysis and figure generation, 2) Scripts contains all code necessary to reproduce figures and analyses for assessing changes in zooplankton community structure over six summer stratified periods, and 3) Inputs contains the meteorological driver file necessary to calculate reservoir discharge (see `TMWB_inflow.R` in Scripts file).

# Instructions to reproduce figures and analyses

1.  Run `01_env.R` to summarize environmental data across months and years that correspond to zooplankton sampling days
2.  Run `02_zoop_dens.R` to summarize zooplankton data and generate zooplankton density figures (Figures 1, S1, and S2)
3.  Run `03_second_stage_NMDS.R` to reproduce the first and second-stage Non-metric MultiDimensional Scaling analysis and figures (Figures 2, 3, S3, S4, S5, and S6; Table S1)
4.  Run `04_zoop_RDA.R` to reprodice redundnacy analyses and figure (Figure 4; Tables 1-4, S2, and S3)
5.  Run `05_indicator_species_analysis.R` to reproduce the indicator species analysis and figure (Figure 5)

Note that `download_NLDAS.R`, `EcoHydRology_functions.R`, and `TMWB_inflow.R` are only necessary to generate the inflow data file used in `01_env.R` and are not necessary to reproduce figures and analyses in this study.
