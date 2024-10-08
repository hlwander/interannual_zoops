# interannual_zoops
Code for the analysis of monthly summer zooplankton in Beaverdam Reservoir over six years

# Repo structure
This repository contains three main folders: 1) Output contains zooplankton and environmental data files used for analysis and figure generation, 2) Scripts contains all code necessary to reproduce figures and analyses for assessing changes in zooplankton community structure over six summer stratified periods, and 3) Inputs contains the meteorological driver file necessary to calculate reservoir discharge (see TMWB_inflow.R in Scripts file).

# Instructions to reproduce figures and analyses
1. Run 01_env.R to summarize environemtnal data across months and years that correspond to zooplankton sampling days
2. Run 02_zoop_dens.R to summarize zooplankton data and generate zooplankton density figures
3. Run 03_second_stage_NMDS.R to reproduce the Non-metric MultiDimensional Scaling analysis
4. Run 04_nmds_drivers.R to visualize differences in environmental drivers across years

Note that download_NLDAS.R, EcoHydRology_functions.R, and TMWB_inflow.R are only necessary to generate the inflow data file used in the 01_env.R and are not necessary to reproduce figures and analyses for this study.
