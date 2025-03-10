########################################################################################################
####          Script to plot the growth of vacuolar KO mutants across all metal depletion media     ####
########################################################################################################
#`---
#`  Title: "metdepKO_growth_vacuolarmutants.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 3 Oct 2024
#`  Description: extract all vacuolar mutants from growth data of KO mutants
#`---


#############################################
### source paths functions and libraries  ###
#############################################

source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")

# general

source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/database_identifier_conversion_functions.R"))
source(paste0(code_dir,"/common_code/input_processed_databases_publisheddatasets.R"))

# specific paths and scripts


plot_dir <- paste0(metdep_KOgrowth_dir,"/output/plots/")
dir.create(plot_dir,recursive = T)

output_tables_dir <- paste0(metdep_KOgrowth_dir,"/output/tables/")
dir.create(output_tables_dir,recursive = T)
