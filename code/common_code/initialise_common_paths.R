#################################################
### Initalise paths common to all experiments ###
#################################################

# main directory containing all code and data for the metallica project
proj_dir <- "/Users/aulakhs/Documents/RalserLab/metallica"
# directory for common functions and graphic params
code_dir <- paste0(proj_dir,"/code")
# directory for all publically available databases used
db_dir <- paste0(proj_dir,"/databases")
# directory for all previously published experimental data
published_dataset_dir <- paste0(proj_dir,"/published_datasets")
# metadata dir
metadata_dir <- paste0(proj_dir,"/metadata")

# experimental data and resuls directories
metdep_hdPCA_dir <- paste0(proj_dir,"/experiment_data/metdep_hdPCA")
metdep_KOgrowth_dir <- paste0(proj_dir,"/experiment_data/metdep_KOgrowth")
metpert_WTgrowth_dir <- paste0(proj_dir,"/experiment_data/metpert_WTgrowth")
metpert_WTmetallomics_dir <- paste0(proj_dir,"/experiment_data/metpert_WTmetallomics")
metpert_WTproteomics_dir <- paste0(proj_dir,"/experiment_data/metpert_WTproteomics")
metpert_ecYeast8simulation_dir <- paste0(proj_dir,"/experiment_data/metpert_Yeast8simulation")
metpert_sim_vs_exp_comparison_dir <- paste0(proj_dir,"/experiment_data/simulation_vs_experiment_comparison")

datasetcomparison_dir <- paste0(proj_dir,"/dataset_comparison")
#
acc_file_dir <-paste0(proj_dir,"/accessory_files")

# Layout directory
lo_dir <- paste0(metadata_dir,"/layouts_96wellplates")

# web app files
metallica_app_dir <- paste0(proj_dir,"/metallica_shiny_app")


