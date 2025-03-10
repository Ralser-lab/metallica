
#`---
#`  Title: "vacuolar_mutants_growth_Zndepletion.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 18 September 2024
#`  Description: check whether vacuolar mutants grow worse in Zn depletion
#`---


#############################################
### source paths functions and libraries  ###
#############################################


# general
source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,'/common_code/graphics_parameters.R'))

source(paste0(code_dir,'/common_code/database_identifier_conversion_functions.R'))
#source(paste0(code_dir,'/common_code/input_processed_databases_publisheddatasets.R'))

plot_dir = paste0( metdep_KOgrowth_dir, "/output/plots")
dir.create(plot_dir,recursive = T)
output_tables_dir = paste0( metdep_KOgrowth_dir, "/output/tables")
dir.create(output_tables_dir,recursive = T)

#################################################################
### Get metdepKOgrowth specific functions and libraries from  ###
#################################################################

source(paste0(code_dir,'/metdep_KOgrowth/metdepKOgrowth_0_libraries_functions.R'))


## GO CC data

GO_CC_vacuole <- filter(GO_gset_CC,
                        grepl("vacuol",term))%>%
                 unique()

## KO data

KOgrowth_data_Zn_vac <- read.csv(paste0(output_tables_dir,"/metdepKOgrowth_statresults_annotated.csv"),
                                 stringsAsFactors = F)%>%
                        filter(grepl("Zn",condition))%>%
                        mutate(vacuolar_mutant = ifelse(ORF %in% unique(GO_CC_vacuole$ORF)))
                                              

KOgrowth_data_vac <- read.csv(paste0(output_tables_dir,"/metdepKOgrowth_statresults_annotated.csv"),
                              stringsAsFactors = F)%>%
                      filter(ORF %in% unique(GO_CC_vacuole$ORF))


ggplot(KOgrowth_data_vac,
       aes(x = mean_effect_size_log2,
           y = -log10(p_Welch_BH),
           colour = metal))+
  geom_point(aes(alpha = sig_in_metal))+
  theme_metallica()

KOgrowth_data_vac_sig <- KOgrowth_data_vac%>%
                         filter(p_Welch_BH < 0.05 & mean_effect_size_log2 < -log(1.2))%>%
                         dplyr::select(ORF, metal)%>%
                         unique()


KOgrowth_data_vac_sig_smry <- KOgrowth_data_vac_sig%>%
                              group_by(metal)%>%
                              summarize(count = n())


