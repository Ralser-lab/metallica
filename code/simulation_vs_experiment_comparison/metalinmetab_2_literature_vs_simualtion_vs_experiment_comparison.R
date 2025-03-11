#`---
#`  Title: "analyze_enzyme_changes_per_metal_v2.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 6 June 2021 
#`  Description: Script to analyze enzyme category wise changes in protein abundance + other datasets
#`---


library(tidyr)
library(RColorBrewer)
library(plotly)
library(corrplot)
library(readr)
#############################################
### source paths functions and libraries  ###
#############################################

source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")

# general

source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/database_identifier_conversion_functions.R"))
source(paste0(code_dir,"/common_code/input_processed_databases_publisheddatasets.R"))


## specific
source(paste0(code_dir,"/metpert_ecYeast8simulation/metpert_ecYeast8simulation_0_libraries_functions.R"))


plot_dir <- paste0(metpert_sim_vs_exp_comparison_dir,"/output/plots")
dir.create(plot_dir,recursive = T)

output_tables_dir <- paste0(metpert_sim_vs_exp_comparison_dir,"/output/tables")
dir.create(output_tables_dir,recursive = T)

####################################################
### Overlap between metabolism related databases ###
####################################################

metabolism_related_dbs <- list(
  ## genome scale metabolic model
  unique(Sc_GEM$ORF),
  # Enzyme Classificatio numbers
  unique(EC_df$ORF),
  # metal binding or transport
  #uGO_MF_all_metal_related_anno
  # KEGG
  unique(KEGG$ORF))

names(metabolism_related_dbs) <- c( "ScGEM",
                                    "EC",
                                    #   "metal_b+t",
                                    "KEGG")

pdf(paste0(plot_dir,"/venn_ScGEM_EC_KEGG_overlap.pdf"), width = 10,height = 10)
venn::venn(metabolism_related_dbs,zcolor = "style",box = F, ilcs=2.5) 
dev.off()


pdf(paste0(plot_dir,"/metabolism_related_datasets_overlap.pdf"),width = 9, height = 6)

upset(fromList(metabolism_related_dbs), 
      keep.order = T,
      order.by = "freq",
      cutoff = 100,
      nintersects = NA,
      main.bar.color ="#CD919E"  ,
      sets.bar.color = "#BCD2EE",
      text.scale = 2,
      point.size = 3,
      matrix.color = "#2b0b57",
      mb.ratio = c(0.5, 0.5)
      
)
dev.off()

########################################################################################
### Fraction of each metabolism related database that has a metal binding annotation ###
########################################################################################

## scGEM

ScGEM_metalfrac <- merge(Sc_GEM ,GO_MF_all_metal_related_anno, by = "ORF", all.x = T)%>%
  dplyr::select(-term.x)%>%
  unique()%>%
  group_by(ORF)%>%
  mutate(metal_anno = ifelse(!is.na(term.y),1,0))%>%
  ungroup()%>%
  dplyr::select(-term.y)%>%
  unique()%>%
  summarise(total = length(ORF),
            num_metanno = sum(metal_anno),
            frac_metanno = num_metanno/total)%>%
  mutate(database = "ScGEM")

## KEGG

KEGG_metalfrac <- merge(KEGG , GO_MF_all_metal_related_anno, by = "ORF", all.x = T)%>%
  dplyr::select(-term.x)%>%
  unique()%>%
  group_by(ORF)%>%
  mutate(metal_anno = ifelse(!is.na(term.y),1,0))%>%
  ungroup()%>%
  dplyr::select(-term.y)%>%
  unique()%>%
  summarise(total = length(ORF),
            num_metanno = sum(metal_anno),
            frac_metanno = num_metanno/total)%>%
  mutate(database = "KEGG")

KEGGnum_metalfrac <- merge(KEGG , GO_MF_all_metal_related_anno, by = "ORF", all.x = T)%>%
  unique()%>%
  group_by(term.x)%>%
  mutate(metal_anno = ifelse(any(!is.na(term.y)),1,0))%>%
  ungroup()%>%
  dplyr::select(metal_anno,term.x)%>%
  unique()%>%
  summarise(total = length(term.x),
            num_metanno = sum(metal_anno),
            frac_metanno = num_metanno/total)%>%
  mutate(database = "KEGGnum")


## EC_df

EC_metalfrac <- merge(EC_df[,c("ORF","EC_level1_name")] , GO_MF_all_metal_related_anno, by = "ORF", all.x = T)%>%
  dplyr::select(-EC_level1_name)%>%
  unique()%>%
  group_by(ORF)%>%
  mutate(metal_anno = ifelse(!is.na(term),1,0))%>%
  ungroup()%>%
  dplyr::select(-term)%>%
  unique()%>%
  summarise(total = length(ORF),
            num_metanno = sum(metal_anno),
            frac_metanno = num_metanno/total)%>%
  mutate(database = "EC")


ECnum_metalfrac <- merge(EC_df[,c("ORF","EC_level1_name","EC_level_2")] , GO_MF_all_metal_related_anno, by = "ORF", all.x = T)%>%
  unique()%>%
  group_by(EC_level1_name,EC_level_2)%>%
  mutate(metal_anno = ifelse(any(!is.na(term)),1,0))%>%
  ungroup()%>%
  dplyr::select(metal_anno,EC_level1_name,EC_level_2)%>%
  unique()%>%
  summarise(total = length(EC_level1_name),
            num_metanno = sum(metal_anno),
            frac_metanno = num_metanno/total)%>%
  mutate(database = "ECnum_level_1_2")



metabolism_metal_frac <- rbind(ScGEM_metalfrac, KEGG_metalfrac,EC_metalfrac)

pdf(paste0(plot_dir,"/percentage_metdatabaseORFs_metalrelated.pdf"),width = 7,height = 6)
ggplot(metabolism_metal_frac,
       aes(x = reorder(database,-frac_metanno)))+
  geom_bar(stat = "identity",
           aes(y = total, color = database), width = 0.5,fill = NA, linewidth = 0.5)+
  geom_bar(stat = "identity",
           aes(y = num_metanno, color = database, fill = database), width = 0.5,)+
  geom_text(aes(y = num_metanno + 50, label = paste0(round(frac_metanno,2)*100,"%")), size =4 )+
  theme_metallica()+
  scale_color_brewer(palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  labs(fill = "", color = "", y = "ORFs")
dev.off()

metabolism_metal_frac_grouped <- rbind(KEGGnum_metalfrac,ECnum_metalfrac)

pdf(paste0(plot_dir,"/percentage_metdatabaseORFs_metalrelated_pathwaywise.pdf"),width = 4,height = 6)

ggplot(metabolism_metal_frac_grouped,
       aes(x = reorder(database,-frac_metanno)))+
  geom_bar(stat = "identity",
           aes(y = total, color = database), width = 0.5,fill = NA, linewidth = 0.5)+
  geom_bar(stat = "identity",
           aes(y = num_metanno, color = database, fill = database), width = 0.5,)+
  geom_text(aes(y = num_metanno + 2, label = paste0(round(frac_metanno,2)*100,"%")),size = 4)+
  theme_metallica()+
  scale_color_brewer(palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  labs(fill = "", color = "")+
  theme(legend.position = "none")
dev.off()


EC_level1_2_metalfrac <- merge(EC_df[,c("ORF","EC_level1_name","EC_level_2")] , GO_MF_all_metal_related_anno, by = "ORF", all.x = T)%>%
  unique()%>%
  group_by(ORF,EC_level1_name, EC_level_2)%>%
  mutate(metal_anno = ifelse(!is.na(term),1,0))%>%
  ungroup()%>%
  dplyr::select(-term)%>%
  unique()%>%
  group_by(EC_level1_name,EC_level_2)%>%
  summarise(total_ORFs = length(ORF),
            num_metanno = sum(metal_anno),
            frac_metanno = num_metanno/total_ORFs)%>%
  ungroup()


pdf(paste0(plot_dir,"/EC_level_1_2_metaldependence.pdf"),width = 6,height = 11)

ggplot(EC_level1_2_metalfrac,
       aes(x = reorder(paste(EC_level1_name,EC_level_2),frac_metanno)))+
  geom_bar(stat = "identity",
           aes(y = total_ORFs), width = 0.5,fill = NA, linewidth = 0.5, color = brewer.pal(3, name = "Set2")[[1]])+
  geom_bar(stat = "identity",
           aes(y = num_metanno), width = 0.5, fill = brewer.pal(3, name = "Set2")[[1]])+
  geom_text(size = 3,aes(y = total_ORFs + 15, label = paste0(round(frac_metanno,2)*100,"%")))+
  theme_metallica()+
  scale_color_brewer(palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  labs(fill = "", color = "")+
  theme(axis.text.y= element_text(size = 8))+
  coord_flip()
dev.off()



KEGG_metalfrac <- merge(KEGG, GO_MF_all_metal_related_anno, by = "ORF", all.x = T)%>%
  unique()%>%
  group_by(ORF,term.x)%>%
  mutate(metal_anno = ifelse(!is.na(term.y),1,0))%>%
  ungroup()%>%
  dplyr::select(-term.y)%>%
  unique()%>%
  group_by(term.x)%>%
  summarise(total_ORFs = length(ORF),
            num_metanno = sum(metal_anno),
            frac_metanno = num_metanno/total_ORFs)%>%
  filter(total_ORFs > 4 & total_ORFs < 200)


pdf(paste0(plot_dir,"/KEGG_pathway_metaldependence.pdf"),width = 10,height = 20)

ggplot(KEGG_metalfrac,
       aes(x = reorder(term.x,frac_metanno)))+
  geom_bar(stat = "identity",
           aes(y = total_ORFs), width = 0.5,fill = NA, linewidth = 0.5, color = brewer.pal(3, name = "Set2")[[2]])+
  geom_bar(stat = "identity",
           aes(y = num_metanno), width = 0.5, fill = brewer.pal(3, name = "Set2")[[2]])+
  geom_text(aes(y = total_ORFs + 5, label = paste0(round(frac_metanno,2)*100,"%")))+
  theme_metallica()+
  labs(fill = "", color = "")+
  coord_flip()+
  theme(axis.text.y= element_text(size = 10))
dev.off()



########################
### based on fluxes ###
#######################

GEMsim_0pt9grwth_fluxes <- read.csv(paste0(metpert_ecYeast8simulation_dir,"/output/tables/cofactorYeast8_flux_change_summary.csv"),stringsAsFactors = F)%>%
  
  ## filter out metal transport fluxes -- these were the ones that were minimized or maximized to achieve simulation 
  filter( !Reaction.Name %in% c("Ca(2+) exchange","Ca(2+) transport", ## Ca  
                                "Uptake of Cu(+)","Cu2(+) transport", "Cu2(+) exchange",
                                "iron (II) transport","	iron(3+) exchange","iron(2+) exchange",
                                "Mg(2+) exchange",
                                "Mn(2+) exchange",
                                "potassium exchange","potassium transport",
                                "sodium exchange",
                                "Zn(2+) exchange"))%>%
  mutate(corrected_flux_change = ifelse(is.na(fold_change_flux)| !is.finite(fold_change_flux), 
                                        ifelse(control_condition_flux == 0 & flux != 0, flux/0.001,
                                               ifelse(control_condition_flux !=0 & flux == 0 , 0.001/control_condition_flux,
                                                      ifelse(control_condition_flux ==0 & flux == 0, 1,NA))),fold_change_flux)) 



GEMsim_0pt9grwth_fluxes_separated<- GEMsim_0pt9grwth_fluxes %>%
  separate_rows(grRules, sep = " and | or ")


# Counting the number of NA values in the 'flux' column
na_count_flux <- sum(GEMsim_0pt9grwth_fluxes_separated$flux ==0)

# Counting the number of NA values in the 'fold_change_flux' column
na_count_fold_change_flux <- sum(is.na(GEMsim_0pt9grwth_fluxes_separated$fold_change_flux))
# Printing the counts
print(paste("Number of NA values in 'flux' column: ", na_count_flux))
print(paste("Number of NA values in 'fold_change_flux' column: ", na_count_fold_change_flux))

# plot all fluxes 

ggplot(GEMsim_0pt9grwth_fluxes_separated,
       aes(x = metal_perturbation,
           y = grRules,
           fill = log2(flux)))+
  geom_tile()


                                   
GEMsim_0pt9grwth_fluxes_separated <- GEMsim_0pt9grwth_fluxes_separated%>%
                                      group_by(grRules, metal_perturbation) %>%
                                      summarize(median_cfc_flux = median(corrected_flux_change, na.rm = TRUE))%>%
                                      ungroup()%>%
                                      dplyr::select(grRules, metal_perturbation,median_cfc_flux)%>%
                                      unique()

            
GEMsim_0pt9grwth_fluxes_separated_ORFsummary <- GEMsim_0pt9grwth_fluxes_separated%>%
                                                      mutate(Significant = abs(log2(median_cfc_flux)) > log2(1.5))%>%
                                                      group_by(grRules)%>%
                                                      mutate(Significant = any(Significant))%>%
                                                      ungroup()%>%
                                                      dplyr::select(grRules,Significant)%>%
                                                      unique()
                                                                                                  
####################################################################
### what % of ScGEM fluxes are changing ? - counter per ORF wise ### ANNOTATED FOR METAL BINDING or no binding
####################################################################

GEMsim_fluxes_fracsig_metrel <- GEMsim_0pt9grwth_fluxes_separated_ORFsummary%>%
  mutate(metal_anno = ifelse(grRules %in% GO_MF_all_metal_related_anno$ORF, TRUE,FALSE))%>%
  group_by(metal_anno)%>%
  summarise( 
    total_meas = length(grRules),
    num_sig = sum(Significant),
    frac_sig = num_sig/total_meas)%>%
    ungroup()%>%
  mutate(type = "flux")

############################################################################################
### Fraction of each metabolism related database that responds to env and int metal conc ###
############################################################################################

## Environmental metal conc

FCvsAE_allproteins_env<- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/extracellularconc_vs_proteinabundance/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),
                                  stringsAsFactors = F)%>%
  dplyr::select(Genes, ORF,Element,Significant,Log2FC_vs_AE,Element.Concentration)%>%
  unique()

FCvsAE_allproteins_env_sig <- FCvsAE_allproteins_env%>%
  dplyr::select(ORF,Significant)%>%
  unique()%>%
  group_by(ORF)%>%
  mutate(Significant = any(Significant))



## Total cellular metal conc

FCvsAE_allproteins_cell<- read.csv(paste0(metpert_WTproteomics_dir,
                                          "/output/tables/intracellularconc_vs_proteinabundance/lmfit_intracellular_metalconc_vs_log2PQ_with_SigNotSig_AdjPVthresh_0.05_FCthresh_1.5.csv"),
                                   stringsAsFactors = F)%>%
  dplyr::select(Genes, ORF,Element, median_relative_intracellular_concentration,median_log2_foldchangevsAE,Significant)%>%
  unique()

FCvsAE_allproteins_cell_sig <- FCvsAE_allproteins_cell%>%
  dplyr::select(ORF,Significant)%>%
  unique()%>%
  group_by(ORF)%>%
  mutate(Significant = any(Significant))


ScGEM_envsig <- merge(data.frame(ORF =metabolism_related_dbs[["ScGEM"]]),FCvsAE_allproteins_env_sig, by = "ORF")%>%
  mutate(type = "env_conc")

ScGEM_cellsig <- merge(data.frame(ORF =metabolism_related_dbs[["ScGEM"]]),FCvsAE_allproteins_cell_sig, by = "ORF")%>%
  mutate(type = "cell_conc")

ScGEM_envcellsig <- rbind(ScGEM_envsig,ScGEM_cellsig)%>%
  mutate(metal_anno = ifelse(ORF %in% GO_MF_all_metal_related_anno$ORF, TRUE,FALSE))%>%
  group_by(metal_anno,type)%>%
  summarize(total_meas = length(ORF),
            num_sig = sum(Significant),
            frac_sig = num_sig/total_meas)%>%
  ungroup()


ScGEM_fluxenvcell_sig <- rbind(GEMsim_fluxes_fracsig_metrel,
                               ScGEM_envcellsig[colnames(GEMsim_fluxes_fracsig_metrel)])%>%
                            mutate(type = factor(type, levels = c("flux","env_conc","cell_conc")))

pdf(paste0(plot_dir,"/percent_fluxenvcell_ScGEM.pdf"),width = 6,height = 4.5)
ggplot(ScGEM_fluxenvcell_sig)+
  geom_bar(aes(x = metal_anno,
               y = total_meas,
               colour = type),
           fill = NA, stat = "identity", position = "dodge",
           width = 0.5,
           linewidth = 0.25)+
  geom_bar(aes(x = metal_anno,
                       y = num_sig,
                       colour = type,
                       fill = type), 
                   alpha = 0.6,
                   linewidth = 0.25,
                   stat = "identity", position = "dodge2",
                   width = 0.5)+
  geom_text(aes(x = metal_anno,
                y = num_sig+10,
                group = type,
                label = paste0(100*round(frac_sig,2),"%") ),
            size = 5)+
  theme_metallica()+
  scale_color_brewer(palette = "Paired", direction = -1)+
  scale_fill_brewer(palette = "Paired",direction = -1)+
  labs(x= "")
dev.off()


metabolism_related_df_envcell_metnomet <- merge(rbind(metabolism_related_df_envsig,metabolism_related_df_cellsig),
                                       GO_MF_all_metal_related_anno, by = "ORF",all.x=T)%>%
  mutate(metal_related = ifelse(!is.na(term),"metal related","unknown"))%>%
  group_by(metal_related,database,type)%>%
  summarize(total_meas = length(ORF),
            num_sig = sum(Significant),
            frac_sig = num_sig/total_meas)%>%
  ungroup()%>%
  dplyr::select(metal_related, total_meas, num_sig, frac_sig, type, database)


metabolism_related_df_envcellflux_metnomet <- rbind(metabolism_related_df_envcell_metnomet,
                                   GEMsim_fluxes_fracsig_allmetaltogether)%>%
  mutate(type = factor(type, levels = c("flux","env_conc","cell_conc")))

pdf(paste0(plot_dir,"/percent_flux_DA_envcell_ScGEM_overall.pdf"),width = 6,height = 5)
ggplot(filter(metab_related_df_flux_exp,frac_sig >0))+
  geom_bar(aes(x = paste(metal_related),
               y = total_meas,
               colour = type),
           fill = NA, stat = "identity", position = "dodge",
           width = 0.5,
           linewidth = 0.25)+
  geom_bar(aes(x = paste(metal_related),
               y = num_sig,
               colour = type,
               fill = type), 
           alpha = 0.6,
           linewidth = 0.25,
           stat = "identity", position = "dodge2",
           width = 0.5)+
  geom_text(aes(x = paste(metal_related),
                y = num_sig+10,
                group = type,
                label = paste0(100*round(frac_sig,2),"%"),
  ),
  size = 3)+
  theme_metallica()+
  scale_color_brewer(palette = "Paired", direction = -1)+
  scale_fill_brewer(palette = "Paired",direction = -1)+
  labs(x= "", title = "% scGEM affected")
dev.off()

write.csv(metabolism_related_df_envcell,paste0(output_tables_dir,"/percentage_flux_envDA_cellDA_metabolism_dbs_overall.csv"),row.names =F)

#################################################
### plot distribution of all simulated fluxes ###
#################################################

pdf(paste0(plot_dir,"/log2_fluxchanges_density.pdf"),width = 7,height = 4.5)
ggplot(GEMsim_0pt9grwth_fluxes,
       aes(x = log2_foldchange_flux,
           y = after_stat(scaled),
           color = metal_perturbation,
           fill = metal_perturbation))+
  geom_density(alpha = 0.2, linewidth = 0.25)+
  scale_fill_manual(values = colkey_EleDir)+  
  scale_color_manual(values = colkey_EleDir)+
  theme_metallica()+
  labs(fill ="", colour = "")+
  scale_x_continuous(limits = c(-log2(10),log2(10)))+
  theme(legend.text = element_text(size = 7))

ggplot(GEMsim_0pt9grwth_fluxes,
       aes(x = log2(corrected_flux_change),
           y = after_stat(scaled),
           color = metal_perturbation,
           fill = metal_perturbation))+
  geom_density(alpha = 0.2, linewidth = 0.25)+
  scale_fill_manual(values = colkey_EleDir)+  
  scale_color_manual(values = colkey_EleDir)+
  theme_metallica()+
  labs(fill ="", colour = "")+
  scale_x_continuous(limits = c(-log2(10),log2(10)))+
  scale_y_continuous(limits = c(0,0.1))+
  theme(legend.text = element_text(size = 7))

ggplot(GEMsim_0pt9grwth_fluxes,
       aes(x = log2(corrected_flux_change),
           y = after_stat(scaled),
           color = metal_perturbation,
           fill = metal_perturbation))+
  geom_density(alpha = 0.2, linewidth = 0.25)+
  scale_fill_manual(values = colkey_EleDir)+  
  scale_color_manual(values = colkey_EleDir)+
  theme_metallica()+
  labs(fill ="", colour = "")+
  scale_x_continuous(limits = c(-log2(2),log2(2)))+
  scale_y_continuous(limits = c(0,0.10))+
  theme(legend.text = element_text(size = 7))


ggplot(filter(GEMsim_0pt9grwth_fluxes, grepl("Zn",metal_perturbation) & 
                abs(log2(corrected_flux_change)) > log2(1.2)),
       aes(x = log2(corrected_flux_change),
           y = after_stat(scaled),
           color = metal_perturbation,
           fill = metal_perturbation))+
  geom_density(alpha = 0.2, linewidth = 0.25)+
  scale_fill_manual(values = colkey_EleDir)+  
  scale_color_manual(values = colkey_EleDir)+
  theme_metallica()+
  labs(fill ="", colour = "")+
  scale_x_continuous(limits = c(-log2(10),log2(10)))+
  theme(legend.text = element_text(size = 7))

ggplot(filter(GEMsim_0pt9grwth_fluxes, grepl("Fe",metal_perturbation) & 
                abs(log2(corrected_flux_change)) > log2(1.2)),
       aes(x = log2(corrected_flux_change),
           y = after_stat(scaled),
           color = metal_perturbation,
           fill = metal_perturbation))+
  geom_density(alpha = 0.2, linewidth = 0.25)+
  scale_fill_manual(values = colkey_EleDir)+  
  scale_color_manual(values = colkey_EleDir)+
  theme_metallica()+
  labs(fill ="", colour = "")+
  scale_x_continuous(limits = c(-log2(10),log2(10)))+
  theme(legend.text = element_text(size = 7))

ggplot(filter(GEMsim_0pt9grwth_fluxes, grepl("Cu",metal_perturbation) & 
                abs(log2(corrected_flux_change)) > log2(1.2)),
       aes(x = log2(corrected_flux_change),
           y = after_stat(scaled),
           color = metal_perturbation,
           fill = metal_perturbation))+
  geom_density(alpha = 0.2, linewidth = 0.25)+
  scale_fill_manual(values = colkey_EleDir)+  
  scale_color_manual(values = colkey_EleDir)+
  theme_metallica()+
  labs(fill ="", colour = "")+
  scale_x_continuous(limits = c(-log2(10),log2(10)))+
  theme(legend.text = element_text(size = 7))

ggplot(filter(GEMsim_0pt9grwth_fluxes, grepl("Ca",metal_perturbation) & 
                abs(log2(corrected_flux_change)) > log2(1.2)),
       aes(x = log2(corrected_flux_change),
           y = after_stat(scaled),
           color = metal_perturbation,
           fill = metal_perturbation))+
  geom_density(alpha = 0.2, linewidth = 0.25)+
  scale_fill_manual(values = colkey_EleDir)+  
  scale_color_manual(values = colkey_EleDir)+
  theme_metallica()+
  labs(fill ="", colour = "")+
  scale_x_continuous(limits = c(-log2(10),log2(10)))+
  theme(legend.text = element_text(size = 7))

ggplot(filter(GEMsim_0pt9grwth_fluxes, grepl("Mg",metal_perturbation) & 
                abs(log2(corrected_flux_change)) > log2(1.2)),
       aes(x = log2(corrected_flux_change),
           y = after_stat(scaled),
           color = metal_perturbation,
           fill = metal_perturbation))+
  geom_density(alpha = 0.2, linewidth = 0.25)+
  scale_fill_manual(values = colkey_EleDir)+  
  scale_color_manual(values = colkey_EleDir)+
  theme_metallica()+
  labs(fill ="", colour = "")+
  scale_x_continuous(limits = c(-log2(10),log2(10)))+
  theme(legend.text = element_text(size = 7))

dev.off()


############################################################################################################
### prepare data for iPATH - metal annotations in GO-MF across metabolic network - visualised with iPATH ###
############################################################################################################


GO_MF_metalanno_for_iPATH <- GO_MF_all_metal_related_anno%>%
  group_by(ORF)%>%
  mutate(specific_anno = ifelse(term == "unspecific", 0,length(unique(term))))%>%
  dplyr::select(ORF, specific_anno)%>%
  unique()%>%
  group_by(ORF)%>%
  summarise(num_specific_anno = sum(specific_anno))%>%
  ungroup()


iPath_devperORF <- GO_MF_metalanno_for_iPATH %>%
  mutate(Uniprot_ID = as.character(lapply(ORF,convert_ORF2Uniprot)))%>%
  filter(num_specific_anno != 0)%>%
  mutate(colour = as.character(lapply(num_specific_anno,convert_num_literature_metalanno_to_color)),
         width = "W12",
         opacity = 1,
         UNIPROT.ID = paste0("UNIPROT:",Uniprot_ID))%>%
  dplyr::select(UNIPROT.ID,colour,width,opacity)

write_delim(iPath_devperORF, paste0(output_tables_dir,"/foriPATHmedian_metalbinding_annots_perORF.tsv"),col_names = F)


## binary - metal or no metal anno

GO_MF_metalanno_for_iPATH <- GO_MF_all_metal_related_anno %>%
  group_by(ORF) %>%
  mutate(specific_anno = ifelse(term == "unspecific", 0, length(unique(term)))) %>%
  dplyr::select(ORF, specific_anno) %>%
  unique() %>%
  group_by(ORF) %>%
  summarise(num_specific_anno = sum(specific_anno)) %>%
  ungroup()

iPath_binary_color <- GO_MF_metalanno_for_iPATH %>%
  mutate(Uniprot_ID = as.character(lapply(ORF, convert_ORF2Uniprot))) %>%
  mutate(colour = ifelse(num_specific_anno > 0, viridis(5)[2], NA), # binary color assignment
         width = "W12",
         opacity = 1,
         UNIPROT.ID = paste0("UNIPROT:", Uniprot_ID)) %>%
  na.omit()%>%
  dplyr::select(UNIPROT.ID, colour, width, opacity)

write_delim(iPath_binary_color, paste0(output_tables_dir, "/foriPATHbinary_metalrelated_annos_binary.tsv"), col_names = FALSE)


## make iPATH df for the flux and protein expression data

all_ORFs <- unique(c(FCvsAE_allproteins_cell_sig$ORF,FCvsAE_allproteins_env_sig$ORF,GEMsim_0pt9grwth_fluxes_separated_ORFsummary$grRules))
env_sig_orfs <- unique(filter(FCvsAE_allproteins_env, Significant == 1)$ORF)
cell_sig_orfs <- unique(filter(FCvsAE_allproteins_cell, Significant== 1)$ORF)
flux_sig_orfs <- unique(filter(GEMsim_0pt9grwth_fluxes_separated_ORFsummary,Significant)$grRules)

# Combine into one dataframe
combined_significance <- data.frame(ORF = all_ORFs) %>%
  mutate(
    Significant_FCvsAE_env = ORF %in% env_sig_orfs,
    Significant_FCvsAE_cell = ORF %in% cell_sig_orfs,
    Significant_GEM = ORF %in% flux_sig_orfs
  )


# Assign colors based on significance conditions
iPath_colors <- combined_significance %>%
  mutate(Color = case_when(
    Significant_FCvsAE_env & Significant_FCvsAE_cell & Significant_GEM ~ viridis(8)[1],
    (Significant_FCvsAE_env | Significant_FCvsAE_cell) & !Significant_GEM ~ viridis(8)[4],
    !Significant_FCvsAE_env & !Significant_FCvsAE_cell & Significant_GEM ~ viridis(8)[7],
    TRUE ~ NA_character_
  ),
  Color = factor(Color, levels = c(viridis(8)[1],viridis(8)[4],viridis(8)[7]))) %>%
  na.omit()%>%
  select(ORF, Color)

iPath_for_output <- iPath_colors %>%
  filter(!is.na(Color)) %>%
  mutate(Uniprot_ID = sapply(ORF, convert_ORF2Uniprot), # Use sapply for vectorized operation
         width = "W12",
         opacity = 1,
         UNIPROT.ID = paste0("UNIPROT:", Uniprot_ID)) %>%
  select(UNIPROT.ID, Color, width, opacity)%>%
  arrange(Color)

write_delim(iPath_for_output, paste0(output_tables_dir, "/foriPATHbinary_flux_expproteomics.tsv"), col_names = FALSE)

