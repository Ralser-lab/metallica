#`---
#`  Title: "compare_ecYeast8simulation_with_measuredproteinabundance.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 6 June 2021 
#`  Description: Script to compare results of simulation of metal depletion and excess with experimentally obtained protein quantity results
#`---

library(tidyr)
#############################################
### source paths functions and libraries  ###
#############################################


source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")

# general

source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/database_identifier_conversion_functions.R"))
source(paste0(code_dir,"/common_code/input_processed_databases_publisheddatasets.R"))



source(paste0(code_dir,"/metpert_ecYeast8simulation/metpert_ecYeast8simulation_0_libraries_functions.R"))


plot_dir <- paste0(metpert_ecYeast8simulation_dir,"/output/plots")
dir.create(plot_dir,recursive = T)

output_tables_dir <- paste0(metpert_ecYeast8simulation_dir,"/output/tables")
dir.create(output_tables_dir,recursive = T)

################################################################################
### Read in experimentally measured protein abundance and simulation results ###
################################################################################

protabun_data<- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/extracellularconc_vs_proteinabundance/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),
                         stringsAsFactors = F)


protabun_data_forcomparison <- protabun_data%>%
  dplyr::select(Genes,ORF,Significant,Element,Element.Concentration,Log2FC_vs_AE)%>%
  mutate(Perturbation = ifelse(Element.Concentration > 1.2,"Excess",
                               ifelse(Element.Concentration < 0.8,"Depletion",
                                      "Control")))%>%
  na.omit()%>%
  mutate(ORF_ElePert = paste(ORF,Element,Perturbation))%>%
  group_by(Genes, Perturbation,Element)%>%
  mutate(MaxMag_log2_experimental_protein_fc = ifelse(abs(max(Log2FC_vs_AE,na.rm=T)) > abs(min(Log2FC_vs_AE,na.rm=T)),
                                                      max(Log2FC_vs_AE,na.rm = T),
                                                      min(Log2FC_vs_AE,na.rm = T)),
         MaxMag_experimental_protein_fc = 2^(MaxMag_log2_experimental_protein_fc))%>%
  ungroup()%>%
  dplyr::select(ORF_ElePert,Element,MaxMag_experimental_protein_fc,Significant)%>%
  unique()

##################################################
### Compare flux to measured protein abundance ###
##################################################

GEMsim_0pt9grwth_fluxes <- data.frame(read_excel(paste0(metpert_ecYeast8simulation_dir,
                                                        "/matlab_simulations/cofactorYeast8_GEM_simulations_from_YuChen_08062023.xlsx"),
                                                 sheet = 3),stringsAsFactors = F)

colnames(GEMsim_0pt9grwth_fluxes) <- GEMsim_0pt9grwth_fluxes[1,]

GEMsim_0pt9grwth_fluxes <- GEMsim_0pt9grwth_fluxes[-1,]%>%
  reshape2::melt(id.vars = c("Reaction ID","Reaction Name","Formula","grRules","Control","KEGG ID","EC number"),
                 variable.name = "metal_perturbation",
                 value.name = "flux")%>%
  mutate(metal_perturbation = gsub("_"," ",metal_perturbation),
         metal = metal_perturbation,
         metal = gsub(" Excess","",metal),
         metal = gsub(" Depletion","",metal),
         flux = as.numeric(flux),
         control_condition_flux = as.numeric(Control),
         fold_change_flux = flux/control_condition_flux,
         log2_foldchange_flux = log2(fold_change_flux),
         log2_foldchange_flux = ifelse(flux == 0 & control_condition_flux == 0 & is.infinite(log2_foldchange_flux),NA,
                                       ifelse(flux == 0 & control_condition_flux != 0 & is.infinite(log2_foldchange_flux),-sign(control_condition_flux),
                                              ifelse(flux !=0 & control_condition_flux == 0 & is.infinite(log2_foldchange_flux),sign(flux),
                                                     log2_foldchange_flux))),
         flux_change_detected = ifelse(abs(log2_foldchange_flux) > log2(1.5),T,F ))

write.csv(GEMsim_0pt9grwth_fluxes, paste0(output_tables_dir,"cofactorYeast8_flux_change_summary.csv"),row.names = F)

GEMsim_0pt9grwth_fluxes <- GEMsim_0pt9grwth_fluxes%>%
  
  ## filter out metal transport fluxes -- these were the ones that were minimized or maximized to achieve simulation 
  filter( !`Reaction Name` %in% c("Ca(2+) exchange","Ca(2+) transport", ## Ca  
                                "Uptake of Cu(+)","Cu2(+) transport", "Cu2(+) exchange",
                                "iron (II) transport","	iron(3+) exchange","iron(2+) exchange",
                                "Mg(2+) exchange",
                                "Mn(2+) exchange",
                                "potassium exchange","potassium transport",
                                "sodium exchange",
                                "Zn(2+) exchange"))%>%
  ## filter out very high flux fold changes 
  filter(abs(log2(fold_change_flux)) < log2 (5) )


GEMsim_0pt9grwth_fluxes_separated<- GEMsim_0pt9grwth_fluxes %>%
  separate_rows(grRules, sep = " and | or ")

#### visualizeflux distriutions across all simulations

pdf(paste0(plot_dir,"/log2_fluxchanges_density.pdf"),width = 7,height = 4.5)
ggplot(GEMsim_0pt9grwth_fluxes,
       aes(x = log2_foldchange_flux,
           y = after_stat(scaled),
           color = metal_perturbation,
           fill = metal_perturbation))+
  geom_density(alpha = 0.4, linewidth = 0.25)+
  scale_fill_manual(values = colkey_EleDir)+  
  scale_color_manual(values = colkey_EleDir)+
  theme_metallica()+
  labs(fill ="", colour = "")+
  scale_x_continuous(limits = c(-1.5,1.5))+
  theme(legend.text = element_text(size = 7))

dev.off()

# get rid of Inf
finite_rows <- is.finite(GEMsim_0pt9grwth_fluxes_separated$flux) & is.finite(GEMsim_0pt9grwth_fluxes_separated$fold_change_flux)
GEMsim_0pt9grwth_fluxes_separated <- GEMsim_0pt9grwth_fluxes_separated[finite_rows, ]

## replace all Inf, 0 and NA values with 1 in flux change
GEMsim_0pt9grwth_fluxes_separated$fold_change_flux[is.na(GEMsim_0pt9grwth_fluxes_separated$fold_change_flux)] <- 1

GEMsim_0pt9grwth_fluxes_separated$ORF_ElePert = paste(GEMsim_0pt9grwth_fluxes_separated$grRules, GEMsim_0pt9grwth_fluxes_separated$metal_perturbation)

##################################################################################
### Summarise and write to disk with significant or not significant annotation ###
##################################################################################

View(GEMsim_0pt9grwth_fluxes_separated)

GEM_flux_res_summary_per_metal <- GEMsim_0pt9grwth_fluxes_separated%>%
                                  group_by(grRules, metal)%>%
                                  summarise(flux_change_detected = any(flux_change_detected,na.rm=T))%>%
                                  ungroup()

write.csv(GEM_flux_res_summary_per_metal, paste0(output_tables_dir,"/flux_change_detected_summary_per_metal.csv"),row.names = F)

GEM_flux_res_summary_per_metalperturbation <- GEMsim_0pt9grwth_fluxes_separated%>%
  group_by(grRules, metal_perturbation)%>%
  summarise(flux_change_detected = any(flux_change_detected,na.rm=T))%>%
  ungroup()

write.csv(GEM_flux_res_summary_per_metal, paste0(output_tables_dir,"/flux_change_detected_summary_per_metalperturbation.csv"),row.names = F)

###############################################################
### merge flux simulation and protein abundanc fold change  ###
###############################################################

simulation_vs_experiment <- merge(GEMsim_0pt9grwth_fluxes_separated[,c("ORF_ElePert","metal_perturbation","grRules","fold_change_flux")],
                                  protabun_data_forcomparison,
                                  by = "ORF_ElePert")
                           ## filter out fc = 1 flux changes
                          #  filter(fold_change_flux != 1 & abs(log2(fold_change_flux)) > log2(0.02))

## visualise raw data
ggplot(simulation_vs_experiment,
       aes(x = log2(fold_change_flux),
           y = log2(MaxMag_experimental_protein_fc)))+
  geom_point(aes(color = factor(Significant)), alpha = 0.7)+
  facet_wrap("metal_perturbation",scales = "free")+
  geom_smooth(method = "lm",color = "black")+
  scale_color_manual(values = c("lightgray","darkred"))+
  theme_minimal()

# Creating a new data frame to store the results
sim_vs_exp_pearsoncor <- data.frame(
  metal_perturbation = unique(simulation_vs_experiment$metal_perturbation),
  correlation = NA,
  p_value = NA
)

# Looping through each metal perturbation to calculate the correlation and p-value
for (i in seq_along(sim_vs_exp_pearsoncor$metal_perturbation)) {
  current_group <- simulation_vs_experiment %>%
    filter(metal_perturbation == sim_vs_exp_pearsoncor$metal_perturbation[i])
  
  cor_test <- cor.test(
    log2(current_group$fold_change_flux),
    log2(current_group$MaxMag_experimental_protein_fc),
    method = "pearson"
  )
  
  sim_vs_exp_pearsoncor$correlation[i] <- cor_test$estimate
  sim_vs_exp_pearsoncor$p_value[i] <- cor_test$p.value
}

pdf(paste0(plot_dir,"/simulated_fluxchanges_vs_expproteinchanges.pdf"),width = 4.5,height = 4)
# Visualizing the results using a scatter plot
ggplot(sim_vs_exp_pearsoncor, 
       aes(x = correlation,
           y = -log10(p_value), 
           label = metal_perturbation,
           color = metal_perturbation)) +
  geom_point(aes(size = -log10(p_value))) +
 # geom_hline(yintercept = -log10(0.05),color = "black",linewidth = 0.25)+
  geom_text_repel(size = 2.5) +
  labs(
    x = "pearson correlation coefficient",
    y = "P-value of cor.test"
  ) +
  scale_color_manual(values = colkey_EleDir)+
  theme_metallica()+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8))
dev.off()


###################################################################
### pathway wise correlation between simulated flux and exp fc  ###
###################################################################

simulation_vs_experiment_KEGG <- merge(KEGG, simulation_vs_experiment,by.x= "ORF", by.y = "grRules")

simulation_vs_experiment_KEGG_group_list <- simulation_vs_experiment_KEGG %>%
  group_by(term) %>%
  filter(n() >4 & length(unique(fold_change_flux)) > 4) %>%
  group_split()

# Step 2: Initialize a list to store the results
sim_vs_exp_KEGG_pearsoncor <- list()

# Step 3: Loop through each term group and compute correlation and p-value
for (i in seq_along(simulation_vs_experiment_KEGG_group_list)) {
  data <- simulation_vs_experiment_KEGG_group_list[[i]]
  
  cor_test <- cor.test(
    log2(data$fold_change_flux),
    log2(data$MaxMag_experimental_protein_fc),
    method = "pearson"
  )
  
  sim_vs_exp_KEGG_pearsoncor[[i]] <- data.frame(
    term = unique(data$term),
    correlation = cor_test$estimate,
    p_value = cor_test$p.value,
    num_ORFs = length(unique(data$ORF))
  )
}

# Step 4: Bind all results into a single data frame
sim_vs_exp_KEGG_pearsoncor <- bind_rows(sim_vs_exp_KEGG_pearsoncor)


pdf(paste0(plot_dir, "/simulated_fluxchanges_vs_expproteinchanges_KEGGpathwaywise.pdf"), width = 9, height = 5.5)

ggplot(sim_vs_exp_KEGG_pearsoncor, 
       aes(x = correlation, 
           y = -log10(p_value), 
           label = term, 
           color = log2(num_ORFs))) +
  geom_point(aes(size = -log10(p_value))) +
  geom_text_repel(size = 2.5) +
  labs(
    x = "pearson correlation coefficient",
    y = "-log10(P-value of cor.test)"
  ) +
  theme_metallica() +
  scale_color_viridis(option = "C",end = 0.8)

dev.off()


###################################################################
### pathway wise correlation between simulated flux and exp fc  ###
###################################################################

simulation_vs_experiment_KEGG <- merge(KEGG, simulation_vs_experiment,by.x= "ORF", by.y = "grRules")

simulation_vs_experiment_KEGG_group_list <- simulation_vs_experiment_KEGG %>%
  group_by(term,metal_perturbation) %>%
  filter(n() >4 & length(unique(fold_change_flux)) > 4) %>%
  group_split()

# Step 2: Initialize a list to store the results
sim_vs_exp_KEGG_pearsoncor <- list()

# Step 3: Loop through each term group and compute correlation and p-value
for (i in seq_along(simulation_vs_experiment_KEGG_group_list)) {
  data <- simulation_vs_experiment_KEGG_group_list[[i]]
  
  cor_test <- cor.test(
    log2(data$fold_change_flux),
    log2(data$MaxMag_experimental_protein_fc),
    method = "pearson"
  )
  
  sim_vs_exp_KEGG_pearsoncor[[i]] <- data.frame(
    term = unique(data$term),
    metal_perturbation = unique(data$metal_perturbation),
    correlation = cor_test$estimate,
    p_value = cor_test$p.value,
    num_ORFs = length(unique(data$ORF))
  )
}

# Step 4: Bind all results into a single data frame
sim_vs_exp_KEGG_pearsoncor <- bind_rows(sim_vs_exp_KEGG_pearsoncor)%>%
                              group_by(term)%>%
                              mutate(pathway_score = correlation * -log10(p_value))%>%
                              ungroup()%>%
                              group_by(metal_perturbation)%>%
                              mutate(metpert_score = correlation * -log10(p_value))%>%
                              ungroup()

pdf(paste0(plot_dir,"/exp_vs_sim_KEGGmetpert_summary.pdf"),width = 6,height = 8)
ggplot(sim_vs_exp_KEGG_pearsoncor, 
        aes(x = reorder(metal_perturbation,-metpert_score), 
            y = reorder(term,pathway_score),
            fill = correlation)) +
   geom_point(aes(size = -log10(p_value)), shape =21, color = "black") +
   labs(
     x = "",
     y = ""
   ) +
   theme_metallica() +
   scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", 
                         midpoint = 0, limits = c(-1, 1)) + 
   theme(legend.position = "bottom",
         axis.text.x = element_text(angle = 90, size = 8),
         axis.text.y = element_text(size = 8),
         legend.title = element_text(size = 8),
         legend.text = element_text(size = 8))
dev.off() 

