#`---
#`  Title: "compare_ecYeast8simulation_with_measuredproteinabundance.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 6 June 2021 
#`  Description: Script to compare results of simulation of metal depletion and excess with experimentally obtained protein quantity results
#`---

library(readr)
#############################################
### source paths functions and libraries  ###
#############################################

source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")

# general

source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/database_identifier_conversion_functions.R"))
source(paste0(code_dir,"/common_code/input_processed_databases_publisheddatasets.R"))

# specific

source(paste0(code_dir,"/metpert_ecYeast8simulation/metpert_ecYeast8simulation_0_libraries_functions.R"))


plot_dir <- paste0(metpert_ecYeast8simulation_dir,"/output/plots/")
dir.create(plot_dir,recursive = T)

output_tables_dir <- paste0(metpert_ecYeast8simulation_dir,"/output/tables/")
dir.create(output_tables_dir,recursive = T)

################################################################################
### Read in experimentally measured protein abundance and simulation results ###
################################################################################

all_datasets <- read.csv(paste0(datasetcomparison_dir,"/output/tables/all_datasets_combined_significant_or_not.csv"),stringsAsFactors = F)

GEMsim_0pt9grwth_fluxes <- data.frame(read_excel(paste0(metpert_ecYeast8simulation_dir,"/matlab_simulations/cofactorYeast8_GEM_simulations_from_YuChen_08062023.xlsx"),sheet = 3),stringsAsFactors = F)

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
                                  log2_foldchange_flux = ifelse(flux == 0 & control_condition_flux == 0 & is.infinite(log2_foldchange_flux),na,
                                                                ifelse(flux == 0 & control_condition_flux != 0 & is.infinite(log2_foldchange_flux),1,
                                                                       ifelse(flux !=0 & control_condition_flux == 0 & is.infinite(log2_foldchange_flux),1,log2_foldchange_flux))),
                                  flux_change_detected = ifelse(abs(log2_foldchange_flux) > log2(1.5),T,F ))

write.csv(GEMsim_0pt9grwth_fluxes, paste0(output_tables_dir,"cofactorYeast8_flux_change_summary.csv"),row.names = F)

write.csv(GEMsim_0pt9grwth_fluxes[,c("")])                           
## density plot of log2 flux changes

ggplot(GEMsim_0pt9grwth_fluxes,
       aes(x = log2_foldchange_flux,
           colour = metal_perturbation,
           fill = metal_perturbation))+
  geom_density()+
  facet_wrap("metal_perturbation",scales = "free")+
  scale_colour_manual(values = colkey_EleDir)+
  scale_fill_manual(values = colkey_EleDir)+
  theme_metallica()



flux_changes_summary_per_metal <- GEMsim_0pt9grwth_fluxes%>%
                                  group_by(metal_perturbation)%>%
                                  summarize(num_fluxes_changed = sum(flux_change_detected, na.rm = T))

pdf(paste0(plot_dir,"/num_fluxes_changed_metalperturabtion_fcvscntrl1pt5cutoff.pdf"),width = 7,height = 7)
ggplot(flux_changes_summary_per_metal,
       aes(x = metal_perturbation,
           y = num_fluxes_changed,
           fill = metal_perturbation))+
  geom_bar(stat = "identity", width = 0.5, colour = "white")+
  scale_fill_manual(values = colkey_EleDir)+
  theme_metallica()+
  theme(axis.text.x = element_text(angle = 90))+
  labs(x = "", fill = "")
dev.off()


KEGG_rxnID_to_pathway <- GEMsim_0pt9grwth_fluxes%>%
                         dplyr::select(`KEGG ID`)%>%
                         na.omit()%>%
                         unique()

colnames(KEGG_rxnID_to_pathway) <- "reactionID"

## loop over everything to extract all pathway annotations then add as many rows for each reaction ID as pathway annos -- multiple pathways per rxn therefore not using lapply

iters <- nrow(KEGG_rxnID_to_pathway)

for(i in 1:iters){
  
  KEGG_pathways = convert_KEGGrxnID2pathwaylist(KEGG_rxnID_to_pathway[i,"reactionID"])
  
  if(length(KEGG_pathways) == 1){
    
    KEGG_rxnID_to_pathway[i,"KEGG_pathwayName"] = KEGG_pathways[[1]]
    KEGG_rxnID_to_pathway[i,"KEGG_pathwayID"] = names(KEGG_pathways[1])
 
  }
  else if(length(KEGG_pathways) > 1){
    
    KEGG_rxnID_to_pathway[i,"KEGG_pathwayName"] = KEGG_pathways[[1]]
    KEGG_rxnID_to_pathway[i,"KEGG_pathwayID"] = names(KEGG_pathways[1])
    
    for(p in 2:length(KEGG_pathways)){
      
      row2bind <- cbind(reactionID = KEGG_rxnID_to_pathway[i,"reactionID"],
                        KEGG_pathwayName = KEGG_pathways[[p]],
                        KEGG_pathwayID = names(KEGG_pathways[2]))
      
      KEGG_rxnID_to_pathway <- rbind(KEGG_rxnID_to_pathway, row2bind) 
    }
  }else{
    
    KEGG_rxnID_to_pathway[i,"KEGG_pathwayName"] = NA
    KEGG_rxnID_to_pathway[i,"KEGG_pathwayID"] = NA
    
  }
}


## loop in batches of 10 to convert kegg reaction ids

#for(l in seq(1,round(nrow(KEGG_rxnID_to_pathway)/10,0)*10,by = 10)) {
  
 # KEGG_rxnID_to_pathway_list <- c(KEGG_rxnID_to_pathway_list,keggGet(KEGG_rxnID_to_pathway$reactionID[seq(l:l+9)]))
  
#}

## add teh last 4

#KEGG_rxnID_to_pathway_list <- c(KEGG_rxnID_to_pathway_list, 
 #                               keggGet(KEGG_rxnID_to_pathway$reactionID[seq(round(nrow(KEGG_rxnID_to_pathway)/10,0)*10,nrow(KEGG_rxnID_to_pathway))]))

## convert the nested list into dataframe

# create empty dataframe
#KEGGrxn2pathway_mapped <- data.frame("reactionID" = character(), "pathwayID" = character(), "pathwayName" = character(), stringsAsFactors = FALSE)

# suppose your data is stored in a list variable named data_list
#for (list_item in KEGG_rxnID_to_pathway_list) {
 # reaction_id <- list_item$ENTRY
  
  #if (!is.null(list_item$PATHWAY)) {
   # pathway_ids <- names(list_item$PATHWAY)
    #pathway_names <- list_item$PATHWAY
    
    #for (i in seq_along(pathway_ids)) {
     # KEGGrxn2pathway_mapped <- rbind(df, data.frame("reactionID" = reaction_id, 
      #                           "pathwayID" = pathway_ids[i], 
       #                          "pathwayName" = pathway_names[i], 
        #                         stringsAsFactors = FALSE))
    #}
  #}
#}

GEMsim_0pt9grwth_fluxes <- merge(GEMsim_0pt9grwth_fluxes,KEGG_rxnID_to_pathway, by.x = "KEGG ID", by.y = "reactionID",all.x = T)

write.csv(GEMsim_0pt9grwth_fluxes,paste0(output_tables_dir,"/GEM_flux_results_with_KEGGpathwaynames.csv"),row.names =F)

################################################################
### Merge GO MF metalwise and KEGG to get pathway wise metal binding 
################################################################


KEGG_metal <- KEGG%>%
              mutate(metal_related = ifelse(ORF %in% unique(GO_gset_MF_metalbinding$ORF), T, F),
                     metal_related_metwiseGOMF = ifelse(ORF %in% unique(GO_gset_MF_binding_metalwise$ORF), T, F))%>%
              group_by(term)%>%
              mutate(metal_related_interm = any(metal_related),
                     metal_related_metwiseGOMF_interm = any(metal_related_metwiseGOMF),
                     KEGG_pathwayName = term)%>%
              ungroup()%>%
              dplyr::select(KEGG_pathwayName, metal_related_interm, metal_related_metwiseGOMF_interm)%>%
              unique()


flux_change_pathwaywise_metalwise_smry <- merge(
                                      unique(na.omit(GEMsim_0pt9grwth_fluxes[,c("KEGG ID","KEGG_pathwayName", "metal","flux_change_detected")])),
                                      KEGG_metal, by = "KEGG_pathwayName")%>%
                                      mutate(metal_anno = ifelse(metal_related_metwiseGOMF_interm,
                                                                 "metal specific anno",
                                                                 "unspecific/unknown"))%>%
                                      group_by(KEGG_pathwayName, metal)%>%
                                      mutate(flux_change_detected_inpathway = any(flux_change_detected))%>%
                                      ungroup()



# plot
simflux_plot <- ggplot(unique(flux_change_pathwaywise_metalwise_smry[,c("KEGG_pathwayName",
                                                        "metal","flux_change_detected_inpathway",
                                                        "metal_anno")]), 
                           aes(x = metal,
                               y = KEGG_pathwayName,
                               fill = flux_change_detected_inpathway)) + 
                      facet_wrap("metal_anno")+
                      geom_tile(colour ="black",linewidth = 0.05, alpha = 0.7)+
                      scale_fill_manual(values = c(viridis(10,option = "C")[[9]],viridis(5,option = "C")[[1]]))+
                      theme_metallica()+
                      theme(axis.text.y = element_text(size = 13))+
                      labs(x="", y = "", fill = "")

all_datasets_KEGG_metalanno <- merge(all_datasets,KEGG, by= "ORF")%>%
                               filter(term %in% unique(flux_change_pathwaywise_metalwise_smry$KEGG_pathwayName))%>% 
                               mutate(metal_anno = ifelse(ORF %in% unique(GO_gset_MF_binding_metalwise$ORF),T,F))%>%
                               group_by(term)%>%
                               mutate(metal_anno = ifelse(any(metal_anno), "metal specific anno", "unspecific/unknown"))%>%
                               ungroup()%>%
                               group_by(term,metal,dataset)%>%
                               mutate(sig_in_pathway = any(Significant))%>%
                               ungroup()%>%
                               dplyr::select(term,metal_anno,metal,dataset,sig_in_pathway)%>%
                               unique()

all_datasets_KEGG_metalanno_smry <- all_datasets_KEGG_metalanno%>%
                               group_by(term,metal,metal_anno)%>%
                               summarise(num_ds_sig_in = sum(sig_in_pathway,na.rm=T))%>%
                               ungroup()

## get number of pathways that respond in flux simulations 

unique(flux_change_pathwaywise_metalwise_smry$KEGG_pathwayName)
length(unique(filter(flux_change_pathwaywise_metalwise_smry, flux_change_detected_inpathway)$KEGG_pathwayName))

## plot num pathways that flux change was detected in per metal
num_pathways_w_fluxchange <- flux_change_pathwaywise_metalwise_smry %>%
  filter(flux_change_detected_inpathway == TRUE) %>%
  group_by(metal) %>%
  summarise(num_pathways = n_distinct(KEGG_pathwayName))

# Plot with ggplot
np_fluxch_plot <- ggplot(num_pathways_w_fluxchange, 
       aes(x = metal,
           y = num_pathways,
           fill = metal)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(y = num_pathways + 1,
                label = num_pathways)) +
  labs(x = "Metal", y = "Number of Pathways") +
  scale_fill_manual(values = colkey_Ele)+
  theme_metallica()+
  labs(x = "", y = "number of pathways with flux change")

num_pathways_w_fluxchange$dataset = "GEM_simulation"
num_pathways_w_fluxchange$num_pathway_any_dataset = num_pathways_w_fluxchange$num_pathways

## get number of pathways that respond in  atleast 1 dataset
length(unique(all_datasets_KEGG_metalanno$term))
length(unique(filter(all_datasets_KEGG_metalanno, num_ds_sig_in > 0 )$term))

## plot number of responsive pathways vs number of datasets they respond in 
# Group by 'metal' and 'num_ds_sig_in' and count the number of unique 'term'
ds_hit_KEGGpathways <- all_datasets_KEGG_metalanno %>%
                        filter(sig_in_pathway == TRUE) %>%
                        group_by(metal)%>%
                        mutate(num_pathway_any_dataset = n_distinct(term))%>%
                        ungroup()%>%
                        group_by(metal,dataset,num_pathway_any_dataset) %>%
                        summarise(num_pathways = n_distinct(term))%>%
                        ungroup()%>%
                        mutate(dataset = factor(dataset,
                                                levels = c("metdepKOgrowth",
                                                           "KOmetallomics",
                                                           "OEmetallomics",
                                                           "y5kmetalspecific",
                                                           "metpertWTproteomics")))

# Plot with ggplot
np_dshit_plot <-ggplot(filter(ds_hit_KEGGpathways,metal != "Mo"), 
       aes(x = metal,
           y = num_pathways,
           group = dataset,
           fill = dataset)) +
  geom_bar(aes(colour  = dataset),
           stat = "identity", 
           position = "dodge",
           width = 0.6) +
  geom_text(aes(y = 49,
                label = num_pathway_any_dataset)) +
  labs(x = "Metal", y = "Number of Pathways") +
  scale_fill_manual(values = colkey_dataset)+
  scale_colour_manual(values = colkey_dataset)+
  theme_metallica()+
  labs(x = "", y = "number of pathways with hits")

## with flux hits in same plot 

flux_KEGG_hits <- rbind(ds_hit_KEGGpathways,num_pathways_w_fluxchange)

np_dshit_fluxplot <-ggplot(filter(flux_KEGG_hits,metal != "Mo"), 
                       aes(x = metal,
                           y = num_pathways,
                           group = dataset,
                           fill = dataset)) +
  geom_bar(
           stat = "identity", 
           position = "dodge",
           colour = "white",
           linewidth = 0.001,
           width = 0.6) +
  geom_text(aes(y = num_pathway_any_dataset,
                label = num_pathway_any_dataset)) +
  labs(x = "Metal", y = "Number of Pathways") +
  scale_fill_manual(values = colkey_dataset)+
  theme_metallica()+
  labs(x = "", y = "number of pathways with hits")

pdf(paste0(plot_dir,"/num_KEGGpathways_changing_fluxsim_dataset.pdf"), width = 8,height = 6)
np_fluxch_plot
np_dshit_plot
np_dshit_fluxplot
dev.off()

## make venn of responsive pathways based on flux simulations and based on datasets - just for Fe

Fe_kp_ds_hits <- unique(filter(all_datasets_KEGG_metalanno, metal == "Fe", sig_in_pathway)$term)

Fe_kp_flux_hits <- unique(filter(flux_change_pathwaywise_metalwise_smry, 
              metal == "Fe", flux_change_detected_inpathway)$KEGG_pathwayName)

Fe_kp_ds_hits[which(!Fe_kp_ds_hits %in% Fe_kp_flux_hits)]
## 

# plot
hit_inds_plot <- ggplot(all_datasets_KEGG_metalanno_smry, 
                           aes(x = metal,
                               y = term,
                               fill = num_ds_sig_in)) + 
                      facet_wrap("metal_anno")+
                      geom_tile(colour ="black",linewidth = 0.05, alpha = 0.75)+
                      scale_fill_viridis(direction = -1,option = "C")+
                      theme_metallica()+
                      theme(axis.text.y = element_text(size = 13))+
                      labs(x="", y = "", fill = "")

pdf(paste0(plot_dir,"/KEGGpathwaywise_flux_datasethits_comparison.pdf"),width = 12,height = 13)
simflux_plot
hit_inds_plot
dev.off()

############################################################################################################################################
### Make a list of UniprotIDs that correspond to reactions that were identified as responsive in our dataset but not in flux simulations ###
############################################################################################################################################

## get grRules from the GEM simulations

ORF_fluxchange_map <- GEMsim_0pt9grwth_fluxes%>%
                      dplyr::select(metal,grRules,flux_change_detected)%>%
                      unique()

# Split rows where grRules contains 'and'
ORF_fluxchange_map <- ORF_fluxchange_map %>%
  separate_rows(grRules, sep = " and ")

# Remove leading and trailing whitespaces
ORF_fluxchange_map$grRules <- trimws(ORF_fluxchange_map$grRules)

colnames(ORF_fluxchange_map)[which(colnames(ORF_fluxchange_map)=="grRules")] <- "ORF"
ORF_fluxchange_map <- unique(ORF_fluxchange_map)  

ORF_fluxchange_binary_map <- ORF_fluxchange_map%>%
                      group_by(metal,ORF)%>%
                      mutate(flux_change_detected = any(flux_change_detected,na.rm = T))%>%
                      ungroup()%>%
                      unique()
####################################################################
### save ORF flux change binary map for visualisation with iPATH ###
####################################################################

ORF_fluxchange_map_binary_foriPATH <- ORF_fluxchange_binary_map%>%
                                      mutate(Uniprot_ID = as.character(lapply(ORF,convert_ORF2Uniprot)))%>%
                                      mutate(colour = viridis(1),
                                             width = "W12",
                                             opacity = 1,
                                             UNIPROT.ID = paste0("UNIPROT:",Uniprot_ID))%>%
                                      dplyr::select(UNIPROT.ID,colour,width,opacity)

write_delim(ORF_fluxchange_map_binary_foriPATH,
            paste0(output_tables_dir,"/foriPATH_sim_fluxchangedetected_binary.tsv"),col_names = F)

##################################################################################################
### save ORF flux change number of metals as width and colour map for visualisation with iPATH ###
##################################################################################################

ORF_fluxchange_nummetal_map <- ORF_fluxchange_map%>%
  group_by(ORF)%>%
  summarise(num_metals_fluxchange = sum(flux_change_detected,na.rm=T))%>%
  ungroup()
  
## Plot new rxn changes detected using iPATH

convert_numeles2colour <- function(x){
  if(x>0){
  col = viridis(8,option = "C",end = 0.9 ,direction = -1)[x]
  }else{
      col = "#808080"
  }
  return(col)
}

ORF_fluxchange_map_nummetal_foriPATH <- ORF_fluxchange_nummetal_map%>%
                                        mutate(Uniprot_ID = as.character(lapply(ORF,convert_ORF2Uniprot)))%>%
                                        mutate(colour = as.character(lapply(num_metals_fluxchange,convert_numeles2colour)),
                                               width = paste0("W",num_metals_fluxchange+6),
                                               opacity = 1,
                                               UNIPROT.ID = paste0("UNIPROT:",Uniprot_ID))%>%
                                        dplyr::select(UNIPROT.ID,colour,width,opacity)

write_delim(ORF_fluxchange_map_nummetal_foriPATH,
            paste0(output_tables_dir,"/foriPATH_sim_fluxchangedetected_nummetals.tsv"),col_names = F)

#############################################################################
### save ORF flux changes detected in any experimental dataset with iPATH ###
#############################################################################

ORF_dshits_map <- all_datasets%>%
  filter(metal != "Mo")%>%
  dplyr::select(metal,ORF,Significant)%>%
  group_by(metal,ORF)%>%
  mutate(Significant = any(Significant,na.rm=T))%>%
  ungroup()


ORF_changedetect_nummetals_map_foriPATH <- ORF_dshits_map%>%
  unique()%>%
  group_by(ORF)%>%
  summarise(num_metals = sum(Significant,na.rm=T))%>%
  ungroup()%>%
  filter(num_metals != 0)%>%
  mutate(Uniprot_ID = as.character(lapply(ORF,convert_ORF2Uniprot)))%>%
  mutate(colour = as.character(lapply(num_metals,convert_numeles2colour)),
         width = paste0("W",num_metals+6),
         opacity = 1,
         UNIPROT.ID = paste0("UNIPROT:",Uniprot_ID))%>%
  dplyr::select(UNIPROT.ID,colour,width,opacity)

write_delim(ORF_changedetect_nummetals_map_foriPATH,
            paste0(output_tables_dir,"/foriPATH_exp_changedetected_nummetals.tsv"),col_names = F)

                  

#######################################################################################################
### save ORF flux changes detected in metpertWTproteomics data only experimental dataset with iPATH ###
#######################################################################################################

ORF_dshits_map <- all_datasets%>%
  filter(metal != "Mo" & dataset == "metpertWTproteomics")%>%
  dplyr::select(metal,ORF,Significant)%>%
  group_by(metal,ORF)%>%
  mutate(Significant = any(Significant,na.rm=T))%>%
  ungroup()


ORF_changedetect_nummetals_map_foriPATH <- ORF_dshits_map%>%
  unique()%>%
  group_by(ORF)%>%
  summarise(num_metals = sum(Significant,na.rm=T))%>%
  ungroup()%>%
  filter(num_metals != 0)%>%
  mutate(Uniprot_ID = as.character(lapply(ORF,convert_ORF2Uniprot)))%>%
  mutate(colour = as.character(lapply(num_metals,convert_numeles2colour)),
         width = paste0("W",num_metals+6),
         opacity = 1,
         UNIPROT.ID = paste0("UNIPROT:",Uniprot_ID))%>%
  dplyr::select(UNIPROT.ID,colour,width,opacity)

write_delim(ORF_changedetect_nummetals_map_foriPATH,
            paste0(output_tables_dir,"/foriPATH_expmetperWTproteomics_changedetected_nummetals.tsv"),col_names = F)


sim_vs_exp_binary_diff_metalwise <- merge(ORF_fluxchange_map,ORF_dshits_map,
                                    by = c("metal","ORF"), all.y = T)%>%
                                    # filter to keep only those ORFs that are present in the GEM
                                    filter(ORF %in% ORF_fluxchange_map$ORF)%>%
                                    mutate( flux_change_detected = ifelse(is.na(flux_change_detected),F,flux_change_detected),
                                            new_in_ds = ifelse(Significant & !flux_change_detected,T,F))%>%
                                    dplyr::select(metal,ORF,new_in_ds)%>%
                                    unique()

sim_vs_exp_binary_diff_anymetal <- sim_vs_exp_binary_diff_metalwise%>%
                                   group_by(ORF)%>%
                                   mutate(num_metals = sum(new_in_ds))%>%
                                   ungroup()%>%
                                   dplyr::select(ORF, num_metals)%>%
                                   filter(num_metals > 0)%>%
                                   unique()



iPath_flux_ds_diff_colnummet <- sim_vs_exp_binary_diff_anymetal%>%
                        mutate(Uniprot_ID = as.character(lapply(ORF,convert_ORF2Uniprot)))%>%
                        mutate(colour = as.character(lapply(num_metals,convert_numeles2colour)),
                               width = "W12",
                               opacity = 1,
                               UNIPROT.ID = paste0("UNIPROT:",Uniprot_ID))%>%
                        dplyr::select(UNIPROT.ID,colour,width,opacity)

write_delim(iPath_flux_ds_diff_colnummet,paste0(output_tables_dir,"/foriPATH_flux_vs_dsets_bindiff_col_by_nummetals.tsv"),col_names = F)


df <- data.frame(num_eles_hit_in = 1:8, 
                 color = sapply(1:8, convert_numeles2colour))

# Plot color key with ggplot2
pdf(paste0(plot_dir,"/colour_key_iPath_numeles_hit_in.pdf"),width = 4,height = 1)
ggplot(df, aes(x = num_eles_hit_in, y = 1, fill = color)) +
  geom_tile(color = NA) +
  scale_fill_identity() +
  labs(x = NULL, y = NULL) +
  theme_metallica()+
  theme(panel.grid = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "none") +
  labs(title = "colour key iPATH - number of metals ORF is hit in")
dev.off()

############################################################################
### prepare dataset with numds each ORF is hit in for plotting with PATH ###  
############################################################################

convert_numdshits2colour <- function(x){
  viridis(6,option = "C",direction = -1)[x+1]
}

iPath_numds_hitinORF <- all_datasets%>%
                        mutate(hit_in_metal_dataset = ifelse(is.na(Significant),F,Significant))%>%
                        group_by(ORF, dataset)%>%
                        mutate(hit_in_dataset = any(hit_in_metal_dataset))%>%
                        ungroup()%>%
                        dplyr::select(ORF,hit_in_dataset,dataset)%>%
                        unique()%>%
                        group_by(ORF)%>%
                        summarise(num_ds_hit_in = sum(hit_in_dataset))%>%
                        ungroup()%>%
                        mutate(Uniprot_ID = as.character(lapply(ORF,convert_ORF2Uniprot)))%>%
                        mutate(colour = as.character(lapply(num_ds_hit_in,convert_numdshits2colour)),
                               width = ifelse(num_ds_hit_in == 0, "W5","W12"),
                               opacity = 1,
                               UNIPROT.ID = paste0("UNIPROT:",Uniprot_ID))%>%
                        dplyr::select(UNIPROT.ID,colour,width,opacity)

write_delim(iPath_numds_hitinORF, paste0(output_tables_dir,"/foriPATH_num_ds_ORF_is_hit_in.tsv"),col_names = F)

## metal wise

metals = unique(all_datasets$metal)
for(i in 1:length(metals)){
iPath_numds_hitinORF <- all_datasets%>%
  mutate(hit_in_metal_dataset = ifelse(is.na(Significant),F,Significant))%>%
  group_by(ORF, dataset)%>%
  mutate(hit_in_dataset = any(hit_in_metal_dataset))%>%
  ungroup()%>%
  dplyr::select(ORF,hit_in_dataset,dataset)%>%
  unique()%>%
  group_by(ORF)%>%
  summarise(num_ds_hit_in = sum(hit_in_dataset))%>%
  ungroup()%>%
  mutate(Uniprot_ID = as.character(lapply(ORF,convert_ORF2Uniprot)))%>%
  mutate(colour = as.character(lapply(num_ds_hit_in,convert_numdshits2colour)),
         width = ifelse(num_ds_hit_in == 0, "W5","W12"),
         opacity = 1,
         UNIPROT.ID = paste0("UNIPROT:",Uniprot_ID))%>%
  dplyr::select(UNIPROT.ID,colour,width,opacity)

}
write_delim(iPath_numds_hitinORF, paste0(output_tables_dir,"/foriPATH_num_ds_ORF_is_hit_in.tsv"),col_names = F)



# save colour key 


df <- data.frame(num_da_hit_in = 0:5, 
                 color = sapply(0:5, convert_numdshits2colour))

# Plot color key with ggplot2
pdf(paste0(plot_dir,"/colour_key_iPath_numds_hit_in.pdf"),width = 3,height = 1)
ggplot(df, aes(x = num_da_hit_in, y = 1, fill = color)) +
  geom_tile(color = NA) +
  scale_fill_identity() +
  labs(x = NULL, y = NULL) +
  theme_metallica()+
  theme(panel.grid = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "none") +
  labs(title = "colour key iPATH - number of dataset ORF is hit in")
dev.off()

