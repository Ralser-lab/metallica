
#`---
#`  Title: "y5kmetspecific_analyze_y5k_DA_in_metrelatedKOs.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 18 April 2023
#`  Description: summarises results of the growth screen of y5k KO library on metal depletion media
#`---


#############################################
### source paths functions and libraries  ###
#############################################

# general
source("/Users/aulakhs/Documents/Ralser Lab/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,'/common_code/graphics_parameters.R'))

source(paste0(code_dir,'/common_code/database_identifier_conversion_functions.R'))
source(paste0(code_dir,'/common_code/input_processed_databases_publisheddatasets.R'))

#################################################################
### Get y5kmetspecificDA specific functions and libraries from  ###
#################################################################

source(paste0(code_dir,'/published_datasets/publisheddatasets_0_libraries_functions.R'))

#################
### Set Paths ###
################# 

plot_dir <- paste0(metpert_WTproteomics_dir,"/output/plots/y5k_metalspecific_DA")
dir.create(plot_dir, recursive = T)

output_tables_dir <- paste0(metpert_WTproteomics_dir,"/output/tables/y5k_metalspecific_DA")
dir.create(output_tables_dir, recursive = T)


############################################################
### read in  data and filter for metal related knockouts ###
############################################################

GO_gset_MF_binding_metalwise_smry <- GO_gset_MF_binding_metalwise%>%
                                     group_by(term)%>%
                                     summarize(num_ORFs = length(ORF))%>%
                                     ungroup()
                                     

pdf(paste0(plot_dir,"/num_ORFs_in_KOs_permetal.pdf"), width = 7, height = 5)
ggplot(GO_gset_MF_binding_metalwise_smry,
       aes(x = term,
           y = num_ORFs,
           fill = term
          ))+
  geom_bar(stat ="identity",
           width = 0.6)+
  scale_fill_manual(values = colkey_Ele)+
  theme_metallica()+
  labs(x ="",
       y = "number of ORFs annotated for \nmetal binding in GOMF",
       fill = "")
dev.off()

########################################
### Metal specific metal-binding KOs ### 
########################################

                          
y5k_protquants_WT <-  fread(paste0(published_dataset_dir,"/y5k_proteomics/yeast5k_noimpute_wide-2023-01-26_Messner2023.csv.gz"),stringsAsFactors = F)%>%
                                       as.data.frame()%>%
                                       reshape2::melt(id.vars = "Protein.Group",variable.name = "KO",value.name = "Protein.Quantity")%>%
                                       mutate(KO = as.character(KO))%>%
                                       separate(col=KO, into = c(NA,NA,NA,NA,"KO",NA,NA),sep="_",remove=F)%>%
                                       filter(KO == "YOR202W")%>%
                                       group_by(Protein.Group)%>%
                                       summarize(Protein.Quantity_WT = mean(Protein.Quantity,na.rm = T))
                                      


## read in and filter p values and protein quantities of ORFs with metal specific metal annotations

y5k_metalrelated_metalspecific_KOs_pvalues <- read.csv(paste0(published_dataset_dir,"/y5k_proteomics/yeast5k_stat_DE_bySTRAIN_Messner2023.csv"),stringsAsFactors = F)%>%
                                              as.data.frame()%>%
                                              reshape2::melt(id.vars = "Protein.Group",variable.name = "KO",value.name = "p_value")%>%
                                              mutate(KO = as.character(KO))%>%
                                              separate(col=KO, into = c(NA,NA,NA,NA,"KO",NA,NA),sep="_",remove=F)%>%
                                              filter(KO != "qc")%>%
                                              filter(KO %in% c(unique(GO_gset_MF_binding_metalwise$ORF)))%>%
                                              unique()%>%
                                              na.omit()%>%
                                              group_by(KO, Protein.Group)%>%
                                              summarize(p_value = mean(p_value))%>%
                                              ungroup()


y5k_metalrelated_metalspecific_KOs_protquants <-  fread(paste0(published_dataset_dir,"/y5k_proteomics/yeast5k_noimpute_wide-2023-01-26_Messner2023.csv.gz"),stringsAsFactors = F)%>%
                                                    as.data.frame()%>%
                                                    reshape2::melt(id.vars = "Protein.Group",variable.name = "KO",value.name = "Protein.Quantity")%>%
                                                    mutate(KO = as.character(KO))%>%
                                                    separate(col=KO, into = c(NA,NA,NA,NA,"KO",NA,NA),sep="_",remove=F)%>%
                                                    filter(KO != "qc")%>%
                                                    filter(KO %in% c(unique(GO_gset_MF_binding_metalwise$ORF)))%>%
                                                    na.omit()%>%
                                                    group_by(KO, Protein.Group)%>%
                                                    summarize(Protein.Quantity = mean(Protein.Quantity))%>%
                                                    ungroup()


y5k_metalrelated_metalspecific_KOs_protquants <- merge(y5k_metalrelated_metalspecific_KOs_protquants, y5k_protquants_WT, by = "Protein.Group")%>%
                                                 mutate(log2FC = log2(Protein.Quantity/Protein.Quantity_WT))

y5k_metalrelated_metalspecific_KOs <- merge(y5k_metalrelated_metalspecific_KOs_pvalues, 
                                            y5k_metalrelated_metalspecific_KOs_protquants, by = c("KO","Protein.Group"))


### Convert all KO and Protein.Group Uniprot IDs to ProteinNames

y5k_KO_ProtName_map <-  data.frame(KO = unique(y5k_metalrelated_metalspecific_KOs[,c("KO")]))%>%
                        mutate( KO = as.character(KO),
                                deleted_protein_name = as.character(lapply(KO,convert_ORF2SingleGeneName)))

y5k_prmeas_ProtName_map <-data.frame(Protein.Group = unique(y5k_metalrelated_metalspecific_KOs[,c("Protein.Group")]))%>%
                          mutate(Protein.Group = as.character(Protein.Group),
                                 measured_protein_name = as.character(lapply(Protein.Group, convert_Uniprot2SingleGeneName)))

## Merge the y5k data with the protein names for deleted and measured genes

# KO gene name
y5k_metalrelated_metalspecific_KOs <- merge(y5k_metalrelated_metalspecific_KOs,y5k_KO_ProtName_map, by = "KO")

# measured protein name
y5k_metalrelated_metalspecific_KOs <- merge(y5k_metalrelated_metalspecific_KOs,y5k_prmeas_ProtName_map, by = "Protein.Group")

# merge with metal specific annotation
y5k_metalrelated_metalspecific_KOs<- merge(y5k_metalrelated_metalspecific_KOs,GO_gset_MF_binding_metalwise,by.x = "KO", by.y = "ORF")

write.csv(y5k_metalrelated_metalspecific_KOs, paste0(published_dataset_dir,"/y5k_proteomics/y5k_metalrelatedKOs_metalspecific.csv"),row.names = F)


#############################
### Visualise num of hits ###
#############################

y5k_metalrelated_metalspecific_KOs_numhits_smry <- y5k_metalrelated_metalspecific_KOs%>%
                                                   filter(p_value < 0.01 & abs(log2FC) > log2(1.5))%>%
                                                   dplyr::select(Protein.Group,term)%>%
                                                   unique()%>%
                                                   group_by(term)%>%
                                                   summarize(num_significant = length(Protein.Group))
write.csv(y5k_metalrelated_metalspecific_KOs_numhits_smry, paste0(output_tables_dir,"/y5kmetspecific_numhits_summary.csv"),row.names = F)

y5k_metsp_numhits_normalised <- merge(GO_gset_MF_binding_metalwise_smry,y5k_metalrelated_metalspecific_KOs_numhits_smry, by = "term")%>%
                                mutate(num_sig_norm = num_significant/num_ORFs)


pdf(paste0(plot_dir,"/y5k_metrelatedmetspecific_DA.pdf"),width = 7, height = 5)

ggplot(y5k_metalrelated_metalspecific_KOs_numhits_smry,
       aes(x = term,
           y = num_significant,
           colour = term,
           fill = term))+
  geom_bar(stat = "identity",
           width = 0.6)+
  scale_color_manual(values = colkey_Ele)+
  scale_fill_manual(values = colkey_Ele)+
  theme_metallica()

ggplot(y5k_metsp_numhits_normalised,
       aes(x = term,
           y = num_sig_norm,
           colour = term,
           fill = term))+
  geom_bar(stat = "identity",
           width = 0.6)+
  scale_color_manual(values = colkey_Ele)+
  scale_fill_manual(values = colkey_Ele)+
  theme_metallica()+
  labs(x = "",
       y = "num sig DA")
dev.off()

##################################
### Light df for metallica app ###
##################################


for_metallica_app <- y5k_metalrelated_metalspecific_KOs %>% 
                     mutate(KO_gene = paste(deleted_protein_name, KO,sep = "|"),
                            measured_protein = paste(measured_protein_name, Protein.Group,sep = "|"),
                            metal_connec2KO = term,
                            colour = ifelse(p_value < 0.05 & abs(log2FC) > log2(1.5),
                                            as.character(colkey_Ele[metal_connec2KO]),"lightgray"),
                            label = ifelse(p_value < 0.05 & abs(log2FC) > log2(1.5),
                                           measured_protein_name, NA))%>%
                     dplyr::select(KO_gene, measured_protein, p_value, log2FC, metal_connec2KO )%>%
                     unique()


write.csv(for_metallica_app,paste0(metallica_app_dir,"/metallica_app_y5kmetspecificKOs.csv"), row.names = F)


