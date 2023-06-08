
#`---
#`  Title: "metallica_dataset_comparison.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 25 May 2023
#`  Description: scompares results of metpertWTproteomics, metdepKOgrowth, KOmetallomics, OEmetallomics and y5kmetspecific
#`---


#############################################
### source paths functions and libraries  ###
#############################################

library(ggplot2)
library(tidyr)
library(dplyr)
library(UpSetR)

# general
source("/Users/aulakhs/Documents/Ralser Lab/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,'/common_code/graphics_parameters.R'))

source(paste0(code_dir,'/common_code/database_identifier_conversion_functions.R'))
source(paste0(code_dir,'/common_code/input_processed_databases_publisheddatasets.R'))

plot_dir = paste0( datasetcomparison_dir, "/output/plots")
dir.create(plot_dir,recursive = T)
output_tables_dir = paste0( datasetcomparison_dir, "/output/tables")
dir.create(output_tables_dir,recursive = T)

###################################################################
### Set Fold Change and Adj PValue threshold for entire script ###
###################################################################

fcthresh <- 1.5  # Magnitude of change that will be considered significant for DE
adjpvthresh <- 0.05

##########################
### Input all datasets ###
##########################

## metal perturbation proteomics results  

PQ_lm_res_df <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/lmfit DE res with PQ and SigNotSig AdjPVthresh ",
                               adjpvthresh," FCthresh ",fcthresh,".csv"),
                        stringsAsFactors = F)[,c("Element","Genes","ORF","Significant")]

PQ_lm_res_sigstat_df <- PQ_lm_res_df%>%
                    filter(Significant == 1)%>%
                    dplyr::select(Element, ORF)%>%
                    unique()
                    
colnames(PQ_lm_res_sigstat_df) <- c("metal","ORF")

metpert_measured <- unique(PQ_lm_res_df$ORF)

met_prot_corr_df <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/all_metal_protein_correlations.csv"),stringsAsFactors = F)

met_prot_corr_sigcor_df <- filter(met_prot_corr_df,
                       pvalue_pearson < 0.05, abs(r2_pearson) > 0.5)%>%
  dplyr::select(Element_Measured, ORF)%>%
  unique()

colnames(met_prot_corr_sigcor_df) <- c("metal","ORF")

### hits from met pert proteomics to validate

methits_metpertprot <- rbind( PQ_lm_res_sigstat_df,met_prot_corr_sigcor_df)%>%
                        unique()%>%
                        mutate(Protein = as.character(lapply(ORF,convert_ORF2Uniprot)))

methits_metpertprot_smry <- methits_metpertprot%>%
                            group_by(metal)%>%
                            summarise(num_hits = length(Protein))
                            
#################################################################
### Input all genome wide datasets to be used for validation  ###
#################################################################

y5k_measured <- as.character(lapply(unique(read.csv(paste0(published_dataset_dir,"/y5k_proteomics/y5k_metalrelatedKOs_metalspecific.csv"),stringsAsFactors = F)$Protein.Group),
                       convert_Uniprot2singleORF))
                


y5k_da_in_metrelKOs_metalspecific <- read.csv(paste0(published_dataset_dir,"/y5k_proteomics/y5k_metalrelatedKOs_metalspecific.csv"),stringsAsFactors = F)%>%
                                     mutate(Significant_y5kmetspecific = ifelse(p_value < 0.05 & abs(log2FC) > log2(1.5),T,F),
                                             p_values_y5k = p_value,
                                             log2prot_y5k = log2(Protein.Quantity),
                                             KO_y5k = KO,
                                             log2FC_y5k = log2FC,
                                              metal = term)%>%
                                      dplyr::select(-c(p_value,Protein.Quantity,KO,log2FC,term))%>%
                                      na.omit()
                     
y5k_da_in_metrelKOs_metalspecific_smry <- y5k_da_in_metrelKOs_metalspecific%>%
                                          dplyr::select(Protein.Group,metal,Significant_y5kmetspecific)%>%
                                          group_by(Protein.Group,metal)%>%
                                          mutate(Significant_y5kmetspecific = ifelse(any(Significant_y5kmetspecific),T,F))%>%
                                          ungroup()%>%
                                          unique()

  
y5kmetspecific_numhits_summary <- read.csv(paste0(published_dataset_dir,"/y5k_proteomics/y5kmetspecific_numhits_summary.csv"),stringsAsFactors = F)
colnames(y5kmetspecific_numhits_summary) <- c("metal","sig_in_metal")

## phenotypic ( growth ) screen on metal depletion conditions of the BY4741 + pHLUM KO library

KOgrowth_measured <- unique(read.csv(paste0(metdep_KOgrowth_dir,"/output/tables/metdepKOgrowth_statresults_annotated.csv"),stringsAsFactors = F)$ORF)

KOgrowthscreen <- read.csv(paste0(metdep_KOgrowth_dir,"/output/tables/metdepKOgrowth_statresults_annotated.csv"),stringsAsFactors = F)
                 

KOgrowthscreen_smry <- KOgrowthscreen%>%
                       mutate(Significant_KOgrowth = sig_in_metal)%>%
                       dplyr::select(ORF,metal,Significant_KOgrowth)%>%
                       unique()
                      

metdepKOgrowth_numhits_summary <- read.csv(paste0(metdep_KOgrowth_dir,"/output/tables/metdepKOgrowth_numhits_summary.csv"),stringsAsFactors = F)                  

## metallomics screen on metal depletion conditions of the BY4741 KO library

KOmetallomics_measured <- unique(read.csv(paste0(published_dataset_dir,"/metallomics/KOmetallomics_zscores_Iacovacci2020.csv"),stringsAsFactors = F,strip.white = T)$ORF)

KOmetallomics <- read.csv(paste0(published_dataset_dir,"/metallomics/KOmetallomics_zscores_Iacovacci2020.csv"),stringsAsFactors = F,strip.white = T)%>%
                 reshape2::melt(id.vars = "ORF", variable.name = "metal", value.name = "Zscore_KOmetallomics")%>%
                 filter(metal %in% c("Ca","Cu","Fe","K","Mg","Mn","Mo","Na","Zn"))%>%
                 mutate(Significant_KOmetallomics = ifelse(abs(Zscore_KOmetallomics) > 1.959, T,F),
                        ORF_metal = paste(ORF, metal))
                 

KOmetallomics_numhits_summary <- KOmetallomics%>%
                                dplyr::select(ORF,metal,Significant_KOmetallomics)%>%
                                unique()%>%
                                group_by(metal)%>%
                                summarize(sig_in_metal = sum(Significant_KOmetallomics, na.rm=T))

## metallomics screen on metal depletion conditions of the OE library ( check the exact name ? whats the BY number)

OEmetallomics_measured <- unique(read.csv(paste0(published_dataset_dir,"/metallomics/OEmetallomics_zscores_Iacovacci2020.csv"),stringsAsFactors = F,strip.white = T)$ORF)

OEmetallomics <- read.csv(paste0(published_dataset_dir,"/metallomics/OEmetallomics_zscores_Iacovacci2020.csv"),stringsAsFactors = F,strip.white = T)%>%
                 reshape2::melt(id.vars = "ORF", variable.name = "metal", value.name = "Zscore_OEmetallomics")%>%
                 filter(metal %in% c("Ca","Cu","Fe","K","Mg","Mn","Mo","Na","Zn"))%>%
                 mutate(Significant_OEmetallomics = ifelse(abs(Zscore_OEmetallomics) > 1.959, T,F),
                        ORF_metal = paste(ORF, metal))

OEmetallomics_numhits_summary <- OEmetallomics%>%
                                dplyr::select(ORF,metal,Significant_OEmetallomics)%>%
                                unique()%>%
                                group_by(metal)%>%
                                summarize(sig_in_metal = sum(Significant_OEmetallomics, na.rm=T))


###########################################################
### Plot number of hits across all datasets in one plot ###
###########################################################

metdepKOgrowth_numhits_summary$dataset <- "metdepKOgrowth"
KOmetallomics_numhits_summary$dataset <- "KOmetallomics"
OEmetallomics_numhits_summary$dataset <- "OEmetallomics"
y5kmetspecific_numhits_summary$dataset <- "y5kmetspecific"

num_hits_all_datasets <- rbind(metdepKOgrowth_numhits_summary,KOmetallomics_numhits_summary,OEmetallomics_numhits_summary,y5kmetspecific_numhits_summary)
num_hits_all_datasets$dataset <- factor( num_hits_all_datasets$dataset,levels = c("metdepKOgrowth","KOmetallomics","OEmetallomics","y5kmetspecific"))

pdf(paste0(plot_dir,"/num_hits_across_all_val_datasets.pdf"),width = 8,height = 5)
ggplot(filter(num_hits_all_datasets, sig_in_metal > 0),
       aes(x = metal,
           y = sig_in_metal,
           colour = metal,
           fill = metal)) +
  geom_point(aes(shape = dataset), size = 2.8, alpha = 0.8, position = position_jitterdodge(jitter.width = 0.8)) +
  scale_y_log10() +
  scale_colour_manual(values = colkey_Ele) +
  scale_fill_manual(values = colkey_Ele) +
  scale_shape_manual(values = c(17, 23, 21, 8)) +
  labs(x = "", y = "number signficant",colour = "", fill = "",shape = "")+
  theme_metallica()
dev.off()


#############################################
### merge all datasets into one dataframe ### -- only keeping hits of metpertprot
#############################################

metprot_validation_data <- merge(methits_metpertprot,y5k_da_in_metrelKOs_metalspecific,by.x = c("Protein","metal"), 
                           by.y = c("Protein.Group","metal"), all.x =T)

metprot_validation_data<- merge(metprot_validation_data,KOgrowthscreen_smry, by = c("ORF","metal"), all.x = T)

metprot_validation_data<- merge(metprot_validation_data,KOmetallomics, by = c("ORF","metal"), all.x = T)

metprot_validation_data<- merge(metprot_validation_data,OEmetallomics, by = c("ORF","metal"), all.x = T)


##############################################################################
### Note down which datasets are the metprotproteomics hits significant in ###
##############################################################################

metprot_validation_forsummary <- merge(methits_metpertprot,y5k_da_in_metrelKOs_metalspecific_smry,by.x = c("Protein","metal"), 
                            by.y =c("Protein.Group","metal"), all.x =T)

metprot_validation_forsummary<- merge(metprot_validation_forsummary,KOgrowthscreen_smry, by = c("ORF","metal"), all.x = T)

metprot_validation_forsummary<- merge(metprot_validation_forsummary,KOmetallomics[,c("ORF","metal","Significant_KOmetallomics")], by = c("ORF","metal"), all.x = T)

metprot_validation_forsummary<- merge(metprot_validation_forsummary,OEmetallomics[,c("ORF","metal","Significant_OEmetallomics")], by = c("ORF","metal"), all.x = T)

### Count number of datasets that validate each hit ###

metprot_validation_forsummary[is.na(metprot_validation_forsummary)] <- 0

metprot_validation_forsummary$Num_ds_validated <- metprot_validation_forsummary$Significant_y5kmetspecific+
                                                  metprot_validation_forsummary$Significant_KOgrowth+
                                                  metprot_validation_forsummary$Significant_KOmetallomics
                                                  

metprot_validation_summary <- metprot_validation_forsummary%>%
                              group_by(metal)%>%
                              summarise(
                                        y5kmetspec = sum(Significant_y5kmetspecific,na.rm=T),
                                        KOgrowth = sum(Significant_KOgrowth,na.rm=T),
                                        KOmetallomics = sum(Significant_KOmetallomics,na.rm=T),
                                        OEmetallomics = sum(Significant_OEmetallomics,na.rm=T)
                                        )%>%
                              ungroup()


metprot_validation_valbyany <- metprot_validation_forsummary%>%
                               reshape2::melt(id.vars= c("ORF","metal","Protein"),
                                              variable.name = "validation_dataset",
                                              value.name = "validated")%>%
                               na.omit()%>%
                               group_by(ORF,metal,Protein)%>%
                               summarise(validated_by_any = ifelse(any(validated),T,F))%>%
                               ungroup()

metprot_validation_valbyany_smry <- metprot_validation_valbyany%>%
                                    group_by(metal)%>%
                                    summarise(num_validated = sum(validated_by_any),
                                              total_metprothits = length(validated_by_any),
                                              fraction_validated_by_any = num_validated/total_metprothits)

pdf(paste0(crossds_validation_dir,"/fraction_hits_validated_by_any.pdf"),width = 10,height = 6)
ggplot(metprot_validation_valbyany_smry)+
  geom_bar(aes(x = metal,
               y = total_metprothits),
           stat = "identity",
           colour = "#191970",
           width = 0.5,
           alpha = 0)+
  geom_bar(aes(x = metal,
               y = num_validated),
           stat = "identity",
           width = 0.5,
           colour = "#191970",
           fill = "#191970")+
  geom_text(aes(x = metal,
                y = total_metprothits+20,
                label = paste0(round(fraction_validated_by_any,2)*100,"%")),
            size = 6)+
  theme_metallica()+
  labs(x = "",
       y = "border:total hits\nfill:validated hits",
       title = "differential abundance analysis and correlation analysis")
    
dev.off()


##########################
### Venn / Upset plots ###
##########################

## all measured 

measured_across_ds <- list(metpert_measured,
                       y5k_measured,
                       KOgrowth_measured,
                       KOmetallomics_measured,
                       OEmetallomics_measured)

names(measured_across_ds) <- c( "metal_pert",
                                "y5k",
                                "KOgrowth",
                                "KOmetallomics",
                                "OEmetallomics")

pdf(paste0(plot_dir,"/venn_measured_acrooss_datasets.pdf"),width = 10, height = 10)
venn::venn(measured_across_ds,zcolor = "style",box = F, ilcs=2.5) 
dev.off()


pdf(paste0(plot_dir,"/upset_all_measured_across_datasets.pdf"),width = 8, height = 4)

upset(fromList(measured_across_ds), nsets = 6, order.by = "freq",
      cutoff = 5,
      matrix.color = "black", 
      main.bar.color = "#191970",
      sets.bar.color = "#191970")


dev.off()

## all hits

valhits <- list(
                unique(filter(metprot_validation_forsummary,as.logical(Significant_y5kmetspecific))$ORF),
                unique(filter(metprot_validation_forsummary,as.logical(Significant_KOgrowth))$ORF),
                unique(filter(metprot_validation_forsummary,as.logical(Significant_KOmetallomics))$ORF),
                unique(filter(metprot_validation_forsummary,as.logical(Significant_OEmetallomics))$ORF)
)
names(valhits) <- c( 
                    "y5k_metspecific",
                    "KOgrowth",
                    "KOmetallomics",
                    "OEmetallomics")

pdf(paste0(crossds_validation_dir,"/venn_cross_dataset_validation.pdf"),width = 10, height = 10)
venn::venn(valhits,zcolor = "style",box = F, ilcs=2.5) 
dev.off()

intersect(intersect(valhits$y5k_metspecific,valhits$KOgrowth),valhits$KOmetallomics)

valhits_with_metpert <-  list(unique(methits_metpertprot$ORF),
                              unique(filter(metprot_validation_forsummary,as.logical(Significant_y5kmetspecific))$ORF),
                              unique(filter(metprot_validation_forsummary,as.logical(Significant_KOgrowth))$ORF),
                              unique(filter(metprot_validation_forsummary,as.logical(Significant_KOmetallomics))$ORF),
                              unique(filter(metprot_validation_forsummary,as.logical(Significant_OEmetallomics))$ORF)
                            )
names(valhits_with_metpert) <- c( 
                    "metpert",
                    "y5k_metspecific",
                    "KOgrowth",
                    "KOmetallomics",
                    "OEmetallomics")

pdf(paste0(plot_dir,"/upset_all_significant_across_datasets.pdf"),width = 8, height = 4)

upset(fromList(valhits_with_metpert), nsets = 6, order.by = "freq",
      cutoff = 6,
      matrix.color = "black", 
      main.bar.color = "#191970",
      sets.bar.color = "#191970")

dev.off()

#############################################
### metal wise upset plots for validation ###
#############################################

eles = unique(metprot_validation_forsummary$metal)

pdf(paste0(plot_dir,"/upset_all_significant_across_datasets_metalwise.pdf"),width = 8, height = 4)

for(e in 1:length(eles)){
  
  methits_metpertprot_ele = filter(methits_metpertprot, metal == eles[e])
  metprot_validation_forsummary_ele = filter(metprot_validation_forsummary, metal == eles[e])
  
  valhits_with_metpert_ele <-  list(unique(methits_metpertprot_ele$ORF),
                                unique(filter(metprot_validation_forsummary_ele,as.logical(Significant_y5kmetspecific))$ORF),
                                unique(filter(metprot_validation_forsummary_ele,as.logical(Significant_KOgrowth))$ORF),
                                unique(filter(metprot_validation_forsummary_ele,as.logical(Significant_KOmetallomics))$ORF),
                                unique(filter(metprot_validation_forsummary_ele,as.logical(Significant_OEmetallomics))$ORF)
  )
  names(valhits_with_metpert_ele) <- paste(eles[e],c( 
                                    "metpert",
                                    "y5k_metspecific",
                                    "KOgrowth",
                                    "KOmetallomics",
                                    "OEmetallomics"))
  
 
  
print(
  upset(fromList(valhits_with_metpert_ele), nsets = 6, order.by = "freq",
        cutoff = 6,
        matrix.color = "black", 
        main.bar.color = colkey_Ele[eles[e]],
        sets.bar.color =  colkey_Ele[eles[e]])
)
  
}
dev.off()


################################
### Validated by > 1 dataset ###
################################

val_morethan1ds <- filter(val_morethan1ds, Num_ds_validated >1)

val_morethan1ds <- merge(val_morethan1ds,GenProt_SGD[,c("ORF","Gene.Name","Uniprot.Annotation.Score")], 
                         by = c("ORF"))
write.csv(val_morethan1ds,
          paste0(output_tables_dir,"/metperhits_validated_by_morethan1ds.csv"),row.names = F)



# summarise validated by > 1 ds

validated_bymorethan1_smry <- filter(metprot_validation_forsummary, Num_ds_validated>1)%>%
                              group_by(metal)%>%
                              summarize(num_orfs_validated = length(ORF))%>%
                              ungroup()


write.csv(validated_bymorethan1_smry,
          paste0(output_tables_dir,"/metperthits_validated_by_morethan1ds_metalwisesummary.csv"),row.names = F)


#############################################
### cdf plots all dataset hits ### --  keeping hits of all studies
#############################################

y5k_da_in_metrelKOs_metalspecific_pgwise <- y5k_da_in_metrelKOs_metalspecific%>%
                                            group_by(metal,Protein.Group,measured_protein_name)%>%
                                            summarise(Significant_y5kmetspecific = ifelse(any(Significant_y5kmetspecific),T,F))%>%
                                            ungroup()%>%
                                            mutate(Measured_ORF = as.character(lapply(Protein.Group, convert_Uniprot2singleORF)))

methits_metpertprot$Significant_metpertWTproteomics = T


meas_in_any_ds <- unique(unlist(measured_across_ds))## ORFs that were measured measured in at least one screen 

## note what % of metal binders were measured in atleast one screen 

length(intersect(unique(met_specific_metal_binders$ORF), meas_in_any_ds)) / length(unique(met_specific_metal_binders$ORF))

metprot_validation_data_alldshits <- merge(filter(met_specific_metal_binders,
                                                  ORF %in% meas_in_any_ds), methits_metpertprot, 
                                           by. = c("ORF","term") ,
                                           by.y = c("ORF","metal"),
                                           all.x = T)

metprot_validation_data_alldshits <- merge(metprot_validation_data_alldshits,y5k_da_in_metrelKOs_metalspecific_pgwise[c("Measured_ORF","metal","Significant_y5kmetspecific")],
                                           by.x = c("ORF","term"), 
                                           by.y = c("Measured_ORF","metal"), all.x =T)



metprot_validation_data_alldshits <- merge(metprot_validation_data_alldshits,KOgrowthscreen_smry, 
                                           by.x = c("ORF","term"),
                                          by.y = c("ORF","metal"),
                                           all.x = T)

metprot_validation_data_alldshits <- merge(metprot_validation_data_alldshits,KOmetallomics[,c("ORF","metal","Significant_KOmetallomics")],
                                           by.x = c("ORF","term"),
                                           by.y = c("ORF","metal"),
                                           all.x = T)

metprot_validation_data_alldshits <- merge(metprot_validation_data_alldshits,OEmetallomics[,c("ORF","metal","Significant_OEmetallomics")], 
                                           by.x = c("ORF","term"),
                                           by.y = c("ORF","metal"),
                                           all.x = T)

recovery_anno <- metprot_validation_data_alldshits%>%
            reshape2::melt(id.vars = c("ORF","term","Protein"),
                           variable.name = "dataset", 
                           value.name = "recovered")%>%
            mutate(metal = term)%>%
            group_by(ORF, metal)%>%
            mutate(recovered = ifelse(is.na(recovered), F, recovered),
                   recovered_by_any = any(recovered),
                   recovered_by_any = ifelse(is.na(recovered_by_any), F, recovered_by_any))
## all dataset recovery

recovery_smry <- recovery_anno%>%
                 dplyr::select(ORF, metal, recovered_by_any)%>%
                 unique()%>%
                 group_by(metal)%>%
                 summarize(total_metbinders = length(ORF),
                           num_recovered = sum(recovered_by_any),
                           frac_recovered = num_recovered/total_metbinders)

## ds wise recovery 

recovery_smry_dswise <- recovery_anno%>%
                        mutate(dataset = gsub("Significant_","",dataset),
                               dataset = factor(dataset,
                                                   levels = c("KOgrowth","KOmetallomics","OEmetallomics",
                                                              "metpertWTproteomics","y5kmetspecific")))%>%
                        group_by(metal, dataset)%>%
                        summarize(total_metbinders = length(ORF),
                                  num_recovered = sum(recovered),
                                  frac_recovered = num_recovered/total_metbinders)%>%
                        ungroup()
                        
pdf(paste0(plot_dir,"/metbinders_recovery_all_datasets.pdf"),width = 10,height = 6)
ggplot(recovery_smry_dswise,
       aes(x = dataset,
           y = frac_recovered,
           fill = metal))+
  geom_bar(stat = "identity", position = "dodge", width = 0.5,colour = "white")+
  scale_fill_manual(values = colkey_Ele)+
  scale_colour_manual(values = colkey_Ele)+
  theme_metallica()+
  labs(x = "",
       y = "fraction of metal binders recovered")
dev.off()
## number datasets 
 # 0 - KOgrowth
 # 1 - KO metallomics
 # 2 - OE metallomics
 # 3 - metpertWTproteomics
 # 4 - y5kKOproteomics_metspecific

## count how many known metal binders are recovered cumulatively by each dataset for each metal 

recovery_anno$ORF_metal = paste(recovery_anno$ORF, recovery_anno$metal)


cuml_rec_KOgrowth <- filter(recovery_anno,dataset == "Significant_KOgrowth" & recovered)%>%
                     group_by(metal)%>%
                     mutate(num_cuml_recovered = length(ORF))

cuml_rec_KOmetallomcis <- filter(recovery_anno,dataset == "Significant_KOmetallomics" & recovered &
                                   !ORF_metal %in% unique(cuml_rec_KOgrowth$ORF_metal))%>%
                          group_by(metal)%>%
                          mutate(num_cuml_recovered = length(ORF))
                          
                                  
cuml_rec_OEmetallomcis <- filter(recovery_anno,dataset == "Significant_OEmetallomics" & recovered & 
                                   !ORF_metal %in% c(unique(cuml_rec_KOgrowth$ORF_metal,
                                                    cuml_rec_KOmetallomcis$ORF_metal)))%>%
                          group_by(metal)%>%
                          mutate(num_cuml_recovered = length(ORF))

cuml_rec_metpertWTproteomics <- filter(recovery_anno,dataset == "Significant_metpertWTproteomics" & recovered &
                                         !ORF_metal %in% unique(c(cuml_rec_KOgrowth$ORF_metal,
                                                                  cuml_rec_KOmetallomcis$ORF_metal,
                                                                  cuml_rec_OEmetallomcis$ORF_metal)))%>%
                                  group_by(metal)%>%
                                  mutate(num_cuml_recovered = length(ORF))
  
cuml_rec_y5kmetspecific <- filter(recovery_anno,dataset == "Significant_y5kmetspecific" & recovered &
                                    !ORF_metal %in% unique(c( cuml_rec_KOgrowth$ORF_metal,
                                                              cuml_rec_KOmetallomcis$ORF_metal,
                                                              cuml_rec_OEmetallomcis$ORF_metal,
                                                              cuml_rec_metpertWTproteomics$ORF_metal)))%>%
                                  group_by(metal)%>%
                                  mutate(num_cuml_recovered = length(ORF))



for_ecdf <- merge(rbind(  unique(cuml_rec_KOgrowth[,c("metal","dataset","num_cuml_recovered")]),
                          unique(cuml_rec_KOmetallomcis[,c("metal","dataset","num_cuml_recovered")]),
                          unique(cuml_rec_OEmetallomcis[,c("metal","dataset","num_cuml_recovered")]),
                          unique(cuml_rec_metpertWTproteomics[,c("metal","dataset","num_cuml_recovered")]),
                          unique(cuml_rec_y5kmetspecific[,c("metal","dataset","num_cuml_recovered")])
                        
            ), recovery_smry[,c("metal","total_metbinders")], by = "metal")%>%
            mutate(frac_cuml_recovered = num_cuml_recovered/total_metbinders)%>%
            mutate(dataset = gsub("Significant_","", dataset),
                   dataset = factor(dataset, levels = c("KOgrowth","KOmetallomics","OEmetallomics",
                                                      "metpertWTproteomics","y5kmetspecific")))%>%
            group_by(metal)%>%
            arrange(dataset)%>%
            mutate(cumsum = cumsum(frac_cuml_recovered))

pdf(paste0(plot_dir,"/metbinders_cumfracrecovery_all_datasets.pdf"),width = 6,height = 6)
ggplot(for_ecdf,
       aes(x = dataset,
           y = cumsum,
           colour = factor(metal),
           group = metal))+
  geom_point( alpha = 0.6)+
  geom_line(linewidth = 1.5, alpha = 0.5)+
  scale_colour_manual(values = colkey_Ele)+
  theme_metallica()+
  theme(axis.text.x = element_text(angle = 90))+
  labs(x = "",
       y = "cumulative % recovery of metal binders",
       colour = "")
dev.off()


