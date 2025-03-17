#`---
#`  Title: "metallica_dataset_comparison_metalbinders_metaltransporters_acrossdatasets.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 8 Nov 2023
#`  Description: looks at how metal binders and metal transporters change across datasets
#`---

#############################################
### source paths functions and libraries  ###
#############################################

library(tidyr)
library(dplyr)

# general
source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,'/common_code/graphics_parameters.R'))

source(paste0(code_dir,'/common_code/database_identifier_conversion_functions.R'))
source(paste0(code_dir,'/common_code/input_processed_databases_publisheddatasets.R'))
source(paste0(code_dir,"/metpert_WTproteomics/metpertWTproteomics_0_libraries_functions_genesetenrichments.R"))

plot_dir = paste0( datasetcomparison_dir, "/output/plots")
dir.create(plot_dir,recursive = T)
dir.create(paste0(plot_dir,"/HyperGSA"))
output_tables_dir = paste0( datasetcomparison_dir, "/output/tables")
dir.create(output_tables_dir,recursive = T)

##########################
### Input all datasets ###
##########################

## metal perturbation proteomics results  

## lms with extracellular metal concentration

metpertWTproteomics_extracell <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/extracellularconc_vs_proteinabundance/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),
                                          stringsAsFactors = F)


# make a list of all  proteins that were quantified in metpertWTproteomics
metpertWTproteomics_measured <- unique(metpertWTproteomics_extracell$ORF)

## lms with intracellular metal concentration

metpertWTproteomics_intracell <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/intracellularconc_vs_proteinabundance/lmfit_intracellular_metalconc_vs_log2PQ_with_SigNotSig_AdjPVthresh_0.05_FCthresh_1.5.csv"),
                                          stringsAsFactors = F)


###############################################################################
### Input metdepKOgrowth, KOmettalomics, OEmetallomice and y5kmetalspecific ###
###############################################################################

y5kmetalspecific <- read.csv(paste0(published_dataset_dir,"/y5k_proteomics/y5k_metalrelatedKOs_metalspecific.csv"),
                             stringsAsFactors = F)
#  mutate(Significant = ifelse(p_value < 0.05 & abs(log2FC) > log2(1.5),T,F))
   #      metal = term)%>%
  #mutate(ORF = as.character(lapply(Protein.Group, convert_Uniprot2singleORF)),
   #      dataset = "y5kmetalspecific",
    #     dataset_type = "proteomics")

## phenotypic ( growth ) screen on metal depletion conditions of the BY4741 + pHLUM KO library

metdepKOgrowth <- read.csv(paste0(metdep_KOgrowth_dir,"/output/tables/metdepKOgrowth_statresults_annotated.csv"),stringsAsFactors = F)%>%
  mutate(Significant = sig_in_metal,
         dataset = "metdepKOgrowth",
         dataset_type = "KOgrowth")

## metallomics screen on metal depletion conditions of the BY4741 KO library

KOmetallomics <- read.csv(paste0(published_dataset_dir,"/metallomics/KOmetallomics_zscores_Iacovacci2020.csv"),stringsAsFactors = F,strip.white = T)%>%
  reshape2::melt(id.vars = "ORF", variable.name = "metal", value.name = "Zscore_KOmetallomics")%>%
  filter(metal %in% c("Ca","Cu","Fe","K","Mg","Mn","Mo","Na","Zn"))%>%
  mutate(Significant = ifelse(abs(Zscore_KOmetallomics) > 1.959, T,F),
         dataset = "KOmetallomics",
         dataset_type = "metallomics")

## metallomics screen on metal depletion conditions of the OE library ( check the exact name ? whats the BY number)

OEmetallomics <- read.csv(paste0(published_dataset_dir,"/metallomics/OEmetallomics_zscores_Iacovacci2020.csv"),stringsAsFactors = F,strip.white = T)%>%
  reshape2::melt(id.vars = "ORF", variable.name = "metal", value.name = "Zscore_OEmetallomics")%>%
  filter(metal %in% c("Ca","Cu","Fe","K","Mg","Mn","Mo","Na","Zn"))%>%
  mutate(Significant = ifelse(abs(Zscore_OEmetallomics) > 1.959, T,F),
         dataset = "OEmetallomics",
         dataset_type = "metallomics")


############################################################################
### merge metal binder and metal transport annotations with each dataset ###
############################################################################

##################################################
### Along environmental concentration of metal ### 
##################################################

# filter data to keep only metal binding annotated proteins
metb_data <- merge(metpertWTproteomics_extracell,GO_gset_MF_binding_metalwise, by = "ORF")%>%
  mutate(Metal_Match = ifelse(Element == term, "Same Metal", "Different Metal"))%>%
  filter(!Element %in% c("Na","Mo"))

label_data <- metb_data%>%
  group_by(Metal_Match, Element)%>%
  mutate(N_group = length(unique(ORF)),
         label_yloc_group = mean(Log2FC_vs_AE,na.rm=T),
         label_yloc_group = ifelse(Metal_Match == "Same Metal", label_yloc_group-0.2,
                                   label_yloc_group+0.2)
        )%>%
  ungroup()%>%
  dplyr::select(Element,Metal_Match,N_group,label_yloc_group)%>%
  unique()

pdf(paste0(plot_dir,"/metalbinders_along_envmetalconc.pdf"),width = 21.5,height = 3.5)
ggplot(metb_data, aes(x = log2(Element.Concentration),
                          y = Log2FC_vs_AE, 
                      color = Element)) +
  geom_smooth(aes(linetype = Metal_Match), method = "loess",alpha = 0.2) +
  geom_text(data = label_data,
            aes(x = 0, y = label_yloc_group, label = N_group,
                fontface = ifelse(Metal_Match == "Same Metal", "plain", "italic")))+
  theme_minimal() +
  facet_wrap("Element",scales = "free", nrow = 1)+
  geom_hline(yintercept = 0, linewidth = 0.1)+
  geom_vline(xintercept = 0, linewidth = 0.1)+
  labs(color = "", linetype = "")+
  scale_color_manual(values = colkey_Ele)+
  scale_linetype_manual(values = c("Same Metal" = "solid", "Different Metal" = "dashed")) +  
  theme_metallica()
dev.off()

## metal transporters

# filter data to keep only metal binding annotated proteins
mett_data <- merge(metpertWTproteomics_extracell,unique(rbind(GO_gset_MF_transporter_metalwise,
                                                              philpott_metal_transporter_df)), by = "ORF")%>%
  mutate(Metal_Match = ifelse(Element == term, "Same Metal", "Different Metal"))%>%
  filter(!Element %in% c("Mg","Mo"))


label_data <- mett_data%>%
  group_by(Metal_Match, Element)%>%
  mutate(N_group = length(unique(ORF)),
         label_yloc_group = mean(Log2FC_vs_AE,na.rm=T),
         label_yloc_group = ifelse(Metal_Match == "Same Metal", label_yloc_group-0.5,
                                   label_yloc_group+0.5)
  )%>%
  ungroup()%>%
  dplyr::select(Element,Metal_Match,N_group,label_yloc_group)%>%
  unique()

pdf(paste0(plot_dir,"/metaltransporters_along_envmetalconc.pdf"),width = 17.5,height = 3.5)
ggplot(mett_data, aes(x = log2(Element.Concentration),
                      y = Log2FC_vs_AE, 
                      color = Element)) +
  geom_smooth(aes(linetype = Metal_Match), method = "loess",alpha = 0.2) +
  geom_text(data = label_data,
            aes(x = 0, 
                y = label_yloc_group,
                label = N_group,
                fontface = ifelse(Metal_Match == "Same Metal", "plain", "italic")))+
  facet_wrap("Element",scales = "free", nrow = 1)+
  labs(x = "Log2(Environmental Metal Concentration)",
       y = "Log2 Fold Change vs AE", color = "", linetype = "")+
  geom_hline(yintercept = 0, linewidth = 0.1)+
  geom_vline(xintercept = 0, linewidth = 0.1)+
  scale_color_manual(values = colkey_Ele)+
  scale_linetype_manual(values = c("Same Metal" = "solid", "Different Metal" = "dashed")) +  
  theme_metallica()
dev.off()

#############################################
### Along cellular concentration of metal ### 
#############################################

# filter data to keep only metal binding annotated proteins
metb_data <- merge(metpertWTproteomics_intracell,GO_gset_MF_binding_metalwise, by = "ORF")%>%
  mutate(Metal_Match = ifelse(Element == term, "Same Metal", "Different Metal"))%>%
  filter(!Element %in% c("Na","Mo"))

label_data <- metb_data%>%
  group_by(Metal_Match, Element)%>%
  mutate(N_group = length(unique(ORF)),
         label_yloc_group = mean(median_log2_foldchangevsAE,na.rm=T),
         label_yloc_group = ifelse(Metal_Match == "Same Metal", label_yloc_group-0.2,
                                   label_yloc_group+0.2)
  )%>%
  ungroup()%>%
  dplyr::select(Element,Metal_Match,N_group,label_yloc_group)%>%
  unique()

pdf(paste0(plot_dir,"/metalbinders_along_cellmetalconc.pdf"),width = 25,height = 3.5)
ggplot(metb_data, aes(x = log2(median_relative_intracellular_concentration),
                      y = median_log2_foldchangevsAE, 
                      color = Element)) +
  geom_smooth(aes(linetype = Metal_Match), method = "loess",alpha = 0.2) +
  geom_text(data = label_data,
            aes(x = 0, y = label_yloc_group, label = N_group,
                fontface = ifelse(Metal_Match == "Same Metal", "plain", "italic")))+
  theme_minimal() +
  facet_wrap("Element",scales = "free", nrow = 1)+
  geom_hline(yintercept = 0, linewidth = 0.1)+
  geom_vline(xintercept = 0, linewidth = 0.1)+
  labs(color = "", linetype = "")+
  scale_color_manual(values = colkey_Ele)+
  scale_linetype_manual(values = c("Same Metal" = "solid", "Different Metal" = "dashed")) +  
  theme_metallica()
dev.off()

## metal transporters

# filter data to keep only metal binding annotated proteins
mett_data <- merge(metpertWTproteomics_intracell,unique(rbind(GO_gset_MF_transporter_metalwise,
                                                              philpott_metal_transporter_df)), by = "ORF")%>%
  mutate(Metal_Match = ifelse(Element == term, "Same Metal", "Different Metal"))%>%
  filter(!Element %in% c("Cu","Mo","Mg","Na"))


label_data <- mett_data%>%
  group_by(Metal_Match, Element)%>%
  mutate(N_group = length(unique(ORF)),
         label_yloc_group = mean(median_log2_foldchangevsAE,na.rm=T),
         label_yloc_group = ifelse(Metal_Match == "Same Metal", label_yloc_group-0.5,
                                   label_yloc_group+0.5)
  )%>%
  ungroup()%>%
  dplyr::select(Element,Metal_Match,N_group,label_yloc_group)%>%
  unique()

pdf(paste0(plot_dir,"/metaltransporters_along_cellmetalconc.pdf"),width = 21,height = 3.5)
ggplot(mett_data, aes(x = log2(median_relative_intracellular_concentration),
                      y = median_log2_foldchangevsAE, 
                      color = Element)) +
  geom_smooth(aes(linetype = Metal_Match), method = "loess",alpha = 0.2) +
  geom_text(data = label_data,
            aes(x = 0, 
                y = label_yloc_group,
                label = N_group,
                fontface = ifelse(Metal_Match == "Same Metal", "plain", "italic")))+
  facet_wrap("Element",scales = "free", nrow = 1)+
  labs(x = "Log2(Cellular Metal Concentration)",
       y = "Log2 Fold Change vs AE", color = "", linetype = "")+
  geom_hline(yintercept = 0, linewidth = 0.1)+
  geom_vline(xintercept = 0, linewidth = 0.1)+
  scale_color_manual(values = colkey_Ele)+
  scale_linetype_manual(values = c("Same Metal" = "solid", "Different Metal" = "dashed")) +  
  theme_metallica()
dev.off()

############################
### metallomics datasets ###
############################


### KO dataset 

metb_KO_data <- merge(KOmetallomics,GO_gset_MF_binding_metalwise, by = "ORF")%>%
             filter(metal == term)%>%
             filter(metal != "K")%>%
             mutate(Gene_Name = as.character(lapply(ORF,convert_ORF2SingleGeneName)),
                    dataset = "KOm")



mett_KO_data <- merge(KOmetallomics,unique(rbind(GO_gset_MF_transporter_metalwise,
                                              philpott_metal_transporter_df)), by = "ORF")%>%
              filter(metal == term )%>%
              mutate(Gene_Name = as.character(lapply(ORF,convert_ORF2SingleGeneName)),
                     dataset = "KOm")


### OE dataset 

metb_OE_data <- merge(OEmetallomics,GO_gset_MF_binding_metalwise, by = "ORF")%>%
  filter(metal == term )%>%
  filter(metal != "K")%>%
  mutate(Gene_Name = as.character(lapply(ORF,convert_ORF2SingleGeneName)),
         dataset = "OEm")

mett_OE_data <- merge(OEmetallomics,unique(rbind(GO_gset_MF_transporter_metalwise,
                                              philpott_metal_transporter_df)), by = "ORF")%>%
              filter(metal == term)%>%
              mutate(Gene_Name = as.character(lapply(ORF,convert_ORF2SingleGeneName)),
                     dataset = "OEm")

################################################
### Plot KO and OE metallomics data together ###
################################################

# metal_binders

colnames(metb_KO_data)[[3]] <- "Zscore"
colnames(metb_OE_data)[[3]] <- "Zscore"

metb_metallomics <- rbind(metb_KO_data,metb_OE_data)

pdf(paste0(plot_dir,"/KO_OE_metallomics_of_metalbindingproteins.pdf"),width = 12 ,height = 7)

ggplot(metb_metallomics,
       aes(x = paste(metal,dataset),
           y = Zscore,
           colour = metal,
           fill = metal))+
  geom_point( aes( shape = dataset,
                   alpha = dataset), 
              size = 3, alpha = 0.5,
              position = position_jitterdodge(jitter.width = 1))+
  geom_boxplot(alpha = 0.1,na.rm =T,linewidth = 0.2, width = 0.3, 
               colour = "black")+
  scale_colour_manual(values = colkey_Ele)+
  scale_fill_manual(values = colkey_Ele)+
  geom_text_repel(aes(label = Gene_Name,alpha = dataset))+
  scale_alpha_manual(values = c(0.6, 0.9)) +  
  theme_metallica()+
  labs(y = "Z-score metal of metal binding KO or OE mutant",
       x = "")+
  theme( legend.position = "none",
         axis.text.x = element_text(angle = 90))
dev.off()

## t test between zscores of KO and OE mutants

mets = unique(metb_metallomics$metal)

ttest_res_KOmOEm <- data.frame()

for(m in 1:length(unique(mets))){
  
  
  KOdata = filter(metb_metallomics, metal == mets[m] & dataset == "KOm")$Zscore
  
  OEdata = filter(metb_metallomics, metal == mets[m] & dataset == "OEm")$Zscore
  
  if(length(KOdata) > 3 & length(OEdata) > 3){
  ttres = t.test(KOdata,OEdata)
  
  ttest_res_KOmOEm <- rbind(ttest_res_KOmOEm,
                             cbind(p_value = as.numeric(ttres$p.value), 
                                   mean_KOm = as.numeric(ttres$estimate[[1]]),
                                   mean_OEm = as.numeric(ttres$estimate[[2]]),
                                  metal = as.character(mets[m]))
                            )
  }
}


## metal transporters

colnames(mett_KO_data)[[3]] <- "Zscore"
colnames(mett_OE_data)[[3]] <- "Zscore"

mett_metallomics <- rbind(mett_KO_data,mett_OE_data)

pdf(paste0(plot_dir,"/KO_OE_metallomics_of_metaltransportproteins.pdf"),width = 12 ,height = 7)

ggplot(mett_metallomics,
       aes(x = paste(metal,dataset),
           y = Zscore,
           colour = metal))+
  geom_point( aes( shape = dataset,
                   alpha = dataset), 
              size = 3, position = position_jitterdodge(jitter.width = 1))+
  scale_colour_manual(values = colkey_Ele)+
  geom_text_repel(aes(label = Gene_Name,alpha = dataset))+
  scale_alpha_manual(values = c(0.6, 0.9)) +  
  theme_metallica()+
  labs(y = "Z-score metal of metal transport KO or OE mutant",
       x = "")+
  theme( legend.position = "none",
         axis.text.x = element_text(angle = 90))
dev.off()


################
### KOgrowth ###
################


metb_data <- merge(metdepKOgrowth,GO_gset_MF_binding_metalwise, by = "ORF")%>%
  filter(metal == term)%>%
  group_by(metal, ORF)%>%
  mutate(mean_effect_size_log2 = mean(mean_effect_size_log2,na.rm = T))%>%
  ungroup()%>%
  dplyr::select(metal, ORF, mean_effect_size_log2)%>%
  unique()%>%
  mutate(Gene_Name = as.character(lapply(ORF,convert_ORF2SingleGeneName)))

pdf(paste0(plot_dir,"/metdepKOgrowth_of_metalbindingproteins.pdf"),width = 10 ,height = 6)
ggplot(metb_data,
       aes(x = metal,
           y = mean_effect_size_log2,
           group = metal,
           colour = metal))+
   geom_point( alpha = 0.6,  size = 3, position = position_jitterdodge(jitter.width = 0.4))+
  geom_text_repel(aes(label = Gene_Name))+
  scale_color_manual(values = colkey_Ele)+
  theme_metallica()+
  labs(y = "log2 mean effect size (growth) metal binding KOs",
       x = "")+
  theme(legend.position = "none")
dev.off()

mett_data <- merge(metdepKOgrowth,unique(rbind(GO_gset_MF_transporter_metalwise,
                                              philpott_metal_transporter_df)), by = "ORF")%>%
  filter(metal == term)%>%
  group_by(metal, ORF)%>%
  mutate(mean_effect_size_log2 = mean(mean_effect_size_log2,na.rm = T))%>%
  ungroup()%>%
  dplyr::select(metal, ORF, mean_effect_size_log2)%>%
  unique()%>%
  mutate(Gene_Name = as.character(lapply(ORF,convert_ORF2SingleGeneName)))

pdf(paste0(plot_dir,"/metdepKOgrowth_of_metaltransportproteins.pdf"),width = 10 ,height = 6)
ggplot(mett_data,
       aes(x = metal,
           y = mean_effect_size_log2,
           group = metal,
           colour = metal))+
   geom_point( alpha = 0.6,  size = 3, position = position_jitterdodge(jitter.width = 0.4))+
  geom_text_repel(aes(label = Gene_Name))+
  scale_color_manual(values = colkey_Ele)+
  theme_metallica()+
  labs(y = "log2 mean effect size (growth) metal transport KOs",
       x = "")+
  theme(legend.position = "none")
dev.off()

###########
### y5k ###
###########

metb_uniprot = GO_gset_MF_binding_metalwise%>%
               mutate(Uniprot_ID = as.character(lapply(ORF,convert_ORF2Uniprot)))
    

metb_data <- merge(y5kmetalspecific,metb_uniprot, by.x = "Protein.Group", by.y = "Uniprot_ID")%>% ## adding metal annotation of measured protein
  filter(KO %in% unique(GO_gset_MF_binding_metalwise$ORF))%>% # keep only metal binding KOs in df
  filter(term.x != "K")%>%
  filter(term.x == term.y)%>%
  filter(KO != ORF)%>%
  dplyr::select(KO,ORF,term.x,term.y,log2FC,p_value)%>%
  unique()%>%
  mutate(Gene_Name_KO = as.character(lapply(KO, convert_ORF2GeneName)),
         Gene_Name_meas = stringr::str_to_title(as.character(lapply(ORF, convert_ORF2GeneName))),
         label = ifelse(abs(log2FC) > log2(1.5), paste(Gene_Name_KO,Gene_Name_meas,sep="-"), NA))

pdf(paste0(plot_dir,"/y5kmetspec_log2FC_metalbindingproteins.pdf"),width = 9 ,height = 7)
ggplot(metb_data, aes(x = term.x, y = log2FC, colour = term.y)) +
  geom_point(aes(size = abs(log2FC), alpha = abs(log2FC)), alpha = 0.6, position = position_jitterdodge()) +
  geom_hline(aes(yintercept = log2(1.5)), colour = "red", linewidth = 0.2) +
  geom_hline(aes(yintercept = -log2(1.5)), colour = "red", linewidth = 0.2) +
  geom_text_repel(aes(label = label), size = 4) +
  scale_size_continuous(range = c(0.2, 4)) +  # Adjust the size range here
  scale_color_manual(values = colkey_Ele) +
  theme_metallica() +
  labs(
    x = "metal annotation of KO and measured proteins",
    y = "log2FC vs control of measured metal binding protein"
  ) +
  theme(legend.position = "none")

dev.off()

# summary 

y5k_metb_changes_per_KO = metb_data %>%
                mutate(Significant = ifelse(abs(log2FC) > log2(1.5) & p_value < 0.05, TRUE, FALSE))%>%
                filter(Significant)%>%
                group_by(term.x, KO,Gene_Name_KO)%>%
                summarise(num_sig = n())%>%
                ungroup()%>%
                arrange(num_sig)

y5k_metb_changes_per_KO_geneunique <- y5k_metb_changes_per_KO %>%
                                      dplyr::select(KO,num_sig,Gene_Name_KO)%>%
                                      unique()

pdf(paste0(plot_dir, "/y5k_numsigprots_numKOs_metalbinding.pdf"),width = 8,height = 6)
ggplot(y5k_metb_changes_per_KO,
       aes(x = factor(num_sig),
           fill = term.x,
           colour = term.x))+
  geom_histogram(alpha = 0.8,stat = "count", position = "dodge2")+
  geom_text_repel(aes(x = factor(num_sig), y = 5,label = ifelse(num_sig >5, Gene_Name_KO, NA)))+
  scale_colour_manual(values = colkey_Ele)+
  scale_fill_manual(values = colkey_Ele)+
  theme_metallica()+
  labs(x = "number of significant DA metal binders in KO",
       y = "number of KOs of metal binders",
       fill = "",
       colour = "")
dev.off()

########################################################################################
###  check if number of proteins changing per KO is related to number of DA proteins ###  
########################################################################################

###  read in growth rate data from Messner et al 

growth_rates_y5k <- read.csv(paste0(db_dir,"/yeast5k_growthrates_byORF.csv"), stringsAsFactors = F)

gr_rate_vs_numchanges <- merge(growth_rates_y5k, y5k_metb_changes_per_KO, by.x = "orf", by.y = "KO")

ggplot(gr_rate_vs_numchanges,
       aes(x = num_sig,
           y = SM))+
  geom_point()+
  geom_text_repel(aes(label = Gene_Name_KO))+
  geom_smooth(method = 'lm')+
  theme_metallica()


### transporters 

mett_uniprot = rbind(GO_gset_MF_transporter_metalwise,
                     philpott_metal_transporter_df)%>%
  mutate(Uniprot_ID = as.character(lapply(ORF,convert_ORF2Uniprot)))

mett_data <- merge(y5kmetalspecific,mett_uniprot, by.x = "Protein.Group", by.y = "Uniprot_ID")%>% ## adding metal annotation of measured protein
  filter(KO %in% unique(rbind(GO_gset_MF_transporter_metalwise,
                              philpott_metal_transporter_df)$ORF))%>% # keep only metal transporter KOs in df
  filter(term.x == term.y)%>%
  filter(KO != ORF)%>%
  dplyr::select(KO,ORF,term.x,term.y,log2FC,p_value)%>%
  unique()%>%
  mutate(Gene_Name_KO = as.character(lapply(KO, convert_ORF2GeneName)),
         Gene_Name_meas = as.character(lapply(ORF, convert_ORF2GeneName)),
         label = ifelse(abs(log2FC) > log2(1.5) & p_value < 0.05, paste(Gene_Name_KO,Gene_Name_meas,sep="-"), NA))


pdf(paste0(plot_dir,"/y5kmetspec_log2FC_metaltransportproteins.pdf"),width = 7 ,height = 7)
ggplot(mett_data,
       aes(x = term.x,
           y = log2FC,
           colour = term.y))+
  geom_point( alpha = 0.6,  size = 3,position = position_jitterdodge())+
  geom_hline(aes(yintercept = log2(1.5)), colour = "red", linewidth = 0.2)+
  geom_hline(aes(yintercept = -log2(1.5)), colour = "red", linewidth = 0.2)+
  geom_text_repel(aes(label = label), size = 2)+
  scale_color_manual(values = colkey_Ele)+
  theme_metallica()+
  labs(
    x = "metal annotation of KO and measured proteins",
    y = "log2FC vs control of measured metal transport protein")+
  theme(legend.position = "none")

dev.off()

##########################################################################################################
### Number of metal binders and transporters and other metal related proteins assessed by all datasets ###
##########################################################################################################

all_ds_sig_or_not <- read.csv(paste0(datasetcomparison_dir,"/output/tables/all_datasets_combined_significant_or_not.csv"),stringsAsFactors = F)

length(which(unique(GO_MF_all_metal_related_anno$ORF) %in% unique(all_ds_sig_or_not$ORF)))
length(unique(GO_MF_all_metal_related_anno$ORF))

#significant 
length(which(unique(GO_MF_all_metal_related_anno$ORF) %in% unique(filter(all_ds_sig_or_not, Significant)$ORF)))


#########################################################################
### Make a df of which metal binders were recovered in which datasets ###
#########################################################################

ds = c("metdepKOgrowth",
       "KOmetallomics",
       "OEmetallomics",
       "y5kmetalspecific",
       "metpertWTproteomics"
)

all_ds_recovery_metrel <- data.frame()

for(i in 1 :nrow(GO_MF_all_metal_related_anno)){

  ORF_i =   GO_MF_all_metal_related_anno[i,"ORF"]
  metal_anno =   GO_MF_all_metal_related_anno[i,"term"]
  
  if(metal_anno != "unspecific"){
  for(d in 1:length(ds)){
    
    sm_sig = filter(all_ds_sig_or_not,
             dataset == ds[d] &
             ORF == ORF_i &
             metal == metal_anno)$Significant
    
    if(length(sm_sig) == 0 ){
      sm_sig = NA
    }
    
    om_sig = filter(all_ds_sig_or_not,
                    dataset == ds[d] &
                      ORF == ORF_i &
                      metal != metal_anno)$Significant
    if(length(om_sig)>0){
      om_sig = any(om_sig)
      }else if(length(om_sig) == 0 ){
      om_sig = NA
      }
    
    all_ds_recovery_metrel <- rbind(all_ds_recovery_metrel,cbind(
                              dataset = ds[d],
                              ORF = ORF_i,
                              metal = metal_anno,
                              same_metal_sig = sm_sig,
                              other_metal_sig = om_sig
                            ))
  }
  }
}

## summarise

all_ds_recovery_metrel_smry <- all_ds_recovery_metrel%>%
                               group_by(metal,dataset)%>%
                               summarise(num_measured_sm = length(!is.na(same_metal_sig)),
                                         num_measured_om = length(!is.na(other_metal_sig)),
                                         
                                         num_sig_sm = sum(same_metal_sig == "TRUE",na.rm=T),
                                         num_sig_om = sum(other_metal_sig == "TRUE",na.rm=T),
                                         
                                         perc_sig_sm = 100*(num_sig_sm / num_measured_sm),
                                         perc_sig_om = 100* (num_sig_om / num_measured_om))%>%
                                ungroup()

# visualise
pdf(paste0(plot_dir,"/all_ds_metalrelated_recovery_samemetal.pdf"),width = 9,height = 3.5)
ggplot(unique(all_ds_recovery_metrel_smry[,c("dataset","metal","perc_sig_sm","num_measured_sm")]),
       aes(x = metal,
           y = dataset,
           fill = perc_sig_sm))+
  geom_tile()+
  geom_text(aes(label = num_measured_sm), size = 4 , colour = "white")+
  scale_fill_viridis_c(begin = 0.05) +
  theme_metallica()+
  labs(x = "", y = "")

ggplot(unique(all_ds_recovery_metrel_smry[,c("dataset","metal","perc_sig_om","num_measured_om")]),
       aes(x = metal,
           y = dataset,
           fill = perc_sig_om))+
  geom_tile()+
  geom_text(aes(label = num_measured_om), size = 4 , colour = "white")+
  scale_fill_viridis_c(begin = 0.05, limits = c(0,70)) +
  theme_metallica()+
  labs(x = "", y = "")

dev.off()

######################################################################################################
### Number of orphan metal binders and metal transporters with a phenotype in one or more datasets ###
######################################################################################################


metb_sig_acr_ds <- merge(GO_gset_MF_metalbinding_unspecific,
                         all_ds_sig_or_not, by = "ORF", all.x = T)


num_orphan_metal_binders <- length(unique(GO_gset_MF_metalbinding_unspecific$ORF))
num_orphan_metal_binders_measured_anyds <- length(unique(filter(metb_sig_acr_ds,
                                                  !is.na(dataset))$ORF))
num_orphan_metal_binders_measured_perds <- metb_sig_acr_ds%>%
                                           filter(!is.na(dataset))%>%
                                           group_by(dataset)%>%
                                           summarise(num_meas = length(unique(ORF)),
                                                     num_sig = sum(Significant))

orphan_metal_binders_nummetalsigin <- metb_sig_acr_ds%>%
                                      na.omit()%>%
                                      group_by(ORF,metal)%>%
                                      mutate(Sig_in_any = any(Significant))%>%
                                      ungroup()%>%
                                      dplyr::select(ORF,metal,Sig_in_any)%>%
                                      unique()%>%
                                      group_by(ORF)%>%
                                      mutate(num_metals_sig_in = sum(Sig_in_any))%>%
                                      ungroup()


numORF_orphan_metal_binders_nummetalsigin <- orphan_metal_binders_nummetalsigin%>%
                                             dplyr::select(ORF,num_metals_sig_in)%>%
                                             unique()%>%
                                             group_by(num_metals_sig_in)%>%
                                             summarise(num_orphan_metalbinders = length(ORF))%>%
                                             ungroup()%>%
                                             mutate(total_num_orphan_metb_meas = num_orphan_metal_binders_measured_anyds,
                                                    percentage_of_measured = 100*round(num_orphan_metalbinders/
                                                      total_num_orphan_metb_meas,3))

pdf(paste0(plot_dir,"/orphan_metal_binders_nummetals_sig_in.pdf"), width = 6, height = 4)
ggplot(numORF_orphan_metal_binders_nummetalsigin,
       aes(x = factor(num_metals_sig_in),
           y = num_orphan_metalbinders))+
  geom_bar(stat = "identity", width = 0.6,
           fill = viridis(5)[1])+
  geom_text(aes(y = num_orphan_metalbinders+3,
                label = paste0(percentage_of_measured,"%")))+
  theme_metallica()
dev.off()


orphan_metb_sig_only1metal <- orphan_metal_binders_nummetalsigin%>%
                              dplyr::select(ORF,num_metals_sig_in)%>%
                              unique()%>%
                              filter(num_metals_sig_in ==1)%>%
                              pull(ORF)%>%
                              unique()

orphan_metb_sig_only1metal_df <- filter(all_ds_sig_or_not,
                                     ORF %in% orphan_metb_sig_only1metal & Significant)%>%
                                 group_by(ORF, metal)%>%
                                 mutate(num_datasets_sig_in = length(unique(dataset)))
                                


orphan_ORF_metal_unq_comb <- unique(paste(filter(orphan_metb_sig_only1metal_df,
                                                 dataset == "metpertWTproteomics")$ORF,
                                          filter(orphan_metb_sig_only1metal_df,
                                                 dataset == "metpertWTproteomics")$metal))

metpertWTproteomics_extracell$ORF_metal <- paste(metpertWTproteomics_extracell$ORF,
                                                 metpertWTproteomics_extracell$Element)

mpt_filtered <- filter(metpertWTproteomics_extracell, 
                       ORF_metal %in% orphan_ORF_metal_unq_comb)

label_df <-  mpt_filtered%>%
             dplyr::select(ORF, Genes, Element)%>%
             unique()

ggplot(mpt_filtered,
       aes(x = log2(Element.Concentration),
           y = Log2FC_vs_AE,
           colour = Genes))+
  geom_smooth()+
  geom_text_repel(data = label_df,
                  aes(x = 0,
                      y = 0,
                      label = Genes))+
  facet_wrap("Element",scales = "free")+
  theme_metallica()


## OE metallomics example

Ca_OE_df <- filter(OEmetallomics,
                           metal == "Ca")
YDR182W_Ca_OE_df <- filter(Ca_OE_df, ORF == "YDR182W")


ggplot(Ca_OE_df,
       aes(x = Zscore_OEmetallomics,
           fill = metal,
           colour = metal))+
  geom_density(alpha = 0.3)+
  geom_vline(xintercept = YDR182W_Ca_OE_df$Zscore_OEmetallomics)+
  geom_text(data = YDR182W_Ca_OE_df,
            aes(x =  Zscore_OEmetallomics -1.2,
                y = 0.8,
                label = "YDR182W"))+
  scale_colour_manual(values = colkey_Ele)+
  scale_fill_manual(values = colkey_Ele)+
  theme_metallica()


