####################################################################################
####          Script to plot results in all screens. for YMR196W and YBR287W    ####
####################################################################################

#`---
#`  Title: "metallica_datasetcomparison_hypothesisgeneration.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 10 August 2023
#`  Description: script to plot results from all screens for the 2 chosen proteins Ymr196w and Ybr287w
#`---

#############################################
### source paths functions and libraries  ###
#############################################

source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")

# general

source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/database_identifier_conversion_functions.R"))
source(paste0(code_dir,"/common_code/input_processed_databases_publisheddatasets.R"))


plot_dir = paste0( datasetcomparison_dir, "/output/plots")
dir.create(plot_dir,recursive = T)

output_tables_dir = paste0( datasetcomparison_dir, "/output/tables")
dir.create(output_tables_dir,recursive = T)


######################
### ORFs to verify ###
######################

chosen_ORFs = c("YMR196W","YBR287W")

######################
### Input datasets ###
######################

metpertWTproteomics_extracell <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/extracellularconc_vs_proteinabundance/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),
                                          stringsAsFactors = F)%>%
                                 filter(ORF %in% chosen_ORFs)

## lms with intracellular metal concentration

metpertWTproteomics_intracell <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/intracellularconc_vs_proteinabundance/lmfit_intracellular_metalconc_vs_log2PQ_with_SigNotSig_AdjPVthresh_0.05_FCthresh_1.5.csv"),
                                          stringsAsFactors = F)%>%
                                  filter(ORF %in% chosen_ORFs)

y5kmetalspecific <- read.csv(paste0(published_dataset_dir,"/y5k_proteomics/y5k_metalrelatedKOs_metalspecific.csv"),
                             stringsAsFactors = F)%>%
                    filter(measured_protein_name %in% as.character(lapply(chosen_ORFs, convert_ORF2SingleGeneName))|
                             KO %in% chosen_ORFs)
                            
### y5k for KOs of YBR287W and YMR196W
library(data.table)

y5k_protquants_WT <-  fread(paste0(published_dataset_dir,"/y5k_proteomics/yeast5k_noimpute_wide-2023-01-26_Messner2023.csv.gz"),stringsAsFactors = F)%>%
  as.data.frame()%>%
  reshape2::melt(id.vars = "Protein.Group",variable.name = "KO",value.name = "Protein.Quantity")%>%
  mutate(KO = as.character(KO))%>%
  tidyr::separate(col=KO, into = c(NA,NA,NA,NA,"KO",NA,NA),sep="_",remove=F)%>%
  filter(KO == "YOR202W")%>%
  group_by(Protein.Group)%>%
  summarize(Protein.Quantity_WT = mean(Protein.Quantity,na.rm = T))



## read in and filter p values and protein quantities of ORFs with metal specific metal annotations

y5k_chosenORFs_KOs_pvalues <- read.csv(paste0(published_dataset_dir,"/y5k_proteomics/yeast5k_stat_DE_bySTRAIN_Messner2023.csv"),stringsAsFactors = F)%>%
  as.data.frame()%>%
  reshape2::melt(id.vars = "Protein.Group",variable.name = "KO",value.name = "p_value")%>%
  mutate(KO = as.character(KO))%>%
  tidyr::separate(col=KO, into = c(NA,NA,NA,NA,"KO",NA,NA),sep="_",remove=F)%>%
  filter(KO %in% chosen_ORFs)%>%
  unique()%>%
  na.omit()%>%
  group_by(KO, Protein.Group)%>%
  summarize(p_value = mean(p_value))%>%
  ungroup()


y5k_chosenORFs_KOs_protquants <-  fread(paste0(published_dataset_dir,"/y5k_proteomics/yeast5k_noimpute_wide-2023-01-26_Messner2023.csv.gz"),stringsAsFactors = F)%>%
  as.data.frame()%>%
  reshape2::melt(id.vars = "Protein.Group",variable.name = "KO",value.name = "Protein.Quantity")%>%
  mutate(KO = as.character(KO))%>%
  tidyr::separate(col=KO, into = c(NA,NA,NA,NA,"KO",NA,NA),sep="_",remove=F)%>%
  filter(KO %in% chosen_ORFs)%>%
  unique()%>%
  na.omit()%>%
  group_by(KO, Protein.Group)%>%
  summarize(Protein.Quantity = mean(Protein.Quantity))%>%
  ungroup()


y5k_chosenORFs_KOs_protquants <- merge(y5k_chosenORFs_KOs_protquants, y5k_protquants_WT, by = "Protein.Group")%>%
                                 mutate(log2FC = log2(Protein.Quantity/Protein.Quantity_WT))

y5k_chosenORFs_KOs <- merge(y5k_chosenORFs_KOs_pvalues, 
                            y5k_chosenORFs_KOs_protquants, by = c("KO","Protein.Group"))
# add protein names of measured proteins

y5k_chosenORFs_KOs <- y5k_chosenORFs_KOs%>%
                      mutate(measured_protein_name = as.character(lapply(Protein.Group, convert_Uniprot2SingleGeneName)),
                             measured_protein_ORF = as.character(lapply(Protein.Group, convert_Uniprot2singleORF)))


#rm(y5k_chosenORFs_KOs_protquants,y5k_protquants_WT,y5k_chosenORFs_KOs_pvalues)


y5k_chosenORFs_KOs<- merge(y5k_chosenORFs_KOs,rbind(GO_gset_MF_binding_metalwise,
                                                    GO_gset_MF_transporter_metalwise,
                                                    philpott_metal_transporter_df),
                           by.x = "measured_protein_ORF", by.y = "ORF")%>%
  unique()
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
         dataset_type = "metallomics",
         Zscore = Zscore_KOmetallomics)%>%
  dplyr::select(-Zscore_KOmetallomics)


## metallomics screen on metal depletion conditions of the OE library ( check the exact name ? whats the BY number)

OEmetallomics <- read.csv(paste0(published_dataset_dir,"/metallomics/OEmetallomics_zscores_Iacovacci2020.csv"),stringsAsFactors = F,strip.white = T)%>%
  reshape2::melt(id.vars = "ORF", variable.name = "metal", value.name = "Zscore_OEmetallomics")%>%
  filter(metal %in% c("Ca","Cu","Fe","K","Mg","Mn","Mo","Na","Zn"))%>%
  mutate(Significant = ifelse(abs(Zscore_OEmetallomics) > 1.959, T,F),
         dataset = "OEmetallomics",
         dataset_type = "metallomics",
         Zscore = Zscore_OEmetallomics)%>%
  dplyr::select(-Zscore_OEmetallomics)

KO_OE_metallomics = rbind(KOmetallomics,OEmetallomics)

#################################
### plot data from each layer ###  YMR196W -- only in relation to Fe
#################################


# metpertWTproteomics_extracell  
pdf(paste0(plot_dir,"/YMR196W_metperWTproteomics.pdf"), width = 4,height = 4)

ggplot(filter(metpertWTproteomics_extracell, ORF == "YMR196W" & Element == "Fe"),
       aes(x = log2(Element.Concentration),
           y = Log2FC_vs_AE,
           colour = BioSpecID))+
  geom_smooth(aes(group = Element),
              method = "lm",
              formula = y ~ poly(x, 2),
              alpha = 0.1,
              colour = "black") + 
  geom_point(size = 4)+
  facet_wrap("Element",scales = "free")+
  geom_hline(yintercept = 0, linewidth = 0.1)+
  geom_vline(xintercept = 0, linewidth = 0.1)+
  scale_colour_manual(values = colkey_BioSpecID)+
  theme_metallica()+
  labs(x = "log2(environmental Fe concentration",
       y = 'log2(FC vs AE Ymr196w protein abundance)')+
  theme(legend.position = "none")

dev.off()

# KOgrowth 

pdf(paste0(plot_dir,"/YMR196W_KOgrowth.pdf"))

ggplot(metdepKOgrowth,
       aes(x = metal,
           y = mean_effect_size_log2,
           colour = metal))+
  geom_violin()+
  geom_point(data = filter(metdepKOgrowth, ORF == "YMR196W"),
             aes(x = metal,
                     y = mean_effect_size_log2,
                     colour = metal),
             size = 10, shape = "*")+
  geom_hline(yintercept = 0, linewidth = 0.1)+
  scale_colour_manual(values = colkey_Ele)+
  theme_metallica()+
  ylim(-1.5,1.5)+
  labs(x = " ",
       y = 'log2(FC vs AE Ymr196w protein abundance)')+
  theme(legend.position = "none")
dev.off()


pdf(paste0(plot_dir,"/YMR196W_KOgrowth_onlyFe.pdf"))
ggplot(filter(metdepKOgrowth,metal == "Fe"),
       aes(x = mean_effect_size_log2,
           colour = metal,
           fill = metal))+
  geom_density(alpha = 0.3)+
  geom_vline(xintercept = mean(filter(metdepKOgrowth, 
                                      ORF == "YMR196W" & metal == "Fe")$mean_effect_size_log2))+  scale_colour_manual(values = colkey_Ele)+
  theme_metallica()+
  labs(x = " ",
       y = 'log2(FC vs AE Ymr196w protein abundance)')+
  theme(legend.position = "none")
dev.off()


## ko and oe metallomics 
pdf(paste0(plot_dir,"/YMR196W_KOOEmetallomics.pdf"), width = 8,height = 5)

ggplot(filter(KO_OE_metallomics, ORF == "YMR196W"),
       aes(x = metal,
           y = Zscore,
           colour = metal,
           fill = metal))+
  geom_point(aes(shape = dataset),
             size = 4)+
  geom_line(aes(group = dataset,
                linetype = dataset), colour = "black")+
  geom_hline(yintercept = 0, linewidth = 0.1)+  
  geom_hline(yintercept = 1.96, linewidth = 0.1)+
  geom_hline(yintercept = 1.96, linewidth = 0.1)+
  scale_colour_manual(values = colkey_Ele)+ 
  scale_fill_manual(values = colkey_Ele)+
  scale_shape_manual(values = c(21,22))+
  theme_metallica()+
  labs(x = " ",
       y = 'Zscore metal concentration')

dev.off()


label_data = filter(KO_OE_metallomics,
                  ORF == "YMR196W" & metal == "Fe")

pdf(paste0(plot_dir,"/YMR196W_KOOEmetallomics_onlyFe.pdf"), width = 6,height = 4)
ggplot(filter(KO_OE_metallomics , metal == "Fe"),
       aes(x = Zscore,
           colour = metal,
           fill = metal))+
  geom_density(alpha = 0.3)+
  geom_vline(data = label_data,
             aes(xintercept = Zscore))+ 
  geom_text(data = label_data,
             aes(x = Zscore +1,
                 y = 0.8,
                 label = round(Zscore,2)), colour = "black")+ 
  scale_colour_manual(values = colkey_Ele)+
  theme_metallica()+
  facet_wrap("dataset", scales = "free")+
  labs(y = "density",
       x = 'Zscore')+
  theme(legend.position = "none")
dev.off()




## y5kproteome

pdf(paste0(plot_dir,"/YMR196W_y5k_measprot.pdf"), width = 5, height = 5)
ggplot(filter(y5kmetalspecific , measured_protein_name == "YMR196W" & 
                !(deleted_protein_name %in% chosen_ORFs) & term == "Fe"),
       aes(x = log2FC,
           y = -log10(p_value),
           colour = term,
           fill = term))+
  geom_point(aes(size = -log10(p_value)),alpha = 0.55)+
  geom_text_repel(aes(label = deleted_protein_name))+
  geom_hline(yintercept = -log10(0.05), linewidth = 0.1)+  
  geom_vline(xintercept = -log2(1.5), linewidth = 0.1)+
  geom_vline(xintercept = log2(1.5), linewidth = 0.1)+
    scale_colour_manual(values = colkey_Ele)+ 
  scale_fill_manual(values = colkey_Ele)+
  scale_shape_manual(values = c(21,22))+
  theme_metallica()+
  labs(x = "log2(FC vs AE) of Ymr196w protein abundance",
       y = '-log10(p-value)')+
  theme(legend.position = "none")
dev.off()

## what happens to metal binding or metal transport proteins when you delete YMR196W

ggplot(filter(y5k_chosenORFs_KOs , KO == "YMR196W"),
       aes(x = log2FC,
           y = -log10(p_value),
           colour = term,
           fill = term))+
  geom_point(aes(size = -log10(p_value)),alpha = 0.55)+
  geom_text_repel(aes(label = stringr::str_to_title(measured_protein_name)))+
  geom_hline(yintercept = -log10(0.05), linewidth = 0.1)+  
  geom_vline(xintercept = -log2(1.5), linewidth = 0.1)+
  geom_vline(xintercept = log2(1.5), linewidth = 0.1)+
  scale_colour_manual(values = colkey_Ele)+ 
  scale_fill_manual(values = colkey_Ele)+
  scale_shape_manual(values = c(21,22))+
  facet_wrap("term",scales = "free")+
  theme_metallica()+
  labs(x = " ",
       y = 'log2(FC vs AE metal related protein abundance)')+
  theme(legend.position = "none")


#################################
### plot data from each layer ###  YBR287W -- only in relation to Fe
#################################

# metpertWTproteomics_extracell  
pdf(paste0(plot_dir,"/YBR287W_metperWTproteomics.pdf"),width = 6,height = 6)

ggplot(filter(metpertWTproteomics_extracell, Significant == 1 &
              ORF == "YBR287W"),
       aes(x = log2(Element.Concentration),
           y = Log2FC_vs_AE,
           colour = BioSpecID))+
  geom_smooth(aes(group = Element),
              method = "lm",
              formula = y ~ poly(x, 2),
              alpha = 0.1,
              colour = "black") + 
  geom_point(size = 3)+
  facet_wrap("Element",scales = "free")+
  geom_hline(yintercept = 0, linewidth = 0.1)+
  geom_vline(xintercept = 0, linewidth = 0.1)+
  scale_colour_manual(values = colkey_BioSpecID)+
  theme_metallica()+
  labs(x = "log2(environmental Fe concentration",
       y = 'log2(FC vs AE Ybr287w protein abundance)')+
  theme(legend.position = "none")

dev.off()

# KOgrowth 

pdf(paste0(plot_dir,"/YBR287W_KOgrowth.pdf"), height = 4,width = 4.7)

ggplot(filter(metdepKOgrowth,!grepl("AE",metal)),
       aes(x = metal,
           y = mean_effect_size_log2,
           colour = metal,
           fill = metal))+
  geom_violin(linewidth = 0.1, alpha = 0.1)+
  geom_point(data = filter(metdepKOgrowth, ORF == "YBR287W" & 
                             metal != "AE"),
             aes(x = metal,
                 y = mean_effect_size_log2,
                 colour = metal),
             size = 2.5, alpha = 0.7)+
  geom_hline(yintercept = 0, linewidth = 0.1)+
  scale_colour_manual(values = colkey_Ele)+
  scale_fill_manual(values = colkey_Ele)+
  theme_metallica()+
  ylim(-0.5,0.5)+
  labs(x = " ",
       y = 'log2 effect size')+
  theme(legend.position = "none")
dev.off()



## ko and oe metallomics 
pdf(paste0(plot_dir,"/YBR287W_KOOEmetallomics.pdf"), width = 8.5,height = 5)

ggplot(filter(KO_OE_metallomics , ORF == "YBR287W"),
       aes(x = metal,
           y = Zscore,
           colour = metal,
           fill = metal))+
  geom_point(aes(shape = dataset),
             size = 4)+
  geom_line(aes(group = dataset,
                linetype = dataset), colour = "black")+
  geom_hline(yintercept = 0, linewidth = 0.1)+  
  geom_hline(yintercept = 1.96, linewidth = 0.1)+
  geom_hline(yintercept = 1.96, linewidth = 0.1)+
  scale_colour_manual(values = colkey_Ele)+ 
  scale_fill_manual(values = colkey_Ele)+
  scale_shape_manual(values = c(21,22))+
  theme_metallica()+
  labs(x = " ",
       y = 'Zscore metallomics',
       colour = "",
       fill = "",
       shape = "",
       linetype = "")

dev.off()
## y5kproteome


## ko and oe metallomics 
pdf(paste0(plot_dir,"/YBR287W_y5k_measprot.pdf"), width = 4.2, height = 4)
ggplot(filter(y5kmetalspecific , measured_protein_name == "YBR287W" & 
                !(deleted_protein_name %in% chosen_ORFs) ),
       aes(x = log2FC,
           y = -log10(p_value),
           colour = term,
           fill = term))+
  geom_point(aes(size = -log10(p_value)),alpha = 0.55)+
  geom_text_repel(aes(label = deleted_protein_name))+
  geom_hline(yintercept = -log10(0.05), linewidth = 0.1)+  
  geom_vline(xintercept = -log2(1.5), linewidth = 0.1)+
  geom_vline(xintercept = log2(1.5), linewidth = 0.1)+
  scale_colour_manual(values = colkey_Ele)+ 
  scale_fill_manual(values = colkey_Ele)+
  scale_shape_manual(values = c(21,22))+
  theme_metallica()+
 # facet_wrap("term",scales = "free_x")+
  labs(x = "log2 FC vs AE Ybr287w ",
       y = '-log10(p-value)')+
  theme(legend.position = "none")
dev.off()

pdf(paste0(plot_dir,"/YBR287W_y5k_measprot_facet.pdf"),width = 12,height = 12)
ggplot(filter(y5kmetalspecific , measured_protein_name == "YBR287W" & 
                !(deleted_protein_name %in% chosen_ORFs) ),
       aes(x = log2FC,
           y = -log10(p_value),
           colour = term,
           fill = term))+
  geom_point(aes(size = -log10(p_value)),alpha = 0.55)+
  geom_text_repel(aes(label = deleted_protein_name))+
  geom_hline(yintercept = -log10(0.05), linewidth = 0.1)+  
  geom_vline(xintercept = -log2(1.5), linewidth = 0.1)+
  geom_vline(xintercept = log2(1.5), linewidth = 0.1)+
  scale_colour_manual(values = colkey_Ele)+ 
  scale_fill_manual(values = colkey_Ele)+
  scale_shape_manual(values = c(21,22))+
  theme_metallica()+
  facet_wrap("term",scales = "free_x")+
  labs(x = " ",
       y = 'log2(FC vs AE Ymr196w protein abundance)')+
  theme(legend.position = "none")
dev.off()


pdf(paste0(plot_dir,"/YBR287W_y5k_KO_metrelprots.pdf"),width = 4.5,height = 4)
ggplot(filter(y5k_chosenORFs_KOs , KO == "YBR287W"),
       aes(x = log2FC,
           y = -log10(p_value),
           colour = term,
           fill = term))+
  geom_point(aes(size = -log10(p_value)),alpha = 0.55)+
  geom_text_repel(aes(label = stringr::str_to_title(measured_protein_name)))+
  geom_hline(yintercept = -log10(0.05), linewidth = 0.1)+  
  geom_vline(xintercept = -log2(1.5), linewidth = 0.1)+
  geom_vline(xintercept = log2(1.5), linewidth = 0.1)+
  scale_colour_manual(values = colkey_Ele)+ 
  scale_fill_manual(values = colkey_Ele)+
  scale_shape_manual(values = c(21,22))+
  #facet_wrap("term",scales = "free")+
  theme_metallica()+
  labs(x = "log2FC metal related protein ",
       y = '-log10(p-value)')+
  theme(legend.position = "none")
dev.off()


pdf(paste0(plot_dir,"/YBR287W_y5k_KO_metrelprots_facet.pdf"),width = 12,height = 12)
ggplot(filter(y5k_chosenORFs_KOs , KO == "YBR287W"),
       aes(x = log2FC,
           y = -log10(p_value),
           colour = term,
           fill = term))+
  geom_point(aes(size = -log10(p_value)),alpha = 0.55)+
  geom_text_repel(aes(label = stringr::str_to_title(measured_protein_name)))+
  geom_hline(yintercept = -log10(0.05), linewidth = 0.1)+  
  geom_vline(xintercept = -log2(1.5), linewidth = 0.1)+
  geom_vline(xintercept = log2(1.5), linewidth = 0.1)+
  scale_colour_manual(values = colkey_Ele)+ 
  scale_fill_manual(values = colkey_Ele)+
  scale_shape_manual(values = c(21,22))+
  facet_wrap("term",scales = "free")+
  theme_metallica()+
  labs(x = " ",
       y = 'log2(FC vs AE metal related protein abundance)')+
  theme(legend.position = "none")
dev.off()

