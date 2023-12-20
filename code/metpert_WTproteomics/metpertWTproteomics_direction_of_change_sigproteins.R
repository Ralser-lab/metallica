##`---
#`  Title: "metpertWTproteomics_direction_of_change_sigproteins.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 9 Oct 2023 
#`  Description: Script to check whether the significantly differentially expressed proteins correlate positively or negatively with metal concentraiton

#############################################
### source paths functions and libraries  ###
#############################################

# general
source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/layout_conversion.R"))
source(paste0(code_dir,"/metpert_WTproteomics/metpertWTproteomics_0_libraries_functions_genesetenrichments.R"))

plot_dir <- paste0(metpert_WTproteomics_dir,"/output/plots/direction_of_correlation_sigproteins")
dir.create(plot_dir, recursive = T)

outputtables_dir <- paste0(metpert_WTproteomics_dir,"/output/tables/direction_of_correlation_sigproteins")
dir.create(outputtables_dir, recursive = T)


# Read the data
lmfitres_extracell <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/extracellularconc_vs_proteinabundance/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),stringsAsFactors = F) %>%
  dplyr::select(Genes,ORF, BioSpecID, Element, Element.Concentration,Significant, PValAdj_BH,Log2FC_vs_AE, LeastComplexModel)%>%
  unique()%>%
  mutate(gene_element = paste(Genes, Element))

all_meas_proteins = unique(lmfitres_extracell$ORF)
all_metals = unique(lmfitres_extracell$Element)

lmfitres_intracell <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/intracellularconc_vs_proteinabundance/lmfit_intracellular_metalconc_vs_log2PQ_with_SigNotSig_AdjPVthresh_0.05_fcthresh_1.5.csv"),stringsAsFactors = F) %>%
  dplyr::select(Genes,ORF, Element,Significant, median_relative_intracellular_concentration,PValAdj_BH,median_log2_foldchangevsAE, LeastComplexModel)%>%
  mutate(gene_element = paste(Genes, Element))

###########################################################################################
### Analyse genes DA along extra along one metal and DA along intra along another metal ###
###########################################################################################

# Filter rows where Significant equals to 1
lmfitres_extracell_filtered <- lmfitres_extracell[lmfitres_extracell$Significant == 1,]
lmfitres_intracell_filtered<- lmfitres_intracell[lmfitres_intracell$Significant == 1,]


lmfitres_extracell_corr <- lmfitres_extracell_filtered %>%
  dplyr::group_by(ORF, Element, Genes) %>%
  dplyr::mutate(correlation = cor.test(Element.Concentration, Log2FC_vs_AE, method = "spearman")$estimate,
                   p.value = cor.test(Element.Concentration, Log2FC_vs_AE, method = "spearman")$p.value) %>%
  dplyr::ungroup()%>%
  mutate(correlation_type = ifelse(correlation > 0.25, "positive",
                                   ifelse(correlation < -0.25, "negative",
                                          "other")))

lmfitres_extracell_corr_summary <- lmfitres_extracell_corr%>%
                                   dplyr::select(Element,correlation_type,ORF)%>%
                                   unique()%>%
                                   group_by(Element, correlation_type)%>%
                                   summarize(num = n())%>%
                                   ungroup()

                                   

lmfitres_intracell_corr <- lmfitres_intracell_filtered %>%
  dplyr::group_by(ORF, Element, Genes) %>%
  dplyr::mutate(correlation = cor.test(median_relative_intracellular_concentration, median_log2_foldchangevsAE, method = "spearman")$estimate,
                   p.value = cor.test(median_relative_intracellular_concentration, median_log2_foldchangevsAE, method = "spearman")$p.value) %>%
  dplyr::ungroup()%>%
  mutate(correlation_type = ifelse(correlation > 0.25, "positive",
                                   ifelse(correlation < -0.25, "negative",
                                          "other"))) 

lmfitres_intracell_corr_summary <- lmfitres_intracell_corr %>%
                                    dplyr::select(Element,correlation_type,ORF)%>%
                                    unique()%>%
                                    group_by(Element, correlation_type) %>%
                                    summarize(num = n()) %>%
                                    ungroup()



lmfitres_extracell_corr_summary$type <- "Extracellular"
lmfitres_intracell_corr_summary$type <- "Intracellular"

merged_data <- rbind(lmfitres_extracell_corr_summary, lmfitres_intracell_corr_summary)%>%
               group_by(Element, type )%>%
               mutate(total = sum(num))%>%
               ungroup()%>%
               group_by(Element, type, correlation_type)%>%
               mutate(perc_of_sig = round(num/total,2)*100)%>%
               ungroup()

pdf(paste0(plot_dir,"/direction_of_correlation_signficant_proteins.pdf"),width = 15, height = 7)
ggplot(merged_data, 
       aes(x = Element,
           y = num,
           group = paste(Element, type),
           fill = correlation_type)) + 
  geom_bar(position = "dodge2",
           stat = "identity",
           width = 0.6) +
  facet_wrap("type")+
  labs(title = "",
        x = "",
       fill = "") +
  geom_text(aes(x = Element,
                y = num+10,
                label = paste0(perc_of_sig,"%")),
            size = 3)+
  ylab("number of significant ORFs") +
  scale_fill_manual(values = c("darkblue","darkred","gray"))+
  theme_metallica()
dev.off()


pdf(paste0(plot_dir,"/direction_of_correlation_signficant_proteins_percplot.pdf"),width = 15, height = 7)
ggplot(merged_data, 
       aes(x = Element,
           y = perc_of_sig,
           group = paste(Element, type),
           fill = correlation_type)) + 
  geom_bar(position = "stack",
           stat = "identity",
           width = 0.6) +
  facet_wrap("type")+
  labs(title = "",
       x = "",
       fill = "") +
  geom_text(aes(x = Element,
                y = 102,
                label = total),
            size = 3)+
  ylab("number of significant ORFs") +
  scale_fill_manual(values = c("darkblue","darkred","gray"))+
  theme_metallica()
dev.off()

### plot tile plots ##

pdf(paste0(plot_dir, "/correlation_type_geomsmooth.pdf"),width = 15,height = 12)
ggplot(lmfitres_extracell_corr, aes(x = log2(Element.Concentration),
                                    y = Log2FC_vs_AE,
                                    color = correlation_type,
                                    group = correlation_type)) +
  geom_smooth(method = "loess")+
  labs(title = " extracellular",
       x = "Element",
       y = "Log2FC_vs_AE",
       colour = "") +
  facet_wrap("Element", ncol = 3, scales = "free") +
  theme_metallica()+
  scale_colour_manual(values = c("darkblue","gray","darkred"))


ggplot(lmfitres_intracell_corr, aes(x = log2(median_relative_intracellular_concentration),
                                    y = median_log2_foldchangevsAE,
                                    color = correlation_type,
                                    group = correlation_type)) +
  geom_smooth(method = "loess")+
  labs(title = " intracellular",
       x = "Element",
       y = "Log2FC_vs_AE",
       colour = "") +
  facet_wrap("Element", ncol = 3, scales = "free") +
  theme_metallica()+
  scale_colour_manual(values = c("lightblue","gray","lightpink"))
dev.off()


#####################################################
### look at positive negative metal binding wise  ###
#####################################################

metal_bindingrbind(GO_gset_MF_binding_metalwise, GO_gset_MF_metalbinding_unspecific,
      GO_gset_MF_transporter_metalwise, GO_gset_MF_metaltransporter_unspecific)


lmfitres_extracell_corr_summary_metbinders <- merge(lmfitres_extracell_corr,GO_gset_MF_binding_metalwise, 
                                                    by.x = c("Element","ORF"), by.y = c("term","ORF")) %>%
  dplyr::select(Element,correlation_type,ORF)%>%
  unique()%>%
  group_by(Element, correlation_type)%>%
  summarize(num = n())%>%
  ungroup()%>%
  mutate(type = "Extracellular")



lmfitres_intracell_corr_summary_metbinders <- merge(lmfitres_intracell_corr,GO_gset_MF_binding_metalwise, 
                                                    by.x = c("Element","ORF"), by.y = c("term","ORF")) %>%
  dplyr::select(Element,correlation_type,ORF)%>%
  unique()%>%
  group_by(Element, correlation_type)%>%
  summarize(num = n())%>%
  ungroup()%>%
  mutate(type = "Intracellular")


merged_metbinder_corr <- rbind(lmfitres_extracell_corr_summary_metbinders,lmfitres_intracell_corr_summary_metbinders)%>%
                          group_by(Element, type )%>%
                          mutate(total = sum(num))%>%
                          ungroup()%>%
                          group_by(Element, type, correlation_type)%>%
                          mutate(perc_of_sig = round(num/total,2)*100)%>%
                          ungroup()


pdf(paste0(plot_dir,"/direction_of_correlation_signficant_proteins_percplot_metalbinders.pdf"),width = 15, height = 7)
ggplot(merged_metbinder_corr, 
       aes(x = Element,
           y = perc_of_sig,
           group = paste(Element, type),
           fill = correlation_type)) + 
  geom_bar(position = "stack",
           stat = "identity",
           width = 0.6) +
  facet_wrap("type")+
  labs(title = "",
       x = "",
       fill = "") +
  geom_text(aes(x = Element,
                y = 102,
                label = total),
            size = 3)+
  ylab("number of significant ORFs") +
  scale_fill_manual(values = c("darkblue","darkred","gray"))+
  theme_metallica()
dev.off()


######################################################################################################
### Enrichment analysis to check which proteins are in the positive, negative and other categories ###
######################################################################################################

## create dataframe with all metals and all ORFs
for_enrichment <- expand.grid(ORF = all_meas_proteins, Element = all_metals)%>%
                  
                  mutate(ORF_metal = paste(ORF, Element),
                         Significant = ifelse(ORF_metal %in% paste(
                                                             filter(lmfitres_extracell_corr,
                                                              correlation_type == "ct")$ORF,
                                                             filter(lmfitres_extracell_corr,
                                                                    correlation_type == "ct")$Element),1,0))%>%
                 dplyr::select(Element,Significant,ORF)%>%
                 unique()


extracell_hits_HyperGSAres<- run_HyperGSA(for_enrichment, "extrametVsprotabun")


