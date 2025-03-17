#########################################################################################################################################
####          Script to plot the abundance of proteins involved in creating and maintaining proton motive force across membranes     ####
#########################################################################################################################################

#`---
#`  Title: "metperWTproteomics_PMF_proteins.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 11 Jan 2024
#`  Description: calculate signalling, metabolic pathways, TFs and protein complexes affected
#`---


#############################################
### source paths functions and libraries  ###
#############################################

source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")

# general

source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/database_identifier_conversion_functions.R"))
source(paste0(code_dir,"/common_code/input_processed_databases_publisheddatasets.R"))

# specific paths and scripts
library(RColorBrewer)
source(paste0(code_dir,"/metpert_WTproteomics/metpertWTproteomics_0_libraries_functions.R"))


plot_dir <- paste0(metpert_WTproteomics_dir,"/output/plots/proton_motive_force")
dir.create(plot_dir,recursive = T)

output_tables_dir <- paste0(metpert_WTproteomics_dir,"/output/tables/proton_motive_force")
dir.create(output_tables_dir,recursive = T)

PMF_proteins_BP <- GO_gset_BP%>%
  filter(grepl("proton",term))

PMF_proteins_MF <- GO_gset_MF%>%
  filter(grepl("proton",term))

PMF_proteins_CC <- GO_gset_CC%>%
  filter(grepl("proton",term))

PMF_proteins <- merge(PMF_proteins_BP,PMF_proteins_MF, by = "ORF", all=T)

PMF_proteins <- merge(PMF_proteins, PMF_proteins_CC, by = "ORF", all = T )
## how many do we measure proteins from 

# Read the data
lmfitres_extracell <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/extracellularconc_vs_proteinabundance/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),stringsAsFactors = F)

lmfitres_intracell <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/intracellularconc_vs_proteinabundance/lmfit_intracellular_metalconc_vs_log2PQ_with_SigNotSig_AdjPVthresh_0.05_fcthresh_1.5.csv"),stringsAsFactors = F)



###################################
### Plot PMF protein abundances ###
###################################

###############
## extracell ##
###############

PMF_BP_PQ_extracell <- filter(lmfitres_extracell,
                              ORF %in% PMF_proteins_BP$ORF)

PMF_BP_PQ_extracell <- merge(PMF_BP_PQ_extracell, PMF_proteins_BP, by = "ORF")

pdf(paste0(plot_dir,"/PMF_BP_terms_sig_genes_extracell.pdf"), width = 17, height = 10)

ggplot(filter(PMF_BP_PQ_extracell, Significant == 1), 
       aes(x = log2(Element.Concentration),
           y = Log2FC_vs_AE,
           colour = Element)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth(aes(group = paste(ORF, Element), colour = Element), alpha = 0.2) +
  facet_wrap(c("term"), scales = "free") +
  scale_color_manual(values = colkey_Ele) +
  theme_metallica()
dev.off()

## focusing on proton transmembrane transport

pdf(paste0(plot_dir,"/PMF_BP_terms_sig_genesextracell.pdf"), width = 12, height = 10)
ggplot(filter(PMF_BP_PQ_extracell, Significant == 1 & term == "proton transmembrane transport"), 
       aes(x = log2(Element.Concentration),
           y = Log2FC_vs_AE,
           colour = Genes)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(aes(group = paste(ORF, Element), colour = Genes), size = 0.75, alpha = 0.1) +
  facet_wrap(c( "Element"), scales = "free") +
  theme_metallica()+
  scale_colour_manual(values = c(
    rep(brewer.pal(8,"Dark2")[6],8),
    rep(brewer.pal(8,"Dark2")[3],3),
    rep(brewer.pal(8,"Dark2")[5],4)
  ))
dev.off()

PMF_BP_PQ_extracell_sig <- filter(PMF_BP_PQ_extracell, Significant == 1 & term == "proton transmembrane transport")

PMF_BP_PQ_extracell_sig_labels <- PMF_BP_PQ_extracell_sig %>%
  group_by(Genes, Element) %>%
  filter(Element.Concentration == min(Element.Concentration)) %>%  
  mutate(minlog2_extracell = log2(Element.Concentration)) %>%  
  ungroup() %>%
  dplyr::select(Genes, Element, minlog2_extracell, Mean_of_reps_Log2FC_vs_AE) %>%  
  unique()


pdf(paste0(plot_dir, "/PMF_proteins_significant_sep_by_loc_extracellular.pdf"), width = 7.5, height = 6)

# 1. Mitochondrial Proteins
ggplot(filter(PMF_BP_PQ_extracell_sig, Genes %in% c("ATP1", "COX2", "COX5B", "COX6", "COX8", "COX9", "COX12", "COX13")), 
       aes(x = log2(Element.Concentration),
           y = Log2FC_vs_AE,
           colour = Element)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(aes(group = paste(ORF, Element), colour = Element, fill = Element),
              method = "lm", formula = y ~ poly(x, 2),
              size = 0.75, alpha = 0.05) +
  geom_text(data = filter(PMF_BP_PQ_extracell_sig_labels, Genes %in% c("ATP1", "COX2", "COX5B", "COX6", "COX8", "COX9", "COX12", "COX13")), 
            aes(x = minlog2_extracell, 
                y = Mean_of_reps_Log2FC_vs_AE, 
                label = Genes), 
            size = 3, hjust = 0.5, vjust = 1) +
  theme_metallica() +
  scale_colour_manual(values = colkey_Ele) +
  scale_fill_manual(values = colkey_Ele) +
  labs(title = "Mitochondria")

# 2. Vacuolar Proteins
ggplot(filter(PMF_BP_PQ_extracell_sig, Genes %in% c("VMA1", "VMA5", "VMA8", "VMA10")), 
       aes(x = log2(Element.Concentration),
           y = Log2FC_vs_AE,
           colour = Element)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(aes(group = paste(ORF, Element), colour = Element, fill = Element),
              method = "lm", formula = y ~ poly(x, 2),
              size = 0.75, alpha = 0.05) +
  geom_text(data = filter(PMF_BP_PQ_extracell_sig_labels, Genes %in% c("VMA1", "VMA5", "VMA8", "VMA10")), 
            aes(x = minlog2_extracell, 
                y = Mean_of_reps_Log2FC_vs_AE, 
                label = Genes), 
            size = 3, hjust = 0.5, vjust = 1) +
  theme_metallica() +
  scale_colour_manual(values = colkey_Ele) +
  scale_fill_manual(values = colkey_Ele) +
  labs(title = "Vacuole")

# 3. Plasma Membrane Proteins
ggplot(filter(PMF_BP_PQ_extracell_sig, Genes %in% c("HXT1", "HXT3", "HXT5")), 
       aes(x = log2(Element.Concentration),
           y = Log2FC_vs_AE,
           colour = Element)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(aes(group = paste(ORF, Element), colour = Element, fill = Element),
              method = "lm", formula = y ~ poly(x, 2),
              size = 0.75, alpha = 0.05) +
  geom_text(data = filter(PMF_BP_PQ_extracell_sig_labels, Genes %in% c("HXT1", "HXT3", "HXT5")), 
            aes(x = minlog2_extracell, 
                y = Mean_of_reps_Log2FC_vs_AE, 
                label = Genes), 
            size = 3, hjust = 0.5, vjust = 1) +
  theme_metallica() +
  scale_colour_manual(values = colkey_Ele) +
  scale_fill_manual(values = colkey_Ele) +
  labs(title = "Plasma Membrane")

dev.off()

###############
## intracell ##
###############


PMF_BP_PQ_intracell <- filter(lmfitres_intracell,
                              ORF %in% PMF_proteins_BP$ORF)

PMF_BP_PQ_intracell <- merge(PMF_BP_PQ_intracell, PMF_proteins_BP, by = "ORF")

pdf(paste0(plot_dir,"/PMF_BP_terms_sig_genes_intracell.pdf"), width = 24, height = 10)
ggplot(filter(PMF_BP_PQ_intracell, Significant ==1),
       aes(x = log2(median_relative_intracellular_concentration),
           y = median_log2_foldchangevsAE,
           colour = Element))+
  geom_point(size = 3, alpha = 0.75)+
  geom_smooth(aes(group = paste(ORF,Element),colour = Element))+
  facet_wrap(c("term"), scales = "free", nrow = 2)+
  scale_color_manual(values = colkey_Ele)+
  theme_metallica()
dev.off()



PMF_BP_PQ_intracell_sig <- filter(PMF_BP_PQ_intracell, Significant == 1 & term == "proton transmembrane transport")

### separate by type of protein 
PMF_BP_PQ_intracell_sig_labels <- PMF_BP_PQ_intracell_sig %>%
  group_by(Genes, Element) %>%
  filter(median_relative_intracellular_concentration == min(median_relative_intracellular_concentration)) %>%
  mutate(minlog2_intcell = log2(median_relative_intracellular_concentration)) %>%
  ungroup() %>%
  dplyr::select(Genes, Element, minlog2_intcell, median_log2_foldchangevsAE) %>%
  unique()


pdf(paste0(plot_dir, "/PMF_proteins_significant_sep_by_loc_intracellular.pdf"), width = 7.5, height = 6)
# mitochondrial
ggplot(filter(PMF_BP_PQ_intracell_sig, Genes %in% c("ATP1","ATP2", "COX2","COX5B","COX6","COX8","COX9","COX12","COX13")), 
       aes(x = log2(median_relative_intracellular_concentration),
           y = median_log2_foldchangevsAE,
           colour = Element)) +
  geom_point(size = 3, alpha = 0.2) +
  geom_smooth(aes(group = paste( ORF, Element), colour = Element, fill = Element),
              method = "lm",formula = y ~ poly(x, 2),
              size = 0.75, alpha = 0.05) +
  geom_text(data = filter(PMF_BP_PQ_intracell_sig_labels, Genes %in% c("ATP1","ATP2", "COX2","COX5B","COX6",
                                                                       "COX8","COX9","COX12","COX13")), 
            aes(x = minlog2_intcell, 
                y = median_log2_foldchangevsAE, 
                label = Genes), 
            size = 3, hjust = 0.5, vjust = 1) +  
  #facet_wrap("Element", scales = "free")+
  theme_metallica()+
  scale_colour_manual(values = colkey_Ele)+
  scale_fill_manual(values = colkey_Ele)+
  labs(title = "Mitochondria")



# Vacuole
ggplot(filter(PMF_BP_PQ_intracell_sig, Genes %in% c("STV1", "VMA2", "VMA4", "VMA5", "VMA7",
                                                    "VMA8", "VNX1")), 
       aes(x = log2(median_relative_intracellular_concentration),
           y = median_log2_foldchangevsAE, 
           colour = Element)) +
  geom_point(size = 3, alpha = 0.2) +
  geom_text(data = filter(PMF_BP_PQ_intracell_sig_labels, Genes %in% c("STV1", "VMA2", "VMA4", "VMA5", "VMA7",
                                                                "VMA8", "VNX1")), 
            aes(x = minlog2_intcell, 
                y = median_log2_foldchangevsAE, 
                label = Genes), 
            size = 3, hjust = 0.5, vjust = 1) +  
  geom_smooth(aes(group = paste(ORF, Element), colour = Element, fill = Element), 
              method = "lm", formula = y ~ poly(x, 2), 
              size = 0.75, alpha = 0.05) +
  #facet_wrap(~Element, scales = "free") +  
  theme_metallica()+
  scale_colour_manual(values = colkey_Ele) +
  scale_fill_manual(values = colkey_Ele) +
  labs(title = "Vacuole")


# Plasma Membrane
ggplot(filter(PMF_BP_PQ_intracell_sig, Genes %in% c("HXT3", "HXT4", "HXT5", "PMA1","PMA2")), 
       aes(x = log2(median_relative_intracellular_concentration),
           y = median_log2_foldchangevsAE, 
           colour = Element)) +
  geom_point(size = 3, alpha = 0.2) +
  geom_text(data = filter(PMF_BP_PQ_intracell_sig_labels, Genes %in% c("HXT3", "HXT4", "HXT5", "PMA1","PMA2")), 
            aes(x = minlog2_intcell, 
                y = median_log2_foldchangevsAE, 
                label = Genes), 
            size = 3, hjust = 0.5, vjust = 1) +  
  geom_smooth(aes(group = paste(ORF, Element), colour = Element, fill = Element), 
              method = "lm", formula = y ~ poly(x, 2), 
              size = 0.75, alpha = 0.05) +
  #facet_wrap(~Element, scales = "free") +  
  theme_metallica()+
  scale_colour_manual(values = colkey_Ele) +
  scale_fill_manual(values = colkey_Ele) +
  labs(title = "Plasma Membrane")



# Golgi & ER
ggplot(filter(PMF_BP_PQ_intracell_sig, Genes %in% c("GEF1", "PMR1")), 
       aes(x = log2(median_relative_intracellular_concentration),
           y = median_log2_foldchangevsAE, 
           colour = Element)) +
  geom_point(size = 3, alpha = 0.2) +
  geom_text(data = filter(PMF_BP_PQ_intracell_sig_labels, Genes %in% c("GEF1", "PMR1")), 
            aes(x = minlog2_intcell, 
                y = median_log2_foldchangevsAE, 
                label = Genes), 
            size = 3, hjust = 0.5, vjust = 1) +  
  geom_smooth(aes(group = paste(ORF, Element), colour = Element, fill = Element), 
              method = "lm", formula = y ~ poly(x, 2), 
              size = 0.75, alpha = 0.05) +
  #facet_wrap(~Element, scales = "free") +  
  theme_metallica()+
  scale_colour_manual(values = colkey_Ele) +
  scale_fill_manual(values = colkey_Ele) +
  labs(title = "Golgi and ER")
dev.off()


