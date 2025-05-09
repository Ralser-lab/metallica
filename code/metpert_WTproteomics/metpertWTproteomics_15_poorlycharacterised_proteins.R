######################################################################
####          Script to prepare data for ensemble clustering      ####
######################################################################

#`---
#`  Title: "metperWTproteomics_prepare_data_for_ensembleclustering.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 18 June 2023
#`  Description: prepare data for ensemble clustering by Oliver Lemke (in python)
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

source(paste0(code_dir,"/metpert_WTproteomics/metpertWTproteomics_0_libraries_functions.R"))


plot_dir <- paste0(metpert_WTproteomics_dir,"/output/plots/poorly_characterised_proteins")
dir.create(plot_dir,recursive = T)

output_tables_dir <- paste0(metpert_WTproteomics_dir,"/output/tables/poorly_characterised_proteins")
dir.create(output_tables_dir,recursive = T)

#####################################################################################
### Read in results of lmfits, along intra and extracellular metal concentrations ###
#####################################################################################

# Read the data
lmfitres_extracell <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/extracellularconc_vs_proteinabundance/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),stringsAsFactors = F) %>%
  dplyr::select(Genes,ORF, BioSpecID, Element, Element.Concentration,Significant, PValAdj_BH,Log2FC_vs_AE, LeastComplexModel)%>%
  unique()

lmfitres_intracell <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/intracellularconc_vs_proteinabundance/lmfit_intracellular_metalconc_vs_log2PQ_with_SigNotSig_AdjPVthresh_0.05_fcthresh_1.5.csv"),stringsAsFactors = F) %>%
  dplyr::select(Genes,ORF, Element,Significant, median_relative_intracellular_concentration,PValAdj_BH,median_log2_foldchangevsAE, LeastComplexModel)

### Summarize hits wrt uniprot annotaiton status
intcell2bind <- unique(lmfitres_intracell[,c("ORF","Element","Significant")])%>%
                mutate(type ="int")

extcell2bind <- unique(lmfitres_extracell[,c("ORF","Element","Significant")])%>%
                mutate(type ="ext")

lmfitres_annoscore <- alldata <- rbind(intcell2bind,extcell2bind)%>%
                      na.omit()%>%
                      group_by(ORF,Element)%>%
                      mutate(Significant_exORin = ifelse(any(Significant) ==1, T, F))%>%
                      dplyr::select(-Significant)%>%
                      unique()

lmfitres_annoscore <- merge(lmfitres_annoscore, GenProt_SGD[,c("ORF","Uniprot.Annotation.Score")], by = "ORF")


extcell2bind_annoscore <- merge(extcell2bind, GenProt_SGD[,c("ORF","Uniprot.Annotation.Score")], by = "ORF")
extcell2bind_annoscore_ORFs <- unique(filter(extcell2bind_annoscore, Significant & Uniprot.Annotation.Score < 3)$ORF)
print(paste("num poorly characterised proteins identified along extracellular metal pert:",length(extcell2bind_annoscore_ORFs)))

intcell2bind_annoscore <- merge(intcell2bind, GenProt_SGD[,c("ORF","Uniprot.Annotation.Score")], by = "ORF")
intcell2bind_annoscore_ORFs <- unique(filter(intcell2bind_annoscore, Significant & Uniprot.Annotation.Score < 3)$ORF)
print(paste("num poorly characterised proteins identified along intracellular metal conc:",length(intcell2bind_annoscore_ORFs)))

length(which(!intcell2bind_annoscore_ORFs %in% extcell2bind_annoscore_ORFs))
length(which(!extcell2bind_annoscore_ORFs %in% intcell2bind_annoscore_ORFs))

lmfitres_annoscore_smry <- lmfitres_annoscore%>%
                           group_by(Element)%>%
                           mutate(total_measured = length(ORF),
                                  total_significant = sum(Significant_exORin))%>%
                           ungroup()%>%
                           group_by(Element,Uniprot.Annotation.Score,
                                    total_measured,total_significant)%>%
                           summarise(num_significant = sum(Significant_exORin))%>%
                           ungroup()%>%
                           mutate(poorly_char = ifelse(Uniprot.Annotation.Score < 3, T,F))
        
                           
## plot
pdf(paste0(plot_dir,"/numproteins_DA_by_uniprotscore.pdf"),width = 7,height = 6)
ggplot()+
  geom_bar(data = unique(lmfitres_annoscore_smry[,c("Element","total_measured")]),
           aes(
               x = Element,
               y = total_measured),
               width = 0.6,
               stat = "identity",
               fill = NA, 
               colour = "black",
               linewidth = 0.1)+
  geom_text(data = unique(lmfitres_annoscore_smry[,c("Element","total_significant")]),
           aes(
             x = Element,
             y = total_significant+50,
             label = total_significant),
             colour = "black")+
  geom_bar(data = lmfitres_annoscore_smry,
           aes(
               x = Element,
               y = num_significant,
               fill = Uniprot.Annotation.Score),
               width = 0.6,
               stat = "identity")+
  scale_fill_viridis_d(option = "C", end = 0.95)+
  theme_metallica()+
  labs(x = "", 
       y = "number of proteins",
       fill = "uniprot\nannotation\nscore")
dev.off()


## Count number uncharacterised

unique(filter(lmfitres_annoscore, Uniprot.Annotation.Score < 3)$ORF)

as.character(lapply(unique(filter(lmfitres_annoscore, Uniprot.Annotation.Score < 3)$ORF), convert_ORF2SingleGeneName))

as.character(lapply(unique(filter(lmfitres_annoscore, Uniprot.Annotation.Score < 3)$ORF), convert_ORF2Uniprot))


unique(filter(lmfitres_annoscore, Uniprot.Annotation.Score < 2)$ORF)

as.character(lapply(unique(filter(lmfitres_annoscore, Uniprot.Annotation.Score < 2)$ORF), convert_ORF2SingleGeneName))

