#`---
#`  Title: "metpertWTproteomics_compare_DAproteins_across_metals.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 7 June 2023 
#`  Description: Script to compare results of linear models across metals

#############################################
### source paths functions and libraries  ###
#############################################

#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("krassowski/complex-upset")

library(UpSetR)

# general
source("/Users/aulakhs/Documents/Ralser Lab/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/layout_conversion.R"))


plot_dir <- paste0(metpert_WTproteomics_dir,"/output/plots/comparison_between_metals")
dir.create(plot_dir, recursive = T)

outputtables_dir <- paste0(metpert_WTproteomics_dir,"/output/tables/comparison_between_metals")
dir.create(outputtables_dir, recursive = T)


# Read the data
lmfitres_extracell <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/extracellularconc_vs_proteinabundance/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),stringsAsFactors = F) %>%
  dplyr::select(Genes,ORF, BioSpecID, Element, Element.Concentration,Significant, PValAdj_BH,Log2FC_vs_AE, LeastComplexModel)%>%
  unique()%>%
  mutate(gene_element = paste(Genes, Element))


lmfitres_intracell <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/intracellularconc_vs_proteinabundance/lmfit_intracellular_metalconc_vs_log2PQ_with_SigNotSig_AdjPVthresh_0.05_fcthresh_1.5.csv"),stringsAsFactors = F) %>%
  dplyr::select(Genes,ORF, Element,Significant, median_relative_intracellular_concentration,PValAdj_BH,median_log2_foldchangevsAE, LeastComplexModel)%>%
  mutate(gene_element = paste(Genes, Element))


### Compare Hits

hitlist_eles <- list()
hitlist_eles_ext_int <- list()

eles <- unique(lmfitres_extracell$Element)

for( e in 1: length(eles)){
  
  hitlist_eles_ext_int[[paste(eles[[e]],"ext" )]] <- sort(as.character(unlist(unique(filter(lmfitres_extracell, Significant == 1 & Element == eles[e])$ORF))))
  hitlist_eles_ext_int[[paste(eles[[e]],"int" )]] <- sort(as.character(unlist(unique(filter(lmfitres_intracell, Significant == 1 & Element == eles[e])$ORF))))
  
  hitlist_eles[[eles[e]]] <- unique(c(hitlist_eles_ext_int[[paste(eles[[e]],"ext" )]] ,hitlist_eles_ext_int[[paste(eles[[e]],"int" )]]))
}


forUpsetR <- UpSetR::fromList(hitlist_eles)

if(any(colSums(forUpsetR)==0 )){
  forUpsetR <- forUpsetR[,-which(colSums(forUpsetR)==0)]
}


# Upset Plot
pdf(paste0(plot_dir,"/upsetplot_across_metals_intrextra.pdf"),height=6,width=8)

UpSetR::upset(forUpsetR,order.by = "freq",
                        group.by = "sets",
                        cutoff = 4,
                        nintersects = NA,
                        sets = sort(colnames(forUpsetR)),
                        main.bar.color ="#3c508b"  ,
                        sets.bar.color = "#26818e",
                        text.scale = 2,
                        point.size = 3,
                        matrix.color = "#2b0b57",
                        mb.ratio = c(0.5, 0.5)
              )
dev.off()


##################################################
### Calculate pairwise overlaps between metals ###
##################################################

# Initialize an empty data frame to store the overlap counts

metal_hits_overlap_counts <- data.frame(
  Metal1 = character(),
  Metal2 = character(),
  Overlap = numeric()
)

# Get the names of the metals
metals <- names(hitlist_eles)

# Loop over each pair of metals
for (i in 1:(length(metals)-1)) {
  for (j in (i+1):length(metals)) {
    # Get the names of the two metals
    metal1 <- metals[i]
    metal2 <- metals[j]
    
    # Calculate the overlap between the two metals
    overlap <- intersect(hitlist_eles[[metal1]], hitlist_eles[[metal2]])
    
    # Store the metals and overlap count in the data frame
    metal_hits_overlap_counts <- rbind(metal_hits_overlap_counts, 
                                       data.frame(Metal1 = metal1, Metal2 = metal2, Overlap = length(overlap)))
  }
}

# Print the overlap counts
metal_hits_overlap_counts

## read in metal-metal correlation data 

dir(metpert_WTmetallomics_dir)
