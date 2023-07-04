##`---
#`  Title: "metpertWTproteomics_compare_DAproteins_across_metals_intraextraseparated.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 7 June 2023 
#`  Description: Script to compare results of linear models across metals

#############################################
### source paths functions and libraries  ###
#############################################

#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("krassowski/complex-upset")

library(UpSetR)
library(ggrepel)

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

###########################################################################################
### Analyse genes DA along extra along one metal and DA along intra along another metal ###
###########################################################################################

# Filter rows where Significant equals to 1
lmfitres_extracell_n <- lmfitres_extracell[lmfitres_extracell$Significant == 1,]
lmfitres_intracell_n<- lmfitres_intracell[lmfitres_intracell$Significant == 1,]

# Rename the columns
names(lmfitres_extracell_n)[which(names(lmfitres_extracell_n) == "Element")] <- "extracellular_element_sig_in"
names(lmfitres_intracell_n)[which(names(lmfitres_intracell_n) == "Element")] <- "intracellular_element_sig_in"

# Join datasets by ORF and Genes
extra_intra_cross_sig <- merge(lmfitres_extracell_n, lmfitres_intracell_n, by=c("ORF", "Genes"))

# Select required columns
extra_intra_cross_sig <- extra_intra_cross_sig[,c("ORF", "Genes", "extracellular_element_sig_in", "intracellular_element_sig_in")]

## Count total proteins that were diff abundt alogn env of 1 and intarcellular of another 

tot_cross_metal_extInt_DA <- unique(filter(extra_intra_cross_sig, extracellular_element_sig_in != intracellular_element_sig_in)$ORF)
length(tot_cross_metal_extInt_DA)

# Filter rows where extracellular_element_sig_in is not equal to intracellular_element_sig_in
#extra_intra_cross_sig <- extra_intra_cross_sig[extra_intra_cross_sig$extracellular_element_sig_in != extra_intra_cross_sig$intracellular_element_sig_in,]

# Remove rows with NA values
extra_intra_cross_sig <- unique(na.omit(extra_intra_cross_sig))

# Print the dataframe
print(extra_intra_cross_sig)

extra_intra_cross_sig_smry <- extra_intra_cross_sig%>%
  group_by(extracellular_element_sig_in)%>%
  mutate(total_extracell_sig = length(ORF))%>%
  ungroup()%>%
  filter(extracellular_element_sig_in != intracellular_element_sig_in)%>%
  group_by(extracellular_element_sig_in,intracellular_element_sig_in)%>%
  mutate( num_met_extra_sig = length(ORF),
          frac_tot_extracell_sig = num_met_extra_sig/total_extracell_sig)%>%
  ungroup()%>%
  dplyr::select(extracellular_element_sig_in,intracellular_element_sig_in,
                total_extracell_sig,num_met_extra_sig,frac_tot_extracell_sig)%>%
  unique()

# Create the nested list
extra_intra_cross_sig_list <- list()

for (element in unique(extra_intra_cross_sig$extracellular_element_sig_in)) {
  extra_intra_cross_sig_list[[paste(element, "env")]] <- extra_intra_cross_sig$ORF[extra_intra_cross_sig$extracellular_element_sig_in == element]
}

for (element in unique(extra_intra_cross_sig$intracellular_element_sig_in)) {
  extra_intra_cross_sig_list[[paste(element, "itc")]] <- extra_intra_cross_sig$ORF[extra_intra_cross_sig$intracellular_element_sig_in == element]
}

# Prepare the data for UpSet plot
forUpsetR <- UpSetR::fromList(extra_intra_cross_sig_list)

if(any(colSums(forUpsetR)==0 )){
  forUpsetR <- forUpsetR[,-which(colSums(forUpsetR)==0)]
}


# Upset Plot
pdf(paste0(plot_dir,"/upsetplot_across_metals_intrextra_mettypewise.pdf"),height=6,width=8)

UpSetR::upset(forUpsetR, order.by = "freq",
              cutoff = 4,
              nintersects = 15,
              sets = sort(colnames(forUpsetR)),
              main.bar.color = "#CD919E",
              sets.bar.color = "#BCD2EE",
              text.scale = 2,
              point.size = 3,
              matrix.color = "#2b0b57",
              mb.ratio = c(0.6, 0.4)
)
dev.off()


################################################################
### Compare num of overlaps with correlations between metals ###
################################################################

overlap_counts <- extra_intra_cross_sig%>%
                  filter(extracellular_element_sig_in != intracellular_element_sig_in)%>%
                  group_by(extracellular_element_sig_in, intracellular_element_sig_in)%>%
                  summarize(Overlap = length(ORF))%>%
                  ungroup()%>%
                  unique()

## read in metal-metal correlation data 

metpert_metmeas_corr <- read.csv(paste0(metpert_WTmetallomics_dir,"/output/tables/metalperturbed_metalmeasured_correlations.csv"), stringsAsFactors = F)


met_hits_corrs <- merge(overlap_counts, metpert_metmeas_corr[,c( "element_measured", "element_perturbed","spearman","p_value_spearman")],
                             by.x = c("extracellular_element_sig_in","intracellular_element_sig_in"),
                             by.y = c("element_perturbed","element_measured"))%>%
  mutate(type = "pert-meas")%>%
  mutate(met_met_label = ifelse(p_value_spearman < 0.05 , 
                                paste(extracellular_element_sig_in,intracellular_element_sig_in,"*"),
                                paste(extracellular_element_sig_in,intracellular_element_sig_in)),
         env_metal_perturbed = extracellular_element_sig_in)

pdf(paste0(plot_dir,"/overlap_btw_extraintrasig_vs_metpert_metmeasured_scatterplot.pdf"),width = 8,height = 6)
ggplot(filter(met_hits_corrs),
       aes(y = Overlap,
           x = spearman))+
  geom_point(aes(colour = env_metal_perturbed),
             size = 4, alpha = 0.6)+
  #facet_wrap("metal_perturbed",scales = "free_x")+
  geom_text_repel(aes(label = met_met_label), size = 3, colour = "black")+
  scale_colour_manual(values = colkey_Ele)+
  theme_metallica()+
  labs(y = "size of overlap",
       x = "spearman correlation coefficient",
       colour = "metal \nperturbed")
dev.off()





# Initialize an empty data frame to store the overlap counts
metal_hits_overlap_counts <- data.frame(
  Metal1 = character(),
  Type1 = character(),
  Metal2 = character(),
  Type2 = character(),
  Overlap = numeric()
)

# Get the names of the metals
metals <- unique(c(extra_intra_cross_sig$extracellular_element_sig_in,
                   extra_intra_cross_sig$intracellular_element_sig_in))

# Loop over each pair of metals
for (i in 1:(length(metals))) {
  for (j in 1:length(metals)) {
    # We only compare extracellular vs intracellular, skip others
    if (i == j) next 
    
    # Get the names of the two metals and their types
    metal1 <- metals[i]
    type1 <- "env"  # assuming hitlist_eles contains extracellular metal data
    
    metal2 <- metals[j]
    type2 <- "itc"  # assuming hitlist_eles contains intracellular metal data
    
    # Calculate the overlap between the two metals
    overlap <- intersect(unique(filter(extra_intra_cross_sig,extracellular_element_sig_in == metal1)$ORF),
                         unique(filter(extra_intra_cross_sig,intracellular_element_sig_in == metal2)$ORF))
    
    # Store the metals, their types and overlap count in the data frame
    metal_hits_overlap_counts <- rbind(metal_hits_overlap_counts, 
                                       data.frame(Metal1 = metal1, 
                                                  Type1 = type1, 
                                                  Metal2 = metal2,
                                                  Type2 = type2,
                                                  Overlap = length(overlap)))
  }
}

# Print the overlap counts
print(metal_hits_overlap_counts)

###################################################################################################
### Calcualte how many overlaps are explained by high correlation at each correlation threshold ###
###################################################################################################

results_df <- data.frame(
  extracellular_element = character(), 
  overlapped_ORFs = integer(),
  overlapped_ORFs_explained_by_corr = integer(),
  correlation_threshold = numeric(),
  stringsAsFactors = FALSE
)

# Loop over each extracellular element

for(corr_threshold in seq(0.3, 1, by = 0.05)){
  
  for(extracellular_element in unique(extra_intra_cross_sig$extracellular_element_sig_in)) {
    
    ORFs_overlapping_with_any_metal = length(unique(filter(extra_intra_cross_sig,extracellular_element_sig_in == extracellular_element)$ORF))
    
    ORFs_explained_by_anymetal_corr <- list()
    
    # Loop over each intracellular element
    for(intracellular_element in unique(extra_intra_cross_sig$intracellular_element_sig_in)) {
      
      ORFs_explained_by_metal_im = unique(filter(extra_intra_cross_sig,
                                                        extracellular_element_sig_in == extracellular_element &
                                                        intracellular_element_sig_in == intracellular_element)$ORF)
      
      if(extracellular_element != intracellular_element){
  
        # Find correlation row in metpert_metmeas_corr data frame
        corr_row <- metpert_metmeas_corr %>%
          filter(
            (element_measured == extracellular_element & element_perturbed == intracellular_element) | 
              (element_measured == intracellular_element & element_perturbed == extracellular_element)
          ) %>%
          slice(which.max(abs(spearman)))
          # If Spearman correlation is > corr_threshold and p-value is < 0.05
          if (abs(corr_row$spearman) > corr_threshold & corr_row$p_value_spearman < 0.05) {
            # add ORFsthat are explained by the high correlation to explained overlaps list
            ORFs_explained_by_anymetal_corr <- c(ORFs_explained_by_anymetal_corr,ORFs_explained_by_metal_im)
          }
        }
    }
    
    # Add row to results data frame
    results_df <- rbind(results_df, data.frame(
      extracellular_element = extracellular_element,
      overlapped_ORFs = ORFs_overlapping_with_any_metal, # Total ORFs is the same as overlapped ORFs in this context
      overlapped_ORFs_explained_by_corr = length(unique(ORFs_explained_by_anymetal_corr)),
      correlation_threshold = corr_threshold))
    }
}

#############
### Plot ###
############

# Preparing the data
results_df_forplot <- results_df %>%
  mutate(fraction_explained_by_corr = overlapped_ORFs_explained_by_corr / overlapped_ORFs) %>%
  na.omit()

# Plotting
pdf(paste0(plot_dir,"/fraction_env_int_crosshits_explained_by_metmetcorr_vscorrthresh.pdf"),height = 6,width=8.3)
ggplot(results_df_forplot, aes(x = correlation_threshold, 
                               y = fraction_explained_by_corr, color = extracellular_element)) +
  geom_point(size = 1)+
  geom_line(alpha = 0.8) +
  scale_color_manual(values = colkey_Ele) +
  labs(
    y = "Fraction cross-metal overlaps between env & itc \nexplained by metmeas-metpert correlation",
    x = "Correlation Threshold",
    color = ""
  )+
  theme_metallica()
dev.off()


pdf(paste0(plot_dir,"/fraction_env_int_crosshits_explained_by_metmetcorr_barcorrthresh0pt8.pdf"),height = 5,width=6)
ggplot(data = filter(results_df_forplot, correlation_threshold == 0.8),
       aes(x = extracellular_element,
           colour = extracellular_element,
           fill = extracellular_element))+
  geom_bar(aes(y =  overlapped_ORFs), stat = "identity", width = 0.6,alpha = 0)+
  geom_bar(aes(y =  overlapped_ORFs_explained_by_corr),
               stat = "identity", width = 0.6, alpha = 0.6)+
  geom_text(aes(y =  overlapped_ORFs_explained_by_corr+15,
                label = paste0(100*round(fraction_explained_by_corr,2),"%")),
                colour = "black")+
  theme_metallica()+
  theme(legend.position = "none")+
  scale_fill_manual(values = colkey_Ele)+
  scale_colour_manual(values = colkey_Ele)+
  labs(x= "",y = "border:total overlapping ORFs\n fill:overlaps explained by metal correlations")
dev.off()
  
