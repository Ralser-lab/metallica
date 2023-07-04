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


### Compare Hits

hitlist_eles <- list()
hitlist_eles_ext_int <- list()

eles <- c("Ca","Cu","Fe","K","Mg","Mn","Mo","Na","Zn")

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
                        #group.by = "sets",
                        cutoff = 4,
                        nintersects = 20,
                        sets = sort(colnames(forUpsetR)),
                        main.bar.color ="#CD919E"  ,
                        sets.bar.color = "#BCD2EE",
                        text.scale = 2,
                        point.size = 3,
                        matrix.color = "#2b0b57",
                        mb.ratio = c(0.6, 0.4)
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

metpert_metmeas_corr <- read.csv(paste0(metpert_WTmetallomics_dir,"/output/tables/metalperturbed_metalmeasured_correlations.csv"), stringsAsFactors = F)


met_hits_corrs_dir1 <- merge(metal_hits_overlap_counts, metpert_metmeas_corr[,c( "element_measured", "element_perturbed","spearman","p_value_spearman")],
                        by.x = c("Metal1","Metal2"),
                        by.y = c("element_perturbed","element_measured"))%>%
                        mutate(type = "pert-meas")%>%
                        mutate(met_met_label = ifelse(p_value_spearman < 0.05, paste(Metal1,Metal2),NA),
                               metal_perturbed = Metal1)

met_hits_corrs_dir2 <- merge(metal_hits_overlap_counts, metpert_metmeas_corr[,c( "element_measured", "element_perturbed","spearman","p_value_spearman")],
                             by.x = c("Metal1","Metal2"),
                             by.y = c("element_measured","element_perturbed"))%>%
                             mutate(type = "meas-pert")%>%
                             mutate(met_met_label = ifelse(p_value_spearman < 0.05, paste(Metal2,Metal1),NA),
                                    metal_perturbed = Metal2)

met_hits_corrs <- rbind(met_hits_corrs_dir1,met_hits_corrs_dir2)%>%
                  filter(Overlap > 0)
                  
pdf(paste0(plot_dir,"/overlap_btw_metals_vs_metpert_metmeasured_scatterplot.pdf"),width = 6.1,height = 4)
ggplot(filter(met_hits_corrs),
       aes(y = Overlap,
           x = spearman))+
  geom_point(aes(colour = metal_perturbed),
             size = 4, alpha = 0.6)+
  #facet_wrap("metal_perturbed",scales = "free_x")+
  geom_text_repel(aes(label = met_met_label), size = 2, colour = "black")+
  scale_colour_manual(values = colkey_Ele)+
  theme_metallica()+
  labs(y = "size of overlap",
       x = "spearman correlation coefficient",
       colour = "metal \nperturbed\nin correlation \ncalculation")
dev.off()

#####################################################################################################################################################################
### Calculate what % of total proteins DA in a metal along intracellular or environmental metal concentration are also DA in a metal with high + or - correlation ###
#####################################################################################################################################################################

# Initialize a data frame to hold results
results_df <- data.frame(
  metal = character(), 
  total_ORFs = integer(), 
  unique_ORFs = integer(),
  shared_ORFs_not_corr = integer(),
  shared_ORFs_corr_gt = integer(), 
  correlation_threshold = numeric(),
  stringsAsFactors = FALSE
)

# Loop over each metal
for (metal in names(hitlist_eles)) {
  
  # List of ORFs for the current metal
  metal_orfs <- hitlist_eles[[metal]]
  total_ORFs <- length(metal_orfs)
  
  # Loop over each threshold
  for (corr_threshold in seq(0.3, 1, by = 0.05)) {
    
    # Initialize vectors to hold ORFs
    shared_ORFs_corr_gt <- c()
    shared_ORFs_not_corr <- c()
    
    # Loop over all other metals
    for (other_metal in names(hitlist_eles)) {
      if (other_metal != metal) {
        
        # List of ORFs for the other metal
        other_metal_orfs <- hitlist_eles[[other_metal]]
        
        # Find shared ORFs
        shared_orfs <- intersect(metal_orfs, other_metal_orfs)
        
        # Find correlation row in either order
        corr_row <- metpert_metmeas_corr %>%
          filter(
            (element_measured == metal & element_perturbed == other_metal) | 
              (element_measured == other_metal & element_perturbed == metal)
          ) %>%
          slice(which.max(abs(spearman)))
        
        # If Spearman correlation is > corr_threshold and p-value is < 0.05
        if (abs(corr_row$spearman) > corr_threshold & corr_row$p_value_spearman < 0.05) {
          shared_ORFs_corr_gt <- union(shared_ORFs_corr_gt, shared_orfs)
        }
      }
    }
    
    # Add remaining shared ORFs that do not meet criteria and are not already in shared_ORFs_corr_gt
    for (other_metal in names(hitlist_eles)) {
      if (other_metal != metal) {
        other_metal_orfs <- hitlist_eles[[other_metal]]
        shared_orfs <- intersect(metal_orfs, other_metal_orfs)
        shared_orfs_not_in_gt <- setdiff(shared_orfs, shared_ORFs_corr_gt)
        shared_ORFs_not_corr <- union(shared_ORFs_not_corr, shared_orfs_not_in_gt)
      }
    }
    
    # Count ORFs
    shared_ORFs_corr_gt_count <- length(shared_ORFs_corr_gt)
    shared_ORFs_not_corr_count <- length(shared_ORFs_not_corr)
    
    # Unique ORFs are the total ORFs minus the ORFs that are in either shared set
    unique_ORFs <- setdiff(metal_orfs, union(shared_ORFs_corr_gt, shared_ORFs_not_corr))
    unique_ORFs_count <- length(unique_ORFs)
    
    # Add row to results data frame
    results_df <- rbind(results_df, data.frame(
      metal = metal,
      total_ORFs = total_ORFs,
      unique_ORFs = unique_ORFs_count,
      shared_ORFs_not_corr = shared_ORFs_not_corr_count,
      shared_ORFs_corr_gt = shared_ORFs_corr_gt_count,
      correlation_threshold = corr_threshold
    ))
  }
}

# Show results data frame
print(results_df)

### Visualise

results_df_forplot <- results_df%>%
              mutate(fraction_unique = unique_ORFs/total_ORFs,
                     fraction_shared_with_atleast1corr_metal = shared_ORFs_corr_gt/total_ORFs,
                     fraction_shared_but_nocorr_metal = shared_ORFs_not_corr/total_ORFs,
                     fraction_explbymetmetcorr_of_all_shared = shared_ORFs_corr_gt/(shared_ORFs_corr_gt+shared_ORFs_not_corr))%>%
              na.omit()
                      
# Create the plot
pdf(paste0(plot_dir,"/fraction_sharedDAproteins_explained_by_metmetcorr_vscorrthresh.pdf"),height = 5,width=7.3)
ggplot(results_df_forplot, aes(x = correlation_threshold, 
                               y = fraction_explbymetmetcorr_of_all_shared, color = metal)) +
  geom_point(size = 1)+
  geom_line(alpha = 0.8) +
  scale_color_manual(values = colkey_Ele) +
  labs(
    y = "fraction of DA proteins explained \nby atleast 1 metal-metal correlation",
    x = "Correlation Threshold",
    color = ""
  )+
  theme_metallica()
dev.off()


pdf(paste0(plot_dir,"/metalwise_unique_sharedexpl_sharedunexpl_hits_barplot.pdf"),height = 5,width=7)
ggplot(filter(results_df_forplot,correlation_threshold == 0.8),
       aes(x = metal))+
  # add total DA proteins
  geom_bar(aes(y = total_ORFs,
               colour = metal),stat = "identity",
           width = 0.6,alpha = 0)+
  # add number that are unique 
  geom_bar(aes(y = unique_ORFs,
               fill = metal),stat = "identity",width = 0.6,alpha = 1)+
  # add % label of unique
  #geom_label(aes(y = unique_ORFs/2,label = fraction_unique))+
  # add number that are shared and the metal-metal correlation is > 0.8 and pvalue < 0.05
  geom_bar(aes(y = unique_ORFs+shared_ORFs_corr_gt,
               fill = metal),stat = "identity", width = 0.6,alpha = 0.4)+
  # add % label of shared and explained by metal-metal correlation
  geom_label(aes(y = unique_ORFs + shared_ORFs_corr_gt/2,
                label = paste0(round(100*fraction_explbymetmetcorr_of_all_shared,0),"%")))+
  # add number that are shared and the metal-metal correlation is < 0.8 and pvalue > 0.05 i.e. not explained by metal-metal correlation
  geom_bar(aes(y = unique_ORFs+shared_ORFs_corr_gt+shared_ORFs_not_corr,
               fill = metal),stat = "identity",width = 0.6,alpha = 0.1)+
  # add % label of shared but not explained by metal-metal correlation
  #geom_label(aes(y = unique_ORFs + shared_ORFs_corr_gt +shared_ORFs_not_corr +50,label = fraction_shared_with_atleast1corr_metal))+
  scale_colour_manual(values = colkey_Ele)+
  scale_fill_manual(values = colkey_Ele)+
  theme_metallica()+
  labs(x= '', colour = "", fill = "", y = "number of proteins")
dev.off()



