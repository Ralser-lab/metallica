#`---
#`  Title: "metpertWTproteomics_compare_lmfits_intracellVsprotein_extracellVsprotein_hits.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 30 May 2023 
#`  Description: Script to fit linear models to protein abundance data along intracellular metal concentraion

#############################################
### source paths functions and libraries  ###
#############################################
# general
source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/layout_conversion.R"))


plot_dir <- paste0(metpert_WTproteomics_dir,"/output/plots/comparison_extra_vs_intra")
dir.create(plot_dir, recursive = T)

outputtables_dir <- paste0(metpert_WTproteomics_dir,"/output/tables/comparison_extra_vs_intra")
dir.create(outputtables_dir, recursive = T)


# Read the data
lmfitres_extracell <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/extracellularconc_vs_proteinabundance/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),stringsAsFactors = F) %>%
  dplyr::select(Genes,ORF, BioSpecID, Element, Element.Concentration,Significant, PValAdj_BH,Log2FC_vs_AE, LeastComplexModel)%>%
  unique()%>%
  mutate(gene_element = paste(Genes, Element))


lmfitres_intracell <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/intracellularconc_vs_proteinabundance/lmfit_intracellular_metalconc_vs_log2PQ_with_SigNotSig_AdjPVthresh_0.05_fcthresh_1.5.csv"),stringsAsFactors = F) %>%
  dplyr::select(Genes,ORF, Element,Significant, median_relative_intracellular_concentration,PValAdj_BH,median_log2_foldchangevsAE, LeastComplexModel)%>%
  mutate(gene_element = paste(Genes, Element))

# Filter rows where Significant == 1 to create filtered dataframes
lmfitres_extracell_filtered <- unique(dplyr::filter(lmfitres_extracell[,c("Genes","ORF", "Element","Significant", "PValAdj_BH","Log2FC_vs_AE","gene_element")], Significant == 1))
lmfitres_intracell_filtered <- unique(dplyr::filter(lmfitres_intracell[,c("Genes","ORF", "Element","Significant", "PValAdj_BH","median_log2_foldchangevsAE","gene_element")], Significant == 1))

# Assign a source label to each filtered dataframe
lmfitres_extracell_filtered$Source <- "extracell"
lmfitres_intracell_filtered$Source <- "intracell"

# Combine the two filtered dataframes into one
combined_df <- dplyr::bind_rows(lmfitres_extracell_filtered, lmfitres_intracell_filtered)

# Determine unique ORFs per Source and Element
unique_orfs_df <- combined_df %>% 
  dplyr::group_by(Element, Source) %>% 
  dplyr::summarise(unique_orfs = length(unique(ORF)), .groups = "drop")

# Spread the Source column into multiple columns
spread_df <- unique_orfs_df %>% 
  tidyr::spread(Source, unique_orfs, fill = 0)

# Calculate the total unique ORFs and common ORFs
overlap_df <- combined_df %>% 
  dplyr::group_by(Element, ORF) %>% 
  dplyr::summarise(n_sources = n_distinct(Source), .groups = "drop") %>% 
  dplyr::group_by(Element) %>% 
  dplyr::summarise(
    total_unique_significant_orfs = n_distinct(ORF),
    common_orfs = sum(n_sources > 1),
    .groups = "drop"
  )

# Merge the spread and overlap data frames
overlap_summary_df <- merge(spread_df, overlap_df, by = "Element") %>% 
  dplyr::mutate(
    extracell_only = extracell - common_orfs,
    intracell_only = intracell - common_orfs
  )
# Calculate total proteins per element in both dataframes
total_proteins_extracell <- aggregate(ORF ~ Element, lmfitres_extracell, function(x) length(unique(x)))
total_proteins_intracell <- aggregate(ORF ~ Element, lmfitres_intracell, function(x) length(unique(x)))

# Rename the second column to "total_proteins"
names(total_proteins_extracell)[2] <- "total_proteins_extracell"
names(total_proteins_intracell)[2] <- "total_proteins_intracell"

# Merge total_proteins_extracell and total_proteins_intracell
total_proteins <- data.frame(merge(total_proteins_extracell, total_proteins_intracell, by = "Element", all = T))
total_proteins[is.na(total_proteins)] <-0
# Add a new column "total_proteins" which is the maximum of total_proteins_extracell and total_proteins_intracell
total_proteins$total_proteins <- apply(total_proteins[,c("total_proteins_extracell","total_proteins_intracell")], 1, max)

# Merge total_proteins with overlap_summary_df
overlap_summary_df <- merge(overlap_summary_df, total_proteins[,c("Element","total_proteins")], by = "Element")

# Convert to long format for plotting
overlap_summary_df_long <- tidyr::pivot_longer(
  overlap_summary_df,
  cols = c(common_orfs, extracell_only, intracell_only, total_proteins),
  names_to = "category",
  values_to = "count"
)

#########################
### Visualise results ###
#########################

# Create the plot

pdf(paste0(plot_dir,"/number_of_significant_proteins_intracellular_extracellular_lmfits.pdf"), width = 7.2,height = 5)
ggplot(overlap_summary_df_long, aes(x = Element, y = count, fill = category)) +
  geom_bar(data = subset(overlap_summary_df_long, category != "total_proteins"), 
           stat = "identity", position = "stack", width = 0.4, colour = "black", alpha = 0.7, size = 0.25) +
  geom_bar(data = subset(overlap_summary_df_long, category == "total_proteins"),
           stat = "identity", width = 0.4, fill = NA, color = "black", size = 0.25) +
  labs(x = " ", y = "number of proteins", fill = "linear model type",
       title = "number of significant proteins linear model fits extracellular and intracellular metal conc") +
  scale_fill_manual(values = c("#CD919E", "#CDC1C5","#BCD2EE"))+
theme_metallica()
dev.off()

overlap_summary_df <- overlap_summary_df%>%
                      mutate(percentage_hits_any = round(total_unique_significant_orfs/total_proteins,3)*100,
                             percentage_hits_intra_only = round(intracell_only/total_proteins,3)*100,
                             percentage_hits_extra_only = round(extracell_only/total_proteins,3)*100,
                             percentage_hits_both = round(common_orfs/total_proteins,3)*100,
                             percentage_hits_both_of_totalhits = round(common_orfs/total_unique_significant_orfs,3)*100
                             )

write.csv(overlap_summary_df,paste0(outputtables_dir,"/number_significant_lmresults_intracellular_extracellular.csv"),row.names = F)

########################################################################################
### Plot top 3 proteins (with lowest pvalues) that are unique in each set and common ###
########################################################################################

extracell_df <- dplyr::filter(combined_df, Source == "extracell")
intracell_df <- dplyr::filter(combined_df, Source == "intracell")

extracell_only <- dplyr::anti_join(unique(extracell_df[,c("ORF","Genes","Element","PValAdj_BH","Source","gene_element")]), 
                                   unique(intracell_df[,c("ORF","Genes","Element","PValAdj_BH","Source","gene_element")]), by = c("ORF", "Genes","Element"))
intracell_only <- dplyr::anti_join(unique(intracell_df[,c("ORF","Genes","Element","PValAdj_BH","Source","gene_element")]),
                                   unique(extracell_df[,c("ORF","Genes","Element","PValAdj_BH","Source","gene_element")]), by = c("ORF","Genes", "Element"))

da_in_both_df <- dplyr::inner_join(unique(intracell_df[,c("ORF","Genes","Element","Source","gene_element")]),
                               unique(extracell_df[,c("ORF","Genes","Element","PValAdj_BH","Source","gene_element")]),by = c("ORF", "Genes", "Element","gene_element"))

# For unique hits in extracellular vs protein abundance

top3_extracellonly <- extracell_only %>%
                      group_by(Element)%>%
                      top_n(-3,PValAdj_BH)
                    
top3_intracellonly <- intracell_only %>%
                      group_by(Element)%>%
                      top_n(-3,PValAdj_BH)

top3_both <- da_in_both_df %>%
             group_by(Element)%>%
             top_n(-3,PValAdj_BH)

pdf(paste0(plot_dir, "/trajectories_top3_perelement_da_extracell_only_linear.pdf"), width = 26, height = 4)
print(ggplot(filter(lmfitres_extracell, 
                    gene_element %in% unique(top3_extracellonly$gene_element) &LeastComplexModel =="linear"),
             aes(x = log2(Element.Concentration),
                 y = Log2FC_vs_AE,
                 color = BioSpecID)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", formula = y ~ x, color = "black") +
        facet_wrap(Element~Genes, scales = "free", nrow = 1) +
        scale_colour_manual(values = colkey_BioSpecID) +
        geom_hline(yintercept = 0 , size = 0.2)+
        geom_vline(xintercept = 0, size = 0.2)+
        theme_metallica() +
        theme(legend.position = "none") +
        labs(title = "Top 3 Trajectories DA in Extrcellular only - linear only", 
             x = "Extrcellular Metal Concentration", y = "Log2 Fold Change vs control"))
dev.off()

pdf(paste0(plot_dir, "/trajectories_top3_perelement_da_extracell_only_quadratic.pdf"), width = 12, height = 4)
print(ggplot(filter(lmfitres_extracell, 
                    gene_element %in% unique(top3_extracellonly$gene_element) & LeastComplexModel == "quadratic"),
             aes(x = log2(Element.Concentration),
                 y = Log2FC_vs_AE,
                 color = BioSpecID)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", formula = y ~ poly(x,2), color = "black") +
        facet_wrap(Element~Genes, scales = "free", nrow = 1) +
        scale_colour_manual(values = colkey_BioSpecID) +
        geom_hline(yintercept = 0 , size = 0.2)+
        geom_vline(xintercept = 0, size = 0.2)+
        theme_metallica() +
        theme(legend.position = "none") +
        labs(title = "Top 3 Trajectories DA in Extrcellular only - quadratic ", 
             x = "Extrcellular Metal Concentration", y = "Log2 Fold Change vs control"))
dev.off()

pdf(paste0(plot_dir, "/trajectories_top3_perelement_da_extracell_only_cubic.pdf"), width = 6, height = 4)
print(ggplot(filter(lmfitres_extracell, 
                    gene_element %in% unique(top3_extracellonly$gene_element) & LeastComplexModel == "cubic"),
             aes(x = log2(Element.Concentration),
                 y = Log2FC_vs_AE,
                 color = BioSpecID)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", formula = y ~ poly(x,3), color = "black") +
        facet_wrap(Element~Genes, scales = "free", ncol = 6) +
        scale_colour_manual(values = colkey_BioSpecID) +
        theme_metallica() +
        theme(legend.position = "none") +
        geom_hline(yintercept = 0 , size = 0.2)+
        geom_vline(xintercept = 0, size = 0.2)+
        labs(title = "Top 3 Trajectories DA in Extrcellular only - cubic", 
             x = "Extrcellular Metal Concentration", y = "Log2 Fold Change vs control"))

dev.off()

##########################
### Intracellular only ###
##########################

pdf(paste0(plot_dir, "/trajectories_top3_perelement_da_intracell_only_linear.pdf"), width = 26.5, height = 4)
print(ggplot(filter(lmfitres_intracell, 
                    gene_element %in% unique(top3_intracellonly$gene_element) & LeastComplexModel == "linear"),
             aes(x = log2(median_relative_intracellular_concentration),
                 y = median_log2_foldchangevsAE,
                 color = Element)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", formula = y ~ x, color = "black") +
        facet_wrap(Element~Genes, scales = "free", nrow = 1) +
        scale_colour_manual(values = colkey_Ele) +
        theme_metallica() +
        theme(legend.position = "none") +
        geom_hline(yintercept = 0 , size = 0.2)+
        geom_vline(xintercept = 0, size = 0.2)+
        labs(title = "Top 3 Trajectories DA in Intracellular only - linear", 
             x = "intracellular metal concentration", y = "Log2 Fold Change vs control"))
dev.off()


pdf(paste0(plot_dir, "/trajectories_top3_perelement_da_intracell_only_quadratic.pdf"), width = 23, height = 4)
print(ggplot(filter(lmfitres_intracell, 
                    gene_element %in% unique(top3_intracellonly$gene_element) & LeastComplexModel == "quadratic"),
             aes(x = log2(median_relative_intracellular_concentration),
                 y = median_log2_foldchangevsAE,
                 color = Element)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", formula = y ~ poly(x,2), color = "black") +
        facet_wrap(Element~Genes, scales = "free", nrow = 1) +
        scale_colour_manual(values = colkey_Ele) +
        theme_metallica() +
        theme(legend.position = "none") +
        geom_hline(yintercept = 0 , size = 0.2)+
        geom_vline(xintercept = 0, size = 0.2)+
        labs(title = "Top 3 Trajectories DA in Intracellular only - quadratic", 
             x = "intracellular metal concentration", y = "Log2 Fold Change vs control"))
dev.off()

pdf(paste0(plot_dir, "/trajectories_top3_perelement_da_intracell_only_cubic.pdf"), width = 29, height = 4)
print(ggplot(filter(lmfitres_intracell, 
                    gene_element %in% unique(top3_intracellonly$gene_element) & LeastComplexModel == "cubic"),
             aes(x = log2(median_relative_intracellular_concentration),
                 y = median_log2_foldchangevsAE,
                 color = Element)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", formula = y ~ poly(x,3), color = "black") +
        facet_wrap(Element~Genes, scales = "free", nrow = 1) +
        scale_colour_manual(values = colkey_Ele) +
        theme_metallica() +
        theme(legend.position = "none") +
        geom_hline(yintercept = 0 , size = 0.2)+
        geom_vline(xintercept = 0, size = 0.2)+
        labs(title = "Top 3 Trajectories DA in Intracellular only - cubic", 
             x = "intracellular metal concentration", y = "Log2 Fold Change vs control"))
dev.off()


###########################################
### Both extracelluar and intracellular ###
###########################################

pdf(paste0(plot_dir, "/trajectories_top3_perelement_da_both_ext_linear.pdf"), width = 18, height = 4)
print(ggplot(filter(lmfitres_extracell, 
                    gene_element %in% unique(top3_both$gene_element) & LeastComplexModel == "linear"),
             aes(x = log2(Element.Concentration),
                 y = Log2FC_vs_AE,
                 color = BioSpecID)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", formula = y ~ x, color = "black") +
        facet_wrap(Element~Genes, scales = "free", nrow = 1) +
        scale_colour_manual(values = colkey_BioSpecID) +
        theme_metallica() +
        theme(legend.position = "none") +
        geom_hline(yintercept = 0 , size = 0.2)+
        geom_vline(xintercept = 0, size = 0.2)+
        labs(title = "Top 3 Trajectories DA in both - linaer", 
             x = "extracellular metal concentration", y = "Log2 Fold Change vs control"))
dev.off()

pdf(paste0(plot_dir, "/trajectories_top3_perelement_da_both_int_linear.pdf"), width = 18, height = 4)
print(ggplot(filter(lmfitres_intracell, 
                    gene_element %in% unique(top3_both$gene_element) & LeastComplexModel == "linear"),
             aes(x = log2(median_relative_intracellular_concentration),
                 y = median_log2_foldchangevsAE,
                 color = Element)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", formula = y ~ x, color = "black") +
        facet_wrap(Element~Genes, scales = "free", nrow = 1) +
        scale_colour_manual(values = colkey_Ele) +
        theme_metallica() +
        theme(legend.position = "none") +
        geom_hline(yintercept = 0 , size = 0.2)+
        geom_vline(xintercept = 0, size = 0.2)+
        labs(title = "Top 3 Trajectories DA in both - linear", 
             x = "intracellular metal concentration", y = "Log2 Fold Change vs control"))
dev.off()

## quadratic

# Plot for extracellular concentration with quadratic model
pdf(paste0(plot_dir, "/trajectories_top3_perelement_da_both_ext_quadratic.pdf"), width = 20, height = 4)
print(ggplot(filter(lmfitres_extracell, 
                    gene_element %in% unique(top3_both$gene_element) & LeastComplexModel == "quadratic"),
             aes(x = log2(Element.Concentration),
                 y = Log2FC_vs_AE,
                 color = BioSpecID)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "black") +
        facet_wrap(Element~Genes, scales = "free", nrow = 1) +
        scale_colour_manual(values = colkey_BioSpecID) +
        theme_metallica() +
        theme(legend.position = "none") +
        geom_hline(yintercept = 0 , size = 0.2)+
        geom_vline(xintercept = 0, size = 0.2)+
        labs(title = "Top 3 Trajectories DA in both - quadratic", 
             x = "extracellular metal concentration", y = "Log2 Fold Change vs control"))
dev.off()

# Plot for intracellular concentration with quadratic model
pdf(paste0(plot_dir, "/trajectories_top3_perelement_da_both_int_quadratic.pdf"), width = 19.5, height = 4)
print(ggplot(filter(lmfitres_intracell, 
                    gene_element %in% unique(top3_both$gene_element) & LeastComplexModel == "quadratic"),
             aes(x = log2(median_relative_intracellular_concentration),
                 y = median_log2_foldchangevsAE,
                 color = Element)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "black") +
        facet_wrap(Element~Genes, scales = "free", nrow = 1) +
        scale_colour_manual(values = colkey_Ele) +
        theme_metallica() +
        geom_hline(yintercept = 0 , size = 0.2)+
        geom_vline(xintercept = 0, size = 0.2)+
        theme(legend.position = "none") +
        labs(title = "Top 3 Trajectories DA in both - quadratic", 
             x = "intracellular metal concentration", y = "Log2 Fold Change vs control"))
dev.off()

## Cubic 

# Plot for extracellular concentration with cubic model
pdf(paste0(plot_dir, "/trajectories_top3_perelement_da_both_ext_cubic.pdf"), width = 8, height = 4)
print(ggplot(filter(lmfitres_extracell, 
                    gene_element %in% unique(top3_both$gene_element) & LeastComplexModel == "cubic"),
             aes(x = log2(Element.Concentration),
                 y = Log2FC_vs_AE,
                 color = BioSpecID)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", formula = y ~ poly(x, 3), color = "black") +
        facet_wrap(Element~Genes, scales = "free", nrow = 1) +
        scale_colour_manual(values = colkey_BioSpecID) +
        theme_metallica() +
        theme(legend.position = "none") +
        geom_hline(yintercept = 0 , size = 0.2)+
        geom_vline(xintercept = 0, size = 0.2)+
        labs(title = "Top 3 Trajectories DA in both - cubic", 
             x = "extracellular metal concentration", y = "Log2 Fold Change vs control"))
dev.off()

# Plot for intracellular concentration with cubic model
pdf(paste0(plot_dir, "/trajectories_top3_perelement_da_both_int_cubic.pdf"), width = 8, height = 4)
print(ggplot(filter(lmfitres_intracell, 
                    gene_element %in% unique(top3_both$gene_element) & LeastComplexModel == "cubic"),
             aes(x = log2(median_relative_intracellular_concentration),
                 y = median_log2_foldchangevsAE,
                 color = Element)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", formula = y ~ poly(x, 3), color = "black") +
        facet_wrap(Element~Genes, scales = "free", nrow = 1) +
        scale_colour_manual(values = colkey_Ele) +
        theme_metallica() +
        theme(legend.position = "none") +
        geom_hline(yintercept = 0 , size = 0.2)+
        geom_vline(xintercept = 0, size = 0.2)+
        labs(title = "Top 3 Trajectories DA in both - cubic", 
             x = "intracellular metal concentration", y = "Log2 Fold Change vs control"))
dev.off()


#############################################################################
### Export hitlist and protein quantity dataframe for ensemble clustering ###
#############################################################################

hitlist_for_ensemble_clustering <- combined_df%>%
                                   dplyr::select(ORF,Genes,Element,gene_element, Source)%>%
                                   unique()

write.csv(hitlist_for_ensemble_clustering, paste0(metpert_WTproteomics_dir,"/output/tables/ensemble_clustering/for_ensclustering_significant_along_extracell_or_intracell.csv"),row.names = F)

alldata_for_ensemble_clustering <- lmfitres_extracell%>%
                                   dplyr::select(BioSpecID, ORF, Genes, Element, Log2FC_vs_AE)
write.csv(alldata_for_ensemble_clustering, paste0(metpert_WTproteomics_dir,"/output/tables/ensemble_clustering/for_ensclustering_alldata_log2FCvsAE.csv"),row.names = F)


length(unique(filter(hitlist_for_ensemble_clustering,  Element == "Fe")$ORF))
length(unique(filter(hitlist_for_ensemble_clustering,  Element == "Fe" & Source == "extracell")$ORF))
length(unique(filter(hitlist_for_ensemble_clustering,  Element == "Fe" & Source == "intracell")$ORF))

