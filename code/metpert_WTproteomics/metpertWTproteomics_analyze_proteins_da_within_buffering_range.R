#`---
#`  Title: "metpertWTproteomics_20_analyze_proteins_da_within_buffering_range.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 18 May April 2023 
#`  Description: Script to determine which differentially abundant proteins showed changes in expression level before intracellular metal concentrations changed
#`---

#############################################
### source paths functions and libraries  ###
#############################################

# general
source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/layout_conversion.R"))

# experiment specific


#################
### set paths ###
#################

# inputs

# outputs

plot_dir <- paste0(metpert_WTproteomics_dir,"/output/plots/metallomics_proteomics_coanalysis")
dir.create(plot_dir, recursive = T)

proteomics_data_env <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/extracellularconc_vs_proteinabundance/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),stringsAsFactors = F)%>%
  dplyr::select(BioSpecID,Element, Genes,ORF,PValAdj_BH, Mean_of_reps_Log2FC_vs_AE, Significant)%>%
  unique()%>%
  filter(Element != "Mn")


total_da_env <- proteomics_data_env %>%
  dplyr::select(Element,Genes,Significant)%>%
  group_by(Element,Significant)%>%
  unique()%>%
  group_by(Element)%>%
  summarise(num_da = sum(Significant))


metallomicsbuffering_data <- read.csv(paste0(metpert_WTmetallomics_dir,"/output/tables/metallomics_buffering_summary.csv"),stringsAsFactors = F)%>%
  filter(element_measured %in% filter(total_da, num_da >1)$Element & element_measured != "Cu")

max_change_in_relmetquant = 2^max(abs(log2(metallomicsbuffering_data$mean_Ratio_to_AEngperwell))) + 0.1
#mean_Ratio_to_AEngperwell is the name of the relative change in intracellular conc

buff_thresh <- c(seq(1.0,max_change_in_relmetquant,by = 0.1))


frac_inside_buffrange <- vector()

for(b in 1:length(buff_thresh)){
  
  metprot_buffrange_bt <- merge(proteomics_data, metallomicsbuffering_data, by = "BioSpecID")%>%
    mutate(da_inside_buffrange = ifelse(abs(log2(mean_Ratio_to_AEngperwell)) <= log2(buff_thresh[b]) &
                                          Significant == 1 &
                                          abs(Mean_of_reps_Log2FC_vs_AE) >= log2(1.10), T, F))
  
  
  metprot_buffrange_smry_bt <- metprot_buffrange_bt%>%
    dplyr::select(Element,Genes,da_inside_buffrange)%>%
    unique()%>%
    group_by(Element)%>%
    summarise(num_da_inside_buffrange = sum(da_inside_buffrange, na.rm = T))
  
  
  frac_inside_buffrange_bt <- merge(metprot_buffrange_smry_bt,total_da, by = "Element")%>%
    mutate(fraction_inside_buffrange = num_da_inside_buffrange/num_da,
           buff_threshold = abs(1-buff_thresh[b]))
  
  frac_inside_buffrange <- rbind(frac_inside_buffrange,frac_inside_buffrange_bt)
  
}

pdf(paste0(plot_dir,"/fraction_da_vs_bufferingthreshold.pdf"),width = 7,height = 5)
ggplot(na.omit(frac_inside_buffrange),
       aes(x = 100*buff_threshold,
           y = fraction_inside_buffrange,
           colour = Element))+
  geom_point()+
  geom_line(alpha = 0.7)+
  theme_metallica()+
  scale_color_manual(values = colkey_Ele)+
  labs(x = "percentage change threshold for relative change \nin intracellular metal",
       y = "fraction proteins \ndifferentially abundant")
  

dev.off()


pdf(paste0(plot_dir,"/fraction_da_vs_bufferingthreshold_xlim.pdf"),width = 7,height = 5)
ggplot(na.omit(frac_inside_buffrange),
       aes(x = 100*buff_threshold,
           y = fraction_inside_buffrange,
           colour = Element))+
  geom_point()+
  geom_line(alpha = 0.7)+ 
  theme_metallica()+
  scale_color_manual(values = colkey_Ele)+
  xlim(0, 200)+
  labs(x = "percentage change \nin intracellular metal",
       y = "fraction proteins \ndifferentially abundant")
dev.off()

#########################################################
### based on DA proteins along cellular concentration ###
#########################################################

proteomics_data_cell <- read.csv(paste0(metpert_WTproteomics_dir,
                                        "/output/tables/intracellularconc_vs_proteinabundance/lmfit_intracellular_metalconc_vs_log2PQ_with_SigNotSig_AdjPVthresh_0.05_FCthresh_1.5.csv"),
                                 stringsAsFactors = F)%>%
  dplyr::select(Element, Genes,ORF,PValAdj_BH, Significant, median_relative_intracellular_concentration,median_log2_foldchangevsAE)%>%
  unique()


total_da_cell <- proteomics_data_cell %>%
  dplyr::select(Element,Genes,Significant)%>%
  group_by(Element,Significant)%>%
  unique()%>%
  group_by(Element)%>%
  summarise(num_da = sum(Significant))


max_change_in_relmetquant = 2^max(abs(log2(proteomics_data_cell$median_relative_intracellular_concentration))) + 0.1

buff_thresh <- c(seq(1.0,max_change_in_relmetquant,by = 0.025))
# Now, perform similar calculations for cell data
frac_inside_buffrange_cell <- vector()

for(b in 1:length(buff_thresh)){
  
  metprot_buffrange_cell <- proteomics_data_cell %>%
    na.omit()%>%
    mutate(da_inside_buffrange = ifelse(abs(log2(median_relative_intracellular_concentration)) <= log2(buff_thresh[b]) &
                                          Significant == 1 &
                                          abs(median_log2_foldchangevsAE) >= log2(1.20), TRUE, FALSE))
  
  metprot_buffrange_smry_cell <- metprot_buffrange_cell %>%
    dplyr::select(Element, Genes, da_inside_buffrange) %>%
    unique() %>%
    group_by(Element) %>%
    summarise(num_da_inside_buffrange = sum(da_inside_buffrange, na.rm = TRUE))
  
  frac_inside_buffrange_cell_bt <- merge(metprot_buffrange_smry_cell, total_da_cell, by = "Element") %>%
    mutate(fraction_inside_buffrange = num_da_inside_buffrange / num_da,
           buff_threshold = abs(1 - buff_thresh[b]))
  
  frac_inside_buffrange_cell <- rbind(frac_inside_buffrange_cell, frac_inside_buffrange_cell_bt)
  
}

# Assuming 'plot_dir' and 'colkey_Ele' are already defined
pdf(paste0(plot_dir, "/fraction_da_vs_bufferingthreshold_cell.pdf"), width = 7, height = 5)
ggplot(na.omit(frac_inside_buffrange_cell),
       aes(x = 100 * buff_threshold,
           y = fraction_inside_buffrange,
           colour = Element)) +
  geom_point() +
  geom_line(alpha = 0.7) +
  theme_metallica() +
  scale_x_continuous(limits = c(0,100))+
  scale_color_manual(values = colkey_Ele) +
  labs(x = "Percentage change threshold for relative change \nin intracellular metal",
       y = "Fraction proteins \ndifferentially abundant (Cellular Data)")
dev.off()








