
#`---
#`  Title: "metpertWTproteomics_run_elementwise_GSA.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 1 June 2023
#`  Description: runs gene set analysis per element - i.e. no categorisation of depletion or excess
#`---

#################################################################################################
###  Gene Set Enrichment Analysis ### Element wise - no categorisation of excess or depletion ###
#################################################################################################

#### Script to make fit linar models along concentration gradient for each protein-element combination


source("/Users/aulakhs/Documents/Ralser Lab/metallica/code/common_code/initialise_common_paths.R")

# general

source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/layout_conversion.R"))
source(paste0(code_dir,"/common_code/input_processed_databases_publisheddatasets.R"))


# script specific

source(paste0(code_dir,"/metpert_WTproteomics/metpertWTproteomics_0_libraries_functions_genesetenrichments.R"))

#################
### Set Paths ###
#################

outputtables_dir <-paste0(metpert_WTproteomics_dir,"/output/tables/GSA")
dir.create(outputtables_dir, recursive = T)

plot_dir <- paste0(metpert_WTproteomics_dir,"/output/plots/comparison_extra_vs_intra/lmfits/gsea_piano")
dir.create(plot_dir, recursive = T)

###########################################################################################################################
### Read in results of statistical model fits - both intracellular and extracellular concentration vs protein abundance ###
###########################################################################################################################

lmfit_envconcVsprotabun_DAres <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/extracellularconc_vs_proteinabundance/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),
                                          stringsAsFactors = F)

lmfit_intrametconcVsprotabun_DAres <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/intracellularconc_vs_proteinabundance/lmfit_intracellular_metalconc_vs_log2PQ_with_SigNotSig_AdjPVthresh_0.05_FCthresh_1.5.csv"),
                                               stringsAsFactors = F)


# define Gene list that is significant per element either along intracellular or extracellular metal concentration
extra_hits_df <- unique(lmfit_envconcVsprotabun_DAres[,c("Element", "Significant","ORF","Genes")])
intra_hits_df <- unique(lmfit_intrametconcVsprotabun_DAres[,c("Element","Significant","ORF","Genes")])

#########################################################################################################################
### Hypergeometric ### Elewise ### using hits along  intracellular and extracellular metal concentration - separately ###
#########################################################################################################################

dir.create(paste0(plot_dir,"/HyperGSA/elewise"), recursive = T)

run_elewise_hyperGSA <- function(hits_df,type_of_lmresults){
  
  ele_HyperGSA_res_df <- data.frame()
  
  eles <- hits_df$Element
  
  for(e in 1:length(eles)){
    
    df2HyperGSA <- filter(hits_df, Element == eles[e])%>%
      dplyr::select(ORF,Significant)
    
    hres <- get_plot_HyperGSA(df2HyperGSA,
                              EnrichPV_thresh=0.05)
    
    if(length(hres)>0){
      
      hres$Element = eles[e]
      
      ele_HyperGSA_res_df <- rbind(ele_HyperGSA_res_df,hres)
    }
  }
  
  write.csv(ele_HyperGSA_res_df,paste0(outputtables_dir,"/HyperGSA_on_lmfit_",type_of_lmresults,"_elewise.csv"), row.names = F)
  
  ### Visualize ele wise Hypergeometric test results ### Sankey Plot ### 2 column ###
  
  # Get unique Gset.Type 
  gsetnames <- unique(ele_HyperGSA_res_df$Gset.Type)
  # Loop over unique Gset.Type
  for(gs in gsetnames) {
    
    sn <- plot_twocolumn_Sankey(ele_HyperGSA_res_df, gs)
    
    htmlwidgets::saveWidget(as_widget(sn), paste0(plot_dir,"/HyperGSA/elewise/HyperGSA_ele_",type_of_lmresults,"_Sankey",gs,".html"))
  }
  
  return(ele_HyperGSA_res_df)
}

intracell_hits_HyperGSAres <- run_elewise_hyperGSA(intra_hits_df,"intrametVsprotabun")
extracell_hits_HyperGSAres<- run_elewise_hyperGSA(extra_hits_df,"extrametVsprotabun")

##############################################################################################################################
### Make a combined 3 column Sankey plot of enrichments in proteins DA along intracellular and extracellular conc gradients###
##############################################################################################################################

all_hits_HyperGSAres <- rbind(intracell_hits_HyperGSAres,extracell_hits_HyperGSAres)%>%
                        filter(!Gset.Term.Enriched %in% c("biological_process"))


enriched_gsnames <- unique(all_hits_HyperGSAres$Gset.Type)

# Loop over unique Gset.Type

for(gs in 1:length(enriched_gsnames)){
  
  if(enriched_gsnames[gs] %in% gsetnames){
      sn <- plot_threecolumn_Sankey(all_hits_HyperGSAres, enriched_gsnames[gs])
      htmlwidgets::saveWidget(as_widget(sn), paste0(plot_dir,"/HyperGSA/elewise/HyperGSA_ele_extint_Sankey",enriched_gsnames[gs],".html"))
      plotly::export(sn, file = paste0(plot_dir,"/HyperGSA/elewise/HyperGSA_ele_extint_Sankey",enriched_gsnames[gs],".pdf"))
      
  }

}
