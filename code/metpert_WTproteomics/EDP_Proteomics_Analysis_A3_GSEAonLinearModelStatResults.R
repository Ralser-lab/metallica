################################################################################################
### EDP Script - A3 - Gene Set Enrichment Analysis ### hypergeometric test and Stouffer GSEA ###
################################################################################################

#### Script to make fit linar models along concentration gradient for each protein-element combination

## Make Sure you have run Scripts "EPP Analysis" 0, 1A and 2A before running this 
# This script will do the following :
# 1. 
# 2. 
# 3. 

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

plot_dir <- paste0(metpert_WTproteomics_dir,"/output/plots/extracellularconc_vs_proteinabundance/lmfits/gsea_piano")
dir.create(plot_dir, recursive = T)

###########################################################################################################################
### Read in results of statistical model fits - both intracellular and extracellular concentration vs protein abundance ###
###########################################################################################################################

lmfit_envconcVsprotabun_DAres <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),
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
  
  eles <- unique(hits_df$Element)
  
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
    
    sn <- make_twocolumn_Sankey_plot(ele_HyperGSA_res_df, gs)
    
    htmlwidgets::saveWidget(as_widget(sn), paste0(plot_dir,"/HyperGSA/elewise/HyperGSA_ele_",type_of_lmresults,"_Sankey",gs,".html"))
  }
}

run_elewise_hyperGSA(intra_hits_df,"intrametVsprotabun")
run_elewise_hyperGSA(extra_hits_df,"extrametVsprotabun")

###############################################################################################################################################
#################  Enrichments on whole ranked vector per element -ie.e which proteins go up vs down in each metal ############################
###############################################################################################################################################




##################################################################################################
#################  Enrichments on element-perturbation eg "Fe Excess" ############################
##################################################################################################

## define excess and depletion and what the fold change in excess vs fold change in depletion is 

Excess_Depletion_df <- lmfit_envconcVsprotabun_DAres %>%
                        mutate(Perturbation = ifelse(Element.Concentration > 1.2, "Excess",
                                                     ifelse(Element.Concentration < (1/1.2),"Depletion","Control")),
                               Ele_Pert = paste(Element, Perturbation, sep = " ")) %>%
                        group_by(Element, Perturbation, Genes) %>%
                        mutate(log2fc_elepert_direction = ifelse(all(is.na(Log2FC_vs_AE)), NA,
                                                                 ifelse(abs(max(Log2FC_vs_AE, na.rm=T)) > abs(min(Log2FC_vs_AE, na.rm=T)),
                                                                        max(Log2FC_vs_AE, na.rm=T),
                                                                        min(Log2FC_vs_AE, na.rm=T)))) %>%
                        ungroup() %>%
                        dplyr::select(Element, Genes, Perturbation, Ele_Pert, PValAdj_BH, log2fc_elepert_direction, Significant) %>%
                        unique() %>%
                        mutate(log2fc_elepert_direction = ifelse(is.infinite(log2fc_elepert_direction), NA, log2fc_elepert_direction)) %>%
                        na.omit()


# all unique element perturbations tested
ele_perts <- unique(Excess_Depletion_df$Ele_Pert)


######################################################################################
### Hypergeometric test on genes differentially expressed in element-perturbations ###
######################################################################################

# Define an empty data frame to store results
elepert_HyperGSA_res_df <- data.frame()

# Loop over each element perturbation
for(ep in 1:length(ele_perts)){
  
  # Split data for the current element perturbation into upregulated and downregulated genes
  df2HyperGSA_up <- Excess_Depletion_df %>%
                    filter(Ele_Pert == ele_perts[ep], log2fc_elepert_direction > 0) %>%
                    mutate(ORF = as.character(lapply(Genes,convert_GeneName2ORF))) %>%
                    dplyr::select(ORF,Significant)%>%
                    filter(ORF!="character(0)")
  
  df2HyperGSA_down <- Excess_Depletion_df %>%
                      filter(Ele_Pert == ele_perts[ep], log2fc_elepert_direction < 0) %>%
                      mutate(ORF = as.character(lapply(Genes,convert_GeneName2ORF))) %>%
                      dplyr::select(ORF,Significant)%>%
                      filter(ORF!="character(0)")
  
  # Run hypergeometric tests for the current element perturbation
  hres_up <- get_plot_HyperGSA(df2HyperGSA_up, paste0("HyperGSA_overexpressed_", ele_perts[ep], ".pdf"), EnrichPV_thresh = 0.05)
  hres_down <- get_plot_HyperGSA(df2HyperGSA_down, paste0("HyperGSA_underexpressed_", ele_perts[ep], ".pdf"), EnrichPV_thresh = 0.05)
  
  # Add the results to the overall results data frame
  if(length(hres_up) > 0){
    hres_up$Ele_Pert <- ele_perts[ep]
    hres_up$Direction <- "OE"
    elepert_HyperGSA_res_df <- rbind(elepert_HyperGSA_res_df, hres_up)
  }
  
  if(length(hres_down) > 0){
    hres_down$Ele_Pert <- ele_perts[ep]
    hres_down$Direction <- "UE"
    elepert_HyperGSA_res_df <- rbind(elepert_HyperGSA_res_df, hres_down)
  }
}

write.csv(elepert_HyperGSA_res_df, paste0(outputtables_dir,"/HyperGSA_on_lmfit_envconcVsprotabun_elepertwise.csv"), row.names = F)

###################################################################################
### StoufferGSA test on genes differentially expressed in element-perturbations ###
###################################################################################

### Run on Excess and Depletion separately ###

StoufferGSA_res_df <- data.frame()

for(ep in 1:length(ele_perts)){
  
    forStoufferGSA <- Excess_Depletion_df%>%
                      data.frame()%>%
                      filter(Ele_Pert == ele_perts[ep])%>%
                      mutate(ORF = as.character(lapply(Genes,convert_GeneName2ORF)))%>%
                      filter(ORF!="character(0)")%>%
                      dplyr::select(ORF,PValAdj_BH,log2fc_elepert_direction)%>%
                      na.omit()%>%
                      as.data.frame()
    
    sres <- get_plot_StoufferGSA(forStoufferGSA,paste("StoufferGSA",unlist(strsplit(ele_perts[ep]," "))[[1]],
                                                      unlist(strsplit(ele_perts[ep]," "))[[2]],"lmfit DEres"),EnrichPV_thresh=0.05)
    
    if(length(sres) > 0){
      
      sres$Element <- unlist(strsplit(ele_perts[ep]," "))[[1]]
      sres$Perturbation <- unlist(strsplit(ele_perts[ep]," "))[[2]]
      
      StoufferGSA_res_df <- rbind(StoufferGSA_res_df,sres)
    }
}
write.csv(StoufferGSA_res_df, paste0(outputtables_dir,"/StoufferGSA on lmfit DEres elewise.csv"))


###################################################################
### Make Sankey plots ### - for results of Stouffer enrichments ###
###################################################################

for(gs in 1:length(gsets)){
  
  sn <- plot_Sankey_DE(StoufferGSA_res_df,
                       GeneSetType = gsetnames[[gs]],
                       minGsetSize = 0,
                       maxGsetSize = Inf)
  
  
  
  htmlwidgets::saveWidget(as_widget(sn), paste0(plot_dir,"/EPP_Stouffer_Sankey",gsetnames[[gs]],".html"))
}

###############################
### fgsea ### Ele Pert wise ###
###############################

ele_perts <- unique(Excess_Depletion_df$Ele_Pert)

fgseaGSA_res_df <- data.frame()

for(ep in 1:length(ele_perts)){
  
  forfgseaGSA <- Excess_Depletion_df%>%
                    data.frame()%>%
                    filter(Ele_Pert == ele_perts[ep])%>%
                    mutate(ORF = as.character(lapply(Genes,convert_GeneName2ORF)))%>%
                    filter(ORF!="character(0)")%>%
                    dplyr::select(ORF,log2fc_elepert_direction,PValAdj_BH)%>%
                    na.omit()%>%
                    as.data.frame()
  
  sres <- get_plot_fgseaGSA(forfgseaGSA,paste("fgseaGSA",unlist(strsplit(ele_perts[ep]," "))[[1]],
                                                    unlist(strsplit(ele_perts[ep]," "))[[2]],"lmfit DEres"),
                           EnrichPV_thresh=0.05)
  
  if(length(sres) > 0){
    
    sres$Element <- unlist(strsplit(ele_perts[ep]," "))[[1]]
    sres$Perturbation <- unlist(strsplit(ele_perts[ep]," "))[[2]]
    
    fgseaGSA_res_df <- rbind(fgseaGSA_res_df,sres)
  }
}
write.csv(fgseaGSA_res_df, paste0(outputtables_dir,"/fgseaGSA on lmfit DEres elewise.csv"))


###################################################################
### Make Sankey plots ### - for results of fgsea enrichments ###
###################################################################

dir.create(paste0(plot_dir,"/fgsea"))

for(gs in 1:length(gsets)){
  
  sn <- plot_Sankey_DE(fgseaGSA_res_df,
                       GeneSetType = gsetnames[[gs]],
                       minGsetSize = 0,
                       maxGsetSize = Inf)
  
  
  
  htmlwidgets::saveWidget(as_widget(sn), paste0(plot_dir,"/fgsea/fgsea_Sankey",gsetnames[[gs]],".html"))
}






##################################################
### Hypergeometric ### least complex model wise ###
##################################################

mts <- unique(lmfit_envconcVsprotabun_DAres$LeastComplexModel)
mts <- mts[-which(mts == "null")]

HyperGSA_res_df <- data.frame()

for(e in 1:length(eles)){
  
  for(m in 1:length(mts)){
    
    df2HyperGSA <- filter(lmfit_envconcVsprotabun_DAres,Element ==eles[e])%>%
      dplyr::select(Genes,LeastComplexModel,Significant)%>%
      mutate(Significant = ifelse(LeastComplexModel == mts[m] & Significant ==1,1,0))%>%
      unique()%>%
      mutate(ORF = as.character(lapply(Genes,convert_GeneName2ORF)))%>%
      dplyr::select(ORF,Significant)%>%
      filter(ORF!="character(0)")
    
    if(sum(df2HyperGSA$Significant) > 5){
      hres<-get_plot_HyperGSA(df2HyperGSA,paste("HyperGSA",eles[e],mts[m],"lmfitDEres.pdf"),
                              EnrichPV_thresh=0.05)
      
      if(length(hres)>0){
        
        hres$Element = eles[e]
        hres$LeastComplexModel = mts[m]
        
        HyperGSA_res_df<-rbind(HyperGSA_res_df,hres)
      }
    }
  }
  
}
write.csv(HyperGSA_res_df,"HyperGSA on lmfit DEres elewise modelwise.csv")


