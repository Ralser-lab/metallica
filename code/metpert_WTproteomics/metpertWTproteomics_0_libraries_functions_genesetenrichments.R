
# Load required libraries
library(dplyr)
library(plotly)
library(doParallel)
library(piano)
library(visNetwork)


run_HyperGSA<-function(forHyperGSA,EnrichPV_thresh=0.05,
                            gsLim_low=3,gsLim_high=400,
                            gsLim_low_reacGEM = 3,
                            gsLim_high_reacGEM = 100
){ 
  
  # Function to carry out Hypergeometric Tests on Genesets of all databases,
  # report results in a df 
  # Input needs to be a df with 2 columns :
  # Column 1 : List of all ORFs measured # Background for HyperGSA
  # Column 2 : Boolean indicating which ORFs are significant / list to use for the hyperGSA
  
  colnames(forHyperGSA) <-  c("ORF","Use_for_HyperGSA")
  ORF_meas_universe = unique(forHyperGSA$ORF)
  
  All_gsets_HyperGSAres <- vector()
  
  for(gs in 1:length(gsets)){
    
    # load each gene set one by one
    gsc.forEnrich <- loadGSC(gsets[[gs]])
    
    # Check that there are more than 3 terms in the ORFs to use for HyperGSA
    
    if(sum(forHyperGSA$Use_for_HyperGSA)>3){
      
      pert_genes <- filter(forHyperGSA,as.logical(Use_for_HyperGSA))$ORF
      
      if(gsetnames[[gs]] %in%c("Sc_GEM","Reactome")){
        
        
        res.HyperGSA<- runGSAhyper(genes = unique(pert_genes), 
                                   universe=ORF_meas_universe,
                                   gsc=gsc.forEnrich,
                                   gsSizeLim = c(gsLim_low_reacGEM,gsLim_high_reacGEM),
                                   adjMethod = "BH")
        
      }else{
        res.HyperGSA<- runGSAhyper(genes = unique(pert_genes), 
                                   universe=ORF_meas_universe,
                                   gsc=gsc.forEnrich,
                                   gsSizeLim = c(gsLim_low,gsLim_high),
                                   adjMethod = "BH")
      }
      
      
      # filter and store HyperGSA enrichment results
      
      enrch.gsets<-data.frame(res.HyperGSA$resTab[,c(2:4)])
      enrch.gsets$Gset.name<-rownames(enrch.gsets)
      colnames(enrch.gsets)<-c("Adj.pval","Sig.in.geneset","NonSig.in.geneset","Gset.Term.Enriched")
      enrch.gsets<-enrch.gsets%>%
        mutate(Gene.Set.Size = Sig.in.geneset
               +NonSig.in.geneset)%>%
        mutate(Gset.Type=gsetnames[gs],
               Enrch.pval.threshold = EnrichPV_thresh)%>%
        filter(Adj.pval< EnrichPV_thresh)
      
      if(nrow(enrch.gsets) > 1 ){
        All_gsets_HyperGSAres<-rbind(All_gsets_HyperGSAres,enrch.gsets)
      }
    }
  }
  return(All_gsets_HyperGSAres)
}

########################################
### Helper functions for Sankey plot ###
########################################

hex2colWopacity <- function(hx) {
  if (hx == "gray") {
    x = rbind(red = 0, green = 0 , blue = 0 )
    x = rgb(x[1], x[2], x[3], max = 255, alpha = 1)
    return(x)
  } else {
    c <- col2rgb(hx)
    c <- rgb(c[1], c[2], c[3], max = 255, alpha = 150)
    return(c)
  }
}

###########################################################################################
### Sankey Plot for results of hyperGSA element wise ### 2 columns only ## no direction ###
###########################################################################################

plot_twocolumn_Sankey <- function(df, gset_type) {
  # Filter data frame based on Gset.Type
  df_type <- df %>% filter(Gset.Type == gset_type)
  
  # Add ID fields to your dataframe
  df_type <- df_type %>%
    mutate(Element_id = as.numeric(factor(Element)), 
           Gset_id = as.numeric(factor(Gset.Term.Enriched)) + max(Element_id))
  
  # Create a list with unique elements for the nodes and labels
  nodes <- data.frame(id = c(unique(df_type$Element_id), 
                             unique(df_type$Gset_id)),
                      name = c(as.character(unique(df_type$Element)), 
                               as.character(unique(df_type$Gset.Term.Enriched))))
  
  # Create links data frame
  links <- data.frame(source = df_type$Element_id - 1, # -1 because of zero-based index
                      target = df_type$Gset_id - 1, 
                      value = df_type$Gene.Set.Size,
                      color = sapply(df_type$Element, function(x) paste0(colkey_Ele[[x]], "80"))) # assign color from colkey_Ele with 0.8 opacity
  
  # Color for nodes
  node_colors <- c(sapply(unique(df_type$Element), function(x) colkey_Ele[[x]]), rep("rgba(0,0,0,0)", length(unique(df_type$Gset.Term.Enriched)))) # assign color from colkey_Ele for Elements and transparent for Gset.Term.Enriched
  
  # Generate Sankey plot using plotly
  p <- plot_ly(type = "sankey",
               domain = list(
                 x =  c(0,1),
                 y =  c(0,1)
               ),
               orientation = "h",
               
               node = list(
                 label =  nodes$name,
                 color =  node_colors,
                 pad = 15,
                 thickness = 20,
                 line = list(
                   color =  "black",
                   width =  0.5
                 )),
               
               link = list(
                 source =  links$source,
                 target =  links$target,
                 value =  links$value,
                 color =  links$color)
  ) %>% 
    layout(
      title =  paste("Sankey Diagram for", gset_type),
      font =  list(
        size = 15,
        family="Times"
      )
    )
  
  return(p)
}

##########################################################################
### 3 column Sankey plot for extrac vs prot and intrac vs protein hits ###
##########################################################################

plot_threecolumn_Sankey <- function(df,
                                    GeneSetType){
  
  ## create dataframe that links terms enriched along extracellular metal conc vs protein to each metal
  
  extra_hits_sankey <- filter(df, type_of_lmresults == "extracellular" & Gset.Type == GeneSetType)%>%
                       mutate(source = paste(Gset.Term.Enriched,"_ext_"),
                              target = paste(Element,"_metal_"),
                              value = 1,
                              link_colour = sapply(Element, function(x) paste0(colkey_Ele[[x]], "80")))
  
  ## create dataframe that links terms enriched along intracellular metal conc vs protein to each metal
  
  intra_hits_sankey <- filter(df, type_of_lmresults == "intracellular" & Gset.Type == GeneSetType)%>%
                      mutate(target = paste(Gset.Term.Enriched,"_int_"),
                             source = paste(Element,"_metal_"),
                             value = 1,
                             link_colour = sapply(Element, function(x) paste0(colkey_Ele[[x]], "80")))
                    
  # combine the two dataframe to generate nodes and links later 
  extr_intr_hits_df_for_sankey <- rbind(extra_hits_sankey,intra_hits_sankey)
  
  ## generate all nodes
  buffer <- 0.001 # adjust this value to get the desired spacing between y-coordinates
  
  sankey_nodes <- data.frame(name = c(as.character(extr_intr_hits_df_for_sankey$source), 
                                      as.character(extr_intr_hits_df_for_sankey$target)) %>% 
                               unique())%>%
    # Add node positions 
    # x coord - if its extracellular put it on the left
    mutate(x_coord = ifelse(grepl("_metal_",name),0.5,
                            ifelse(grepl("_ext_",name),0.1,0.9)))%>%
    mutate(n=1)%>%
    # y coord
    group_by(x_coord)%>%
    mutate(y_cood = length(x_coord),
           counter = cumsum(n)
    )%>%
    ungroup()%>%
    mutate(y_cood = (counter + buffer) / (y_cood + 2 * buffer) )%>% # add buffer to both the counter and y_cood to maintain the relative positioning
    dplyr::select(-n,-counter)
  
  
  extr_intr_hits_df_for_sankey$IDsource = match(extr_intr_hits_df_for_sankey$source, sankey_nodes$name)-1 
  extr_intr_hits_df_for_sankey$IDtarget = match(extr_intr_hits_df_for_sankey$target, sankey_nodes$name)-1
  
  node_color_vector <- ifelse(sankey_nodes$name == "Ca _metal_",colkey_Ele[["Ca"]],
                              ifelse(sankey_nodes$name == "Cu _metal_",colkey_Ele[["Cu"]],
                                     ifelse(sankey_nodes$name == "Fe _metal_",colkey_Ele[["Fe"]],
                                            ifelse(sankey_nodes$name == "K _metal_",colkey_Ele[["K"]],
                                                   ifelse(sankey_nodes$name == "Mg _metal_",colkey_Ele[["Mg"]],
                                                          ifelse(sankey_nodes$name == "Mn _metal_",colkey_Ele[["Mn"]],
                                                                 ifelse(sankey_nodes$name == "Mo _metal_",colkey_Ele[["Mo"]],
                                                                        ifelse(sankey_nodes$name == "Na _metal_",colkey_Ele[["Na"]],
                                                                               ifelse(sankey_nodes$name == "Zn _metal_",colkey_Ele[["Zn"]],
                                                                                      "gray")))))))))
  
  node_color_vector <- as.character(lapply(node_color_vector,hex2colWopacity))
  
  sankey_nodes <- sankey_nodes%>%
                  mutate(name = gsub(" _metal_","",name),
                         name = gsub(" _ext_","",name),
                         name = gsub(" _int_","",name))
  
  
  fig <- plotly::plot_ly(
    type = "sankey",
    node = list(
      label = sankey_nodes$name,
      x = sankey_nodes$x_coord,
      y = sankey_nodes$y_cood,
      color = node_color_vector,
      pad = 20,
      thickness = 20,
      line = list(
        color = "black",
        width = 0.5
      )
    ),
    link = list(
      source = extr_intr_hits_df_for_sankey$IDsource,
      target = extr_intr_hits_df_for_sankey$IDtarget,
      value =  extr_intr_hits_df_for_sankey$value,
      color = extr_intr_hits_df_for_sankey$link_colour
    ),
  )
  
  
  
  fig <- fig %>% layout(
    title = paste(GeneSetType),
    font = list(
      size = 15,
      family = "Times"
    )
  )
  
  fig
  
  return(fig)
  
}

###########################################################################################################
### 3 column Sankey plot for dataset to metal DA in dataset to gset term enriched in dataset-metal pair ###
###########################################################################################################

plot_threecolumn_dataset2metal2genesetterm_Sankey <- function(df,GeneSetType){
  
  ## create dataframe that links datasets to metals
  
  datasets_to_metals_sankey_df <- filter(df, Gset.Type == GeneSetType)%>%
                                  group_by(dataset,metal)%>%
                                  mutate(num_gsterms_enriched = length(Gset.Term.Enriched))%>%
                                  ungroup()%>%
                                  dplyr::select(dataset,metal,num_gsterms_enriched)%>%
                                  unique()%>%
                                  mutate(source = paste0(dataset,"_dataset_"),
                                          target = paste(dataset,metal,"dsmet_",sep="_"),
                                          value = num_gsterms_enriched,
                                          link_colour = sapply(metal, function(x) paste0(colkey_Ele[[x]], "80")))%>%
                                  dplyr::select(-dataset,-metal,-num_gsterms_enriched)
                                  
  ## create dataframe that links metal-dataset values to gset terms enriched
  
  metals_to_enrichedterms_sankey_df <- filter(df, Gset.Type == GeneSetType)%>%
                                       dplyr::select(dataset,metal,Gset.Term.Enriched)%>%
                                       unique()%>%
                                       mutate(source = paste(dataset,metal,"dsmet_",sep="_"),
                                               target = Gset.Term.Enriched,
                                               value = 1,
                                               link_colour = sapply(metal, function(x) paste0(colkey_Ele[[x]], "80")))%>%
                                       dplyr::select(-dataset, -metal, -Gset.Term.Enriched)
  
  # combine the two dataframe to generate nodes and links later 
  datasets2metals2gseterms_for_sankey <- rbind(datasets_to_metals_sankey_df,metals_to_enrichedterms_sankey_df)
  
  ## generate all nodes
  buffer <- 0.001 # adjust this value to get the desired spacing between y-coordinates
  
  sankey_nodes <- data.frame(name = c(as.character(datasets2metals2gseterms_for_sankey$source), 
                                      as.character(datasets2metals2gseterms_for_sankey$target)) %>% 
                  unique())%>%
    # Add node positions 
    # x coord - if its extracellular put it on the left
    mutate(x_coord = ifelse(grepl("_dsmet_",name),0.5,
                            ifelse(grepl("_dataset_",name),0.1,0.9)))%>%
    mutate(n=1)%>%
    # y coord
    group_by(x_coord)%>%
    mutate(y_cood = length(x_coord),
           counter = cumsum(n)
    )%>%
    ungroup()%>%
    mutate(y_cood = (counter + buffer) / (y_cood + 2 * buffer) )%>% # add buffer to both the counter and y_cood to maintain the relative positioning
    dplyr::select(-n,-counter)
  
  
  datasets2metals2gseterms_for_sankey$IDsource = match(datasets2metals2gseterms_for_sankey$source, sankey_nodes$name)-1 
  datasets2metals2gseterms_for_sankey$IDtarget = match(datasets2metals2gseterms_for_sankey$target, sankey_nodes$name)-1
  
  convert_nodenames_to_hexcolour <- function(x){
    
    colour = ifelse(grepl("Ca_dsmet_",x),colkey_Ele[["Ca"]],
                    ifelse(grepl("Cu_dsmet_",x),colkey_Ele[["Cu"]],
                           ifelse(grepl("Fe_dsmet_",x),colkey_Ele[["Fe"]],
                                  ifelse(grepl("K_dsmet_",x),colkey_Ele[["K"]],
                                         ifelse(grepl("Mg_dsmet_",x),colkey_Ele[["Mg"]],
                                                ifelse(grepl("Mn_dsmet_",x),colkey_Ele[["Mn"]],
                                                       ifelse(grepl("Mo_dsmet_",x),colkey_Ele[["Mo"]],
                                                              ifelse(grepl("Na_dsmet_",x),colkey_Ele[["Na"]],
                                                                     ifelse(grepl("Zn_dsmet_",x),colkey_Ele[["Zn"]],
                                                                            "gray")))))))))
    return(colour)
  }
  node_color_vector <- as.character(lapply(sankey_nodes$name,convert_nodenames_to_hexcolour))
  
  node_color_vector <- as.character(lapply(node_color_vector,hex2colWopacity))
  
  sankey_nodes <- sankey_nodes%>%
                  mutate(name = gsub("_dataset_","",name),
                         name = gsub("_dsmet_","",name),
                         name = gsub("metdepKOgrowth_","",name),
                         name = gsub("KOmetallomics_","",name),
                         name = gsub("OEmetallomics_","",name),
                         name = gsub("y5kmetalspecific_","",name),
                         name = gsub("metpertWTproteomics_","",name),
                         )
  
  
  fig <- plotly::plot_ly(
    type = "sankey",
    node = list(
      label = sankey_nodes$name,
      x = sankey_nodes$x_coord,
      y = sankey_nodes$y_cood,
      color = node_color_vector,
      pad = 20,
      thickness = 20,
      line = list(
        color = "black",
        width = 0.5
      )
    ),
    link = list(
      source = datasets2metals2gseterms_for_sankey$IDsource,
      target = datasets2metals2gseterms_for_sankey$IDtarget,
      value =  datasets2metals2gseterms_for_sankey$value,
      color = datasets2metals2gseterms_for_sankey$link_colour
    ),
  )
  
  
  
  fig <- fig %>% layout(
    title = paste(GeneSetType),
    font = list(
      size = 15,
      family = "Times"
    )
  )
  
  fig
  
  return(fig)
}
