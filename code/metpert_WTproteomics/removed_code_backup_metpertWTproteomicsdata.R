
get_plot_StoufferGSA <- function (df_forStouffer,plotname,EnrichPV_thresh=0.05){
  
  registerDoParallel(cores=length(gsets))
  ## Function to run Stouffer method on df with ORFs, qvalues (adjusted pvalues ) and log2fc 
  
  # Input dataframe should have :
  # Column 1 : ORF
  # Column 2 : Qvalues ( adjusted p values or only pvalues if you want to run it without multiple testing correction )
  # Column 3 : log(Fold Change) - vs your control condition or median / mean of all
  
  ###############################################
  # Run GSA with Stouffer method # Paralellized #
  ###############################################
  colnames(df_forStouffer) <- c("ORF","pv.Adj","Log2FC")
  ORF_meas_universe = unique(df_forStouffer$ORF)
  
  # Create the ranking statistics
  df_forStouffer$rank_stat <- -log10(df_forStouffer$pv.Adj) * sign(df_forStouffer$Log2FC)
  
  rank_vector <- df_forStouffer$rank_stat
  names(rank_vector) <- df_forStouffer$ORF
  
  All_gsets_Stoufferres <- foreach(gs = 1:length(gsets), .combine=rbind) %dopar%{
    
    get_Stouffer_enrichments <- function(gset,df_forStouffer){
      
      gsc.forEnrich <- loadGSC(gsets[[gs]])
      
      res.Stouffer <- runGSA(geneLevelStats = rank_vector,
                             geneSetStat = "stouffer",
                             gsc = gsc.forEnrich,
                             adjMethod = "BH")
      
      StouffResTab <- GSAsummaryTable(res.Stouffer)[,c(1,2,5,8)]
      colnames(StouffResTab)=c("Gset.Term.Enriched","Gene.Set.Size",
                               "pAdj_disdir_up","pAdj_disdir_down")
      
      StouffResTab$Gset.Type <- gsetnames[gs]
      
      if(nrow(StouffResTab)==0){
        dfE=cbind("No Significant Results",NA,NA,NA,gsetnames[gs])
        colnames(dfE)=c("Gset.Term.Enriched","Gene.Set.Size",
                        "pAdj_disdir_up","pAdj_disdir_down",
                        "Gset.Type")
      }
      
      return(StouffResTab)
    }
    gset_2func=gsets[[gs]]
    
    data.frame(get_Stouffer_enrichments(gset_2func,df_forStouffer))
    
  }
  return(All_gsets_Stoufferres)
}

### fgsea 
get_plot_fgseaGSA <- function (df_forfgsea,plotname,EnrichPV_thresh=0.05){
  
  registerDoParallel(cores=length(gsets))
  
  # Set column names
  colnames(df_forfgsea) <- c("ORF","Log2FC","p_value")
  
  # Create a named vector from the fold change, with gene IDs as names, and rank them
  rank_vector <- df_forfgsea$Log2FC
  names(rank_vector) <- df_forfgsea$ORF
  
  # Rank the vector in decreasing order
  rank_vector <- rank_vector[order(rank_vector, decreasing = TRUE)]
  
  ORF_meas_universe = unique(df_forfgsea$ORF)
  
  All_gsets_fgseares <- foreach(gs = 1:length(gsets), .combine=rbind) %dopar%{
    
    get_fgsea_enrichments <- function(gset,df_forfgsea){
      
      # load each gene set one by one
      gsc.forEnrich <- loadGSC(gsets[[gs]])
      
      res.fgsea <- runGSA(geneLevelStats = rank_vector,
                          geneSetStat = "fgsea",
                          gsc = gsc.forEnrich,
                          nPerm = 1000,
                          adjMethod = "BH")
      
      
      fgseaResTab <- GSAsummaryTable(res.fgsea)[,c(1,2,5,8)]
      colnames(fgseaResTab)=c("Gset.Term.Enriched","Gene.Set.Size",
                              "pAdj_disdir_up","pAdj_disdir_down")
      
      fgseaResTab$Gset.Type <- gsetnames[gs]
      
      if(nrow(fgseaResTab)==0){
        dfE=cbind("No Significant Results",NA,NA,NA,gsetnames[gs])
        colnames(dfE)=c("Gset.Term.Enriched","Gene.Set.Size",
                        "pAdj_disdir_up","pAdj_disdir_down",
                        "Gset.Type")
      }
      
      return(fgseaResTab)
    }
    gset_2func=gsets[[gs]]
    
    data.frame(get_fgsea_enrichments(gset_2func,df_forfgsea))
  }
  return(All_gsets_fgseares)
}
###################################################################
### 2 column Sankey for dataset to gene set term enriched plots ###
###################################################################

plot_twocolumn_dataset2GSterm_Sankey <- function(df, gset_type) {
  # Filter data frame based on Gset.Type
  df_type <- df %>% filter(Gset.Type == gset_type)
  
  # Add ID fields to your dataframe
  df_type <- df_type %>%
    mutate(dataset_id = as.numeric(factor(dataset)), 
           source = paste(dataset),
           target = paste(Gset.Term.Enriched),
           value = 1,
           Gset_id = as.numeric(factor(Gset.Term.Enriched)) + max(dataset_id),
           color = as.character(lapply(df_type$dataset,function(x) colkey_dataset[[x]] ) ))
  
  # Create a list with unique elements for the nodes and labels
  sankey_nodes <- data.frame(id = c(unique(df_type$dataset_id), 
                                    unique(df_type$Gset_id)),
                             name = c(as.character(unique(df_type$dataset)), 
                                      as.character(unique(df_type$Gset.Term.Enriched))))
  
  # Color for nodes
  node_colors <- c(sapply(unique(df_type$dataset), function(x) colkey_dataset[[x]]), rep("rgba(0,0,0,0)", length(unique(df_type$Gset.Term.Enriched)))) # assign color from colkey_Ele for Elements and transparent for Gset.Term.Enriched
  
  # ad node id number to dataframe 
  df_type$IDsource = match(df_type$source, sankey_nodes$name)-1 
  df_type$IDtarget = match(df_type$target, sankey_nodes$name)-1
  
  # Generate Sankey plot using plotly
  p <- plot_ly(type = "sankey",
               domain = list(
                 x =  c(0,1),
                 y =  c(0,1)
               ),
               orientation = "h",
               
               node = list(
                 label =  sankey_nodes$name,
                 color =  node_colors,
                 pad = 15,
                 thickness = 20,
                 line = list(
                   color =  "black",
                   width =  0.5
                 )),
               
               link = list(
                 source =  df_type$IDsource,
                 target =  df_type$IDtarget,
                 value =  df_type$value,
                 color =  df_type$color)
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

############################################
### Run HyperGSA on hits in each dataset ###
############################################

## function for running hyperGSA per element per dataset

run_dbwise_hyperGSA <- function(hits_df,dataset){
  
  HyperGSA_res_df <- data.frame()
  
  hres <- run_HyperGSA(data.frame(hits_df),
                            EnrichPV_thresh=0.05)
  
  if(length(hres)>0){
    
    hres$dataset = dataset
    
    HyperGSA_res_df <- rbind(HyperGSA_res_df,hres)
  }
  
  write.csv(HyperGSA_res_df,paste0(outputtables_dir,"/HyperGSA_on_lmfit_",dataset,"_datasetwise.csv"), row.names = F)
  
  return(HyperGSA_res_df)
}

hyperGSA_res_all_datasets <- vector()

for(ds in 1:length(ds_list)){
  
  
  for_hypgsa <- ds_list[[ds]]
  sig_column <- grep("Significant", colnames(for_hypgsa))
  colnames(for_hypgsa)[sig_column] <- "Significant"
  
  for_hypgsa <- for_hypgsa[,c("metal","ORF","Significant")]%>%
    group_by(ORF)%>%
    summarise(Signficant = ifelse(any(Significant,na.rm = T), T, F ))%>%
    ungroup()
  
  hyperGSA_res <- run_dbwise_hyperGSA(for_hypgsa,names(ds_list)[ds])
  hyperGSA_res_all_datasets <- rbind(hyperGSA_res_all_datasets, hyperGSA_res)
  
}

write.csv(hyperGSA_res_all_datasets, paste0(output_tables_dir,"/HyperGSAresults_alldatasets.csv"),row.names = F)

#####################################################
### Plot Sankey plot dataset to gene set enriched ###
#####################################################

enriched_gsnames <- unique(hyperGSA_res_all_datasets$Gset.Type)


for(egs in 1:length(enriched_gsnames)){
  sn <- plot_twocolumn_dataset2GSterm_Sankey(hyperGSA_res_all_datasets,enriched_gsnames[egs])
  htmlwidgets::saveWidget(as_widget(sn), paste0(plot_dir,"/HyperGSA/HyperGSA_dswise_Sankey",enriched_gsnames[egs],".html"))
  plotly::export(sn, file = paste0(plot_dir,"/HyperGSA/HyperGSA_dswise_Sankey",enriched_gsnames[egs],".pdf"))
}

########################################
### Helper functions for Sankey plot ###
########################################

pval2colour_OE <- function(x){
  
  col = colorRamp2(c(0.001,0.05), c("darkred","#FFC1C1"),
                   transparency = 0.7)
  
  return(col(x))
  
}

pval2colour_UE <- function(x){
  
  col = colorRamp2(c(0.001,0.05), c("darkblue","#B0E2FF"),
                   transparency = 0.7)
  
  return(col(x))
}


##################################################################
### Sankey Plot for element perturbation up and down regulated ###
##################################################################

plot_Sankey_DE <- function(df,
                           GeneSetType,
                           minGsetSize,
                           maxGsetSize){
  
  
  ## Make source to target df for Gset terms enriched in under expressed 
  
  UE_for_Sankey <- df %>%
    filter(Gset.Type==GeneSetType &
             pAdj_disdir_down < 0.05)%>%
    mutate(source=paste(Gset.Term.Enriched,"UE"),
           target=paste(Element,Perturbation,"pert"),
           value=1,
           pvAdjEnch = pAdj_disdir_down)%>%
    dplyr::select(source,target,value,pvAdjEnch)%>%
    unique()%>%
    # Add colours according to pvAdj of enrichment
    # mutate(link_colour = as.character(lapply(pvAdjEnch,pval2colour_UE)))
    # Add colour according to perturbation 
    mutate(link_colour = ifelse(target =="Ca Depletion pert",
                                grDevices::colorRampPalette(c("#BCD2EE","#1F78B4"))(9)[[1]],
                                ifelse(target =="Ca Excess pert",
                                       grDevices::colorRampPalette(c("#BCD2EE","#1F78B4"))(9)[[9]],
                                       ifelse(target =="Cu Depletion pert",
                                              grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[1]],
                                              ifelse(target =="Cu Excess pert",
                                                     grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[12]],
                                                     ifelse(target =="Fe Depletion pert",
                                                            grDevices::colorRampPalette(c("#F08080","#B0171F"))(13)[[1]],
                                                            ifelse(target =="Fe Excess pert",
                                                                   grDevices::colorRampPalette(c("#F08080","#B0171F"))(13)[[12]],
                                                                   ifelse(target =="K Depletion pert",
                                                                          grDevices::colorRampPalette(c("#CDC9C9","#5B5B5B"))(6)[[1]],
                                                                          ifelse(target =="K Excess pert",
                                                                                 grDevices::colorRampPalette(c("#CDC9C9","#5B5B5B"))(6)[[6]],
                                                                                 ifelse(target =="Mg Depletion pert",
                                                                                        grDevices::colorRampPalette(c("#E3A869","#8B4513"))(9)[[1]],
                                                                                        ifelse(target =="Mg Excess pert",
                                                                                               grDevices::colorRampPalette(c("#E3A869","#8B4513"))(9)[[8]],
                                                                                               ifelse(target =="Mn Depletion pert",
                                                                                                      grDevices::colorRampPalette(c("#FFE1FF","#6A3D9A"))(13)[[1]],
                                                                                                      ifelse(target =="Mn Excess pert",
                                                                                                             grDevices::colorRampPalette(c("#FFE1FF","#6A3D9A"))(13)[[13]],  
                                                                                                             ifelse(target =="Mo Depletion pert",
                                                                                                                    grDevices::colorRampPalette(c("#40E0D0","#008080"))(13)[[1]],
                                                                                                                    ifelse(target =="Mo Excess pert",
                                                                                                                           grDevices::colorRampPalette(c("#40E0D0","#008080"))(13)[[12]], 
                                                                                                                           ifelse(target =="Na Depletion pert",
                                                                                                                                  grDevices::colorRampPalette(c("#FFEC8B","#EE9A00"))(13)[[1]],
                                                                                                                                  ifelse(target =="Na Excess pert",
                                                                                                                                         grDevices::colorRampPalette(c("#FFEC8B","#EE9A00"))(13)[[13]], 
                                                                                                                                         ifelse(target =="Zn Depletion pert",
                                                                                                                                                grDevices::colorRampPalette(c("#CAE1FF","#000080"))(14)[[3]],
                                                                                                                                                ifelse(target =="Zn Excess pert",
                                                                                                                                                       grDevices::colorRampPalette(c("#CAE1FF","#000080"))(14)[[14]],
                                                                                                                                                       "gray"))))))))))))))))))
    )%>%
    mutate(link_colour =as.character(lapply(link_colour, hex2colWopacity)))
  ## Overexpressed Sankey Nodes
  
  OE_for_Sankey <- df %>%
    filter(Gset.Type==GeneSetType &
             pAdj_disdir_up < 0.05  )%>%
    mutate(source=paste(Element,Perturbation,"pert"),
           target=paste(Gset.Term.Enriched,"OE"),
           value=1,
           pvAdjEnch = pAdj_disdir_up)%>%
    dplyr::select(source,target,value,pvAdjEnch)%>%
    unique()%>%
    # Add colours according to pvAdj of enrichment
    # mutate(link_colour = as.character(lapply(pvAdjEnch,pval2colour_OE)))
    # Add colour according to perturbation 
    mutate(link_colour = ifelse(source =="Ca Depletion pert",
                                grDevices::colorRampPalette(c("#BCD2EE","#1F78B4"))(9)[[1]],
                                ifelse(source =="Ca Excess pert",
                                       grDevices::colorRampPalette(c("#BCD2EE","#4F78B4"))(9)[[9]],
                                       ifelse(source =="Cu Depletion pert",
                                              grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[1]],
                                              ifelse(source =="Cu Excess pert",
                                                     grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[12]],
                                                     ifelse(source =="Fe Depletion pert",
                                                            grDevices::colorRampPalette(c("#F08080","#B0171F"))(13)[[1]],
                                                            ifelse(source =="Fe Excess pert",
                                                                   grDevices::colorRampPalette(c("#F08080","#B0171F"))(13)[[12]],
                                                                   ifelse(source =="K Depletion pert",
                                                                          grDevices::colorRampPalette(c("#CDC9C9","#5B5B5B"))(6)[[1]],
                                                                          ifelse(source =="K Excess pert",
                                                                                 grDevices::colorRampPalette(c("#CDC9C9","#5B5B5B"))(6)[[6]],
                                                                                 ifelse(source =="Mg Depletion pert",
                                                                                        grDevices::colorRampPalette(c("#E3A869","#8B4513"))(9)[[1]],
                                                                                        ifelse(source =="Mg Excess pert",
                                                                                               grDevices::colorRampPalette(c("#E3A869","#8B4513"))(9)[[8]],
                                                                                               ifelse(source =="Mn Depletion pert",
                                                                                                      grDevices::colorRampPalette(c("#FFE1FF","#6A3D9A"))(13)[[1]],
                                                                                                      ifelse(source =="Mn Excess pert",
                                                                                                             grDevices::colorRampPalette(c("#FFE1FF","#6A3D9A"))(13)[[13]],  
                                                                                                             ifelse(source =="Mo Depletion pert",
                                                                                                                    grDevices::colorRampPalette(c("#40E0D0","#008080"))(13)[[1]],
                                                                                                                    ifelse(source =="Mo Excess pert",
                                                                                                                           grDevices::colorRampPalette(c("#40E0D0","#008080"))(13)[[12]], 
                                                                                                                           ifelse(source =="Na Depletion pert",
                                                                                                                                  grDevices::colorRampPalette(c("#FFEC8B","#EE9A00"))(13)[[1]],
                                                                                                                                  ifelse(source =="Na Excess pert",
                                                                                                                                         grDevices::colorRampPalette(c("#FFEC8B","#EE9A00"))(13)[[13]], 
                                                                                                                                         ifelse(source =="Zn Depletion pert",
                                                                                                                                                grDevices::colorRampPalette(c("#CAE1FF","#000080"))(14)[[3]],
                                                                                                                                                ifelse(source =="Zn Excess pert",
                                                                                                                                                       grDevices::colorRampPalette(c("#CAE1FF","#000080"))(14)[[14]],
                                                                                                                                                       "gray"))))))))))))))))))
    )%>%
    mutate(link_colour =as.character(lapply(link_colour, hex2colWopacity)))
  ## Combined 
  
  DE_for_Sankey <- rbind(UE_for_Sankey,OE_for_Sankey)
  
  DE_for_Sankey_nodes <- data.frame(name=c(as.character(DE_for_Sankey$source), 
                                           as.character(DE_for_Sankey$target)) %>% 
                                      unique())%>%
    # Add node positions 
    # x coord
    mutate(x_coord = ifelse(grepl("pert",name),0.5,
                            ifelse(grepl("UE",name),0.1,0.9)))%>%
    mutate(n=1)%>%
    # y coord
    group_by(x_coord)%>%
    mutate(y_cood = length(x_coord),
           counter = cumsum(n)
    )%>%
    ungroup()%>%
    mutate(y_cood= counter/y_cood)%>%
    dplyr::select(-n,-counter)
  
  DE_for_Sankey$IDsource=match(DE_for_Sankey$source, DE_for_Sankey_nodes$name)-1 
  DE_for_Sankey$IDtarget=match(DE_for_Sankey$target, DE_for_Sankey_nodes$name)-1
  
  pert_cat <- ifelse(grepl("UE",DE_for_Sankey_nodes$name),"darkblue",
                     ifelse(grepl("OE",DE_for_Sankey_nodes$name),"darkred",
                            ifelse(grepl("pert",DE_for_Sankey_nodes$name),
                                   ifelse(DE_for_Sankey_nodes$name =="Ca Depletion pert",
                                          grDevices::colorRampPalette(c("#BCD2EE","#1F78B4"))(9)[[1]],
                                          ifelse(DE_for_Sankey_nodes$name =="Ca Excess pert",
                                                 grDevices::colorRampPalette(c("#BCD2EE","#1F78B4"))(9)[[9]],
                                                 ifelse(DE_for_Sankey_nodes$name =="Cu Depletion pert",
                                                        grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[1]],
                                                        ifelse(DE_for_Sankey_nodes$name =="Cu Excess pert",
                                                               grDevices::colorRampPalette(c("#90EE90","#228B22"))(13)[[12]],
                                                               ifelse(DE_for_Sankey_nodes$name =="Fe Depletion pert",
                                                                      grDevices::colorRampPalette(c("#F08080","#B0171F"))(13)[[1]],
                                                                      ifelse(DE_for_Sankey_nodes$name =="Fe Excess pert",
                                                                             grDevices::colorRampPalette(c("#F08080","#B0171F"))(13)[[12]],
                                                                             ifelse(DE_for_Sankey_nodes$name =="K Depletion pert",
                                                                                    grDevices::colorRampPalette(c("#CDC9C9","#5B5B5B"))(6)[[1]],
                                                                                    ifelse(DE_for_Sankey_nodes$name =="K Excess pert",
                                                                                           grDevices::colorRampPalette(c("#CDC9C9","#5B5B5B"))(6)[[6]],
                                                                                           ifelse(DE_for_Sankey_nodes$name =="Mg Depletion pert",
                                                                                                  grDevices::colorRampPalette(c("#E3A869","#8B4513"))(9)[[1]],
                                                                                                  ifelse(DE_for_Sankey_nodes$name =="Mg Excess pert",
                                                                                                         grDevices::colorRampPalette(c("#E3A869","#8B4513"))(9)[[8]],
                                                                                                         ifelse(DE_for_Sankey_nodes$name =="Mn Depletion pert",
                                                                                                                grDevices::colorRampPalette(c("#FFE1FF","#6A3D9A"))(13)[[1]],
                                                                                                                ifelse(DE_for_Sankey_nodes$name =="Mn Excess pert",
                                                                                                                       grDevices::colorRampPalette(c("#FFE1FF","#6A3D9A"))(13)[[13]],  
                                                                                                                       ifelse(DE_for_Sankey_nodes$name =="Mo Depletion pert",
                                                                                                                              grDevices::colorRampPalette(c("#40E0D0","#008080"))(13)[[1]],
                                                                                                                              ifelse(DE_for_Sankey_nodes$name =="Mo Excess pert",
                                                                                                                                     grDevices::colorRampPalette(c("#40E0D0","#008080"))(13)[[12]], 
                                                                                                                                     ifelse(DE_for_Sankey_nodes$name =="Na Depletion pert",
                                                                                                                                            grDevices::colorRampPalette(c("#FFEC8B","#EE9A00"))(13)[[1]],
                                                                                                                                            ifelse(DE_for_Sankey_nodes$name =="Na Excess pert",
                                                                                                                                                   grDevices::colorRampPalette(c("#FFEC8B","#EE9A00"))(13)[[13]], 
                                                                                                                                                   ifelse(DE_for_Sankey_nodes$name =="Zn Depletion pert",
                                                                                                                                                          grDevices::colorRampPalette(c("#CAE1FF","#000080"))(14)[[3]],
                                                                                                                                                          ifelse(DE_for_Sankey_nodes$name =="Zn Excess pert",
                                                                                                                                                                 grDevices::colorRampPalette(c("#CAE1FF","#000080"))(14)[[14]],
                                                                                                                                                                 "gray")))))))))))))))))),"gray"))) 
  pert_cat <- as.character(lapply(pert_cat,hex2colWopacity))
  
  DE_for_Sankey_nodes<- DE_for_Sankey_nodes%>%
    mutate(name = gsub("pert","",name),
           name = gsub("UE","",name),
           name = gsub("OE","",name))
  
  
  
  fig <- plotly::plot_ly(
    type = "sankey",
    node = list(
      label = DE_for_Sankey_nodes$name,
      x=DE_for_Sankey_nodes$x_coord,
      y=DE_for_Sankey_nodes$y_cood,
      color = pert_cat,
      pad = 20,
      thickness = 20,
      line = list(
        color = "black",
        width = 0.5
      )
    ),
    
    link = list(
      source = DE_for_Sankey$IDsource,
      target = DE_for_Sankey$IDtarget,
      value =  DE_for_Sankey$value,
      color = DE_for_Sankey$link_colour
    )
  )
  fig <- fig %>% layout(
    title = paste(GeneSetType),
    font = list(
      size = 15,
      family="Times"
    )
  )
  return(fig)
}

