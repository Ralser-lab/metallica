
# Load required libraries
library(dplyr)
library(plotly)
library(doParallel)
library(piano)
library(visNetwork)


get_plot_HyperGSA<-function(forHyperGSA,EnrichPV_thresh=0.05,
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

hex2colWopacity <- function(hx){
  
  c <- col2rgb(hx)
  c <- rgb(c[1],c[2],c[3],max=255,alpha=125)
  
  return(c)
}


###########################################################################################
### Sankey Plot for results of hyperGSA element wise ### 2 columns only ## no direction ###
###########################################################################################

make_twocolumn_Sankey_plot <- function(df, gset_type) {
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
                      color = sapply(df_type$Element, function(x) paste0(colkey_Ele[[x]], "80"))) # assign color from colkey_Ele with 0.5 opacity
  
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
        size = 10
      )
    )
  
  return(p)
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
      family="Arial"
    )
  )
  return(fig)
}

