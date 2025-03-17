
library(diann)

library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(readxl)
library(proBatch)
library(DEP)

plot_DIANNstat_distrb <- function(df,plotname){
  
  df <- reshape2::melt(df,id.vars = "File.Name")
  pdf(paste0("QC stats ",plotname,".pdf"),width=20,height=30)
  print(
    ggplot(df,
           aes(x=value))+
      geom_histogram(bins=75,fill ="#191970")+
      facet_wrap("variable",scales="free",ncol=3)+
      theme_metallica()+
      theme(axis.text.x = element_text(angle=90))
  )
  dev.off()
}

##########################################################
### proteomics data processing and wrangling functions ###
##########################################################

## Function for getting output from diann_maxLFQ() using the column specified in pepquant_cname as precursor Intensity column
get_PGQ_df<-function(df){
  
  PGQ <- as.data.frame(diann_maxlfq(df, 
                                    sample.header = "File.Name",
                                    group.header="Protein.Ids",
                                    id.header = "Precursor.Id",
                                    quantity.header = "Precursor.Normalised"))
  
  PGQ$Protein.Ids=rownames(PGQ)
  
  PGQ <- reshape2::melt( PGQ,id.vars="Protein.Ids",
                         value.name="Protein.Quantity",
                         variable.name="File.Name" )
  
  return(PGQ)
}



### Function to set up data for proBatch
setup_proBatch_annos<-function(df,onlyAEdat=F,perele=F){ 
  
  
  ###############################################################
  ### Clean up and retain only essential columns for ProBatch ###
  ###############################################################
  
  PGcol <<- 'Protein.Ids'
  PGintcol <<- 'Protein.Quantity'
  fname_col <<- 'File.Name'
  PGintcol_AEcorr<<-"AEBatchCorr.PGquant"
  PGintcol_QNDC_AEcorr<<-"AEBatchCorr.PGquant.QNDC"
  PGintcol_qsm_AEcorr<<-"AEBatchCorr.PGquant.qsm"
  PGintcol_qsmBC_AEcorr<<-"AEBatchCorr.PGquant.qsmBC"
  
  ##########################################
  ### Make Experimental annotation table ###
  ##########################################
  
  exp_anno <- unique(df[,c("File.Name","Layout","TimePoint","LO_TP","Digestion_Batch",
                           "Element","Element.Concentration",
                           "BioSpecID","Date",
                           "Time","Sample.Type")])
  
  technical_factors <<- c('LO_TP','Digestion_Batch') 
  batch_col <<- 'LO_TP'
  
  biological_factors <<- c("LO_TP","Element","Sample.Type")
  
  biological_factors_AE <<- c("Sample.Type")
  technical_factors_AE <<- c('LO_TP','Digestion_Batch') 
  
  sample_anno <<- date_to_sample_order(exp_anno, 
                                       time_column = c('Date','Time'),
                                       new_time_column = 'DateTime', 
                                       dateTimeFormat = c('%Y-%m-%d', '%H:%M'), 
                                       new_order_col = 'order', instrument_col = NULL)%>%
    mutate(
      Element.Concentration = ifelse(is.na(Element.Concentration),0,Element.Concentration),
      BioSpecID = gsub("NA","",BioSpecID))
  
  
  #####################
  ### Colour Scheme ###
  #####################
  
  color_list_AE <<- sample_annotation_to_colors(sample_anno,
                                                factor_columns = c("Sample.Type","Digestion_Batch","LO_TP"), 
                                                numeric_columns = c('DateTime','order'),
                                                numeric_palette_type ="viridis")
  color_list_AE$Sample.Type[c(names(color_list_AE$Sample.Type))] <- c("#00008B","#800000","#EE6AA7")
  color_list_AE <<-color_list_AE
  color_list_AE$Digestion_Batch <-c("#EE9A00","#104E8B")
  
  if(!onlyAEdat){
    color_list <<- sample_annotation_to_colors(sample_anno,
                                               factor_columns = c("Element",'LO_TP',
                                                                  "Sample.Type","Digestion_Batch"), 
                                               numeric_columns = c('DateTime','order',
                                                                   "Element.Concentration"),
                                               numeric_palette_type ="viridis")
    color_list$Digestion_Batch <-c("#EE9A00","#104E8B")
    color_list$Element.Concentration <- c("#FFFFFF",viridis(10,begin=0.2))
    if(perele == F){
      color_list$Element[c(names(color_list$Element))] <- c("#8B8682","#FFD700",
                                                            "#B8860B","#B0171F","#AB82FF","#6E8B3D",
                                                            "#EE82EE","#808A87","#B8860B","#104E8B")
    }
    color_list$Sample.Type[c(names(color_list$Sample.Type))] <- c("#00008B","#800000"	,"#008B00")
    color_list<<-color_list
    
  }
}


### Function for Batch correction at precursor level ###

AEmed_correction_Precursor <- function(df){
  
  
  ### Correct for Digestion Batches ###
  
  AEmeds_DigBatch <- df %>%
    filter(Sample.Type == "Control")%>%
    ## Store digestion batch medians
    group_by(Precursor.Id,Digestion_Batch)%>%
    summarize(AE.median.IntraDigBatch = median (Precursor.Normalised,na.rm = T))%>%
    ungroup()%>%
    group_by(Precursor.Id)%>%
    mutate(AE.median.InterDigBatch = median (AE.median.IntraDigBatch,na.rm = T))%>%
    ungroup()%>%
    
    dplyr::select(Precursor.Id,Digestion_Batch,
                  AE.median.IntraDigBatch,AE.median.InterDigBatch
    )%>%
    unique()
  
  DigBatch_corr_df <- merge(df,AEmeds_DigBatch,  by = c("Digestion_Batch", "Precursor.Id"))%>%
    # Normalize values and bring them back to log2 scale
    mutate( DigBatch.corr.Precursor.Quantity = (Precursor.Normalised/AE.median.IntraDigBatch)*AE.median.InterDigBatch,
            File.Name = as.character(File.Name),
            DigBatch.corr.Precursor.Quantity = ifelse(is.infinite(DigBatch.corr.Precursor.Quantity),
                                                      NA,DigBatch.corr.Precursor.Quantity )
    )
  
  ## Store 96-well plate medians
  AEmeds_Batch <- DigBatch_corr_df%>%
    filter(Sample.Type == "Control")%>%
    group_by(Precursor.Id, LO_TP)%>%
    summarize(AE.median.IntraBatch = median (DigBatch.corr.Precursor.Quantity,na.rm = T))%>%
    ungroup()%>%
    group_by(Precursor.Id)%>%
    mutate(AE.median.InterBatch = median (AE.median.IntraBatch,na.rm = T))%>%
    ungroup()%>%
    dplyr::select(Precursor.Id,LO_TP,
                  AE.median.IntraBatch,AE.median.InterBatch
    )%>%
    unique()
  
  
  corr_df <- merge(DigBatch_corr_df,AEmeds_Batch,  by = c("LO_TP", "Precursor.Id"))%>%
    # Normalize values and bring them back to log2 scale
    mutate(AEmed.corr.Precursor.Quantity = (DigBatch.corr.Precursor.Quantity/AE.median.IntraBatch)*AE.median.InterBatch,
           File.Name = as.character(File.Name))%>%
    mutate(Precursor.Normalised = ifelse(is.infinite(AEmed.corr.Precursor.Quantity),
                                         NA,AEmed.corr.Precursor.Quantity ))%>%
    dplyr::select(-c(AE.median.IntraDigBatch,AE.median.InterDigBatch,
                     DigBatch.corr.Precursor.Quantity,AE.median.IntraBatch,AE.median.InterBatch,AEmed.corr.Precursor.Quantity))
  return(corr_df)
  
}


## Function to replace NAs absent in all replicates of a sample with the minimum detected value

fill_in_NAs <- function(df){
  
  df = reshape2::melt(df,value.name = "Protein.Quantity")
  colnames(df) = c("Protein.Ids","File.Name","Protein.Quantity")
  
  minVals = df%>%
    group_by(Protein.Ids)%>%
    summarize(Min.Protein.Quantity = min(Protein.Quantity,na.rm = T))%>%
    unique()%>%
    ungroup()%>%
    mutate(Min.Protein.Quantity = ifelse(is.infinite(Min.Protein.Quantity),NA,Min.Protein.Quantity))%>%
    data.frame()
  rownames(minVals) = minVals$Protein.Ids
  
  AEdat<- merge(df,unique(PG_quants[,c("File.Name","BioSpecID","Sample.Type")]),by="File.Name")%>%
    filter(Sample.Type =="Control")%>%
    group_by(Protein.Ids)%>%
    summarize(Median.Protein.Quant = median(Protein.Quantity,na.rm=T))%>%
    as.data.frame()
  
  rownames(AEdat) = AEdat$Protein.Ids
  
  df_fNA = merge(df,unique(PG_quants[,c("File.Name","BioSpecID","Element","Element.Concentration")]),by="File.Name")%>%
    mutate(Pert= ifelse(Element.Concentration <1 ,"Depletion",
                        ifelse(Element.Concentration >1,"Excess","Control")))%>%
    group_by(Element,Pert,Protein.Ids)%>%
    mutate(Missing.in.Perturbation = all(is.na(Protein.Quantity)))%>%
    ungroup()%>%
    group_by(BioSpecID,Protein.Ids)%>%
    mutate(mPQ = median (Protein.Quantity,na.rm = T))%>%
    mutate(Missing.In.All.Reps = all(is.na(Protein.Quantity)))%>%
    mutate(Protein.Quantity.Imp = ifelse( is.na(Protein.Quantity),
                                          ifelse(Missing.in.Perturbation,minVals[Protein.Ids,"Min.Protein.Quantity"],
                                                 ifelse(Missing.In.All.Reps,AEdat[Protein.Ids,"Median.Protein.Quant"] , mPQ)
                                          ), 
                                          Protein.Quantity 
    )
    )%>%
    reshape2::dcast(Protein.Ids~File.Name,value.var = "Protein.Quantity.Imp")%>%
    na.omit()
  
  rwn<-df_fNA$Protein.Ids
  df_fNA<-df_fNA[,-1]
  df_fNA<-do.call(cbind, df_fNA)
  rownames(df_fNA)<-rwn
  
  return(df_fNA)
}



## Function for getting the intensity values stored in assays(SumExpObject)[[1]] into a df with filenames as colnames and rownames as protein group names

SumExpObj_to_df<-function(SEObj){
  
  df<-SummarizedExperiment::assays(SEObj)[[1]]
  colD<-data.frame(SummarizedExperiment::colData(SEObj))
  colnames(df)<-colD[colnames(df),"label"]
  
  return(df)
  
}

#################################################
### proteomics processing visualisation tools ###
#################################################

### Plots to check intensity profiles ###


plot_protint_dist <- function(df,plotname){
  
  df<- reshape2::melt(df)
  
  colnames(df)<- c("Protein Id", "File Name","Log2 Protein Quantity")
  df <- df%>%
    separate(`File Name`,into = c(NA,"Layout"),sep="_LO_",remove=F)%>%
    separate(Layout, into = c("Layout",NA),sep="_")
  print(
    ggplot(df,
           aes(x=`Log2 Protein Quantity`,
               group = `File Name`,
               colour = Layout))+
      geom_density(alpha=0.2)+
      labs(title = plotname)+
      theme_metallica()
  )
}


### Function for making all QC plots from probatch matrix

make_proBatchQCplots <- function(datMatrix,sel_anno_hclust,onlyAEdat=F,samLimBCAEdat=F,
                                 dsname="",AE_or_samp="AllEle",
                                 sanno=sample_anno){
  ## Filer sample annotation file by filenames in datamatrix
  
  sanno <- filter(sanno,
                  File.Name %in% colnames(datMatrix))
  if(onlyAEdat){
    cl <- color_list_AE
    htmpc <- seq(0.4, 1, by = 0.6/101)
  }else{
    cl <- color_list
    htmpc <- seq(0, 1, by = 1/101)
  }
  
  
  datMatrix_hc <- datMatrix
  datMatrix_hc[is.na(datMatrix_hc)] <- 0
  datMatrix_hc[is.infinite(datMatrix_hc)] <- 0
  
  pdf(paste0(plot_dir,"/",AE_or_samp,"_",dsname,"_hclust.pdf"), width=15,height=6)
  plot_hierarchical_clustering(datMatrix_hc,
                               sample_annotation = sanno,
                               color_list = cl,
                               sample_id_col = fname_col,
                               factors_to_plot = sel_anno_hclust,
                               distance = 'euclidean', 
                               agglomeration = 'complete',
                               label_samples = FALSE,
                               plot_title=dsname) 
  dev.off()
  
  datMatrix_PCA <- datMatrix
  datMatrix_PCA[is.infinite(datMatrix_PCA)] <- NA
  
  datMatrix_PCA<-datMatrix_PCA[complete.cases(datMatrix_PCA), ] 
  
  
  
  PCA_LOTP <- plot_PCA(datMatrix_PCA, sanno, 
                     color_by = 'LO_TP',
                     plot_title = 'Batch',
                     sample_id_col = fname_col,
                     color_scheme = cl[['LO_TP']] )
  
  
  PCA_Element <- plot_PCA(datMatrix_PCA, sanno, 
                        color_by = 'Element',
                        plot_title = 'Element',
                        sample_id_col = fname_col,
                        color_scheme = cl[['Element']] )
  
  
  PCA_SampleType <- plot_PCA(datMatrix_PCA, sanno, 
                           color_by = 'Sample.Type',
                           plot_title = 'Sample Type',
                           sample_id_col = fname_col,
                           color_scheme = cl[['Sample.Type']] )
  PCA_DigestionBatch <- plot_PCA(datMatrix_PCA, sanno, 
                               color_by = 'Digestion_Batch',
                               plot_title = 'Digestion_Batch',
                               sample_id_col = fname_col,
                               color_scheme = cl[['Digestion_Batch']] )
  
  PCA_order <- plot_PCA(datMatrix_PCA, sanno, 
                      color_by = 'order',
                      plot_title = 'Measurement order',
                      sample_id_col = fname_col )
  
  
  if(onlyAEdat){
    pdf(paste0(plot_dir,"/",AE_or_samp,"_",dsname,"_PCA_and_PVCA.pdf"), width=15,height=6)
    
    print((PCA_LOTP+PCA_SampleType)/(PCA_DigestionBatch+plot_spacer()))
    print(
      plot_PVCA(datMatrix_PCA,
                sample_id_col = fname_col,
                feature_id_col = PGcol,
                sample_annotation = sanno,
                biological_factors = biological_factors_AE,
                technical_factors = technical_factors_AE[which(technical_factors_AE %in% sel_annos)],
                plot_title = dsname
      )
      
    )
    dev.off()
  }else{
    pdf(paste0(plot_dir,"/",AE_or_samp,"_",dsname,"_PCA_and_PVCA.pdf"), width=12,height=8)
    
    print((PCA_LOTP+PCA_Element+PCA_SampleType)/(PCA_DigestionBatch+PCA_order))
    
    if(length(technical_factors[which(technical_factors %in% sel_annos)]) > 0){
      print(
        plot_PVCA(datMatrix_PCA,
                  sample_id_col = fname_col,
                  feature_id_col = PGcol,
                  sample_annotation = sanno,
                  biological_factors = biological_factors,
                  technical_factors = technical_factors[which(technical_factors %in% sel_annos)],
                  plot_title = dsname
        )
      )
    }
    dev.off()
    
  }
}

DE_analysis_and_plots_for_AllEleQC <- function(df,dsname,expdes,onlyAEdat=T,pcacompatible=F){
  
  if(onlyAEdat)  {
    bioCntrls <- expdes%>%
      filter(!grepl("Extraction Control",condition))
    
    
    df_fDEa<-data.frame(df)%>%
      mutate(ID = rownames(df),
             name = rownames(df))
    df_fDEa <- df_fDEa[,c(bioCntrls$label,"name","ID")]
    
    se_fDE<- make_se(df_fDEa,seq(1,ncol(df_fDEa)-2),bioCntrls)
    
  }else{
    samples <- expdes%>%
      filter(!grepl("AllEleprepQC",condition))
    
    df_fDEa<-data.frame(df)%>%
      mutate(ID = rownames(df),
             name = rownames(df))
    df_fDEa <- df_fDEa[,c(samples$label,"name","ID")]
    
    se_fDE<- make_se(df_fDEa,seq(1,ncol(df_fDEa)-2),samples)
    
  }
  if(onlyAEdat){
    deres <- test_diff(se_fDE, type = "all")
  }else{
    deres <- test_diff(se_fDE, type = "control",control = "AllEle", design_formula = formula(~0+condition+replicate))
  }
  
  deres <- add_rejections(deres, alpha = 0.05, lfc = log2(1.5))
  
  # Get results in a df format 
  
  deres_df<-get_df_long(deres)
  
  # Find out if there are any significant proteins
  anySig <- sum(colSums(as.matrix(deres_df[grepl("_significant",colnames(deres_df))]))) > 0
  
  
  all_comparisons <- colnames(SummarizedExperiment::elementMetadata(deres))[grepl("_vs_",colnames(elementMetadata(deres)))]
  all_comparisons <- gsub("_CI.L","",all_comparisons)
  all_comparisons <- gsub("_CI.R","",all_comparisons)
  all_comparisons <- gsub("_p.adj","",all_comparisons)
  all_comparisons <- gsub("_p.val","",all_comparisons)
  all_comparisons <- gsub("_significant","",all_comparisons)
  all_comparisons <- unique(gsub("_diff","",all_comparisons))
  
  vol_plots<-list()
  for(ac in 1:length(all_comparisons)){
    vol_plots[[ac]] <-  DEP::plot_volcano(deres,contrast = all_comparisons[ac],
                                          adjusted = T,
                                          plot=T)
  }
  nCol <- floor(sqrt(length(all_comparisons)))
  
  
  if(pcacompatible){
    
    pca <- DEP::plot_pca(deres, x = 1, y = 2, n = 500, point_size = 4)+ ## PCA can't handle NAs yet
      labs(title=paste0("PCA top 500 variables on ",dsname,"data"))
    
    pca <- pca+theme(legend.position = "none")
  }else{
    pca <- plot_spacer()
  }
  
  pvhist <- DEP::plot_p_hist(deres,
                             adjusted = T,
                             wrap =F)
  
  
  if(anySig){
    
    bp1<-DEP::plot_cond_freq(deres)+scale_fill_viridis_d(end=0.95)+
      labs(title="Frequency of significant\n conditions per protein")
  }else{
    bp1<-plot_spacer()
  }
  if(onlyAEdat){
    pdf(paste0("AllEle DE results ",dsname,".pdf"),width = 40,height=40 )
  }else{
    pdf(paste0("Sample DE results ",dsname,".pdf"),width = 40,height=40 )
  }
  #if(onlyAEdat){
  
  
  # DEP::plot_heatmap(deres, type = "centered",
  #                  kmeans = F, 
  #                 col_limit = 2, 
  #                indicate = c("condition", "replicate"),
  #               show_row_names = T,
  #              col_font_size = 8,
  #             row_font_size = 10 )
  # }else{
  #  DEP::plot_heatmap(deres, type = "centered",
  #                     kmeans = F, 
  #                      col_limit = 2, 
  #                      indicate = c("condition"),
  #                     show_row_names = T,
  #                    col_font_size = 8,
  #                   row_font_size = 8 )
  #}
  
  
  
  do.call("grid.arrange", c(vol_plots, ncol=nCol))
  
  
  print((pca+pvhist)/bp1)
  
  dev.off()
  
  return(deres_df)
}










  