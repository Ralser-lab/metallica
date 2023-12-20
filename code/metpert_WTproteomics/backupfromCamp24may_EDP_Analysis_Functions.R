options(bitmapType='cairo')
 #install.packages('ncdf4', configure.args = c("--with-nc-config=/camp/apps/eb/dev/software/netCDF/4.7.1-gompi-2019b/bin/nc-config"))
require(org.Sc.sgd.db)
require(proBatch)
require(dplyr)
require(tibble)
require(tidyverse)
require(ggplot2)
require(raster)
require(viridis)
require(patchwork)
require(readxl)
require(reshape2)
require(diann)
require(quantro)
require(zoo)
require(rlang)
require(DEP)
require(qvalue)
require(gridExtra)
require(SummarizedExperiment)
require(matrixStats)
require(ComplexHeatmap)
require(dendextend)
require(circlize)
require(piano)
require(visNetwork)
#require(clusterProfiler)
require(enrichplot)
require(ggupset)
#require(ipath)
require(IDPmisc)
require(foreach)
require(doParallel)
require(data.table)
require(RColorBrewer)

require(edgeR)
require(splines)
require(UpSetR)
require(venneuler)
require(ggrepel)
require(venn)
require(ggpolypath)
require(Cairo)

#library(pathview)
require(plotly)
require(impute)
require(gprofiler2)

#####################################
### Set Paths used in all scripts ###
#####################################

# proteomics_dir path from camp onDemand

dbdir = '/camp/lab/ralserm/working/Simran Aulakh/Ionome Project/Element Depletion/EleDepChel January 2020/Databases'

datdir = paste0(proteomics_dir,"/DIANNout")
procdatdir = paste0(proteomics_dir,"/ProcData")
Int_datdir =paste0(paste0(proteomics_dir,"/IntermediateDFs"))

# for storing results

FinalOutCSVdir <-paste0(proteomics_dir,"/Final Output Stat Imp")

# for plots

plotdir=paste0(proteomics_dir,"/Plots")
QCplotdir=paste0(plotdir,"/QC_plots")

dir.create(procdatdir)
dir.create(Int_datdir)
dir.create(plotdir)
dir.create(QCplotdir)
dir.create(FinalOutCSVdir)

## Uniprot--SGD orfname -- gene id converter

GenProt_SGD <- read.csv(paste0(dbdir,"/uniprot_allScerevisiae_20230208.tsv"), sep="\t",stringsAsFactors = F)[c("Entry","Gene.Names..ordered.locus.","Gene.Names..primary.","Annotation")]%>%
  unique()

colnames(GenProt_SGD)<- c("Uniprot.ID","ORF","Gene.Name","Uniprot.Annotation.Score")


## Convert multiple ORF names into multiple rows
ORFs2bind<- data.frame()

for(i in 1:nrow(GenProt_SGD)){
  
  orfs <- trimws(unlist(strsplit(GenProt_SGD[i,"ORF"],";")))
  
  if(length(orfs)>1){
    
    for(j in 2:length(orfs)){
      
      ORFs2bind <- rbind(ORFs2bind,
                         cbind(Uniprot.ID =GenProt_SGD[i,"Uniprot.ID"],
                               ORF = orfs[j],
                               Gene.Name = GenProt_SGD[i,"Gene.Name"],
                               Uniprot.Annotation.Score = GenProt_SGD[i,"Uniprot.Annotation.Score"]
                         )
      )
    }
    GenProt_SGD[i,"ORF"] <- orfs[1]
    
  }
}

GenProt_SGD<-rbind(GenProt_SGD,ORFs2bind)

## Convert multiple gene names into multiple rows
genes2bind<- data.frame()

for(i in 1:nrow(GenProt_SGD)){
  
  genes <- unlist(strsplit(GenProt_SGD[i,"Gene.Name"],";"))
  
  if(length(genes)>1){
    
    for(j in 2:length(genes)){
      
      genes2bind <- rbind(genes2bind,
                          cbind(Uniprot.ID =GenProt_SGD[i,"Uniprot.ID"],
                                ORF =  GenProt_SGD[i,"ORF"],
                                Gene.Name = genes[j],
                                Uniprot.Annotation.Score = GenProt_SGD[i,"Uniprot.Annotation.Score"]
                          )
      )
    }
    GenProt_SGD[i,"Gene.Name"] <- genes[1]
    
  }
}
GenProt_SGD<-rbind(GenProt_SGD,genes2bind)

GenProt_SGD<-GenProt_SGD%>%
  mutate(ORF = ifelse(ORF =="", NA, ORF),
         Gene.Name = trimws(ifelse(Gene.Name =="",ORF,Gene.Name)))


#################
#################
### Functions ###
#################
#################

#################################################
### Functions for Quality control of raw data ###
#################################################


plot_DIANNstat_distrb <- function(df,plotname){
  
  df <- melt(df,id.vars = "File.Name")
  pdf(paste0("QC stats ",plotname,".pdf"),width=20,height=30)
  print(
    ggplot(df,
           aes(x=value))+
      geom_histogram(bins=75,fill ="#191970")+
      facet_wrap("variable",scales="free",ncol=3)+
      theme_SKA()+
      theme(axis.text.x = element_text(angle=90))
  )
  dev.off()
}


### Fix function in DEP to annotate heatmap legends with more than 12 annotation columns ###

# Internal function to get ComplexHeatmap::HeatmapAnnotation object
get_annotation_fixedfunc <- function(dep, indicate) {
  assertthat::assert_that(
    inherits(dep, "SummarizedExperiment"),
    is.character(indicate))
  
  # Check indicate columns
  col_data <- colData(dep) %>%
    as.data.frame()
  columns <- colnames(col_data)
  if(all(!indicate %in% columns)) {
    stop("'",
         paste0(indicate, collapse = "' and/or '"),
         "' column(s) is/are not present in ",
         deparse(substitute(dep)),
         ".\nValid columns are: '",
         paste(columns, collapse = "', '"),
         "'.",
         call. = FALSE)
  }
  if(any(!indicate %in% columns)) {
    indicate <- indicate[indicate %in% columns]
    warning("Only used the following indicate column(s): '",
            paste0(indicate, collapse = "', '"),
            "'")
  }
  
  # Get annotation
  anno <- select(col_data, indicate)
  
  # Annotation color
  names <- colnames(anno)
  anno_col <- vector(mode="list", length=length(names))
  names(anno_col) <- names
  for(i in names) {
    var = anno[[i]] %>% unique() %>% sort()
    if(length(var) == 1)
      cols <- c("black")
    if(length(var) == 2)
      cols <- c("orangered", "cornflowerblue")
    if(length(var) < 7 & length(var) > 2)
      cols <- RColorBrewer::brewer.pal(length(var), "Pastel1")
    if(length(var) > 7)
      cols <- RColorBrewer::brewer.pal(length(var), "Set3")
    if(length(var) > 12)
      cols <- viridis(length(var))
    names(cols) <- var
    anno_col[[i]] <-  cols
  }
  
  # HeatmapAnnotation object
  HeatmapAnnotation(df = anno,
                    col = anno_col,
                    show_annotation_name = TRUE)
}

environment(get_annotation_fixedfunc) <- asNamespace('DEP')
assignInNamespace("get_annotation", get_annotation_fixedfunc, ns = "DEP")



#################
### Functions ###
#################


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





## Function for getting the intensity values stored in assays(SumExpObject)[[1]] into a df with filenames as colnames and rownames as protein group names

SumExpObj_to_df<-function(SEObj){
  
  df<-SummarizedExperiment::assays(SEObj)[[1]]
  colD<-data.frame(SummarizedExperiment::colData(SEObj))
  colnames(df)<-colD[colnames(df),"label"]
  
  return(df)
  
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
  
  pdf(paste(AE_or_samp,dsname,"hclust.pdf"), width=20,height=10)
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
  
  
  
  PCA_LOTP<-plot_PCA(datMatrix_PCA, sanno, 
                     color_by = 'LO_TP',
                     plot_title = 'Batch',
                     sample_id_col = fname_col,
                     color_scheme = cl[['LO_TP']] )
  
  
  PCA_Element<-plot_PCA(datMatrix_PCA, sanno, 
                        color_by = 'Element',
                        plot_title = 'Element',
                        sample_id_col = fname_col,
                        color_scheme = cl[['Element']] )
  
  
  PCA_SampleType<-plot_PCA(datMatrix_PCA, sanno, 
                           color_by = 'Sample.Type',
                           plot_title = 'Sample Type',
                           sample_id_col = fname_col,
                           color_scheme = cl[['Sample.Type']] )
  PCA_DigestionBatch<-plot_PCA(datMatrix_PCA, sanno, 
                               color_by = 'Digestion_Batch',
                               plot_title = 'Digestion_Batch',
                               sample_id_col = fname_col,
                               color_scheme = cl[['Digestion_Batch']] )
  
  PCA_order<-plot_PCA(datMatrix_PCA, sanno, 
                      color_by = 'order',
                      plot_title = 'Measurement order',
                      sample_id_col = fname_col )
  
  
  if(onlyAEdat){
    pdf(paste(AE_or_samp,dsname,"PCA and PVCA.pdf"), width=8,height=6)
    
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
    pdf(paste(AE_or_samp,dsname,"PCA and PVCA.pdf"), width=12,height=8)
    
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
  
  # print(plot_sample_corr_heatmap(datMatrix,
  #                                 sample_annotation = sample_anno,
  #                                 color_list = cl,
  #                                 sample_id_col = fname_col,
  #                                 show_rownames = F,
  #                                 show_colnames = F,
  #                                 cluster_rows = T,
  #                                 cluster_cols = T,
  #                                 breaks = htmpc,
  #                                 samples_to_plot = colnames(datMatrix),
  #                                 factors_to_plot= sel_annos,
  #                                 annotation_legend=F,
  #                                 plot_title = paste('Correlation btw samples',
  #                                                    dsname))
  #  )
  
  
}

## Function to return result of ttest between groups specified by a group column using protein quantities specified by a protein quantity column 

get_TTest_Res <- function (df,PQcol,gp_fac,AEdat = F){
  
  ## Function to return BH corrected pvalues from ttest results between groups -- protein wise
  
  gps_totest= unique(as.character(df[,gp_fac]))
  
  if(!AEdat){ ## If dealing with sample data -- take out All Ele conditions
    
    gps_totest = gps_totest[-which(gps_totest %in% c("AllEle","AllEleprepQC"))]
  }
  
  PGs<- unique(df$Genes)
  
  ttres_df<-vector()
  
  for(gp1 in 1:length(gps_totest) ){ 
    
    for(gp2 in gp1:length(gps_totest) ) {
      
      if(gps_totest[gp1] != gps_totest[gp2]) {
        
        gps_ttres <- vector()
        
        for(pg in 1:length(PGs)){
          
          if(AEdat){
            
            PQ_group1 <- na.omit(filter(df,!!sym(gp_fac) == gps_totest[gp1] 
                                        & Genes == PGs[pg])[,PQcol])
            
          }else{
            
            ## If dealing with sample data, compare all conditions against the 
            # all element controls, not amongst each other
            PQ_group1 <- na.omit(filter(df,!!sym(gp_fac) == "AllEle" 
                                        & Genes == PGs[pg])[,PQcol])
            
          }
          PQ_group2 <- na.omit(filter(df,!!sym(gp_fac) == gps_totest[gp2] 
                                      & Genes == PGs[pg])[,PQcol])
          
          if(length(PQ_group1) > 2 & length(PQ_group2) > 2 & var(PQ_group1) !=0 & var(PQ_group2)!=0 ){
            
            ttres = t.test(PQ_group1,PQ_group2,var.equal = F)
            
            gps_ttres = rbind(gps_ttres,
                              cbind(paste(gps_totest[gp1],"vs",gps_totest[gp2],sep="_"),
                                    PGs[pg],
                                    ttres$p.value,
                                    ttres$estimate[[1]]/ttres[[2]] ))
          }
        }
        gps_ttres = cbind(gps_ttres,p.adjust(gps_ttres[,3],method = "BH"))
        ttres_df <- rbind(ttres_df,gps_ttres)
      }
    }
  }
  ttres_df<-data.frame(ttres_df,stringsAsFactors = F)
  colnames(ttres_df)<-c("Group for DE","Genes","Pvalue","FC_gp1_by_gp2","Adj.pValue")
  for(k in 3:5){
    ttres_df[,k] <- as.numeric(ttres_df[,k])
  }
  return(ttres_df)
}


## function to convert protein group to Gene Name -- return uniprot for proteins without gene names -- and multiple gene names for protein groups with more than 1 uniprot id

convert_Uniprot2GeneName <- function(unp_id){
  
  uids <- trimws(unlist(strsplit(unp_id,";")))
  
  allgns <- ""
  
  for( i in 1:length(uids)){
    
    gn <- GenProt_SGD%>%
      filter(Uniprot.ID == uids[i])%>%
      pull(Gene.Name)%>%
      unique()
    
    if(length(gn) <= 0 ){
      
      gns <-unp_id
      
    }else if(length(gn) > 0 ){  
      gns <- ""
      
      for(j in 1:length(gn)){ 
        gns <- paste(gns,gn[j],sep=" ") }
      gns <- trimws(gns)
      allgns<-paste(allgns,gns,sep=" ")
      allgns <- trimws(allgns)}
    
    else{
      allgns <- unp_id
    }
  }
  return(allgns)
}

## function to convert protein group to Gene Name -- return uniprot for proteins without gene names -- and multiple gene names for protein groups with more than 1 uniprot id

convert_Uniprot2SingleGeneName <- function(unp_id){
  
  uid <- trimws(unlist(strsplit(unp_id,";")))[[1]]
  
  gn <- GenProt_SGD%>%
    filter(Uniprot.ID == uid)%>%
    pull(Gene.Name)%>%
    unique()
  
  if(length(gn) > 1){
    gn <- trimws(gn[[1]])
  }else if( length(gn)==0 ){
    gn <- unp_id
  }else if(is.na(gn)){
    gn <- GenProt_SGD%>%
      filter(Uniprot.ID == uid)%>%
      pull(ORF)%>%
      unique()
    gn <- trimws(gn[[1]])
  }
  return(trimws(gn))
}

convert_ORF2GeneName <- function(ORFname){
  
  gn <- unique(filter(GenProt_SGD,ORF == ORFname)[,"Gene.Name"])
  
  if(length(gn)>1){
  }
  else if( length(gn)==0 ){
    gn <- ORFname
  }else if(is.na(gn)){
    gn <- ORFname
  }
  return(trimws(gn))
}

convert_ORF2SingleGeneName<-function(ORFname){
  
  gn <- unique(filter(GenProt_SGD,ORF == ORFname)[,"Gene.Name"])
  
  if(length(gn) > 1){
    gn <- trimws(gn[[1]])
  }else if( length(gn)==0 ){
    gn <- ORFname
  }else if(is.na(gn)){
    gn <- ORFname
  }
  return(trimws(gn))
}

convert_ORF2Uniprot <- function(ORFname){
  
  unp <- unique(filter(GenProt_SGD,ORF == ORFname)[,"Uniprot.ID"])
  
  if( length(unp)==0 ){
    unp <- ORFname
  }else if(is.na(unp)){
    unp <- ORFname
  }
  return(trimws(unp))
}

convert_GeneName2ORF <- function(genename){
  
  orf <- GenProt_SGD%>%
    filter(Gene.Name == genename)%>%
    dplyr::select(Gene.Name,ORF)%>%
    unique()
  
  orf <- unname(trimws(orf[["ORF"]]))
  
  if(length(orf) <= 0 ){
    
    orf <- genename
  }
  return(orf)
}


convert_GeneName2SingleORF <-function(genename){
  
  orf <- GenProt_SGD%>%
    filter(Gene.Name == genename)%>%
    dplyr::select(Gene.Name,ORF)%>%
    unique()
  
  orf <- unname(trimws(orf[["ORF"]]))
  
  if(length(orf) <= 0 ){
    
    orf <- genename
  }else if (length(orf) > 1){
    orf <- orf[[1]]
  }
  return(orf)
  
}


convert_GeneName2UNIPROT <- function(genename){
  
  unp <- GenProt_SGD%>%
    filter(Gene.Name == genename)%>%
    dplyr::select(Gene.Name,Uniprot.ID)%>%
    unique()
  
  
  unp <-unname(trimws(unp[["Uniprot.ID"]]))
  
  if(length(unp) <= 0 ){
    
    unp <- genename
  }
  return(unp)
}

# Function to normalize a df with AE meds and return AE med norm df

AEmed_correction <- function(df){
  
  # Convert incoming dataframe to 2^ as they are at log2 scale
  
  df <- 2^df
  
  df <- matrix_to_long(df,
                       sample_annotation = sample_anno,
                       feature_id_col = PGcol,
                       measure_col = "Protein.Quantity",
                       sample_id_col = fname_col)  
  
  df <- merge(unique(PG_quants[,c("File.Name","Protein.Ids","LO_TP",
                                  "Element","Digestion_Batch",
                                  "Sample.Type")]),
              df,by=c("File.Name","Protein.Ids"))
  
  ### Correct for Digestion Batches ###
  
  AEmeds_DigBatch <- df %>%
    filter(Sample.Type == "Control")%>%
    ## Store digestion batch medians
    group_by(Protein.Ids,Digestion_Batch)%>%
    summarize(AE.median.IntraDigBatch = median (Protein.Quantity,na.rm = T))%>%
    ungroup()%>%
    group_by(Protein.Ids)%>%
    mutate(AE.median.InterDigBatch = median (AE.median.IntraDigBatch,na.rm = T))%>%
    ungroup()%>%
    
    dplyr::select(Protein.Ids,Digestion_Batch,
                  AE.median.IntraDigBatch,AE.median.InterDigBatch
    )%>%
    ungroup()%>%
    unique()
  
  DigBatch_corr_df <- merge(df,AEmeds_DigBatch,  by = c("Digestion_Batch", "Protein.Ids"))%>%
    # Normalize values and bring them back to log2 scale
    mutate( DigBatch.corr.Protein.Quantity = (Protein.Quantity/AE.median.IntraDigBatch)*AE.median.InterDigBatch,
            File.Name = as.character(File.Name),
            DigBatch.corr.Protein.Quantity = ifelse(is.infinite(DigBatch.corr.Protein.Quantity),
                                                    NA,DigBatch.corr.Protein.Quantity )
    )
  
  ## Store 96-well plate medians
  AEmeds_Batch <- DigBatch_corr_df%>%
    filter(Sample.Type == "Control")%>%
    group_by(Protein.Ids, LO_TP)%>%
    summarize(AE.median.IntraBatch = median (DigBatch.corr.Protein.Quantity,na.rm = T))%>%
    ungroup()%>%
    group_by(Protein.Ids)%>%
    mutate(AE.median.InterBatch = median (AE.median.IntraBatch,na.rm = T))%>%
    ungroup()%>%
    dplyr::select(Protein.Ids,LO_TP,
                  AE.median.IntraBatch,AE.median.InterBatch
    )%>%
    ungroup()%>%
    unique()
  
  
  corr_df <- merge(DigBatch_corr_df,AEmeds_Batch,  by = c("LO_TP", "Protein.Ids"))%>%
    # Normalize values and bring them back to log2 scale
    mutate(AEmed.corr.Protein.Quantity = log2((DigBatch.corr.Protein.Quantity/AE.median.IntraBatch)*AE.median.InterBatch),
           File.Name = as.character(File.Name))%>%
    mutate(AEmed.corr.Protein.Quantity = ifelse(is.infinite(AEmed.corr.Protein.Quantity),
                                                NA,AEmed.corr.Protein.Quantity ))%>%
    reshape2::dcast(Protein.Ids~File.Name,value.var = "AEmed.corr.Protein.Quantity")
  
  rwn<-corr_df$Protein.Ids
  corr_df<-corr_df[,-1]
  corr_df<-do.call(cbind, corr_df)
  rownames(corr_df)<-rwn
  
  return(corr_df)
  
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
    ungroup()%>%
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
    ungroup()%>%
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

fill_in_NAs_knn <- function(df){
  
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
                                                 ifelse(Missing.In.All.Reps,NA , mPQ)
                                          ), 
                                          Protein.Quantity )
    )%>%
    reshape2::dcast(Protein.Ids~File.Name,value.var = "Protein.Quantity.Imp")
  
  rownames(df_fNA) = df_fNA$Protein.Ids
  df_fNA = df_fNA[,-1]
  
  df_imputed = impute::impute.knn(as.matrix(df_fNA),k=1)$data
  #eles = unique(df_fNA$Element)
  # eles = eles[-c(grep("AllEle",eles))]
  #perts = c("Excess","Depletion")
  
  #df_imputed = vector()
  
  #for(e in 1:length(eles)){
  
  # for(p in 1:length(perts)){
  
  #  df_fKNN = filter(df_fNA, Element ==eles[e] & Pert == perts[p])%>%
  #   reshape2::dcast(Protein.Ids~File.Name,value.var = "Protein.Quantity.Imp")
  
  #    rownames(df_fKNN) = df_fKNN$Protein.Ids
  #    df_fKNN = df_fKNN[,-1]
  
  #    df_knnimp = impute::impute.knn(as.matrix(df_fKNN),k=2)$data
  
  #    df_imputed = cbind(df_imputed,df_knnimp)
  #  }
  #}
  
  #cntrl_imp = filter(df_fNA, grepl("AllEle",Element))%>%
  #           reshape2::dcast(Protein.Ids~File.Name,value.var = "Protein.Quantity.Imp")
  #rownames(cntrl_imp) = cntrl_imp$Protein.Ids
  #cntrl_imp = cntrl_imp[,-1]
  
  #cntrl_imp = impute::impute.knn(as.matrix(cntrl_imp),k=2)$data
  
  #df_imputed = cbind(df_imputed,cntrl_imp)
  return(df_imputed)
}





## Function to carry out Differential expression analysis using DEP and limma, make DE plots and return result of DEanalysis

DE_analysis_and_plots <- function(df,dsname,expdes,onlyAEdat=T,pcacompatible=F){
  
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

##
return_Condition <-function(x){
  c=unlist(strsplit(x,"_vs_"))[[1]]
  return(c)
}
return_StatType <-function(x){
  s=unlist(strsplit(x,"_"))[[4]]
  return(s)
}

# Function to convert UniprotIDs to ORFs

convert_Uniprot2singleORF <- function(x){
  
  uid = unlist(strsplit(x,";"))[[1]]
  
  singleORF = unlist(filter(GenProt_SGD,Uniprot.ID==uid)$ORF[[1]])
  
  return(singleORF)
  
}


######################################
### Functional Enrichment Analysis ###
######################################
setwd(dbdir)

GO_gset=read.csv("GO2allGSC.csv",stringsAsFactors = F)
GO_gset_MF=unique(filter(GO_gset,ontology=="MF")[,c("ORF","term")])
GO_gset_BP=unique(filter(GO_gset,ontology=="BP")[,c("ORF","term")])
GO_gset_CC=unique(filter(GO_gset,ontology=="CC")[,c("ORF","term")])

GOslim_BP=read.delim("go_slim_mapping_BP.tab")[,c(1,5)]
colnames(GOslim_BP)<-c("ORF","term")
GOslim_CC=read.delim("go_slim_mapping_CC.tab")[,c(1,5)]
colnames(GOslim_CC)<-c("ORF","term")
GOslim_MF=read.delim("go_slim_mapping_MF.tab")[,c(1,5)]
colnames(GOslim_MF)<-c("ORF","term")

rm(GO_gset)

Reactome<-read.csv("ORF2Reactome.csv",stringsAsFactors = F)
colnames(Reactome)<-c("ORF","term")

KEGG<-read.csv("ORF2KEGG.csv",stringsAsFactors = F)[,c("ORF","KEGG.Term.Name")]
colnames(KEGG)<-c("ORF","term")

KEGG$term<-gsub("[ ][-][ ]Saccharomyces cerevisiae[ ][(]budding yeast[)]","", KEGG$term)

KEGG <- filter(KEGG,!term %in%  c("Not Found in SCe KEGG","Metabolic pathways"))

SGDPhen<-read.csv("ORF2SGDPhenotypes.csv",stringsAsFactors = F)[,c("ORF","phenotype")]
colnames(SGDPhen)<-c("ORF","term")

EC<-read.csv("ORF2EC.csv",stringsAsFactors = F)[,c("ORF","EC.ID")]
colnames(EC)<-c("ORF","term")

Sc_GEM<-read.delim("Scerevisiae_GEM.txt",stringsAsFactors = F)
colnames(Sc_GEM)<-c("ORF","term")

MetPDB <- read.csv("ORF2MetPDB.csv",stringsAsFactors = F)

gsets<-list(GOslim_BP,GOslim_CC,GOslim_MF,
            GO_gset_BP,GO_gset_CC,GO_gset_MF,
            Reactome,
            KEGG,
            SGDPhen,
            EC,
            Sc_GEM,
            MetPDB)

gsetnames<-c("GOslim_BP","GOslim_CC","GOslim_MF",
             "GO_BP","GO_CC","GO_MF",
             "Reactome",
             "KEGG",
             "SGDPhen",
             "EC",
             "Sc_GEM",
             "MetalPDB")

GO_gset_MF_metalbinding <- GO_gset_MF %>%
  filter(grepl("binding",term))%>%
  filter(grepl("calcium",term) |
           grepl("copper",term) |
           grepl("iron",term) |
           grepl("magnesium",term) |
           grepl("manganese",term) |
           grepl("molybdenum", term) |
           grepl("potassium", term) |
           grepl("sodium", term) |
           grepl("zinc", term) |
           grepl("metal",term))%>%
  filter(!grepl("calcium-dependent",term))

MetPDB_biosig <- filter(MetPDB,term %in% c("Ca","Cu","Fe","K","Mg","Mn","Mo","Na","Zn"))


GO_gset_MF_binding_metalwise <- GO_gset_MF %>%
  filter(grepl("binding",term))%>%
  filter(grepl("calcium",term) |
           grepl("copper",term) |
           grepl("iron",term) |
           grepl("magnesium",term) |
           grepl("manganese",term) |
           grepl("molybdenum", term) |
           grepl("potassium", term) |
           grepl("sodium", term) |
           grepl("zinc", term) )%>%
  filter(!grepl("calcium-dependent",term))%>%
  mutate(term = ifelse(grepl("calcium",term),"Ca",
                       ifelse(grepl("copper", term),"Cu",
                              ifelse(grepl("iron",term),"Fe",
                                     ifelse(grepl("magnesium",term),"Mg",
                                            ifelse(grepl("manganese",term),"Mn",
                                                   ifelse(grepl("molybdenum",term),"Mo",
                                                          ifelse(grepl("potassium",term),"K",
                                                                 ifelse(grepl("sodium",term),"Na",
                                                                        ifelse(grepl("zinc",term),"Zn",NA))))))))))%>%
  dplyr::select(ORF,term)%>%
  unique()


#######################################################################
### true positive for metal specific metal binding GOMF and MetPDB ####
#######################################################################

met_specific_metal_binders <- unique(rbind(GO_gset_MF_binding_metalwise,MetPDB_biosig))


##########################################
### Metal transporters in S.cerevisiae ###
##########################################

metal_transporter_genes <- list(
  # Na
  Na = c("NHA1","NHX1","VNX1","PHO89"),
  # K
  K =  c("ENA1","ENA2","ENA5","TRK1","TRK2","KHA1",
         "VNX1","TOK1","NHX1","VHC1"),
  # Cant find new name of ENA3 and ENA4 in yeast databases.  
  # MRS7 renamed to YLH47  -- also removed because not a transporter
  # MDM38 removed because its not a transporter itself 
  
  # Ca
  Ca = c("MID1","CCH1","PMC1","YVC1","PMR1","VCX1"),
  # ecm7 removed because not transporter
  # Fe
  Fe = c("FIT1","FIT2","FIT3","FRE1","FRE2","FRE3","FRE4","FET3",
         "FTR1","CCC2","ARN2","SIT1","FRE6","SMF3",
         "FET5","FTH1","COT1","MRS4","FRE5","ARN1" ,"CCC1"),
  # HMX1 removed because its not a trasporter but involved in heme degradation 
  # Zn
  Zn = c("ZRT1", "ZRT2","ZRT3","FET4","PHO84","MSC2","ZRG17","ZRC1"),
  # Cu
  Cu = c("CTR1","CTR3","FET4","CCC2","FRE6","CTR2"),
  # Mn
  Mn = c("SMF1","SMF2","PHO84","PMR1")
  
)

metal_transporter_philpott_ORFs <- unique(as.character(lapply(unlist(metal_transporter_genes), convert_GeneName2ORF) ))
all_transporters_GOMF_ORFs <- unique(filter(GO_gset_MF,grepl("transporter",term))$ORF)
metal_transporters_GOMF <- unique(filter(filter(GO_gset_MF,grepl("transporter",term)),grepl("metal",term)))
metal_transporters_GOMF_ORFs <- unique(filter(filter(GO_gset_MF,grepl("transporter",term)),grepl("metal",term))$ORF)

# Relevant published datasets

StressResp_Brauer2007 <- read.csv("Brauer2007_SupplementaryTable1.csv",stringsAsFactors = F)

# Read in Common Metal Responsive Genes from Jin et Al 2008
CMRgenes_Jin2008 <- read.csv("Jin2008_CMRgenes.csv",stringsAsFactors = F)



get_plot_HyperGSA<-function(forHyperGSA,plotname,EnrichPV_thresh=0.05,
                            gsLim_low=3,gsLim_high=400,
                            gsLim_low_reacGEM = 3,
                            gsLim_high_reacGEM = 100
){ 
  
  # Function to carry out Hypergeometric Tests on Genesets of all databases,
  # report results in a df 
  # plot results in pdf file ( needs to be opened before calling the function )
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
      
      # plot result of HyperGSA
      
      if( sum(res.HyperGSA$p.adj < EnrichPV_thresh) > 2 ){
        
        # network plot
        EAhyper=networkPlot2(res.HyperGSA, 
                             class="non", 
                             adjusted=T, 
                             significance=EnrichPV_thresh,
                             lay = "layout_nicely",
                             ncharLabel=Inf,
                             nodeSize=c(1,50), 
                             edgeWidth=c(2,20), 
                             labelSize=12,
                             overlap=10, 
                             scoreColors=c(viridis(1,begin=0.2),viridis(1,begin=0.8)),
                             main = paste(gsetnames[gs],plotname,"Hypergeometric Enrichment pAdj < ",EnrichPV_thresh)
        )
        visSave(EAhyper,file=paste0(gsetnames[gs],plotname,"HyperGSA",".html"))
        
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
  
  All_gsets_Stoufferres <- foreach(gs = 1:length(gsets), .combine=rbind) %dopar%{
    
    get_Stouffer_enrichments <- function(gset,df_forStouffer){
      
      # load each gene set one by one
      gsc.forEnrich <- loadGSC(gsets[[gs]])
      
      pvl = df_forStouffer$pv.Adj
      names(pvl) = df_forStouffer$ORF
      drs = df_forStouffer$Log2FC
      names(drs) = df_forStouffer$ORF
      
      res.Stouffer <- runGSA(pvl,drs, 
                             geneSetStat = "stouffer",
                             gsc = gsc.forEnrich,
                             adjMethod = "BH")
      
      
      StouffResTab <- GSAsummaryTable(res.Stouffer)[,c(1,2,5,8)]
      colnames(StouffResTab)=c("Gset.Term.Enriched","Gene.Set.Size",
                               "pAdj_disdir_up","pAdj_disdir_down")
      # Visualize Distinct Directional    
      
      if(gs %in% c(2,5) ){
        if((length(which(res.Stouffer$pAdjDistinctDirUp < EnrichPV_thresh)) +  length(which(res.Stouffer$pAdjDistinctDirDn < EnrichPV_thresh)))> 2 ){
          GSEA_stouff=networkPlot2(res.Stouffer, "distinct", "both", adjusted=T, significance=EnrichPV_thresh, ncharLabel=Inf,
                                   nodeSize=c(2,30), edgeWidth=c(1,8), labelSize=12, overlap=10,
                                   scoreColors=c("#FFFFE0","#FFEC8B","#FFC125","#EE9A00","#F0FFF0","#ADD8E6","#87CEFF","#4F94CD")) 
          visSave(GSEA_stouff,paste0(plotname,"_",gsetnames[gs],"_distinctdir_FDR",EnrichPV_thresh,".html")) }
      }else{
        if((length(which(res.Stouffer$pAdjDistinctDirUp < EnrichPV_thresh)) +  length(which(res.Stouffer$pAdjDistinctDirDn < EnrichPV_thresh)))> 2 ){
          GSEA_stouff=networkPlot2(res.Stouffer, "distinct", "both", adjusted=T, significance=EnrichPV_thresh, ncharLabel=Inf,
                                   nodeSize=c(2,30), edgeWidth=c(1,8), labelSize=12, overlap=10,physics = F,
                                   scoreColors=c("#FFFFE0","#FFEC8B","#FFC125","#EE9A00","#F0FFF0","#ADD8E6","#87CEFF","#4F94CD")) 
          visSave(GSEA_stouff,paste0(plotname,"_",gsetnames[gs],"_distinctdir_FDR",EnrichPV_thresh,".html"))
          
        } 
      }
      
      pdf(file=paste0(plotname,"_",gsetnames[gs],"_Heatmap.pdf"),width = 10,height=20)
      GSAheatmap(res.Stouffer, adjusted=T, ncharLabel=Inf, cellnote='pvalue', cex=1.4)
      dev.off()
      
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

#########################
### Cluster Profiler ### 
#########################

run_clusterProfiler_on_df <- function(df_forClusProf,plotname,numscoretype="Log2(FC)"){
  
  ### Make all relevant plots using ClusterProfiler
  
  # Input : dataframe with :
  # Column 1 : ORFs
  # Column 2 : numeric vector ( can be FC log2(FC) or other numeric that tells you the phenotype)
  
  # clusterProfiler requires a sorted vector but this will be done by the function
  # Note : this function expects all measured proteins/ORFs to be supplied, therefor ## DO NOT prefilter the results
  
  ## numeric vector with phenotype
  geneList <- df_forClusProf[,2]
  
  ## Assign names to the numeric vector
  names(geneList) <- as.character(df_forClusProf[,1])
  
  ## feature 3: sort in decreasing order
  geneList <- sort(geneList, decreasing = TRUE)
  
  ## Hypergeometric test
  
  
  forHyperGeo <- names(geneList)[abs(geneList) > 2]
  forHyperGeo.df <- bitr(forHyperGeo, fromType = "ORF",
                         toType = c("ENTREZID","GENENAME"),
                         OrgDb = org.Sc.sgd.db,
                         drop = F)
  
  gsetnames <- c("CC","MF","BP")
  
  for(gs in 1:length(gsetnames)){
    
    gsn=gsetnames[gs]
    
    hypergeo_res <- enrichGO(gene         = forHyperGeo.df$GENENAME,
                             OrgDb         = org.Sc.sgd.db,
                             keyType       = 'GENENAME',
                             ont           = gsn,
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05)
    
    
    ## Visualize hypergeometric test results 
    
    if(length(hypergeo_res$qvalue) > 2 ){
      
      print(
        cnetplot(hypergeo_res, foldChange=geneList)+
          labs(title=paste("Hypergeometric test",gsn,plotname,"qvalue enrichment 0.05") )
      )
      #heatplot(hypergeo_res, foldChange=geneList)
      print(
        emapplot(hypergeo_res)+
          labs(title=paste("Hypergeometric test",gsn,plotname,"qvalue enrichment 0.05") )
      )
      print(
        goplot(hypergeo_res)+
          labs(title=paste("Hypergeometric test",gsn,plotname,"qvalue threshold 0.05") )
      )
    }
    
    ## GSEA using the named numeric vector called genelist
    
    geneList.forGSEA <- geneList
    nms <- bitr(names(geneList.forGSEA), fromType = "ORF",
                toType = c("GENENAME"),
                OrgDb = org.Sc.sgd.db,
                drop = F)
    names(geneList.forGSEA) <- nms$GENENAME
    
    GSEA_res <- gseGO ( geneList     = geneList.forGSEA,
                        OrgDb        = org.Sc.sgd.db,
                        keyType      = "GENENAME",
                        ont          = gsn,
                        nPerm        = 1000,
                        minGSSize    = 40,
                        maxGSSize    = 500,
                        pvalueCutoff = 0.05,
                        verbose      = FALSE)
    ## Visualize GSEA results 
    
    if(length(GSEA_res$p.adjust) > 2 ){
      
      print(
        dotplot(GSEA_res)+ scale_colour_viridis_c()+
          theme_classic(base_size = 25)+
          labs(title = paste("GSEA",gsn,plotname,"pvalue cutoff 0.05"))
      )
      
      print(
        cnetplot(GSEA_res, 
                 foldChange = geneList.forGSEA,
                 colorEdge = TRUE)+
          labs(title = paste("GSEA",gsn,plotname,"pvalue cutoff 0.05"))
      )
      
      print(
        enrichplot::gseaplot2(GSEA_res, 
                              geneSetID = 1:length(GSEA_res$Description),
                              pvalue_table = F,
                              ES_geom ="dot",
                              paste("GSEA",gsn,plotname,"pvalue cutoff 0.05"))
      )
      print(
        heatplot(GSEA_res, foldChange=geneList.forGSEA)+
          scale_fill_gradient2( low = "blue", mid = "white",
                                high = "red", space = "Lab" )+
          theme(axis.text.x = element_text(angle=90))+
          labs(title = paste("GSEA",gsn,plotname,"pvalue cutoff 0.05"),
               fill = numscoretype)
      )
      print(
        upsetplot(GSEA_res,n=5)+
          theme_classic(base_size = 25)+
          labs(title = paste("GSEA",gsn,plotname,"pvalue cutoff 0.05"))
      )
      
    }
    
  }
  
  
  #############
  ###  KEGG ###
  #############
  
  
  geneList.forKEGG <- geneList
  
  nms <- bitr(names(geneList.forKEGG), 
              fromType = "ORF",
              toType = c("UNIPROT"),
              OrgDb = org.Sc.sgd.db,
              drop = F)
  ## keep only one mapping for each ORF 
  if(sum(duplicated(nms$ORF))>0){
    nms<-nms[-which(duplicated(nms$ORF)),]
  }
  
  names(geneList.forKEGG) <- nms$UNIPROT
  
  # KEGG Hypergeometric test
  
  genes4KEGG_hypG <- names(geneList.forKEGG)[abs(geneList.forKEGG) > 2]
  
  HyperGeo_kegg_res<- enrichKEGG( gene         = genes4KEGG_hypG,
                                  organism     = 'sce',
                                  keyType      = "uniprot",
                                  pvalueCutoff = 0.05 )
  # KEGG module Hypergeometric test
  
  HyperGeo_keggMod_res<-enrichMKEGG( gene         = genes4KEGG_hypG,
                                     organism     = 'sce',
                                     keyType      = "uniprot",
                                     pvalueCutoff = 0.05 )
  
  ## Visualize kegg HyperGeo results
  
  # normal
  
  if(length(HyperGeo_kegg_res$qvalue) > 2 ){
    
    print(
      cnetplot(HyperGeo_kegg_res, foldChange=geneList.forKEGG)+
        labs(title=paste("Hypergeometric test kegg",plotname,"qvalue enrichment 0.05") )
    )
    #heatplot(hypergeo_res, foldChange=geneList)
    print(
      emapplot(HyperGeo_kegg_res)+
        labs(title=paste("Hypergeometric test kegg",plotname,"qvalue enrichment 0.05") )
    )
    
  }
  
  # kegg module
  
  if(length(HyperGeo_keggMod_res$qvalue) > 2 ){
    
    print(
      cnetplot(HyperGeo_keggMod_res, foldChange=geneList)+
        labs(title=paste("Hypergeometric test keggmodule",plotname,"qvalue enrichment 0.05") )
    )
    #heatplot(hypergeo_res, foldChange=geneList)
    print(
      emapplot(HyperGeo_keggMod_res)+
        labs(title=paste("Hypergeometric test keggmodule",plotname,"qvalue enrichment 0.05") )
    )
    
  }
  
  # KEGG GSEA
  
  GSEA_kegg_res <- gseKEGG(geneList     = geneList.forKEGG,
                           organism     = 'sce',
                           keyType = "uniprot",
                           nPerm        = 1000,
                           pvalueCutoff = 0.05,
                           verbose      = FALSE)
  
  
  
  GSEA_keggMod_res<-gseMKEGG(geneList     = geneList.forKEGG,
                             organism     = 'sce',
                             keyType = "uniprot",
                             nPerm        = 1000,
                             pvalueCutoff = 0.05,
                             verbose      = FALSE)
  ## Visualize kegg gsea results
  
  
  ## normal kegg
  
  if(length(GSEA_kegg_res$p.adjust) > 2 ){
    
    print(
      dotplot(GSEA_kegg_res)+ scale_colour_viridis_c()+
        theme_classic(base_size = 25)+
        labs(title = paste("GSEA kegg",plotname,"pvalue cutoff 0.05"))
    )
    
    print(
      cnetplot(GSEA_kegg_res, 
               foldChange = geneList.forKEGG,
               colorEdge = TRUE)+
        labs(title = paste("GSEA kegg",plotname,"pvalue cutoff 0.05"))
    )
    
    print(
      enrichplot::gseaplot2(GSEA_kegg_res, 
                            geneSetID = 1:length(GSEA_kegg_res$Description),
                            pvalue_table = F,
                            ES_geom ="dot",
                            paste("GSEA kegg",plotname,"pvalue cutoff 0.05"))
    )
    print(
      heatplot(GSEA_kegg_res, foldChange=geneList.forKEGG)+
        scale_fill_gradient2( low = "blue", mid = "white",
                              high = "red", space = "Lab" )+
        theme(axis.text.x = element_text(angle=90))+
        labs(title = paste("GSEA kegg",plotname,"pvalue cutoff 0.05"),
             fill = numscoretype)
    )
    print(
      upsetplot(GSEA_kegg_res,n=5)+
        theme_classic(base_size = 25)+
        labs(title = paste("GSEA kegg",plotname,"pvalue cutoff 0.05"))
    )
  }
  
  ## kegg module
  
  if(length(GSEA_keggMod_res$p.adjust) > 2 ){
    
    print(
      dotplot(GSEA_keggMod_res)+ scale_colour_viridis_c()+
        theme_classic(base_size = 25)+
        labs(title = paste("GSEA kegg module",plotname,"pvalue cutoff 0.05"))
    )
    
    print(
      cnetplot(GSEA_keggMod_res, 
               foldChange = geneList.forKEGG,
               colorEdge = TRUE)+
        labs(title = paste("GSEA kegg module",plotname,"pvalue cutoff 0.05"))
    )
    
    print(
      enrichplot::gseaplot2(GSEA_keggMod_res, 
                            geneSetID = 1:length(GSEA_keggMod_res$Description),
                            pvalue_table = F,
                            ES_geom ="dot",
                            paste("GSEA kegg module",plotname,"pvalue cutoff 0.05"))
    )
    print(
      heatplot(GSEA_keggMod_res, foldChange=geneList.forKEGG)+
        scale_fill_gradient2( low = "blue", mid = "white",
                              high = "red", space = "Lab" )+
        theme(axis.text.x = element_text(angle=90))+
        labs(title = paste("GSEA kegg module",plotname,"pvalue cutoff 0.05"),
             fill = numscoretype)
    )
    print(
      upsetplot(GSEA_keggMod_res,n=5)+
        theme_classic(base_size = 25)+
        labs(title = paste("GSEA kegg module",plotname,"pvalue cutoff 0.05"))
    )
  }
}




prepare_data_for_iPATH_logFC <- function (df){
  
  ### Function to prepare dataframe for iPATH
  ## Input dataframe should have :
  # Column 1 : UNIPROT IDs
  # Column 2 : numeric column with log2(fc) or Zscore or other numeric vector
  
  colnames(df) <- c("UNIPROT.ID","log2_FC")
  
  df_for_iPATH <- df %>%
    mutate(colour = ifelse(log2_FC < -1.5, "#2D42D5", 
                           ifelse(log2_FC > -1.5 & log2_FC < 1.5,"#CDC1C5",
                                  "#C2002D" ))) %>%
    mutate(width = cut_interval( abs(log2_FC), 
                                 n=10,
                                 labels = paste0("W",seq(11,20,1) )
    )
    )%>%
    mutate(UNIPROT.ID=paste0("UNIPROT:",UNIPROT.ID))%>%
    dplyr::select(UNIPROT.ID,colour,width)
  return(df_for_iPATH)
}


prepare_data_for_iPATH_upDown <- function (df){
  
  ### Function to prepare dataframe for iPATH
  ## Input dataframe should have :
  # Column 1 : UNIPROT IDs
  # Column 2 : Factor  column with  "Up" "Down" or "NoSigChange"
  
  colnames(df) <- c("UNIPROT.ID","Direction")
  
  df_for_iPATH <- df %>%
    mutate(colour = ifelse(Direction == "Down", "#3A5FCD", 
                           ifelse(Direction == "Up","#B0171F",
                                  "#A2B5CD")),
           width = "W13",
           opacity="1",
           UNIPROT.ID=paste0("UNIPROT:",UNIPROT.ID))%>%
    dplyr::select(UNIPROT.ID,colour,width,opacity)%>%
    data.frame()
  return(df_for_iPATH)
}



prepare_data_for_iPATH_pvalue <- function (df,ele){
  
  ### Function to prepare dataframe for iPATH
  ## Input dataframe should have :
  # Column 1 : UNIPROT IDs
  # Column 2 : numeric column with adjusted pvalues
  
  colnames(df) <- c("UNIPROT.ID","adjPVal")
  
  ## Make colour scale based on element
  ele_cols<-c("white","white","white",colkey_EleDir[which(grepl(ele,names(colkey_EleDir)))])
  
  
  ele_cols = colorRampPalette(ele_cols,space="Lab")(101)
  
  # Function to convert pvalues values to colour 
  
  pvAdj2colour_elewise<-function(x){
    
    col_df <- data.frame(palcol=ele_cols,
                         numcol=round(seq(from = 1,to = 0,length.out = 101),2))
    rownames(col_df) <- col_df$numcol
    res_col <- col_df[which(rownames(col_df)== round(x,2)) ,"palcol"]
    
    cols<-col_df$palcol
    names(cols)<-col_df$palcol
    
    pdf(paste(ele,"iPath Colour scale.pdf"),width=4,height=2)
    ggplot(col_df,
           aes(x=numcol,
               y=1,
               fill=palcol))+
      geom_tile()+
      scale_fill_manual(values=cols)+
      theme_minimal()+
      theme(
        legend.position="none",
        axis.text.x = element_text(size=25),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
      labs(x=paste(ele,"Adj Pvalue"),y="")
    dev.off()
    if(x > 0.2){
      res_col <- colkey_EleDir[which(grepl(ele,names(colkey_EleDir)))][[1]]
    }
    return(as.character(res_col))
  }
  
  pvAdj2width<-function(x){
    
    if(x > 0.2){
      res_w <- "W3"
    }else{
      width_df <- data.frame(widths = round(seq(4,13,length.out =21),1),
                             numcol=round(seq(from = 0.2,to = 0,length.out = 21),2))
      rownames(width_df) <- width_df$numcol
      res_w <- paste0("W",width_df[which(rownames(width_df)== round(x,0)) ,"widths"])
    }
    return(as.character(res_w))
  }
  
  pvAdj2opac<-function(x){
    
    if(x > 0.2){
      res_o <- 1
    }else{
      opacity_df <- data.frame(opacity = seq(0.1,0.8,length.out = 21),
                               numcol=round(seq(from = 0.2,to = 0,length.out = 21),2))
      rownames(opacity_df) <- opacity_df$numcol
      res_o <- opacity_df[which(rownames(opacity_df)== round(x,2)) ,"opacity"]
    }
    return(res_o)
  }
  
  df_for_iPATH <- df %>%
    mutate(colour = as.character(lapply(adjPVal,pvAdj2colour_elewise))) %>% 
    mutate(width = as.character(lapply(adjPVal,pvAdj2width)))%>%
    mutate(opacity = as.character(lapply(adjPVal,pvAdj2opac)))%>%
    mutate(UNIPROT.ID=paste0("UNIPROT:",UNIPROT.ID))%>%
    dplyr::select(UNIPROT.ID,colour,width,opacity)
  
  return(df_for_iPATH)
}

##############################################################
### Scale Rows of a dataframe with a given control column ####
##############################################################

# Calculates Z score with control as centre and sd of entire dataframe 
## takes  logged values
scale_rows_vsMeanAEsdAS <- function(df,cntrl_column="AllEle 1",ASstats_df,idtype = "Genes"){
  
  ## Input in in log
  df=2^df
  cntrl = df[,cntrl_column]
  df = sweep(df,MARGIN = 1,STATS = cntrl,FUN = "-")
  
  df_forscaling = df%>%
    mutate( IDs = rownames(df))
  colnames(df_forscaling)[which(colnames(df_forscaling) =="IDs")] <- idtype
  
  df_forscaling = merge(df_forscaling,ASstats_df,by=idtype)
  rownames(df_forscaling)= df_forscaling[,idtype]
  
  row_sds = df_forscaling[,"SD.Protein.Quant"]
  
  df_scaled = df_forscaling[,-c(which(colnames(df_forscaling) %in%
                                        c("Genes","Protein.Ids","SD.Protein.Quant","Mean.Protein.Quant")))]
  
  df_scaled = sweep(df_scaled,MARGIN = 1,STATS = row_sds,FUN = "/")
  
  return(df_scaled)
}

scale_rows_vsMeanASsdAS <- function(df,ASstats_df,idtype="Genes"){
  
  df=2^df
  df= round(df,4)
  ASstats_df[,2:3]= round(ASstats_df[,2:3],4)
  
  df_forscaling = df%>%
    mutate(IDs = rownames(df))
  colnames(df_forscaling)[which(colnames(df_forscaling) =="IDs")] <- idtype
  
  df_forscaling = merge(df_forscaling,ASstats_df,by=idtype)
  rownames(df_forscaling)= df_forscaling[,idtype]
  
  row_means = df_forscaling[,"Mean.Protein.Quant"]
  row_sds = df_forscaling[,"SD.Protein.Quant"]
  
  df_scaled = df_forscaling[,-c(which(colnames(df_forscaling) %in%
                                        c("Genes","Protein.Ids","SD.Protein.Quant",
                                          "Mean.Protein.Quant")))]
  
  df_scaled = sweep(df_scaled,MARGIN = 1,STATS = row_means,FUN = "-")
  df_scaled = sweep(df_scaled,MARGIN = 1,STATS = row_sds,FUN = "/")
  
  return(df_scaled)
}


scale_rows_vs_AllEle_meansdAE <- function(df,AEstats_df){
  
  df= round(df,4)
  AEstats_df[,2:5]= round(AEstats_df[,2:5],4)
  
  df_forscaling = df%>%
    mutate(Protein.Ids = rownames(df))
  
  df_forscaling = merge(df_forscaling,AEstats_df,by="Protein.Ids")
  rownames(df_forscaling)= df_forscaling$Protein.Ids
  
  row_means = df_forscaling[,"Mean.Protein.Quant"]
  row_sds = df_forscaling[,"SD.Protein.Quant"]
  
  df_scaled = df_forscaling[,-c(which(colnames(df_forscaling) %in%
                                        c("Protein.Ids","SD.Protein.Quant",
                                          "Mean.Protein.Quant",
                                          "MAD.Protein.Quant",
                                          "Median.Protein.Quant",
                                          "Cov")))]
  
  df_scaled = sweep(df_scaled,MARGIN = 1,STATS = row_means,FUN = "-")
  df_scaled = sweep(df_scaled,MARGIN = 1,STATS = row_sds,FUN = "/")
  
  return(df_scaled)
}

scale_rows_vs_AllEle_medianmadAE <- function(df,AEstats_df){
  
  df= round(df,4)
  AEstats_df[,2:5]= round(AEstats_df[,2:5],4)
  
  df_forscaling = df%>%
    mutate(Protein.Ids = rownames(df))
  
  df_forscaling = merge(df_forscaling,AEstats_df,by="Protein.Ids")
  rownames(df_forscaling)= df_forscaling$Protein.Ids
  
  row_medians = df_forscaling[,"Median.Protein.Quant"]
  row_mads = df_forscaling[,"MAD.Protein.Quant"]
  
  df_scaled = df_forscaling[,-c(which(colnames(df_forscaling) %in%
                                        c("Protein.Ids","SD.Protein.Quant",
                                          "Mean.Protein.Quant",
                                          "MAD.Protein.Quant",
                                          "Median.Protein.Quant",
                                          "Cov")))]
  
  df_scaled = sweep(df_scaled,MARGIN = 1,STATS = row_medians,FUN = "-")
  df_scaled = sweep(df_scaled,MARGIN = 1,STATS = row_mads,FUN = "/")
  
  return(df_scaled)
}


#####################################################
### Function to plot result of linear regression  ### 
#####################################################

ggplotRegression <- function (fit,label_col) {
  
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], 
                               y = names(fit$model)[1])) + 
    geom_point(size=3,colour="navyblue",alpha=0.5  ) +
    stat_smooth(method = "lm", col = "darkred") +
    geom_text_repel(aes(label=!!sym(label_col)))+
    labs(size=2,title = paste("R2 = ",signif(summary(fit)$r.squared, 3),
                              "Adj R2 = ",signif(summary(fit)$adj.r.squared, 3),
                              "Intercept =",signif(fit$coef[[1]],3 ),"\n",
                              "Slope =",signif(fit$coef[[2]], 3),
                              " P =",signif(summary(fit)$coef[2,4], 3)))+
    theme_SKA(base_size = 20)
}



############################
### Annotation Status db ###
############################

setwd(dbdir)
annostatdb=read.delim("go_slim_mapping.tab",stringsAsFactors=F,sep="\t",header=F)
colnames(annostatdb)=c("ORF","Protein Name","SGD Number","GO Term Type","GO Term Description","GO Term Number","ORF Type")
annostatdb<-unique(annostatdb[,c("ORF","ORF Type")])



############################
### Ipath in R functions ###
############################

to_parameters_iPATH_SKA = function(ipath_data,
                                   default_color = '#bdbdbd',
                                   default_width = 3,
                                   default_radius = 7,
                                   default_opacity = 1,
                                   keep_colors = FALSE,
                                   whole_modules = FALSE,
                                   whole_pathways = FALSE,
                                   query_reactions = FALSE,
                                   tax_filter = 'sce',
                                   export_dpi = 1200){
  
  selection  = paste0(apply(ipath_data, 1, function(x) paste(x, collapse =" ")),
                      collapse = "\n")
  
  ipath_parameters = list(selection = selection,
                          export_type = 'svg',
                          keep_colors = ifelse(keep_colors, 1, 0),
                          include_metabolic = 1,
                          include_secondary = 1,
                          include_antibiotic = 0,
                          include_microbial = 1,
                          whole_modules = ifelse(whole_modules, 1, 0),
                          default_opacity = ifelse(default_opacity, 1, 0),
                          whole_pathways = ifelse(whole_pathways, 1, 0),
                          default_width = default_width,
                          default_color = default_color,
                          default_radius = default_radius,
                          query_reactions = ifelse(query_reactions, 1, 0),
                          tax_filter = tax_filter,
                          export_dpi = export_dpi
  )
  return(ipath_parameters)
}

################################################################################################################################
### Function to plot protein quantities (as point + fit plot of selected proteins along conc gradient with BioSpecID colours ###
################################################################################################################################


plot_selected_protQuantpoints <- function(dfp,selkeywd){
  
  # dfp for this project should be PQ_lm_res_df or some derivative of it that has log2(protein quant), LeasComplexModel
  # Element.Concentration & BioSpecID values
  
  mods = unique(dfp$LeastComplexModel)
  mods = mods[-which(mods=="null")]
  
  for(i in 1:length(mods)){
    
    
    nlin = nrow(unique(filter(dfp,LeastComplexModel == "linear")[,c("Genes","Element")]))
    
    if(nlin > 0){
      if(nlin >= 3){
        pdf(paste0(selkeywd," linear.pdf"),width = round(sqrt(nlin),0)*2.6,
            height = round(sqrt(nlin),0)*2.75)  
        print(
          ggplot(filter(dfp,LeastComplexModel == "linear"),
                 aes(x=log2(Element.Concentration),
                     y=Log2FC_vs_AE,
                     colour=BioSpecID,
                     group = Element))+
            facet_wrap(c("Genes","Element"),scales="free", ncol =  round(sqrt(nlin),0))+
            geom_point(size=3,alpha=0.8)+
            geom_smooth(aes(group=Element),method="lm",formula=y~poly(x,1),colour="black")+
            scale_colour_manual(values=colkey_BioSpecID)+
            theme_SKA()+
            theme(legend.position = "none")
        )
        dev.off()
      }else{
        pdf(paste0(selkeywd," linear.pdf"),width = 7,
            height = 5)  
        print(
          ggplot(filter(dfp,LeastComplexModel == "linear"),
                 aes(x=log2(Element.Concentration),
                     y=Log2FC_vs_AE,
                     colour=BioSpecID,
                     group = Element))+
            facet_wrap(c("Genes","Element"),scales="free", ncol =  2)+
            geom_point(size=3,alpha=0.8)+
            geom_smooth(aes(group=Element),method="lm",formula=y~poly(x,1),colour="black")+
            scale_colour_manual(values=colkey_BioSpecID)+
            theme_SKA()+
            theme(legend.position = "none")
        )
        dev.off()
        
      }
    }
    
    nquad = nrow(unique(filter(dfp,LeastComplexModel == "quadratic")[,c("Genes","Element")]))
    
    if(nquad > 0){
      if(nquad >=3){
        pdf(paste0(selkeywd," quadratic.pdf"),width = round(sqrt(nquad),0)*2.6,
            height = round(sqrt(nquad),0)*2.9)  
        print(
          ggplot(filter(dfp,LeastComplexModel == "quadratic"),
                 aes(x=log2(Element.Concentration),
                     y=Log2FC_vs_AE,
                     colour=BioSpecID,
                     group = Element))+
            facet_wrap(c("Genes","Element"),scales="free", ncol =  round(sqrt(nquad),0))+
            geom_point(size=3,alpha=0.8)+
            geom_smooth(aes(group=Element),method="lm",formula=y~poly(x,2),colour="black")+
            scale_colour_manual(values=colkey_BioSpecID)+
            theme_SKA()+
            theme(legend.position = "none")
        )
        dev.off() 
      }else{
        pdf(paste0(selkeywd," quadratic.pdf"),width = 7,
            height = 5) 
        print(
          ggplot(filter(dfp,LeastComplexModel == "quadratic"),
                 aes(x=log2(Element.Concentration),
                     y=Log2FC_vs_AE,
                     colour=BioSpecID,
                     group = Element))+
            facet_wrap(c("Genes","Element"),scales="free", ncol =  2)+
            geom_point(size=3,alpha=0.8)+
            geom_smooth(aes(group=Element),method="lm",formula=y~poly(x,2),colour="black")+
            scale_colour_manual(values=colkey_BioSpecID)+
            theme_SKA()+
            theme(legend.position = "none")
        )
        dev.off()
        
      }
    }
    
    ncub = nrow(unique(filter(dfp,LeastComplexModel == "cubic")[,c("Genes","Element")]))
    
    if(ncub > 0){
      if(ncub >= 3){
        
        pdf(paste0(selkeywd," cubic.pdf"),width = round(sqrt(ncub),0)*2.6,
            height = round(sqrt(ncub),0)*2.75)  
        print(
          ggplot(filter(dfp,LeastComplexModel == "cubic"),
                 aes(x=log2(Element.Concentration),
                     y=Log2FC_vs_AE,
                     colour=BioSpecID,
                     group = Element))+
            facet_wrap(c("Genes","Element"),scales="free", ncol =  round(sqrt(ncub),0))+
            geom_point(size=3,alpha=0.8)+
            geom_smooth(aes(group=Element),method="lm",formula=y~poly(x,3),colour="black")+
            scale_colour_manual(values=colkey_BioSpecID)+
            theme_SKA()+
            theme(legend.position = "none")
        )
        dev.off()
      }else{
        
        pdf(paste0(selkeywd," cubic.pdf"),width = 7,
            height = 5)  
        print(
          ggplot(filter(dfp,LeastComplexModel == "cubic"),
                 aes(x=log2(Element.Concentration),
                     y=Log2FC_vs_AE,
                     colour=BioSpecID,
                     group = Element))+
            facet_wrap(c("Genes","Element"),scales="free", ncol =  round(sqrt(ncub),0))+
            geom_point(size=3,alpha=0.8)+
            geom_smooth(aes(group=Element),method="lm",formula=y~poly(x,3),colour="black")+
            scale_colour_manual(values=colkey_BioSpecID)+
            theme_SKA()+
            theme(legend.position = "none")
        )
        dev.off()
        
      }
    }
  }
}


########################################################################################################################
### Function to plot protein quantities (as heatmap) of selected proteins along conc gradient with BioSpecID colours ###
########################################################################################################################

plot_selected_protHeatmap <- function (df,
                                       is_rowanno=F,
                                       rowanno_df = NA,
                                       selkeywd ="none"){
  
  ## Column annotation based on BioSpecID colours
  colanno_df <- as.data.frame(cbind(BioSpecID = colnames(df)))
  rownames(colanno_df)<-colanno_df$BioSpecID
  
  col_anno = HeatmapAnnotation(df=colanno_df,
                               col = list ( BioSpecID = colkey_BioSpecID),
                               show_legend= c("col"=F)
  )
  
  if(is_rowanno==T){
    row_anno = rowAnnotation(term = rowanno_df$term)
    
    pdf(paste0("Heatmap ",selkeywd,".pdf"),width = 30,height = nrow(df)*0.7)
    
    draw(ComplexHeatmap::Heatmap(as.matrix(df),
                                 #col=col_fun,
                                 # cluster_columns = dendc,
                                 # cluster_rows = dendr,
                                 use_raster = T,
                                 raster_device = "CairoPNG",
                                 #column_split = kc,
                                 
                                 #row_split = kr,
                                 row_gap = unit(3, "mm"),
                                 row_names_gp = gpar(fontsize=8),   
                                 
                                 column_names_gp = gpar(fontsize = 11),
                                 column_gap = unit(2, "mm"),
                                 column_names_side = "top",
                                 
                                 heatmap_legend_param = list(title = "FC vs AE"),
                                 top_annotation = col_anno,
                                 right_annotation = row_anno
    ),
    merge_legend = TRUE
    
    )
    dev.off()
  }else{
    
    pdf(paste0("Heatmap ",selkeywd,".pdf"),width = 30,height =  nrow(df)*0.7)
    
    draw(ComplexHeatmap::Heatmap(as.matrix(df),
                                 #col=col_fun,
                                 # cluster_columns = dendc,
                                 # cluster_rows = dendr,
                                 use_raster = T,
                                 raster_device = "CairoPNG",
                                 #column_split = kc,
                                 
                                 #row_split = kr,
                                 row_gap = unit(3, "mm"),
                                 row_names_gp = gpar(fontsize=8),   
                                 
                                 column_names_gp = gpar(fontsize = 11),
                                 column_gap = unit(2, "mm"),
                                 column_names_side = "top",
                                 
                                 heatmap_legend_param = list(title = "FC vs AE"),
                                 top_annotation = col_anno,
    ),
    merge_legend = TRUE
    
    )
    dev.off()
    
  }
}

#############################################################
### Plot selected protein heatmap with multiple row annos ###
#############################################################


plot_selected_protHeatmap_multiRowAnno <- function (df,
                                                    RAnno = NA,
                                                    selkeywd ="none"){
  
  ## Column annotation based on BioSpecID colours
  colanno_df <- as.data.frame(cbind(BioSpecID = colnames(df)))
  rownames(colanno_df)<-colanno_df$BioSpecID
  
  col_anno = HeatmapAnnotation(df=colanno_df,
                               col = list ( BioSpecID = colkey_BioSpecID),
                               show_legend= c("col"=F)
  )
  
  
  
  pdf(paste0("Heatmap ",selkeywd,".pdf"),width = 30,height = nrow(df)*0.2)
  
  draw(ComplexHeatmap::Heatmap(as.matrix(df),
                               #col=col_fun,
                               # cluster_columns = dendc,
                               # cluster_rows = dendr,
                               # raster_device = "CairoPNG",
                               #  use_raster = T,
                               #column_split = kc,
                               
                               #row_split = kr,
                               row_gap = unit(3, "mm"),
                               row_names_gp = gpar(fontsize=8),   
                               
                               column_names_gp = gpar(fontsize = 11),
                               column_gap = unit(2, "mm"),
                               column_names_side = "top",
                               
                               heatmap_legend_param = list(title = "Zscore vs AE"),
                               top_annotation = col_anno,
                               right_annotation = RAnno
  ),
  merge_legend = TRUE
  
  )
  dev.off()
  
}

############################################################################
### Function to plot phenotype in KO screen of a group of selected genes ###
############################################################################

plot_selected_KOpointline <- function(df,
                                      selkeywd ="none"){
  
  # df should be a derivative of KO growth results file in allgenomewide datasets folder with ORF names converted to Genes
  # Should have the following columns : "Condition", "Genes" and "Fitness.Relative.to.AE" columns
  # Columns needed are : "Genes","Condition","Fitness.Relative.to.AE",
  #"Fitness.Relative.to.AE.perRep","Significant"   
  # Fitness.Relative.to.AE.perRep = Colony_size_corr_postQC/Median.ColSize.in.AE
  
  # Per protein tiles plots
  gens = unique(df$Genes)
  
  plot_list <- list()
  for( g in 1:length(gens)){
    
    df <- df%>%
      group_by(Genes, Condition)%>%
      mutate(sd_reps = sd(Fitness.Relative.to.AE.perRep, na.rm = T),
             mean_reps = mean(Fitness.Relative.to.AE.perRep, na.rm=T))
    
    
    
    plot_list[[g]]<- ggplot(filter(df,Genes == gens[[g]]),
                            aes(x=Condition,
                                y=mean_reps))+
      geom_errorbar(aes(ymin = mean_reps - sd_reps, ymax = mean_reps+sd_reps), width =0.25)+
      geom_point(aes(colour=factor(Significant)))+
      scale_colour_manual(values = c("gray","darkred"))+
      geom_line(aes(x=Condition,
                    y=Fitness.Relative.to.AE,
                    group=Genes),
                colour="black")+
      labs(x="",y = "",title=gens[[g]])+
      theme_SKA()+
      theme(legend.position="none",
            axis.text.x = element_text(angle=90))
    
  }
  
  if(length(plot_list) < 3){
    pdf(paste0("KOgrowth phen ",selkeywd,".pdf"),
        width=10,height = length(plot_list)*4)
    multiplot(plotlist = plot_list, cols = 1)
    dev.off()
  }else if(length(plot_list) < 40){
    pdf(paste0("KOgrowth phen ",selkeywd,".pdf"),
        width=10,height = length(plot_list)*3)
    multiplot(plotlist = plot_list, cols = 1)
    dev.off()
  }
  
  ## Plot all genes on the same plot # PS -- don't do this for more than 10-15 proteins,
  # the plot becomes hard to read
  pdf(paste0("KOgrowth phen all prot together ",selkeywd,".pdf"),width = 10,height=7)
  
  print(
    ggplot(unique(df[,c("Genes","Condition","mean_reps","sd_reps")]),
           aes(x=Condition, y=mean_reps,
               colour=Genes))+
      geom_errorbar(aes(ymin = mean_reps - sd_reps, ymax = mean_reps+sd_reps, group = paste(Genes,Condition)),
                    width =0.25, position = position_dodge(width =0.9), colour = "black")+
      geom_point(size=3,alpha=0.8,position = position_dodge(width=0.9))+
      
      #geom_line(aes(x=Condition,
      #             y=mean_reps,
      #            group=Genes))+
      labs(x="",y = "Fitness relative to AE")+
      theme_SKA()+
      theme(legend.position="bottom",
            axis.text.x = element_text(angle=90))
  )
  dev.off()
  
}



##################################################################################
### Function to plot phenotype in hdPCA phenotype of a group of selected genes ###
##################################################################################

plot_selected_hdPCApointline <- function(df,  selkeywd ="none"){
  
  # df should be a derivative of hdpca results file (hdPCA CyclicLoessNorm colony sizes w Statistical Test Results.csv) in allgenomewide datasets folder with ORF names converted to Genes
  # Columns needed are : Genes,Condition,Ttest_Welchs.FDR,Mean.of.Replicates.Log2.Colony.size # Convert ORF to Genes with convertORG2SingleGeneName before calling this function
  #"Fitness.Relative.to.AE.perRep","Significant"   
  # Fitness.Relative.to.AE.perRep = Colony_size_corr_postQC/Median.ColSize.in.AE
  
  # Per protein tiles plots
  gens = unique(df$Genes)
  
  plot_list <- list()
  for( g in 1:length(gens)){
    
    plot_list[[g]]<- ggplot(filter(df,Genes == gens[[g]]),
                            aes(x=Condition,
                                y=Rel.Mean.Log2.CS))+
      geom_point(aes(colour=factor(Significance.ttw.FDR.01)))+
      scale_colour_manual(values = c("gray","darkred"))+
      geom_line(aes(x=Condition,
                    y=Rel.Mean.Log2.CS,
                    group=Genes),
                colour="BLACK")+
      labs(x="",y = "",title=gens[[g]])+
      theme_SKA()+
      theme(legend.position="none",
            axis.text.x = element_text(angle=90))
    
  }
  
  if(length(plot_list) < 3){
    pdf(paste0("hdPCA phen ",selkeywd,".pdf"),
        width=10,height = length(plot_list)*4)
    multiplot(plotlist = plot_list, cols = 1)
    dev.off()
  }else if(length(plot_list) < 40){
    pdf(paste0("hdPCA phen ",selkeywd,".pdf"),
        width=10,height = length(plot_list)*3)
    multiplot(plotlist = plot_list, cols = 1)
    dev.off()
  }
  
  ## Plot all genes on the same plot # PS -- don't do this for more than 10-15 proteins,
  # the plot becomes hard to read
  pdf(paste0("hdPCA phen all prot together ",selkeywd,".pdf"),width = 10,height=7)
  
  print(ggplot(df, aes(x=Condition, y=Rel.Mean.Log2.CS,
                       colour=Genes))+
          geom_point(size=3,alpha=0.8)+
          geom_line(aes(x=Condition,
                        y=Rel.Mean.Log2.CS,
                        group=Genes))+
          labs(x="",y = "protein dimerization relative to AE")+
          theme_SKA()+
          theme(legend.position="bottom",
                axis.text.x = element_text(angle=90))
  )
  dev.off()
  
}

#############################################################################
### Function to plot pathview plots of a set of pathview ids and elements ###
#############################################################################

plot_pathview <- function(df,pids,eles4PV,keywd4PV){
  
  
  # df = PQ_lm_res_df or derivative of it .. i.e. filtered on Genes
  # pids - list of pathview ids to plot 
  # eles - list of elements for which to make the pathview plots
  
  for(e in 1:length(eles4PV)){
    
    df4PV <- unique(filter(df,Element == eles4PV[e])[,c("Genes","Element.Concentration","Median.AE.value",
                                                        "Log2.Protein.Quantity.")])%>%
      mutate(ORF = as.character(lapply(Genes,convert_GeneName2ORF)),
             log2_FC_vs_AE = Log2.Protein.Quantity. - Median.AE.value)%>%
      dplyr::select(ORF,Element.Concentration, log2_FC_vs_AE)%>%
      group_by(ORF,Element.Concentration)%>%
      summarize(mean_log2_FC_vs_AE = mean(log2_FC_vs_AE,na.rm = T))%>%
      ungroup()%>%
      reshape2::dcast(ORF~Element.Concentration, value.var = "mean_log2_FC_vs_AE")
    
    write.csv(df4PV,
              paste0(keywd4PV, eles4PV[e]," df for pathView.csv"),row.names = F)
    
    rownames(df4PV) <- df4PV$ORF
    df4PV <- df4PV[,-1]
    ## Plot pathview pathways
    
    for( p in 1:length(pids)){
      
      pv.out <- pathview(gene.data = df4PV, 
                         species = "sce",
                         pathway.id = pids[p],
                         gene.idtype="ORF",
                         kegg.native = T,
                         out.suffix = paste0(eles4PV[e]," ", pids[p]," pathview"), 
                         low = list(gene = "darkblue", cpd = "purple"), 
                         high = list(gene = "darkred", cpd = "yellow"),
                         mid = list(gene = "white", cpd = "white")
      )
      
    }
    
    
  }
}






#########################
### Accessory Version ###
#########################

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#########
### Plots to check intensity profiles ###


plot_protint_dist <- function(df,plotname){
  
  df<- reshape::melt(df)
  
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
      theme_SKA()
    
  )
  
}
plot_protint_dist_facet <- function(df,plotname){
  
  df<- reshape::melt(df)
  
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
      facet_wrap("Layout")+
      labs(title = plotname)+
      theme_SKA()
    
  )
  
}


plot.umap_SKA<-function(umres, labels_col){
  
  ### Function to plot UMAP results
  # Input : umres- umap results object & a dataframe with 2 columns ( BioSpecID for colour and for shape)
  
  x=cbind(umres$layout,labels_col)
  colnames(x)<-c("UMAP_1","UMAP_2","Label")
  
  print(
    ggplot(x,aes(x=as.numeric(UMAP_1),y=as.numeric(UMAP_2)))+
      geom_point(aes(colour=Label),size=5,alpha=0.8)+
      geom_text_repel(aes(label=Label),size=3)+
      scale_colour_manual(values = colkey_BioSpecID)+
      labs(x="UMAP 1",y="UMAP 2")+
      scale_shape_manual(values=c(16,17,15))+
      theme_SKA(base_size=30)+
      theme(legend.position = "none")
  )
  
}

