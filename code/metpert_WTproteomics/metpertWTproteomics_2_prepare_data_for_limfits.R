###################################################################
### EDP Script - 1A - Data preparation for statistical analysis ###
###################################################################

#################
### Set Paths ###
#################

# general
source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/layout_conversion.R"))

# proteomics specific functions 

source(paste0(code_dir,"/metpert_WTproteomics/metpertWTproteomics_0_libraries_functions.R"))


plot_dir <- paste0(metpert_WTproteomics_dir,"/output/plots/forStat")
qc_plot_dir <- paste0(plot_dir,"/qc_plots")

dir.create(qc_plot_dir,recursive = T)
dir.create(paste0(qc_plot_dir,"/alldata_stats"))

intm_data_dir <- paste0(paste0(metpert_WTproteomics_dir,"/output/intermediate_datafiles"))
dir.create(intm_data_dir,recursive = T)

output_tables_dir <- paste0(metpert_WTproteomics_dir,"/output/forStat") 
dir.create(output_tables_dir,recursive = T)

########################################
### All Element control samples only ###
########################################

AE_QCplot_forStat_dir = paste0(qc_plot_dir,"/AE only")
  
dir.create(AE_QCplot_forStat_dir,recursive = T)

###################
### Input Data ###
###################

######################################################################################
######################################################################################
########################## All Element control samples only ##########################
######################################################################################
######################################################################################


## Read in Experimental Design for AE samples

Exp_Design<-read.csv(paste0(intm_data_dir,"/EDP_ExperimentDesign_AEonly.csv"),stringsAsFactors = F)

### Read in Sample data and filter for AE samples only ###

AE_data <- read.csv(paste0(intm_data_dir,"/AllData_prepared_noPrecFilter.csv"),stringsAsFactors = F)%>%
           filter(grepl("AllEle",BioSpecID))

###########################################################################
### Filter out Precursors that have InterBatch CVs > 30 % in AE samples ###
###########################################################################

## Compute CVs Intra Batch and Inter Batch CVs

AE_CVs_preBC <- AE_data%>%
                filter(BioSpecID == "AllEle 1")%>% # Dont use AllEle prep QC samples for CV determination
                mutate(Lo_Prot_PrecID = paste(Layout, Protein.Ids,Precursor.Id,sep="_"),
                       Prot_PrecID=paste(Protein.Ids,Precursor.Id,sep="_"))%>%
                group_by(Lo_Prot_PrecID)%>%
                mutate(CV.AE.Intra.Batch = sd(Precursor.Normalised,na.rm = T)/mean(Precursor.Normalised,na.rm = T))%>%
                ungroup()%>%
                group_by(Prot_PrecID)%>%
                mutate(CV.AE.Inter.Batch = sd(Precursor.Normalised,na.rm = T)/mean(Precursor.Normalised,na.rm = T))%>%
                mutate(QCtype="Pre Batch Correction")

PrecIDList <- unique(filter(AE_CVs_preBC,CV.AE.Intra.Batch <= 30)$Prot_PrecID)

AE_data <- AE_data%>%
            mutate(Prot_PrecID = paste(Protein.Ids,Precursor.Id,sep="_"))%>%
            filter(Prot_PrecID %in% PrecIDList)

#####################################################
### Compute Protein Group quantities ### Raw Data ###
#####################################################

PG_quants<-get_PGQ_df(AE_data)

## Merge PQ_quants with annotation

PG_quants<-merge(PG_quants,
                 unique(AE_data[,c("File.Name","Layout","Element","Element.Concentration",
                                   "Digestion_Batch","TimePoint",
                                   "LO_TP","BioSpecID","Sample.Type","Date","Time")]),
                 by="File.Name")

## Setup annotations for proBatch
setup_proBatch_annos(PG_quants,onlyAEdat = T)
essential_columns <- c(PGcol, PGintcol, fname_col)

### Prepare data for use with DEP ### Raw ###

DEP_PGquants <- reshape2::dcast ( PG_quants,
                                  Protein.Ids~File.Name,
                                  value.var = "Protein.Quantity")%>%
  mutate ( name = Protein.Ids,
           ID = Protein.Ids)

### Create Summarized Experiment Object ##

## Note : this will log2 transform all values

SumExp <- DEP::make_se(DEP_PGquants,seq(2,ncol(DEP_PGquants)-2),Exp_Design)  

############################################
### Batch Correct at the precursor level ###
############################################

AE_data_MedCorr <- AEmed_correction_Precursor(AE_data)

### Convert batch corrected data into protein quantities ###

PG_quants_BC <- get_PGQ_df(AE_data_MedCorr)

## Merge PQ_quants with annotation

PG_quants_BC <- merge(PG_quants_BC,
                    unique(AE_data[,c("File.Name","Layout","Element","Element.Concentration",
                                      "Digestion_Batch","TimePoint",
                                      "LO_TP","BioSpecID","Sample.Type","Date","Time")]),
                    by="File.Name")


DEP_PGquants_BC <- reshape2::dcast ( PG_quants_BC,
                                     Protein.Ids~File.Name,
                                     value.var = "Protein.Quantity")%>%
  mutate ( name = Protein.Ids,
           ID = Protein.Ids)


SumExp_BatchCorr <- DEP::make_se(DEP_PGquants_BC,seq(2,ncol(DEP_PGquants_BC)-2),Exp_Design)  

######################
### Median Scaling ###
######################

PG_quants_BCMS <- PG_quants_BC%>%
                  group_by(File.Name)%>%
                  mutate( MPQ_perFile= median(Protein.Quantity,na.rm=T))%>%
                  ungroup()%>%
                  # Calculate mean of median protein intensities across all files
                  mutate( Median_MPQ = median(MPQ_perFile,na.rm=T),
                          # Scale each protein quantity such that all files will have same Median 
                          Median_Scaled_Protein_Quantity = (Protein.Quantity/MPQ_perFile)*Median_MPQ)
                



DEP_PGquants_BCMS <- reshape2::dcast ( PG_quants_BCMS,
                                       Protein.Ids~File.Name,
                                       value.var = "Median_Scaled_Protein_Quantity")%>%
  mutate ( name = Protein.Ids,
           ID = Protein.Ids)


SumExp_BatchCorrMedSc <- DEP::make_se(DEP_PGquants_BCMS,seq(2,ncol(DEP_PGquants_BCMS)-2),Exp_Design)  

##################################
### Collect all SumExp objects ###
##################################

SEos<-list(
  SumExp,
  SumExp_BatchCorr,
  SumExp_BatchCorrMedSc
)

names(SEos) <- c("Log2 Raw",
                 "BatchCorrected",
                 "BatchCorrected_medianScaled"
                 )

# Get DFs out of se object ( these will have log2 scale values )

#SEdfs<-lapply(SEos, SumExpObj_to_df)

##################################################################
### Plot effect of filtering, AE correction and NA replacement ###
##################################################################


## Number of proteins identified vs num of samples they're identified in

# Raw
idp1 <- plot_frequency(SumExp)+
        labs(title="Protein identification overlap\nRaw data")

## Plot coverage of protein identification between samples

# Raw
covp1 <- plot_coverage(SumExp)+
         labs(title="Protein coverage \n Raw data") +
         ylim(0,2500)


## Number of proteins identified per sample

# Raw
np1 <- plot_numbers(SumExp)+
       scale_fill_manual(values=c(viridis(n=10),c(viridis(n=10))))+
       coord_flip()+
       labs(title="Proteins per sample \n Raw data")+
       theme(legend.position="none",legend.key.size = unit(3,"mm"))

## Dependence of number of missing values on protein intensity

# Raw
pd1 <- plot_detect(SumExp)

## Sample wise protein quantity boxplots

swbp1 <- plot_normalization(SumExp)+  
          scale_fill_manual(values=c(viridis(n=10),viridis(n=10)) )+
          theme(legend.position = "none")
swbp2 <- plot_normalization(SumExp_BatchCorr)+
          scale_fill_manual(values=c(viridis(n=10),viridis(n=10)) )+
          theme(legend.position = "none")
swbp3 <- plot_normalization(SumExp_BatchCorrMedSc)+
          scale_fill_manual(values=c(viridis(n=10),viridis(n=10)) )+
          theme(legend.position = "none")


pdf(paste0(AE_QCplot_forStat_dir,"/effect_of_filtering_batchcorrection_medianscaling_on_data_forlmfit.pdf"),width=20,height=20)
idp1
covp1
np1
grid.arrange(pd1)
dev.off()

pdf(paste0(AE_QCplot_forStat_dir,"/median_protein_quantities_across_data_processing_stages_forlmfit.pdf"),width=40,height=20)

swbp1+swbp2+swbp3+plot_layout(ncol=3)

dev.off()

##############################
### Make proBatch QC plots ###
##############################

### Set up Pro batch annotations

#sel_annos <- c("LO_TP","Sample.Type","Digestion_Batch") ## annotations to plot by

#for( i in 1:length(SEdfs)){
  
  
  
#  make_proBatchQCplots(SEdfs[[i]],
#                       sel_anno_hclust=sel_annos,
#                       onlyAEdat = T,
#                       dsname=names(SEdfs)[[i]],
#                       AE_or_samp = "AllEle")
  
#}

###################################
### Get DE data between batches ###
###################################

## Carry out DE analysis, save plots as pdfs and store results in a list of dfs 

#DE_res_dfs<-list()

#for( i in 1:2){ # Dont do DE analysis on data with NAs replaced
  # DF 3 throws errors with proBatch because of FDRTOOLS new version or something 
#  DE_res_dfs[[i]] <- DE_analysis_and_plots(df = SEdfs[[i]],
#                                           dsname = names(SEdfs)[[i]],
#                                           expdes =  Exp_Design,
#                                           onlyAEdat = T)
#  names(DE_res_dfs)[[i]] <- names(SEdfs)[[i]]
#}


###########################################
### Write results of DE analysis to csv ###
###########################################

# Keep only unique columns specified results of comparisons
#fin_deres_df <- unique(DE_res_dfs[[2]][,c("ID",colnames(DE_res_dfs[[2]])[grep("_vs_",colnames(DE_res_dfs[[2]]))] ) ])%>%
#  unique()%>%
#  reshape2::melt(id.vars = c("ID"))%>%
#  separate(variable, into = c("Condition_1",NA,"Condition_2","Stat.Type"),sep="_")%>%
#  reshape2::dcast(ID+Condition_1+Condition_2~Stat.Type,value.var = "value")%>%
#  na.omit()%>%
#  unique()

#colnames(fin_deres_df)[[1]] <- "Genes"

#Write results of DE analysis on data with NAs to csv
#write.csv(fin_deres_df,"EDP Differential Expression Results_forStat_AllEleOnly.csv",row.names = F)

# Prepare data with NAs for writing to CSV 

Fin_wNA_df <- merge(PG_quants_BCMS,sample_anno[,c(1,which(!colnames(sample_anno) %in% colnames(PG_quants_BCMS)))],by="File.Name")%>%
              mutate(  BioSpecID=trimws(BioSpecID),
                      `Log2(Protein.Quantity)` = log2(Protein.Quantity))%>%
              filter(BioSpecID != "AllEleprepQC 1")

write.csv(Fin_wNA_df,paste0(output_tables_dir,"/protein_quantities_w_NAs_AllEleOnly.csv"),row.names = F)


##########################################
### Plot number of proteins per Layout ###
##########################################

numprot <- Fin_wNA_df %>%
            na.omit()%>%
            dplyr::select(Protein.Ids, Digestion_Batch, Layout,File.Name)%>%
            unique()%>%
            group_by(File.Name)%>%
            mutate(Num_Prot = length(Protein.Ids))%>%
            dplyr::select(Layout,Num_Prot,Digestion_Batch,File.Name)%>%
            unique()

pdf(paste0(AE_QCplot_forStat_dir,"/number_of_proteins_per_layout_AE_control_samples.pdf"),width = 5,height = 4 )
ggplot(numprot,
       aes(x=factor(Layout),y=Num_Prot,
           colour=Digestion_Batch,
           group=Layout))+
  geom_point(size=3,alpha=0.6,position=position_jitterdodge())+
  scale_color_manual(values= c("darkred","darkblue"))+
  geom_boxplot(alpha=0,colour="black")+
  theme_metallica()+
  theme(legend.position = "bottom",
        panel.grid.major = element_line(colour="black",size = 0.1),
        panel.border = element_rect(size=0.15),
        axis.ticks = element_line(size=0.75))+
  labs(x="Layout", y="Number of proteins quantified")
dev.off()


## Plot protein quantities after median scaling 
pdf(paste0(AE_QCplot_forStat_dir,"/median_protein_quantity_after_diann_maxlfq_AEcorr_medianscale.pdf"),width=20,height=10,family = "Times")

ggplot(PG_quants_BCMS,
       aes( x = File.Name,
            y = log2(Median_Scaled_Protein_Quantity),
            group = File.Name))+
  geom_boxplot(alpha=0.3)+
  theme_metallica()+
  theme(axis.text.x = element_text(angle=90,size=5))
dev.off()

PG_quants_BCMS_mean<-PG_quants_BCMS%>%
  group_by(File.Name)%>%
  summarize(MeanPQ = mean(Median_Scaled_Protein_Quantity,na.rm = T))

pdf(paste0(AE_QCplot_forStat_dir,"/mean_protein_quantity_after_diann_maxlfq_AEcorr_medianscale.pdf"),width=20,height=10,family = "Times")

ggplot(PG_quants_BCMS_mean,
       aes( x = File.Name,
            y = MeanPQ))+
  geom_bar(stat="identity",fill="navyblue",colour="white")+
  theme_metallica()+
  theme(axis.text.x = element_text(angle=90,size=5))

dev.off()


PG_quants_BCMS_TPC <- PG_quants_BCMS%>%
                      group_by(File.Name)%>%
                      summarize(Total_Protein_Quantity = sum(Median_Scaled_Protein_Quantity,na.rm = T))

pdf(paste0(AE_QCplot_forStat_dir,"/total_protein_quantity_after_diann_maxlfq_AEcorr_medianscale.pdf"),width=20,height=10,family = "Times")

ggplot(PG_quants_BCMS_TPC,
       aes( x = File.Name,
            y = log2(Total_Protein_Quantity)))+
  geom_bar(stat="identity",fill="navyblue",colour="white")+
  theme_metallica()+
  theme(axis.text.x = element_text(angle=90,size=5))
dev.off()


###################################################################################################
###################################################################################################
###################################################################################################
########################################### All samples ########################################### 
###################################################################################################
###################################################################################################
###################################################################################################

AS_QCplot_forStat_dir = paste0(qc_plot_dir,"/forStat/All Samples")
dir.create(AS_QCplot_forStat_dir,recursive = T)

###################
### Input Data ###
###################

#######################################
### Read in Experiment Design files ###
#######################################

Exp_Design <- read.csv(paste0(intm_data_dir,"/metpertWTproteomics_ExperimentDesign.csv"),stringsAsFactors = F)

Exp_Design_AEonly <- read.csv(paste0(intm_data_dir,"/metpertWTproteomics_ExperimentDesign_AEonly.csv"),stringsAsFactors = F)

###############################
### Read in All Sample data ###
###############################

Sample_data <- read.csv(paste0(intm_data_dir,"/AllData_prepared_noPrecFilter.csv"),stringsAsFactors = F)

Sample_data <- Sample_data%>% 
               mutate(Prot_PrecID = paste(Protein.Ids,Precursor.Id,sep="_"))%>%
               filter(Prot_PrecID %in% PrecIDList)

########################################
### Compute Protein Group quantities ###
########################################

PG_quants <- get_PGQ_df(Sample_data)

## Merge PG_qunats with annotation

PG_quants <- merge(PG_quants,
                 unique(Sample_data[,c("File.Name","Layout","Element",
                                       "Element.Concentration",
                                       "Digestion_Batch",  "TimePoint",
                                       "LO_TP","BioSpecID","Sample.Type","Date","Time")]),
                 by="File.Name")
    

## Setup annotations for proBatch
setup_proBatch_annos(PG_quants,onlyAEdat = F)
essential_columns <- c(PGcol, PGintcol, fname_col)

### Prepare data for use with DEP ###

DEP_PGquants <- reshape2::dcast ( PG_quants,
                                  Protein.Ids~File.Name,
                                  value.var = "Protein.Quantity")%>%
  mutate ( name = Protein.Ids,
           ID = Protein.Ids)


### Create Summarized Experiment Object ##

## Note : this will log2 transform all values

SumExp <- make_se(DEP_PGquants,seq(2,ncol(DEP_PGquants)-2),Exp_Design )

### Batch Correct at the precursor level ###

Samdata_fMedCorr <- AEmed_correction_Precursor(Sample_data)

############################################################
### Convert batch corrected data into protein quantities ###
############################################################

PG_quants_BC <- get_PGQ_df(Samdata_fMedCorr)

## Merge PQ_quants with annotation

PG_quants_BC <- merge(PG_quants_BC,
                    unique(Sample_data[,c("File.Name","Layout","Element","Element.Concentration",
                                          "Digestion_Batch","TimePoint",
                                          "LO_TP","BioSpecID","Sample.Type","Date","Time")]),
                    by="File.Name")


DEP_PGquants_BC <- reshape2::dcast ( PG_quants_BC,
                                     Protein.Ids~File.Name,
                                     value.var = "Protein.Quantity")%>%
                    mutate ( name = Protein.Ids,
                             ID = Protein.Ids)


SumExp_BatchCorr <- make_se(DEP_PGquants_BC,seq(2,ncol(DEP_PGquants_BC)-2),Exp_Design)  

######################
### Median Scaling ###
######################

PG_quants_BCMS <- PG_quants_BC%>%
                  group_by(File.Name)%>%
                  mutate( MPQ_perFile= median(Protein.Quantity,na.rm=T))%>%
                  ungroup()%>%
                  # Calculate mean of median protein intensities across all files
                  mutate( Median_MPQ = median(MPQ_perFile,na.rm=T),
                          # Scale each protein quantity such that all files will have same Median 
                          Median_Scaled_Protein_Quantity = (Protein.Quantity/MPQ_perFile)*Median_MPQ)


DEP_PGquants_BCMS <- reshape2::dcast ( PG_quants_BCMS,
                                       Protein.Ids~File.Name,
                                       value.var = "Median_Scaled_Protein_Quantity")%>%
                      mutate ( name = Protein.Ids,
                               ID = Protein.Ids)


SumExp_BatchCorrMedSc <- make_se(DEP_PGquants_BCMS,seq(2,ncol(DEP_PGquants_BCMS)-2),Exp_Design)  

###################################
### Batch Correction with Limma ###
###################################


Sam_data_fLimma <- Sample_data%>%
                    dplyr::select(File.Name,Precursor.Id,Precursor.Normalised)%>%
                    spread(File.Name, Precursor.Normalised)
rownames(Sam_data_fLimma)<-Sam_data_fLimma$Precursor.Id

Sam_data_fLimma<-Sam_data_fLimma[,-1]

## Digestion Batch Anno

DigBatch_fLimma <- data.frame(File.Name = colnames(Sam_data_fLimma))

DigBatch_fLimma <- t(merge(DigBatch_fLimma,unique(Sample_data[,c("File.Name","Digestion_Batch")]),by = "File.Name")$Digestion_Batch)

## Layout Anno

LOBatch_fLimma <- data.frame(File.Name = colnames(Sam_data_fLimma))

LOBatch_fLimma <- t(merge(LOBatch_fLimma,unique(Sample_data[,c("File.Name","LO_TP")]),by = "File.Name")$LO_TP)


## Design matrix -- to keep differences between biospecIDs

Design_fLimma <- data.frame(File.Name = colnames(Sam_data_fLimma))

Design_fLimma <- merge(Design_fLimma,unique(Sample_data[,c("File.Name","BioSpecID")]),by = "File.Name")%>%
  reshape2::dcast(File.Name ~ BioSpecID)

Design_fLimma <- Design_fLimma[,-1]

Design_fLimma[!is.na(Design_fLimma)] <- 1
Design_fLimma[is.na(Design_fLimma)] <- 0

for(i in 1:ncol(Design_fLimma)){
  Design_fLimma[,i]<-as.numeric(Design_fLimma[,i])
}


Samdat_limmaBC <- limma::removeBatchEffect(log2(1+as.matrix(Sam_data_fLimma)),
                                           batch = DigBatch_fLimma,
                                           babatch2 = LOBatch_fLimma,
                                           design = Design_fLimma)

Samdat_limmaBC <- Samdat_limmaBC %>%
  reshape2::melt(value.name = "Precursor.Normalised")

colnames(Samdat_limmaBC)[1:2] <- c("Precursor.Id","File.Name")

Samdat_limmaBC<- merge(Samdat_limmaBC,Sample_data[,-which(colnames(Sample_data)=="Precursor.Normalised")],
                       by=c("Precursor.Id","File.Name"))%>%
  mutate(Precursor.Normalised = 2^Precursor.Normalised)%>%
  filter(Precursor.Normalised > 0)


##################################################################
### Convert limma batch corrected data into protein quantities ###
##################################################################

PG_quants_limmaBC <- get_PGQ_df(Samdat_limmaBC)

## Merge PQ_quants with annotation

PG_quants_limmaBC <- merge(PG_quants_limmaBC,
                           unique(Sample_data[,c("File.Name","Layout","Element","Element.Concentration",
                                                 "Digestion_Batch","TimePoint",
                                                 "LO_TP","BioSpecID","Sample.Type","Date","Time")]),
                           by="File.Name")


DEP_PGquants_limmaBC <- reshape2::dcast ( PG_quants_limmaBC,
                                          Protein.Ids~File.Name,
                                          value.var = "Protein.Quantity")%>%
  mutate ( name = Protein.Ids,
           ID = Protein.Ids)


SumExp_limmaBatchCorr <- make_se(DEP_PGquants_limmaBC,seq(2,ncol(DEP_PGquants_limmaBC)-2),Exp_Design)  

######################################################
### Median Scale data after limma Batch Correction ###
######################################################

PG_quants_limmaBCMS <- PG_quants_limmaBC%>%
                        group_by(File.Name)%>%
                        mutate( MPQ_perFile= median(Protein.Quantity,na.rm=T))%>%
                        ungroup()%>%
                        # Calculate mean of median protein intensities across all files
                        mutate( Median_MPQ = median(MPQ_perFile,na.rm=T),
                                # Scale each protein quantity such that all files will have same Median 
                                Median_Scaled_Protein_Quantity = (Protein.Quantity/MPQ_perFile)*Median_MPQ)
                      



DEP_PGquants_limmaBCMS <- reshape2::dcast ( PG_quants_limmaBCMS,
                                            Protein.Ids~File.Name,
                                            value.var = "Median_Scaled_Protein_Quantity")%>%
                              mutate ( name = Protein.Ids,
                                       ID = Protein.Ids)


SumExp_limmaBatchCorrMedSc <- make_se(DEP_PGquants_limmaBCMS,seq(2,ncol(DEP_PGquants_limmaBCMS)-2),Exp_Design)  


##################################
### Collect all SumExp objects ### 
##################################

SEos<-list(
  SumExp,
  SumExp_BatchCorr,
  SumExp_BatchCorrMedSc,
  
  SumExp_limmaBatchCorr,
  SumExp_limmaBatchCorrMedSc
)

names(SEos) <- c("Log2 Raw",
                 "BatchCorrected",
                 "BatchCorrected_medianScaled",
                 "limmaBatchCorrected",
                 "limmaBatchCorrected_medianScaled"
                )

# Get DFs out of se object ( these will have log2 scale values )

#SEdfs<-lapply(SEos, SumExpObj_to_df)


#############################################################################
### For limma Corrected - collect SumExp objects and dfs with AllEle only ### 
#############################################################################

## Batch corrected only

PG_quants_limmaBC_AEonly <- PG_quants_limmaBC%>%
                             filter(grepl("AllEle",BioSpecID))


DEP_PGquants_limmaBC_AEonly <- reshape2::dcast ( PG_quants_limmaBC_AEonly,
                                                 Protein.Ids~File.Name,
                                                 value.var = "Protein.Quantity")%>%
                                  mutate ( name = Protein.Ids,
                                           ID = Protein.Ids)

Exp_Design_AEonly <- filter(Exp_Design_AEonly,label %in% colnames(DEP_PGquants_limmaBC_AEonly))


SumExp_limmaBatchCorr_AEonly <- make_se(DEP_PGquants_limmaBC_AEonly,seq(2,ncol(DEP_PGquants_limmaBC_AEonly)-2),Exp_Design_AEonly)  


# Batch corrected and Median scaled

PG_quants_limmaBCMS_AEonly <- PG_quants_limmaBCMS%>%
                              filter(grepl("AllEle",BioSpecID))


DEP_PGquants_limmaBCMS_AEonly <- reshape2::dcast ( PG_quants_limmaBCMS_AEonly,
                                                   Protein.Ids~File.Name,
                                                   value.var = "Median_Scaled_Protein_Quantity")%>%
                                    mutate ( name = Protein.Ids,
                                             ID = Protein.Ids)


SumExp_limmaBatchCorrMedSc_AEonly <- make_se(DEP_PGquants_limmaBCMS_AEonly,seq(2,ncol(DEP_PGquants_limmaBCMS_AEonly)-2),Exp_Design_AEonly)  




SEos_limma_AEonly <- list(SumExp_limmaBatchCorr_AEonly,SumExp_limmaBatchCorrMedSc_AEonly)

names(SEos_limma_AEonly) <- c( "limmaBatchCorrected AEsamp only", "limmaBatchCorrected_medianScaled AEsamp only")

#SEdfs_limma_AEonly <- lapply(SEos_limma_AEonly, SumExpObj_to_df)

##################################################################
### Plot effect of filtering, AE correction and NA replacement ###
##################################################################


## Number of proteins identified vs num of samples they're identified in

# Raw
idp1 <- plot_frequency(SumExp)+
  labs(title="Protein identification overlap\nRaw data")

## Plot coverage of protein identification between samples

# Raw
covp1 <-plot_coverage(SumExp)+
  labs(title="Protein coverage \n Raw data") +
  ylim(0,2500)


## Number of proteins identified per sample

# Raw
np1 <- plot_numbers(SumExp)+
  scale_fill_manual(values=c(viridis(n=111),c(viridis(n=111))))+
  coord_flip()+
  labs(title="Proteins per sample \n Raw data")+
  theme(legend.position="none",legend.key.size = unit(3,"mm"))

## Dependence of number of missing values on protein intensity

# Raw
pd1 <- plot_detect(SumExp)

## Sample wise protein quantity boxplots

swbp1 <- plot_normalization(SumExp)+  
  scale_fill_manual(values=c(viridis(n=111),c(viridis(n=111))))+
  theme(legend.position = "none")
swbp2 <- plot_normalization(SumExp_BatchCorr)+
  scale_fill_manual(values=c(viridis(n=111),c(viridis(n=111))))+
  theme(legend.position = "none")
swbp3 <- plot_normalization(SumExp_BatchCorrMedSc)+
  scale_fill_manual(values=c(viridis(n=111),c(viridis(n=111))))+
  theme(legend.position = "none")

setwd(AS_QCplot_forStat_dir)

pdf("Effect of Filtering_BatchCorrection_MedianScaling on data.pdf",width=20,height=20)

idp1

covp1

np1

grid.arrange(pd1)
dev.off()

pdf("Median protein Quantities across data processing stages.pdf",width=40,height=20)

swbp1+swbp2+swbp3+plot_layout(ncol=3)

dev.off()

#########################################
### Plots to check intensity profiles ###
#########################################

#pdf("Protein Intensity profiles facet.pdf",width=20,height=15)
#for(i in 1:length(SEdfs)){
#  plot_protint_dist_facet(SEdfs[[i]], names(SEdfs)[[i]])
  
#}

#pdf("Protein Intensity profiles allsamp.pdf",width=20,height=15)
#for(i in 1:length(SEdfs)){
#  plot_protint_dist(SEdfs[[i]], names(SEdfs)[[i]])
#}

##############################
### Make proBatch QC plots ###
##############################

### Set up Pro batch annotations

#sel_annos <- c("LO_TP","Sample.Type","Digestion_Batch") ## annotations to plot by


#for( i in 1:length(SEdfs)){
  
  
 # make_proBatchQCplots(SEdfs[[i]],
  #                     sel_anno_hclust=sel_annos,
   #                    onlyAEdat = F,
    #                   dsname=names(SEdfs)[[i]],
     #                  AE_or_samp = "Samples",
      #                 sanno = sample_anno )
  
#}

## Plot AE samples only from the limma Batch correction

#sel_annos_AEonly <- c("LO_TP","Digestion_Batch") ## annotations to plot by

#sanno_AEonly <- filter(sample_anno,grepl("AllEle",BioSpecID))

#for( i in 1:length(SEdfs_limma_AEonly)){
  
  
#  make_proBatchQCplots(SEdfs_limma_AEonly[[i]],
#                       sel_anno_hclust=sel_annos,
#                       onlyAEdat = T,
#                       dsname=names(SEdfs_limma_AEonly)[[i]],
#                       AE_or_samp = "AllEle",
#                       sanno = sanno_AEonly
#  )
  
#}


## Use unimputed data after batch correction and median scaling for this 

#df4rc <- SEdfs[[3]] %>%
 # data.frame()%>%
#  mutate(Protein.Ids=rownames(SEdfs[[3]]))%>%
#  reshape2::melt(id.vars="Protein.Ids",variable.name= "File.Name",
#                 value.name="Log2(Protein.Quantity)")

#df4rc <- merge(df4rc,sample_anno,by="File.Name")%>%
#  mutate(BioSpecID = trimws(BioSpecID))%>%
#  filter(BioSpecID != "AllEleprepQC 1")

#bspcs <- unique(df4rc$BioSpecID)

#cor_check <- vector()

#for( i in 1:length(bspcs)){
  
 # df = filter(df4rc,BioSpecID == bspcs[[i]])%>%
#    reshape2::dcast(Protein.Ids~File.Name,value.var="Log2(Protein.Quantity)")
#  rownames(df) = df$Protein.Ids
  
  
#  if(ncol(df)>2){
    
#    df = df[,-1]
#    cordf=stats::cor(df,use = "pairwise.complete",method="spearman")
    
#    cor_check = rbind(cor_check,c(bspcs[[i]],min(cordf)))
 # }
#}
#colnames(cor_check)<-c("BioSpecID","Replicate_Correlation")

#write.csv(cordf,"Replicate correlations Spearman.csv")

######################################################
######################################################
### AE median Batch Corrected DFs to write to disk ###
######################################################
######################################################

Fin_wNA_df <- merge(PG_quants_BCMS,sample_anno[,c(1,which(!colnames(sample_anno)%in% colnames(PG_quants_BCMS)))],
                    by="File.Name")%>%
  mutate(BioSpecID=trimws(BioSpecID),
         `Log2(Protein.Quantity)` = log2(Protein.Quantity))%>%
  filter(BioSpecID != "AllEleprepQC 1")

write.csv(Fin_wNA_df,paste0(output_tables_dir,"/metperWTproteomics_protein_quantities_w_NAs_alldata_unfiltered_for_stat.csv"),row.names = F)


# Calculate per Gene mean, sd, median, mad for AllEle replicates across dataset and write to CSV

AERep_Stats <- Fin_wNA_df%>%
                filter(BioSpecID =="AllEle 1")%>%
                group_by(Protein.Ids)%>%
                mutate(Num.Non.NAs = sum(!is.na(`Log2(Protein.Quantity)`)))%>%
                ungroup()

## num proteins AE samples

AERep_numprot <- AERep_Stats %>%
                 dplyr::select(File.Name,Protein.Ids,Protein.Quantity)%>%
                 unique()%>%
                 na.omit()%>%
                 group_by(File.Name)%>%
                 summarise(num_prot_quantified = n())%>%
                 ungroup()

mean(AERep_numprot$num_prot_quantified)              

## Plot distribution of fraction of controls that are NAs

ggplot(unique(AERep_Stats[,c("Protein.Ids","Num.Non.NAs")]),
       aes(x=Num.Non.NAs))+
  geom_histogram(bins=25,colour="darkred",
                 fill="darkred")+
  theme_metallica()

AERep_Stats <- AERep_Stats%>%
               filter(Num.Non.NAs > 10)%>%
               group_by(Protein.Ids)%>%
               summarize(SD.Protein.Quant = sd(2^(`Log2(Protein.Quantity)`),na.rm=T),
                          Mean.Protein.Quant = mean(2^(`Log2(Protein.Quantity)`),na.rm=T)
               )%>%
               ungroup()%>%
               mutate(Cov = SD.Protein.Quant/Mean.Protein.Quant)%>%
               as.data.frame()

rownames(AERep_Stats) <- AERep_Stats$Protein.Ids

write.csv(AERep_Stats,paste0(output_tables_dir,"/AllEle_replicates_protein_quantity_stats.csv"),row.names=F)


# Calculate per Gene mean, sd, median, mad for all samples AllEle replicates across dataset and write to CSV

AllSample_Stats <- Fin_wNA_df%>%
                    filter(!grepl("AllEle",BioSpecID))%>%
                    group_by(Protein.Ids)%>%
                    mutate(Num.Non.NAs = sum(!is.na(`Log2(Protein.Quantity)`)))%>%
                    ungroup()


AllSample_numprot <- AllSample_Stats%>%
                      dplyr::select(File.Name,Protein.Ids,Protein.Quantity)%>%
                      unique()%>%
                      na.omit()%>%
                      group_by(File.Name)%>%
                      summarise(num_prot_quantified = n())%>%
                      ungroup()

mean(AllSample_numprot$num_prot_quantified)   


AllSample_signal_CoV <- AllSample_Stats%>%
  group_by(Protein.Ids)%>%
  summarise(signal_cov = sd(Protein.Quantity, na.rm=T)/mean(Protein.Quantity,na.rm=T))%>%
  ungroup()

mean(AllSample_signal_CoV$signal_cov,na.rm=T)

## Plot distribution of fraction of controls that are NAs

ggplot(unique(AllSample_Stats[,c("Protein.Ids","Num.Non.NAs")]),
       aes(x=Num.Non.NAs))+
  geom_histogram(bins=25,colour="darkred",
                 fill="darkred")+
  theme_metallica()

AllSample_Stats <- AllSample_Stats%>%
                    filter(Num.Non.NAs >10)%>%
                    group_by(Protein.Ids,BioSpecID)%>%
                    summarize(SD.Protein.Quant = sd(2^(`Log2(Protein.Quantity)`),na.rm=T),
                              Mean.Protein.Quant = mean(2^(`Log2(Protein.Quantity)`),na.rm=T)
                    )%>%
                    ungroup()%>%
                    mutate(Cov = SD.Protein.Quant/Mean.Protein.Quant)%>%
                    as.data.frame()

write.csv(AllSample_Stats,paste0(output_tables_dir,"/AllSample_protein_quantity_replicate_stats.csv"),row.names=F)

metalwise_replicate_stats <- Fin_wNA_df%>%
                    filter(!grepl("AllEle",BioSpecID))%>%
                    group_by(Protein.Ids,BioSpecID,Element)%>%
                    summarize(SD.Protein.Quant = sd(2^(`Log2(Protein.Quantity)`),na.rm=T),
                              Mean.Protein.Quant = mean(2^(`Log2(Protein.Quantity)`),na.rm=T)
                    )%>%
                    mutate(Cov = SD.Protein.Quant/Mean.Protein.Quant)%>%
                    ungroup()%>%
                    group_by(Protein.Ids, Element)%>%
                    summarize(Mean.COVreps = mean(Cov,na.rm=T))%>%
                    ungroup()%>%
                    as.data.frame()

rownames(metalwise_replicate_stats) <- paste(metalwise_replicate_stats$Protein.Ids,metalwise_replicate_stats$Element)

write.csv(metalwise_replicate_stats,paste0(output_tables_dir,"/metalwise_protein_quantity_replicate_stats.csv"),row.names=F)

num_proteins <- metalwise_replicate_stats%>%
                na.omit()%>%
                group_by(Element)%>%
                summarize(Num_Unique_ProtIDs = length(unique(Protein.Ids)))%>%
                ungroup()
                  
write.csv(num_proteins,paste0(output_tables_dir,"/metalwise_num_proteins_quantified_for_lmfit_stats.csv"),row.names=F)

### Number of proteins ###

## layout-wise
numprot <- Fin_wNA_df %>%
          na.omit()%>%
          dplyr::select(Protein.Ids,Layout,Digestion_Batch,BioSpecID)%>%
          unique()%>%
          group_by(Layout,BioSpecID,Digestion_Batch)%>%
          summarize(Num_Prot = length(Protein.Ids))

pdf(paste0(qc_plot_dir,"/number_proteins_quantified_per_layout_data_for_lmfit_stats.pdf"),width = 5,height=4)
ggplot(numprot,
       aes(x=factor(Layout),y=Num_Prot,
           colour=Digestion_Batch,
           group=Layout))+
  geom_point(size=3,alpha=0.6,position=position_jitterdodge())+
  scale_color_manual(values= c("darkred","darkblue"))+
  geom_boxplot(alpha=0,colour="black")+
  theme_metallica()+
  theme(legend.position = "bottom",
        panel.grid.major = element_line(colour="black",size = 0.1),
        panel.border = element_rect(size=0.15),
        axis.ticks = element_line(size=0.75))+
labs(x="Layout", y="Number of proteins quantified")
dev.off()

# BioSpecID wise

numprot_BioSpecID <- Fin_wNA_df %>%
                    na.omit()%>%
                    dplyr::select(Protein.Ids,BioSpecID)%>%
                    unique()%>%
                    group_by(BioSpecID)%>%
                    summarize(Num_Prot = length(Protein.Ids))

pdf(paste0(qc_plot_dir,"/number_proteins_quantified_per_condition_data_for_lmfit_stats.pdfylim.pdf"),width = 12,height=5)
ggplot(numprot_BioSpecID,
       aes(x=BioSpecID,y=Num_Prot,
           colour=BioSpecID,
       ))+
  geom_point(size=3,alpha=1.0)+
  scale_color_manual(values=colkey_BioSpecID)+
  theme_metallica()+
  theme(
        panel.grid.major = element_line(colour="lightgray",size = 0.1),
        panel.border = element_rect(size=0.15),
        legend.position = "none",
        axis.text.x = element_text(angle=90,size=9.5,hjust = 0.9),
        axis.ticks = element_line(size=0.75))+
  labs(x="", y="Number of proteins quantified")+
  ylim(1300,2500)

dev.off()

