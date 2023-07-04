#######################################################################
### EDP Script - 1B - Data preparation for complete matrix analysis ###
#######################################################################

nathr =0.85 ## Threshold for fraction of samples that have to be not NA to retain a protein

#################
### Set Paths ###
#################

projdir <- "/camp/lab/ralserm/working/Simran Aulakh/Ionome Project/Element Depletion/EleDepChel January 2020"
proteomics_dir = "/camp/lab/ralserm/working/Simran Aulakh/Ionome Project/Element Depletion/EleDepChel January 2020/Proteomics"

protRscript_dir <- paste0(proteomics_dir,"/RScripts")

source(paste0(projdir,"/Common Functions/Initialize_Paths.R"))
source(paste0(projdir,"/Common Functions/Graphic_Parameters.R"))
source(paste0(protRscript_dir,'/EDP_Analysis_Functions.R'))

datdir = paste0(proteomics_dir,"/DIANNout")
plotdir=paste0(proteomics_dir,"/Plots")
QCplotdir=paste0(plotdir,"/QC_plots")
Int_datdir =paste0(paste0(proteomics_dir,"/IntermediateDFs"))

setwd(Int_datdir)
###############################################################################
### Create precursor  filter based on CV of precursors in all ele controls ###
##############################################################################

AE_data<-read.csv("AllData_prepared_noPrecFilter.csv",stringsAsFactors = F)%>%
  filter(grepl("AllEle",Element))

## Compute CVs Intra Batch and Inter Batch CVs

AE_CVs_preBC <- AE_data%>%
  mutate(Lo_Prot_PrecID = paste(Layout, Protein.Ids,Precursor.Id,sep="_"),
         Prot_PrecID=paste(Protein.Ids,Precursor.Id,sep="_"))%>%
  group_by(Lo_Prot_PrecID)%>%
  mutate(CV.AE.Intra.Batch = cv(Precursor.Normalised,na.rm = T))%>%
  ungroup()%>%
  group_by(Prot_PrecID)%>%
  mutate(CV.AE.Inter.Batch = cv(Precursor.Normalised, na.rm = T))%>%
  mutate(QCtype="Pre Batch Correction")

PrecIDList<-unique(filter(AE_CVs_preBC,CV.AE.Intra.Batch <= 30)$Prot_PrecID)
length(PrecIDList)

################################################################################################
### Create protein completeness filter based on detection of proteins in 75 % of all samples ###
################################################################################################

Sample_data<-read.csv("AllData_prepared_noPrecFilter.csv",stringsAsFactors = F)


df_forProtFilt <- Sample_data%>%
  filter(!grepl("AllEle",BioSpecID))%>%
  mutate(Prot_PrecID = paste(Protein.Ids,Precursor.Id,sep="_"))%>%
  filter(Prot_PrecID %in% PrecIDList)


PG_quants_4protfilt<-get_PGQ_df(df_forProtFilt)

## Merge PQ_quants with annotation

PG_quants_4protfilt<-merge(PG_quants_4protfilt,
                 unique(df_forProtFilt[,c("File.Name","Layout","Element","Element.Concentration",
                                   "Digestion_Batch","TimePoint",
                                   "LO_TP","BioSpecID","Sample.Type","Date","Time")]),
                 by="File.Name")

protcompletness_filter <- PG_quants_4protfilt%>%
                        filter(!grepl("AllEle",BioSpecID ))%>%
                        dplyr::select(Protein.Ids,File.Name,Protein.Quantity)%>%
                        reshape2::dcast(Protein.Ids~File.Name)%>%
                        reshape2::melt()%>%
                        group_by(Protein.Ids)%>%
                        summarize(Fraction_not_NA = 1- (sum(is.na(value))/length(value)) )

ggplot(protcompletness_filter,
       aes(x= Fraction_not_NA))+
  geom_histogram(bins=75)

Prot2keep <- unique(filter(protcompletness_filter, Fraction_not_NA > nathr )$Protein.Ids)

length(Prot2keep)

######################################################################################
######################################################################################
########################## All Element control samples only ##########################
######################################################################################
######################################################################################

AE_QCplot_forFullMatrix_dir = paste0(QCplotdir,"/forFullMatrix/AE only")

dir.create(AE_QCplot_forFullMatrix_dir,recursive = T)

##################
### Input Data ###
##################

setwd(Int_datdir)

Exp_Design<-read.csv("EDP_ExperimentDesign_AEonly.csv",stringsAsFactors = F)

###########################################################################
### Filter out Precursors that have InterBatch CVs > 30 % in AE samples  & Proteins that were not detected in 80 % of all samples ###
###########################################################################

AE_data <- AE_data%>%
  mutate(Prot_PrecID = paste(Protein.Ids,Precursor.Id,sep="_"))%>%
  filter(Prot_PrecID %in% PrecIDList)%>%
  filter(Protein.Ids %in% Prot2keep)

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

## Remove samples with < 1000 proteins quantified 

numprot_filter <- PG_quants %>%
  dplyr::select(File.Name,Protein.Ids,Protein.Quantity)%>%
  na.omit()%>%
  group_by(File.Name)%>%
  summarize(numprot_detected = length(Protein.Ids))%>%
  ungroup()%>%
  filter(numprot_detected >=1300)%>%
  pull(File.Name)

PG_quants <- filter(PG_quants, File.Name %in% numprot_filter)  

### Remove File.Names of that have < 1000 proteins quantified from Exp_Design ###

Exp_Design <- filter( Exp_Design, label %in% numprot_filter)

### Prepare data for use with DEP ### Raw ###

DEP_PGquants <- reshape2::dcast ( PG_quants,
                                  Protein.Ids~File.Name,
                                  value.var = "Protein.Quantity")%>%
  mutate ( name = Protein.Ids,
           ID = Protein.Ids)

### Create Summarized Experiment Object ##

## Note : this will log2 transform all values

SumExp <- make_se(DEP_PGquants,seq(2,ncol(DEP_PGquants)-2),Exp_Design )

############################################
### Batch Correct at the precursor level ###
############################################

AE_data_MedCorr <- AEmed_correction_Precursor(AE_data)


### Convert batch corrected data into protein quantities ###

PG_quants_BC<-get_PGQ_df(AE_data_MedCorr)

## Merge PQ_quants with annotation

PG_quants_BC<-merge(PG_quants_BC,
                    unique(AE_data[,c("File.Name","Layout","Element","Element.Concentration",
                                      "Digestion_Batch","TimePoint",
                                      "LO_TP","BioSpecID","Sample.Type","Date","Time")]),
                    by="File.Name")

PG_quants_BC <- filter(PG_quants_BC, File.Name %in% numprot_filter)      

PG_quants_BC <- filter(PG_quants_BC, Protein.Ids %in% Prot2keep)  

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

SumExp_BatchCorrMedSc <- make_se(DEP_PGquants_BCMS,seq(2,ncol(DEP_PGquants_BCMS)-2),
                                 Exp_Design)  

##############################
### Fill in missing values ###
##############################

for_NAfill <- PG_quants_BCMS %>%
  reshape2::dcast(Protein.Ids~File.Name, value.var = "Median_Scaled_Protein_Quantity")

# If missing in all replicates of a condition --> fill in with minimum value measured in dataset
# Else - replace with median of all measured replicates

NAfree.df <- fill_in_NAs(for_NAfill)

NAfree.df.fSE <- data.frame(NAfree.df)%>%
  mutate(ID = rownames(NAfree.df),
         name = rownames(NAfree.df))

SumExp_noNA <- DEP::make_se(NAfree.df.fSE,seq(1:(ncol(NAfree.df.fSE)-2)),Exp_Design)

###################################
### Batch Correction with Limma ###
###################################

AE_data_fLimma <- AE_data%>%
  dplyr::select(File.Name,Precursor.Id,Precursor.Normalised)%>%
  spread(File.Name, Precursor.Normalised)
rownames(AE_data_fLimma)<-AE_data_fLimma$Precursor.Id

AE_data_fLimma<-AE_data_fLimma[,-1]

## Digestion Batch Anno

DigBatch_fLimma <- data.frame(File.Name = colnames(AE_data_fLimma))

DigBatch_fLimma <- t(merge(DigBatch_fLimma,unique(AE_data[,c("File.Name","Digestion_Batch")]),
                           by = "File.Name")$Digestion_Batch)

## Layout Anno

LOBatch_fLimma <- data.frame(File.Name = colnames(AE_data_fLimma))

LOBatch_fLimma <- t(merge(LOBatch_fLimma,unique(AE_data[,c("File.Name","LO_TP")]),
                          by = "File.Name")$LO_TP)

## Design matrix -- to keep differences between BiospecIDs

Design_fLimma <- data.frame(File.Name = colnames(AE_data_fLimma))

Design_fLimma <- merge(Design_fLimma,unique(AE_data[,c("File.Name","BioSpecID")]),
                       by = "File.Name")%>%
  reshape2::dcast(File.Name ~ BioSpecID)

Design_fLimma <- Design_fLimma[,-1]

Design_fLimma[!is.na(Design_fLimma)] <- 1
Design_fLimma[is.na(Design_fLimma)] <- 0

for(i in 1:ncol(Design_fLimma)){
  Design_fLimma[,i]<-as.numeric(Design_fLimma[,i])
}

AEdat_limmaBC <- limma::removeBatchEffect(log2(1+as.matrix(AE_data_fLimma)),
                                          batch = DigBatch_fLimma,
                                          babatch2 = LOBatch_fLimma,
                                          design = Design_fLimma)

AEdat_limmaBC <- AEdat_limmaBC %>%
  reshape2::melt(value.name = "Precursor.Normalised")

colnames(AEdat_limmaBC)[1:2] <- c("Precursor.Id","File.Name")

AEdat_limmaBC<- merge(AEdat_limmaBC,AE_data[,-which(colnames(AE_data)=="Precursor.Normalised")],
                      by=c("Precursor.Id","File.Name"))%>%
  mutate(Precursor.Normalised = 2^Precursor.Normalised)%>%
  filter(Precursor.Normalised > 0)

##################################################################
### Convert limma batch corrected data into protein quantities ###
##################################################################

PG_quants_limmaBC <- get_PGQ_df(AEdat_limmaBC)

## Merge PQ_quants with annotation

PG_quants_limmaBC <- merge(PG_quants_limmaBC,
                           unique(AE_data[,c("File.Name","Layout","Element","Element.Concentration",
                                             "Digestion_Batch","TimePoint",
                                             "LO_TP","BioSpecID","Sample.Type","Date","Time")]),
                           by="File.Name")

PG_quants_limmaBC <- filter(PG_quants_limmaBC, File.Name %in% numprot_filter)  


DEP_PGquants_limmaBC <- reshape2::dcast ( PG_quants_limmaBC,
                                          Protein.Ids~File.Name,
                                          value.var = "Protein.Quantity")%>%
  mutate ( name = Protein.Ids,
           ID = Protein.Ids)

SumExp_LimmaBatchCorr <- make_se(DEP_PGquants_limmaBC,seq(2,ncol(DEP_PGquants_limmaBC)-2),
                                 Exp_Design)  

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

SumExp_limmaBatchCorrMedSc <- make_se(DEP_PGquants_limmaBCMS,seq(2,ncol(DEP_PGquants_limmaBCMS)-2),
                                      Exp_Design)  

############################################
### Fill in NAs for limma corrected data ###
############################################

for_NAfill_limma <- PG_quants_limmaBCMS %>%
  reshape2::dcast(Protein.Ids~File.Name, value.var = "Median_Scaled_Protein_Quantity")

# If missing in all replicates of a condition --> fill in with minimum value measured in dataset
# Else - replace with median of all measured replicates

NAfree.df_limma <- fill_in_NAs(for_NAfill_limma)

NAfree.df.fSE_limma <- data.frame(NAfree.df_limma)%>%
  mutate(ID = rownames(NAfree.df_limma),
         name = rownames(NAfree.df_limma))

SumExp_limma_noNA <- DEP::make_se(NAfree.df.fSE_limma,seq(1:(ncol(NAfree.df.fSE_limma)-2)),
                                  Exp_Design)  

##################################
### Collect all SumExp objects ###
##################################

SEos<-list(
  SumExp,
  SumExp_BatchCorr,
  SumExp_BatchCorrMedSc,
  
  SumExp_LimmaBatchCorr,
  SumExp_limmaBatchCorrMedSc,
  SumExp_noNA,
  SumExp_limma_noNA
)

names(SEos) <- c("Log2 Raw",
                 "BatchCorrected",
                 "BatchCorrected_medianScaled",
                 "limmaBatchCorrected",
                 "limmaBatchCorrected_medianScaled",
                 "BatchCorrected_medianScaled_noNA",
                 "limmaBatchCorrected_medianScaled_noNa")

# Get DFs out of se object ( these will have log2 scale values )

SEdfs<-lapply(SEos, SumExpObj_to_df)

##################################################################
### Plot effect of filtering, AE correction and NA replacement ###
##################################################################


## Number of proteins identified vs num of samples they're identified in

# Raw
idp1 <- plot_frequency(SumExp)+
  labs(title="Protein identification overlap\nRaw data Full Matrix")

## Plot coverage of protein identification between samples

# Raw
covp1 <-plot_coverage(SumExp)+
  labs(title="Protein coverage \n Raw data Full Matrix") +
  ylim(0,2500)


## Number of proteins identified per sample

# Raw
np1 <- plot_numbers(SumExp)+
  scale_fill_manual(values=c(viridis(n=10),c(viridis(n=10))))+
  coord_flip()+
  labs(title="Proteins per sample \n Raw data Full Matrix")+
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
swbp4 <- plot_normalization(SumExp_noNA)+
  scale_fill_manual(values=c(viridis(n=10),viridis(n=10)) )+
  theme(legend.position = "none")


setwd(AE_QCplot_forFullMatrix_dir)

pdf("Effect of Filtering_AEmedCorrection_NAhandling FullMatrix on data.pdf",
    width=20,height=20)

idp1

covp1

np1

grid.arrange(pd1)

dev.off()

pdf("Median protein Quantities across data processing FullMatrix stages.pdf",width=40,height=20)

swbp1+swbp2+swbp3+ swbp4+plot_layout(ncol=4)

dev.off()

##############################
### Make proBatch QC plots ###
##############################

### Set up Pro batch annotations

sel_annos <- c("LO_TP","Sample.Type","Digestion_Batch") ## annotations to plot by

for( i in 1:length(SEdfs)){
  
  
  
  make_proBatchQCplots(SEdfs[[i]],
                       sel_anno_hclust=sel_annos,
                       onlyAEdat = T,
                       dsname=names(SEdfs)[[i]],
                       AE_or_samp = "AllEle")
  
}

###################################
### Get DE data between batches ###
###################################

## Carry out DE analysis, save plots as pdfs and store results in a list of dfs 

DE_res_dfs<-list()

for( i in 1:5){ # Dont do DE analysis on data with NAs replaced
  
  DE_res_dfs[[i]] <- DE_analysis_and_plots(df = SEdfs[[i]],
                                           dsname = names(SEdfs)[[i]],
                                           expdes =  Exp_Design,
                                           onlyAEdat = T)
  names(DE_res_dfs)[[i]] <- names(SEdfs)[[i]]
}


###########################################
### Write results of DE analysis to csv ###
###########################################

setwd(FinalOutCSVdir)

# Prepare data with NAs for writing to CSV 

Fin_wNA_df <-  merge(PG_quants_BCMS,sample_anno[,c(1,which(!colnames(sample_anno) %in% colnames(PG_quants_BCMS)))],by="File.Name")%>%
  data.frame()%>%
  mutate(BioSpecID=trimws(BioSpecID),
         `Log2(Protein.Quantity)` = log2(Protein.Quantity))%>%
  filter(BioSpecID != "AllEleprepQC 1")

write.csv(Fin_wNA_df,"Protein Quantities w NAs_AllEleOnly FullMatrix.csv",row.names = F)

# limma BC
Fin_wNA_df_limma <-  merge(PG_quants_limmaBCMS,sample_anno[,c(1,which(!colnames(sample_anno) %in% colnames(PG_quants_limmaBCMS)))],by="File.Name")%>%
  data.frame()%>%
  mutate(BioSpecID=trimws(BioSpecID),
         `Log2(Protein.Quantity)` = log2(Protein.Quantity))%>%
  filter(BioSpecID != "AllEleprepQC 1")


write.csv(Fin_wNA_df_limma,"Protein Quantities w NAs_AllEleOnly limmaBC FullMatrix.csv",row.names = F)

# Prepare data without NAs for writing to CSV 
# AEmed batch corr
FinnoNA_df <- NAfree.df%>%
  data.frame()%>%
  mutate(Protein.Ids=rownames(NAfree.df))%>%
  reshape2::melt(id.vars="Protein.Ids",variable.name= "File.Name",
                 value.name="Protein.Quantity")%>%
  mutate(`Log2(Protein.Quantity)`=log2(Protein.Quantity))

FinnoNA_df <- merge(FinnoNA_df,sample_anno,by="File.Name")%>%
  filter(Element != "AllEleprepQC")

write.csv(FinnoNA_df,"Protein Quantities wImpValues_AllEleOnly FullMatrix.csv",row.names = F)

#limma BC
FinnoNA_df_limma <-NAfree.df_limma%>%
  data.frame()%>%
  mutate(Protein.Ids=rownames(NAfree.df))%>%
  reshape2::melt(id.vars="Protein.Ids",variable.name= "File.Name",
                 value.name="Protein.Quantity")%>%
  mutate(`Log2(Protein.Quantity)`=log2(Protein.Quantity))

FinnoNA_df_limma <- merge(FinnoNA_df_limma,sample_anno,by="File.Name")%>%
  filter(Element != "AllEleprepQC")

write.csv(FinnoNA_df,"Protein Quantities wImpValues_AllEleOnly limmaBC FullMatrix.csv",row.names = F)


##########################################
### Plot number of proteins per Layout ###
##########################################

numprot <- Fin_wNA_df %>%
  na.omit()%>%
  dplyr::select(Protein.Ids,Layout,File.Name)%>%
  unique()%>%
  group_by(File.Name)%>%
  mutate(Num_Prot = length(Protein.Ids))%>%
  dplyr::select(Layout,Num_Prot,File.Name)%>%
  unique()

ggplot(numprot,
       aes(x=Layout,y=Num_Prot,group=Layout))+
  geom_point(size=3)+
  geom_boxplot(alpha=0)

## Plot protein quantities after median scaling 

pdf("Median Protein Quantity after diann_maxlfq AEcorr medianscale FullMatrix.pdf",width=20,height=10,family = "Times")

ggplot(PG_quants_BCMS,
       aes( x = File.Name,
            y = log2(Median_Scaled_Protein_Quantity),
            group = File.Name))+
  geom_boxplot(alpha=0.3)+
  theme_SKA()+
  theme(axis.text.x = element_text(angle=90,size=5))
dev.off()

PG_quants_BCMS_mean<-PG_quants_BCMS%>%
  group_by(File.Name)%>%
  summarize(MeanPQ = mean(Median_Scaled_Protein_Quantity,na.rm = T))

pdf("Mean Protein Quantity after diann_maxlfq AEcorr medianscale FullMatrix.pdf",width=20,height=10,family = "Times")

ggplot(PG_quants_BCMS_mean,
       aes( x = File.Name,
            y = MeanPQ))+
  geom_bar(stat="identity",fill="navyblue",colour="white")+
  theme_SKA()+
  theme(axis.text.x = element_text(angle=90,size=5))

dev.off()


PG_quants_BCMS_TPC<-PG_quants_BCMS%>%
  group_by(File.Name)%>%
  summarize(Total_Protein_Quantity = sum(Median_Scaled_Protein_Quantity,na.rm = T))

pdf("Total Protein Quantity after diann_maxlfq AEcorr medianscale FullMatrix.pdf",width=20,height=10,family = "Times")

ggplot(PG_quants_BCMS_TPC,
       aes( x = File.Name,
            y = log2(Total_Protein_Quantity)))+
  geom_bar(stat="identity",fill="navyblue",colour="white")+
  theme_SKA()+
  theme(axis.text.x = element_text(angle=90,size=5))
dev.off()


###################################################################################################
###################################################################################################
###################################################################################################
########################################### All samples ########################################### 
###################################################################################################
###################################################################################################
###################################################################################################

AS_QCplot_forFullMatrix_dir = paste0(QCplotdir,"/forFullMatrix/All Samples")

dir.create(AS_QCplot_forFullMatrix_dir,recursive = T)

###################
### Input Data ###
###################

## Read in uniprot - sgd - gene name conversion file to add gene names ###

GenProt_SGD <- read.csv(paste0(dbdir,"SGD_Uniprot.csv"),stringsAsFactors = F)

#######################################
### Read in Experiment Design files ###
#######################################

setwd(Int_datdir)

Exp_Design<-read.csv("EDP_ExperimentDesign.csv",stringsAsFactors = F)

Exp_Design_AEonly<-read.csv("EDP_ExperimentDesign_AEonly.csv",stringsAsFactors = F)

############################################################################################################################
### Filter out Precursors that have InterBatch CVs > 30 % in AE samples  & Proteins that were detected in < 80 % samples ###
############################################################################################################################

Sample_data <- Sample_data%>% 
  mutate(Prot_PrecID = paste(Protein.Ids,Precursor.Id,sep="_"))%>%
  filter(Prot_PrecID %in% PrecIDList)%>%
  filter(Protein.Ids %in% Prot2keep)

### Compute Protein Group quantities ###

PG_quants<-get_PGQ_df(Sample_data)

## Merge PG_qunats with annotation

PG_quants<-merge(PG_quants,
                 unique(Sample_data[,c("File.Name","Layout","Element",
                                       "Element.Concentration",
                                       "Digestion_Batch",  "TimePoint",
                                       "LO_TP","BioSpecID","Sample.Type","Date","Time")]),
                 by="File.Name")

## Remove samples with < 1300 proteins quantified 

numprot_filter <- PG_quants %>%
  dplyr::select(File.Name,Protein.Ids,Protein.Quantity)%>%
  na.omit()%>%
  group_by(File.Name)%>%
  summarize(numprot_detected = length(Protein.Ids))%>%
  ungroup()%>%
  filter(numprot_detected >=1300)%>%
  pull(File.Name)

PG_quants <- filter(PG_quants, File.Name %in% numprot_filter)      

## Setup annotations for proBatch

setup_proBatch_annos(PG_quants,onlyAEdat = F)
essential_columns <- c(PGcol, PGintcol, fname_col)


Exp_Design <- filter( Exp_Design, label %in% numprot_filter)

### Prepare data for use with DEP ###

DEP_PGquants <- reshape2::dcast ( PG_quants,
                                  Protein.Ids~File.Name,
                                  value.var = "Protein.Quantity")%>%
  mutate ( name = Protein.Ids,
           ID = Protein.Ids)

### Create Summarized Experiment Object ##

## Note : this will log2 transform all values

SumExp <- make_se(DEP_PGquants,seq(2,ncol(DEP_PGquants)-2),Exp_Design )

############################################
### Batch Correct at the precursor level ###
############################################

Samdata_fMedCorr <- AEmed_correction_Precursor(Sample_data)

### Convert batch corrected data into protein quantities ###

PG_quants_BC<-get_PGQ_df(Samdata_fMedCorr)

## Merge PQ_quants with annotation

PG_quants_BC<-merge(PG_quants_BC,
                    unique(Sample_data[,c("File.Name","Layout","Element","Element.Concentration",
                                          "Digestion_Batch","TimePoint",
                                          "LO_TP","BioSpecID","Sample.Type","Date","Time")]),
                    by="File.Name")

PG_quants_BC <- filter(PG_quants_BC, File.Name %in% numprot_filter)    

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

PG_quants_BCMS <- filter(PG_quants_BCMS, File.Name %in% numprot_filter)      


DEP_PGquants_BCMS <- reshape2::dcast ( PG_quants_BCMS,
                                       Protein.Ids~File.Name,
                                       value.var = "Median_Scaled_Protein_Quantity")%>%
  mutate ( name = Protein.Ids,
           ID = Protein.Ids)

SumExp_BatchCorrMedSc <- make_se(DEP_PGquants_BCMS,seq(2,ncol(DEP_PGquants_BCMS)-2),Exp_Design)  

##############################
### Fill in missing values ###
##############################

for_NAfill <- PG_quants_BCMS %>%
  reshape2::dcast(Protein.Ids~File.Name, value.var = "Median_Scaled_Protein_Quantity")

# If missing in all replicates of a condition --> fill in with minimum value measured in dataset
# Else - replace with median of all measured replicates

NAfree.df <- fill_in_NAs(for_NAfill)

NAfree.df.fSE <- data.frame(NAfree.df)%>%
  mutate(ID = rownames(NAfree.df),
         name = rownames(NAfree.df))

SumExp_noNA <- DEP::make_se(NAfree.df.fSE,seq(1:(ncol(NAfree.df.fSE)-2)),Exp_Design)

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

PG_quants_limmaBC <- filter(PG_quants_limmaBC, File.Name %in% numprot_filter)    

DEP_PGquants_limmaBC <- reshape2::dcast ( PG_quants_limmaBC,
                                          Protein.Ids~File.Name,
                                          value.var = "Protein.Quantity")%>%
  mutate ( name = Protein.Ids,
           ID = Protein.Ids)


SumExp_LimmaBatchCorr <- make_se(DEP_PGquants_limmaBC,seq(2,ncol(DEP_PGquants_limmaBC)-2),Exp_Design)  

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

############################################
### Fill in NAs for limma corrected data ###
############################################

for_NAfill_limmaBCMS <- PG_quants_limmaBCMS %>%
  reshape2::dcast(Protein.Ids~File.Name, value.var = "Median_Scaled_Protein_Quantity")

# If missing in all replicates of a condition --> fill in with minimum value measured in dataset
# Else - replace with median of all measured replicates

NAfree.df_limmaBCMS <- fill_in_NAs(for_NAfill_limmaBCMS)

NAfree.df.fSE_limmaBCMS <- data.frame(NAfree.df_limmaBCMS)%>%
  mutate(ID = rownames(NAfree.df_limmaBCMS),
         name = rownames(NAfree.df_limmaBCMS))

SumExp_limma_noNA <- DEP::make_se(NAfree.df.fSE_limmaBCMS,seq(1:(ncol(NAfree.df.fSE_limmaBCMS)-2)),Exp_Design)

##################################
### Collect all SumExp objects ### 
##################################

SEos<-list(
  SumExp,
  SumExp_BatchCorr,
  SumExp_BatchCorrMedSc,
  
  SumExp_LimmaBatchCorr,
  SumExp_limmaBatchCorrMedSc,
  SumExp_noNA,
  SumExp_limma_noNA
)

names(SEos) <- c("Log2 Raw",
                 "BatchCorrected",
                 "BatchCorrected_medianScaled",
                 "limmaBatchCorrected",
                 "limmaBatchCorrected_medianScaled",
                 "BatchCorrected_medianScaled_noNA",
                 "limmaBatchCorrected_medianScaled_noNa")

# Get DFs out of se object ( these will have log2 scale values )

SEdfs<-lapply(SEos, SumExpObj_to_df)

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

Exp_Design_AEonly<-filter(Exp_Design_AEonly,label %in% colnames(DEP_PGquants_limmaBC_AEonly))


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


# Batch corrected and Median scaled and noNAs

DEP_PGquants_limma_noNA_AEonly <- NAfree.df.fSE[,c(which(colnames(NAfree.df.fSE)%in% colnames(DEP_PGquants_limmaBC_AEonly) ))]

SumExp_limma_noNA_AEonly <- make_se(DEP_PGquants_limma_noNA_AEonly,seq(1,ncol(DEP_PGquants_limma_noNA_AEonly)-2),Exp_Design_AEonly)  


SEos_limma_AEonly <- list(SumExp_limmaBatchCorr_AEonly,SumExp_limmaBatchCorrMedSc_AEonly,SumExp_limma_noNA_AEonly)

names(SEos_limma_AEonly) <- c(
  "limmaBatchCorrected AEsamp only",
  "limmaBatchCorrected_medianScaled AEsamp only",
  "limmaBatchCorrected_medianScaled_noNa AEsamp only")

SEdfs_limma_AEonly <- lapply(SEos_limma_AEonly, SumExpObj_to_df)

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
swbp4 <- plot_normalization(SumExp_noNA)+
  scale_fill_manual(values=c(viridis(n=111),c(viridis(n=111))))+
  theme(legend.position = "none")

setwd(AS_QCplot_forFullMatrix_dir)
pdf("Effect of Filtering_AEmedCorrection_NAhandling on data.pdf",width=20,height=20)

idp1

covp1

np1

grid.arrange(pd1)
dev.off()

pdf("Median protein Quantities across data processing stages.pdf",width=40,height=20)

swbp1+swbp2+swbp3+ swbp4+plot_layout(ncol=4)

dev.off()

#########################################
### Plots to check intensity profiles ###
#########################################

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

pdf("Protein Intensity profiles facet.pdf",width=20,height=15)
for(i in 1:length(SEdfs)){
  plot_protint_dist_facet(SEdfs[[i]], names(SEdfs)[[i]])
  
}

pdf("Protein Intensity profiles allsamp.pdf",width=20,height=15)
for(i in 1:length(SEdfs)){
  plot_protint_dist(SEdfs[[i]], names(SEdfs)[[i]])
}

## check which sample in LO5 gives the different shape
#LO5fns <- unique(filter(Sdata,Layout ==5)[,'File.Name'])
#dfLO5 <- as.data.frame(SEdfs[[5]])
#dfLO5 <- dfLO5[,c(LO5fns)]%>%
#        reshape2::melt()

#ggplot(dfLO5,
#       aes(x=value,
#           colour=variable))+
#  geom_density()

# Remove C4 and check

#ggplot(filter(dfLO5,variable!="SKA_202008_EDP_LO_5_C4.wiff"),
#      aes(x=value,
#        colour=variable))+
#geom_density()

##############################
### Make proBatch QC plots ###
##############################

### Set up Pro batch annotations

sel_annos <- c("LO_TP","Sample.Type","Digestion_Batch") ## annotations to plot by

for( i in 1:length(SEdfs)){
  
  
  make_proBatchQCplots(SEdfs[[i]],
                       sel_anno_hclust=sel_annos,
                       onlyAEdat = F,
                       dsname=names(SEdfs)[[i]],
                       AE_or_samp = "Samples",
                       sanno = sample_anno )
  
}

## Plot AE samples only from the limma Batch correction

sel_annos_AEonly <- c("LO_TP","Digestion_Batch") ## annotations to plot by

sanno_AEonly <- filter(sample_anno,grepl("AllEle",BioSpecID))

for( i in 1:length(SEdfs_limma_AEonly)){
  
  
  make_proBatchQCplots(SEdfs_limma_AEonly[[i]],
                       sel_anno_hclust=sel_annos,
                       onlyAEdat = T,
                       dsname=names(SEdfs_limma_AEonly)[[i]],
                       AE_or_samp = "AllEle",
                       sanno = sanno_AEonly
  )
  
}



########################################################
### QC replicate correlation  ### - to find outliers ###
########################################################

## Use unimputed data after batch correction and median scaling for this 

#df4rc <- SEdfs[[3]] %>%
#  data.frame()%>%
 # mutate(Protein.Ids=rownames(SEdfs[[3]]))%>%
  #reshape2::melt(id.vars="Protein.Ids",variable.name= "File.Name",
   #              value.name="Log2(Protein.Quantity)")

#df4rc <- merge(df4rc,sample_anno,by="File.Name")%>%
#  mutate(BioSpecID = trimws(BioSpecID))%>%
#  filter(BioSpecID != "AllEleprepQC 1")

#bspcs <- unique(df4rc$BioSpecID)

#cor_check <- vector()

#for( i in 1:length(bspcs)){
  
#  df = filter(df4rc,BioSpecID == bspcs[[i]])%>%
#    reshape2::dcast(Protein.Ids~File.Name,value.var="Log2(Protein.Quantity)")
#  rownames(df) = df$Protein.Ids
  
  
#  if(ncol(df)>2){
    
#    df = df[,-1]
#    cordf=stats::cor(df,use = "pairwise.complete",method="spearman")
    
#    cor_check = rbind(cor_check,c(bspcs[[i]],min(cordf)))
#  }
#}
#colnames(cor_check)<-c("BioSpecID","Replicate_Correlation")
setwd(FinalOutCSVdir)
#write.csv(cordf,"Replicate correlations Spearman FullMatrix.csv")

##################################################
##################################################
### limma Batch Corrected DFs to write to disk ###
##################################################
##################################################

# Prepare data with NAs for writing to CSV 

Fin_wNA_df_limma <- merge(PG_quants_limmaBCMS,sample_anno[,c(1,which(!colnames(sample_anno) %in% colnames(PG_quants_limmaBCMS)))],
                    by="File.Name")%>%
  mutate(BioSpecID=trimws(BioSpecID),
         `Log2(Protein.Quantity)` = log2(Protein.Quantity))%>%
  filter(BioSpecID != "AllEleprepQC 1")


write.csv(Fin_wNA_df_limma,"Protein Quantities w NAs limmaBC FullMatrix.csv",row.names = F)


# Calculate per Gene mean, sd, median, mad for AllEle replicates across dataset and write to CSV

AERep_Stats <- Fin_wNA_df_limma%>%
  filter(BioSpecID == "AllEle 1")%>%
  group_by(Protein.Ids)%>%
  mutate(Num.Non.NAs = sum(!is.na(`Log2(Protein.Quantity)`)))%>%
  ungroup()

## Plot distribution of fraction of controls that are NAs

ggplot(unique(AERep_Stats[,c("Protein.Ids","Num.Non.NAs")]),
       aes(x=Num.Non.NAs))+
  geom_histogram(bins=25,colour="darkred",
                 fill="darkred")+
  theme_SKA()

AERep_Stats<-AERep_Stats%>%
  filter(Num.Non.NAs > 10)%>%
  group_by(Protein.Ids)%>%
  summarize(SD.Protein.Quant = sd(2^(`Log2(Protein.Quantity)`),na.rm=T),
            Mean.Protein.Quant = mean(2^(`Log2(Protein.Quantity)`),na.rm=T)
  )%>%
  ungroup()%>%
  mutate(Cov = SD.Protein.Quant/Mean.Protein.Quant)%>%
  as.data.frame()

rownames(AERep_Stats)<-AERep_Stats$Protein.Ids

write.csv(AERep_Stats,"AllEle Replicates Protein Quantity Stats limmaBC FullMatrix.csv",row.names=F)


# Calculate per Gene mean, sd, median, mad for all samples AllEle replicates across dataset and write to CSV

AllSample_Stats <- Fin_wNA_df_limma%>%
  filter(!grepl("AllEle",BioSpecID))%>%
  group_by(Protein.Ids)%>%
  mutate(Num.Non.NAs = sum(!is.na(`Log2(Protein.Quantity)`)))%>%
  ungroup()

## Plot distribution of fraction of controls that are NAs

ggplot(unique(AllSample_Stats[,c("Protein.Ids","Num.Non.NAs")]),
       aes(x=Num.Non.NAs))+
  geom_histogram(bins=25,colour="darkred",
                 fill="darkred")+
  theme_SKA()

AllSample_Stats<-AllSample_Stats%>%
  filter(Num.Non.NAs >10)%>%
  group_by(Protein.Ids)%>%
  summarize(SD.Protein.Quant = sd(2^(`Log2(Protein.Quantity)`),na.rm=T),
            Mean.Protein.Quant = mean(2^(`Log2(Protein.Quantity)`),na.rm=T)
  )%>%
  ungroup()%>%
  mutate(Cov = SD.Protein.Quant/Mean.Protein.Quant)%>%
  as.data.frame()

rownames(AllSample_Stats)<-AllSample_Stats$Protein.Ids

write.csv(AllSample_Stats,"AllSample Protein Quantity Stats limmaBC FullMatrix.csv",row.names=F)


Elewise_Stats <- Fin_wNA_df_limma%>%
  filter(!grepl("AllEle",BioSpecID))%>%
  group_by(Protein.Ids,Element)%>%
  summarize(SD.Protein.Quant = sd(2^(`Log2(Protein.Quantity)`),na.rm=T),
            Mean.Protein.Quant = mean(2^(`Log2(Protein.Quantity)`),na.rm=T)
  )%>%
  mutate(Cov = SD.Protein.Quant/Mean.Protein.Quant)%>%
  as.data.frame()

rownames(Elewise_Stats)<-paste(Elewise_Stats$Protein.Ids,Elewise_Stats$Element)

write.csv(Elewise_Stats,"Elewise Protein Quantity Stats limmaBC FullMatrix.csv",row.names=F)


###################################################
### Prepare data without NAs for writing to CSV ### limma Batch Corrected
###################################################

### limma Batch Corrected

FinnoNA_df_limma <- NAfree.df_limmaBCMS%>%
  data.frame()%>%
  mutate(Protein.Ids=rownames(NAfree.df))%>%
  reshape2::melt(id.vars="Protein.Ids",variable.name= "File.Name",
                 value.name="Protein.Quantity")%>%
  mutate(`Log2(Protein.Quantity)`=log2(Protein.Quantity))

FinnoNA_df_limma <- merge(FinnoNA_df_limma,sample_anno,by="File.Name")%>%
  mutate(BioSpecID = trimws(BioSpecID))%>%
  filter(BioSpecID != "AllEleprepQC 1")

write.csv(FinnoNA_df_limma,"Protein Quantities wImpValues limmaBC FullMatrix.csv",row.names = F)

###########################################
### Write unimputed data for clustering ### limma Batch Corrected
###########################################

### limma Batch Corrected

FinwNA_df_forclus <- Fin_wNA_df_limma

write.csv(FinwNA_df_forclus,"Protein Quantities noImp for clustering limmaBC FullMatrix.csv",row.names=F)


######################################################
######################################################
### AE median Batch Corrected DFs to write to disk ###
######################################################
######################################################

Fin_wNA_df <-  merge(PG_quants_BCMS,sample_anno[,c(1,which(!colnames(sample_anno) %in% colnames(PG_quants_BCMS)))],by="File.Name")%>%
  data.frame()%>%
  mutate(BioSpecID=trimws(BioSpecID),
         `Log2(Protein.Quantity)` = log2(Protein.Quantity))%>%
  filter(BioSpecID != "AllEleprepQC 1")


write.csv(Fin_wNA_df,"Protein Quantities w NAs FullMatrix.csv",row.names = F)


# Calculate per Gene mean, sd, median, mad for AllEle replicates across dataset and write to CSV

AERep_Stats <- Fin_wNA_df%>%
  filter(BioSpecID =="AllEle 1")%>%
  group_by(Protein.Ids)%>%
  mutate(Num.Non.NAs = sum(!is.na(`Log2(Protein.Quantity)`)))%>%
  ungroup()

## Plot distribution of fraction of controls that are NAs

ggplot(unique(AERep_Stats[,c("Protein.Ids","Num.Non.NAs")]),
       aes(x=Num.Non.NAs))+
  geom_histogram(bins=25,colour="darkred",
                 fill="darkred")+
  theme_SKA()

AERep_Stats<-AERep_Stats%>%
  filter(Num.Non.NAs > 10)%>%
  group_by(Protein.Ids)%>%
  summarize(SD.Protein.Quant = sd(2^(`Log2(Protein.Quantity)`),na.rm=T),
            Mean.Protein.Quant = mean(2^(`Log2(Protein.Quantity)`),na.rm=T)
  )%>%
  ungroup()%>%
  mutate(Cov = SD.Protein.Quant/Mean.Protein.Quant)%>%
  as.data.frame()

rownames(AERep_Stats)<-AERep_Stats$Protein.Ids

write.csv(AERep_Stats,"AllEle Replicates Protein Quantity Stats FullMatrix.csv",row.names=F)


# Calculate per Gene mean, sd, median, mad for all samples AllEle replicates across dataset and write to CSV

AllSample_Stats <- Fin_wNA_df%>%
  filter(!grepl("AllEle",BioSpecID))%>%
  group_by(Protein.Ids)%>%
  mutate(Num.Non.NAs = sum(!is.na(`Log2(Protein.Quantity)`)))%>%
  ungroup()

## Plot distribution of fraction of controls that are NAs

ggplot(unique(AllSample_Stats[,c("Protein.Ids","Num.Non.NAs")]),
       aes(x=Num.Non.NAs))+
  geom_histogram(bins=25,colour="darkred",
                 fill="darkred")+
  theme_SKA()

AllSample_Stats<-AllSample_Stats%>%
  filter(Num.Non.NAs >10)%>%
  group_by(Protein.Ids)%>%
  summarize(SD.Protein.Quant = sd(2^(`Log2(Protein.Quantity)`),na.rm=T),
            Mean.Protein.Quant = mean(2^(`Log2(Protein.Quantity)`),na.rm=T)
  )%>%
  ungroup()%>%
  mutate(Cov = SD.Protein.Quant/Mean.Protein.Quant)%>%
  as.data.frame()

rownames(AllSample_Stats)<-AllSample_Stats$Protein.Ids

write.csv(AllSample_Stats,"AllSample Protein Quantity Stats FullMatrix.csv",row.names=F)


Elewise_Stats <- Fin_wNA_df%>%
  filter(!grepl("AllEle",BioSpecID))%>%
  group_by(Protein.Ids,Element)%>%
  summarize(SD.Protein.Quant = sd(2^(`Log2(Protein.Quantity)`),na.rm=T),
            Mean.Protein.Quant = mean(2^(`Log2(Protein.Quantity)`),na.rm=T)
  )%>%
  mutate(Cov = SD.Protein.Quant/Mean.Protein.Quant)%>%
  as.data.frame()

rownames(Elewise_Stats)<-paste(Elewise_Stats$Protein.Ids,Elewise_Stats$Element)

write.csv(Elewise_Stats,"Elewise Protein Quantity Stats FullMatrix.csv",row.names=F)

#### Replicate stats

Ele_RepStats <- Fin_wNA_df%>%
  filter(!grepl("AllEle",BioSpecID))%>%
  group_by(Protein.Ids,Element,BioSpecID)%>%
  summarize(SD.Protein.Quant = sd(2^(`Log2(Protein.Quantity)`),na.rm=T),
            Mean.Protein.Quant = mean(2^(`Log2(Protein.Quantity)`),na.rm=T)
  )%>%
  mutate(Cov = SD.Protein.Quant/Mean.Protein.Quant)%>%
  ungroup()%>%
  group_by(Protein.Ids, Element)%>%
  summarize(Mean.COVreps = mean(Cov,na.rm=T))%>%
  ungroup()%>%
  as.data.frame()

rownames(Ele_RepStats)<-paste(Ele_RepStats$Protein.Ids,Ele_RepStats$Element)

write.csv(Ele_RepStats,"Protein Quantity Replicate Stats Branch C.csv",row.names=F)

Ele_RepStats_Summary <-na.omit(Ele_RepStats)%>%
  group_by(Element)%>%
  summarize(NumProt = length(unique(Protein.Ids)),
         Mean_COV = mean(Mean.COVreps,na.rm=T))%>%
  ungroup()

write.csv(Ele_RepStats_Summary,"Protein Quantity Replicate Stats Summary Branch C.csv",row.names=F)


###################################################
### Prepare data without NAs for writing to CSV ### AE median Batch Corrected
###################################################

FinnoNA_df <- NAfree.df%>%
  data.frame()%>%
  mutate(Protein.Ids=rownames(NAfree.df))%>%
  reshape2::melt(id.vars="Protein.Ids",variable.name= "File.Name",
                 value.name="Protein.Quantity")%>%
  mutate(`Log2(Protein.Quantity)`=log2(Protein.Quantity))

FinnoNA_df <- merge(FinnoNA_df,sample_anno,by="File.Name")%>%
  mutate(BioSpecID = trimws(BioSpecID))%>%
  filter(BioSpecID != "AllEleprepQC 1")

write.csv(FinnoNA_df,"Protein Quantities wImpValues FullMatrix.csv",row.names = F)


### Write unimputed data for clustering ###

FinwNA_df_forclus <- Fin_wNA_df

write.csv(FinwNA_df_forclus,"Protein Quantities noImp for clustering FullMatrix.csv",row.names=F)


#############################
### QC number of proteins ###
#############################

numprot <- FinwNA_df_forclus %>%
  na.omit()%>%
  dplyr::select(Protein.Ids,Layout,Digestion_Batch,BioSpecID)%>%
  unique()%>%
  group_by(Layout,Digestion_Batch,BioSpecID)%>%
  summarize(Num_Prot = length(Protein.Ids))

pdf("Number of proteins per Layout Branch C.pdf",width = 5,height=4)
ggplot(numprot,
       aes(x=factor(Layout),y=Num_Prot,
           colour=Digestion_Batch,
           group=Layout))+
  geom_point(size=3,alpha=0.6,position=position_jitterdodge())+
  scale_color_manual(values= c("darkred","darkblue"))+
  geom_boxplot(alpha=0,colour="black")+
  theme_SKA()+
  theme(legend.position = "bottom",
        panel.grid.major = element_line(colour="black",size = 0.1),
        panel.border = element_rect(size=0.15),
        axis.ticks = element_line(size=0.75))+
  labs(x="Layout", y="Number of proteins quantified")
dev.off()

numprot_BioSpecID <- FinwNA_df_forclus %>%
  na.omit()%>%
  dplyr::select(Protein.Ids,BioSpecID)%>%
  unique()%>%
  group_by(BioSpecID)%>%
  summarize(Num_Prot = length(Protein.Ids))

pdf("Num protein per condition Dataset C ylim.pdf",width = 13,height=5)
ggplot(numprot,
       aes(x=BioSpecID,y=Num_Prot,
           colour=BioSpecID,
       ))+
  geom_point(size=3,alpha=1.0)+
  scale_color_manual(values=colkey_BioSpecID)+
  theme_SKA()+
  theme(
    panel.grid.major = element_line(colour="lightgray",size = 0.1),
    panel.border = element_rect(size=0.15),
    legend.position = "none",
    axis.text.x = element_text(angle=90,size=9.5,hjust = 0.9),
    axis.ticks = element_line(size=0.75))+
  labs(x="", y="Number of proteins quantified")+
  ylim(1300,2500)
dev.off()
