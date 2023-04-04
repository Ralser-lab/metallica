#`---
#`  Title: "metallomics_3_analyze_intracellular_icpms_data.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 8 April 2021 
#`  Description: Script to analyze all ICP-MS data from intracellular metallomics measurements
#`---

### source paths functions and libraries  ###
#############################################

# general
source("/camp/lab/ralserm/working/Simran Aulakh/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/layout_conversion.R"))

# experiment specific
source(paste0(code_dir,"/metpert_WTmetallomics/metallomics_0_libraries_functions.R"))

#################
### set paths ###
#################

# inputs
icpms_data_dir <- paste0(metpert_WTmetallomics_dir,"/data_from_agilent_masshunter/intracellular_metallomics")
od_dir <- paste0(metpert_WTmetallomics_dir,"/od600_at_sampling")

# outputs
plot_dir <- paste0(metpert_WTmetallomics_dir,"/results/plots/intracellular_metallomics")
outputtables_dir <- paste0(metpert_WTmetallomics_dir,"/results/output_tables")

dir.create(plot_dir, recursive = T)
dir.create(outputtables_dir, recursive = T)

####################################################
### read in data processed in Agilent MassHunter ###
####################################################

metallomics_data = read_excel(paste0(icpms_data_dir,"/data_intracellular_icpms_from_AgilentMassHunter.xlsx"))

#################################################
### Separate data into Blank + QC and Samples ###
#################################################

BlankQC_df <- metallomics_data %>%
              filter(grepl("P1-",SampleName) | grepl("P4-",SampleName))%>%
              mutate(Type = ifelse(grepl("P4-",SampleName),"QC","Blank"))%>%
              dplyr::select(-SampleName)

Blank_df <- metallomics_data %>%
            filter(grepl("P1-",SampleName) | grepl("P4-",SampleName))%>%
            mutate(Type = ifelse(grepl("P4-",SampleName),"QC","Blank"))%>%
            filter(Type=="Blank")%>%
            dplyr::select(-SampleName)

Sample_df <- as.data.frame(metallomics_data) %>%
              filter(grepl("P3-",SampleName))%>%
              mutate(SampleName=gsub("P3-","",SampleName))%>%
              mutate(plate_position = paste(LO,SampleName,sep = "_"))

###################################
### Subtract blank from samples ###
###################################

Blank_medians <- colMedians(as.matrix(Blank_df[,c(2:12)]), na.rm = T)

Sample_df_sb <- Sample_df[c(3:13)]
Sample_df_sb <- sweep(as.matrix(Sample_df_sb), 2, Blank_medians, FUN="-")

Sample_df[,c(3:13)] <- Sample_df_sb

##########################################################
### Normalize all CPS to dilution volume to get CPS/mL ###
##########################################################

for(i in 1:nrow(Sample_df)){
  
  lo=as.numeric(Sample_df[i,"LO"])
  
  if(lo == 1){
    
    Sample_df[i,3:13] = as.numeric(Sample_df[i,3:13])/1.054
    
  }else if (lo ==2){
    
    Sample_df[i,3:13] = as.numeric(Sample_df[i,3:13])/0.844
    
  }else if (lo ==3){
    
    Sample_df[i,3:13] = as.numeric(Sample_df[i,3:13])/0.747
    
  }else if (lo %in% c(4,7)){
    
    Sample_df[i,3:13] = as.numeric(Sample_df[i,3:13])/0.709
    
  }else if (lo %in% c(5,6)){
    
    Sample_df[i,3:13] = as.numeric(Sample_df[i,3:13])/0.805
    
  }else if (lo == 8){
    
    Sample_df[i,3:13] = as.numeric(Sample_df[i,3:13])/0.748
    
  }
}

###############################
### Read in 96 well layouts ###
###############################

lo_fns <- dir(lo_dir)

layouts_allplates <- read.csv(paste0(acc_file_dir,"/all_layouts_combined.csv"),stringsAsFactors = F)%>%
                     mutate(condition = gsub ("Empty","ProcessBlank",condition))%>%
                     mutate(condition = gsub ("AllEle_prepQC","ProcessBlank", condition))%>%
                     mutate(condition = gsub ("AllEle_New","AllEle", condition))%>%
                     mutate(condition = gsub ("0x","0",condition))%>%
                     mutate(condition = gsub ("_x/5000","_0.005",condition))%>%
                     mutate(condition = gsub ("_x/1000","_0.001",condition))%>%
                     mutate(condition = gsub ("_x/100","_0.01",condition))%>%
                     mutate(condition = gsub ("_x/50","_0.02",condition))%>%
                     mutate(condition = gsub ("_x/20","_0.05",condition))%>%
                     mutate(condition = gsub ("_x/10","_0.10",condition))%>%
                     mutate(condition = gsub ("_x/5","_0.20",condition))%>%
                     mutate(condition = gsub ("_x/2","_0.5",condition))%>%
                     mutate(condition = gsub ("_2x","_2",condition))%>%
                     mutate(condition = gsub ("_5x","_5",condition))%>%
                     mutate(condition = gsub("AllEle","AllEle_1",condition))%>%
                     dplyr::select(-c(position))

#####################################
### Add Sample Names to Sample df ###
#####################################

Sample_df$SampleName <- layouts_allplates[match(Sample_df$plate_position,layouts_allplates$plate_position),"condition"]

Sample_df <- Sample_df%>%
              filter(!plate_position %in% c("2_F2","2_C2","1_C6","2_F1",
                                            "2_B2","2_E2",
                                            "4_F8","6_D2",
                                            # Filter out old Fe media in plate 6
                                            "6_F1","6_B2","6_C2","6_E2","6_F2",
                                            "6_G2","6_D3","6_E7","6_F7","6_B11",
                                            "6_C11","6_D11","6_G11","6_B12"
                                            
                                            ))%>%
              ## Filter out all chelator samples - all of them have a "+" in the sample name
              filter(!grepl("[+]",SampleName))%>%
              ## Layout 8 just had chelators , so no need to process the AllEle and QCs from this plate either
              filter(LO != "8")

######################################################
### Merge Process Blanks with run Blanks and QC df ###
######################################################

ProcessBlank_df <- Sample_df%>%
                    filter(SampleName == "ProcessBlank")%>%
                    filter(grepl("H",plate_position) | grepl("A",plate_position))%>%
                    mutate(Type = SampleName)%>%
                    dplyr::select(-c(SampleName,plate_position))

Sample_df <- Sample_df%>%
              filter(SampleName != "ProcessBlank" & !grepl("Ca_10",SampleName))%>%
              filter(!(grepl("C1",plate_position)|grepl("D1",plate_position)|grepl("E1",plate_position)))

BlPrBlQC_df <- rbind(ProcessBlank_df,BlankQC_df)%>%
               reshape2::melt(id.vars=c("Type","LO"),value.name = "CPS",variable.name="element_measured")

###############################
### Stats on QCs and Blanks ###
###############################

BlPrBlQC_Stat <- BlPrBlQC_df %>%
                  na.omit()%>%
                  group_by(Type,LO,element_measured)%>%
                  summarize(Mean.CPS=mean(CPS),
                            sd.CPS=sd(CPS))%>%
                  mutate(LOQ = ifelse(Type == "ProcessBlank",Mean.CPS+(5*sd.CPS),NA))%>%
                  ungroup()

# Plot

pdf(paste0(plot_dir,"/blanks_qc_loq.pdf"),width=12,height=8)
        ggplot(BlPrBlQC_Stat,aes(x=as.factor(LO),
                                 y=log2(Mean.CPS), 
                                 colour=Type, group=Type))+
          geom_point(size=3,
                     alpha=0.8,
                     position=position_dodge(width=0.8))+
          geom_errorbar(aes(ymin=log2(Mean.CPS-sd.CPS), ymax=log2(Mean.CPS+sd.CPS)), width=.5,
                        position=position_dodge(width=0.8),
                        colour="black") +
          geom_point(shape="-",colour="red",size=8,
                     aes(y=log2(LOQ)))+
          facet_wrap("element_measured",scales = "free") +
          scale_colour_manual(values=c("#A4D3EE","#4F94CD","#8B4789"))+
          theme_metallica()+
          theme( text = element_text(family = "Tinos"))+
          theme(strip.background =element_rect(fill="white"))+
          theme(strip.text  =element_text(size=13,colour="darkblue"))+
          labs(x="batch",y="mean(log2(counts per second))")

dev.off()

#####################################################
### Filter Sample data according to LOQ per Batch ###
####################################################

# Save LOQ df from QC stats df

LOQ_df <- BlPrBlQC_Stat %>%
          dplyr::select(LO,element_measured,LOQ)%>%
          na.omit()

Sample_df_long <- Sample_df %>%
                  reshape2::melt(id.vars=c("LO","SampleName","plate_position"),
                                 variable.name="element_measured", value.name = "CPS")

# Print number of samples before LOQ based filtering
paste("Number of samples BEFORE LOQ filter:",length(unique(paste(Sample_df$SampleName,
                                                                 Sample_df$plate_position))))

Sample_df_LOQfilt <- merge(Sample_df_long, LOQ_df,by=c("LO","element_measured"))%>%
                     filter(CPS > LOQ)%>%
                     mutate(Type ="Samples")

paste("Number of samples AFTER LOQ filter:",length(unique(paste(Sample_df_LOQfilt$SampleName,
                                                                Sample_df_LOQfilt$plate_position))))

####################################################################
### Combine Blanks, Process Blanks and QC with Process Blank data ###
####################################################################

BlPrBlQC_df_LOQ <- merge(BlPrBlQC_df,LOQ_df,by=c("LO","element_measured"))%>%
                    # Filter out QCs that are below LOQ
                   filter( Type != "QC" | CPS > LOQ)%>%
                   mutate(SampleName = Type)


AllData_LOQ_filt <- rbind(Sample_df_LOQfilt[,c("LO","element_measured","SampleName","Type","CPS","LOQ")],
                          BlPrBlQC_df_LOQ[,c("LO","element_measured","SampleName","Type","CPS","LOQ")])

###############################################
### Plot all Sample Values after LOQ filter ###
###############################################

pdf(paste0(plot_dir,"/all_metallomics_samples_after_loq_filter.pdf"),width=12,height=8)

ggplot(AllData_LOQ_filt,
       aes(x=LO,
           y=log2(CPS),
           colour=Type))+
        geom_point(alpha=0.6,size=3)+
        facet_wrap("element_measured",scales="free")+
        geom_point(shape="-",colour="red",size=10,
                   aes(y=log2(LOQ)))+
        scale_colour_manual(values=c("#A4D3EE","#4F94CD","#8B4789","black"))+
        theme_metallica()+
        theme(strip.background =element_rect(fill="white"))+
        theme(strip.text  =element_text(size=13,colour="darkblue"))+
        labs(x="Batch",y="log2(counts per second)")

dev.off()

######################
### Normalize to P ###
######################

Sample_df_LOQfilt_cast <- Sample_df_LOQfilt %>%
                          dplyr::select(plate_position,SampleName,element_measured,CPS)%>%
                          reshape2::dcast(plate_position+SampleName ~ element_measured, value.var ="CPS")

rownames(Sample_df_LOQfilt_cast) <- paste(Sample_df_LOQfilt_cast$plate_position,Sample_df_LOQfilt_cast$SampleName)
Sample_df_LOQfilt_cast_for_Norm <- Sample_df_LOQfilt_cast[,-c(1,2)]

Sam_Pnorm <- Sample_df_LOQfilt_cast_for_Norm/Sample_df_LOQfilt_cast_for_Norm$P

###
# Calculate scaling factor -- median P content in AE cells of that batch to bring each metal back to original PPB range

AErows <- rownames(Sam_Pnorm)[grep("AllEle",rownames(Sam_Pnorm))]

scaling <- colMeans(as.matrix(Sample_df_LOQfilt_cast_for_Norm[which(rownames(Sample_df_LOQfilt_cast_for_Norm) %in% AErows),]))["P"]

scaling <- rep(scaling,ncol(Sam_Pnorm))

Sam_Pnorm <- sweep(as.matrix(Sam_Pnorm), 2, scaling, FUN="*")
Sam_Pnorm <- as.data.frame(Sam_Pnorm)
Sam_Pnorm$normalization <- "Phosphorus"
Sam_Pnorm$Sample.Name <- rownames(Sam_Pnorm)

Sample_df_LOQfilt_cast$normalization <- "Unnormalized"
Sample_df_LOQfilt_cast$Sample.Name <- rownames(Sample_df_LOQfilt_cast)

###############
### Cleanup ###
###############

Sample_data_n <- rbind(Sam_Pnorm, Sample_df_LOQfilt_cast[,c(colnames(Sam_Pnorm))])%>%
                  filter(!grepl("_o",Sample.Name))%>%
                  separate(Sample.Name, into = c("plate_position","Sample.Name"), sep = " ")%>%
                  separate(plate_position,into=c("plate",NA),sep="_",remove = F)%>%
                  separate(Sample.Name,into = c("element_perturbed","rel_element_concentration_th"),sep="_")%>%
                  mutate(rel_element_concentration_th = as.numeric(rel_element_concentration_th),
                         rel_element_concentration_th = ifelse(grepl("AllEle",element_perturbed),1,rel_element_concentration_th))%>%
                  reshape2::melt(id.vars=c("normalization","plate_position","plate","element_perturbed","rel_element_concentration_th" ),
                                 variable.name = "element_measured")


#####################
### Merge with OD ###
#####################

od_data <- read.csv(paste0(od_dir,"/OD_at_sampling_metallomics.csv"),stringsAsFactors = F)

##################################
### Convert PPB to ng per cell ###
##################################

ConvFac <- read.csv(paste0(acc_file_dir,"/conversion_factors_ppb_to_ngpercell.csv"),stringsAsFactors = F)

Sample_data_n <- merge(Sample_data_n,od_data, by.x="plate_position",by.y="lop")

Sample_data_n <- merge(Sample_data_n,ConvFac,by.x="plate",by.y="Layout")%>%
                 filter(SpecOD > 0.1)%>%
                 mutate(ng_perwell = (value * Final_Vol_for_ICPMS))%>%
                 mutate(pg_percell = (ng_perwell*1000)/
                                     (SpecOD*1.5 * 10000000 *Volume_perwell_of_culture_transferred_to_filterplate))

#####################################
### Are there any batch effects ? ###
#####################################

pdf(paste0(plot_dir,"/batch_effects_pgpercell_data_BEFORE_correction.pdf"),width=20,height=15)
ggplot(filter(Sample_data_n,element_perturbed=="AllEle" & normalization =="Phosphorus"),
       aes(x = plate, y= pg_percell))+
  geom_point()+
  facet_wrap("element_measured",scales="free")+
  theme_metallica()+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text  =element_text(size=13,colour="darkblue"))+
  labs(x="batch",y="log2(counts per second)")
dev.off()

# -- YES! tiny ones

#############################
### Correct Batch Effects ###
#############################

## All element values for each batch
## Batch correction using median od all ele samples 

AE_for_batchcorr <- Sample_data_n%>%
                    filter(element_perturbed == "AllEle")%>%
                    group_by(plate,normalization,element_measured)%>%
                    summarize(Median.pgpcell.AE.batch = median(pg_percell,na.rm=T),
                              Median.ngpwell.AE.batch = median(ng_perwell,na.rm=T))%>%
                    ungroup()%>%
                    group_by(normalization,element_measured)%>%
                    mutate(Median.pgpcell.Across.batches = median(Median.pgpcell.AE.batch,na.rm=T),
                           Median.ngpwell.Across.batches = median(Median.ngpwell.AE.batch,na.rm=T))

Sample_data_n <- merge(Sample_data_n,AE_for_batchcorr,by=c("plate","normalization","element_measured"))%>%
                  ## Batch correct and re-scale ##
                 mutate(pg_percell_BC = (pg_percell/Median.pgpcell.AE.batch)* Median.pgpcell.Across.batches,
                        ng_perwell_BC = (ng_perwell/Median.ngpwell.AE.batch)* Median.ngpwell.Across.batches)

## Visualize effects after batch correction

pdf(paste0("batch_effects_pgpercell_data_AFTER_correction.pdf"),width=20,height=15)
ggplot(filter(Sample_data_n,element_perturbed=="AllEle" & normalization =="Phosphorus"),
       aes(x = plate,group=plate))+
       geom_point(aes( y= pg_percell),
                   alpha=0.7,size=6,shape="*",
                   colour="black",
                   position = position_jitter(width = 0.2))+
       geom_point(aes( y= pg_percell_BC), 
                   alpha=0.5,size=2,
                   colour="darkred",
                   position = position_jitter(width = 0.2))+
       facet_wrap("element_measured",scales="free")+
       theme_metallica()+
       theme(strip.background =element_rect(fill="white"))+
       theme(strip.text  =element_text(size=13,colour="darkblue"))+
       labs(x="batch",y="pg / cell")
dev.off()

## Visualize effects after batch correction

pdf(paste0(plot_dir,"/batch_effects_ngperwell_data_AFTER_correction.pdf"),width=20,height=15)
ggplot(filter(Sample_data_n,element_perturbed=="AllEle" & normalization =="Phosphorus"),
       aes(x = plate,group=plate))+
  geom_point(aes( y= ng_perwell),
             alpha=0.7,size=6,shape="*",
             colour="black",
             position = position_jitter(width = 0.2))+
  geom_point(aes( y= ng_perwell_BC), 
             alpha=0.5,size=2,
             colour="darkred",
             position = position_jitter(width = 0.2))+
  facet_wrap("element_measured",scales="free")+
  theme_metallica()+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text  =element_text(size=13,colour="darkblue"))+
  labs(x="Batch",y="ng / well (P normalized)")
dev.off()


#########################################################
## AE qc after batch correction ### Unnormalized Data ###
#########################################################

AE_Unnorm_OD<- filter(Sample_data_n,
                      element_perturbed == "AllEle" &
                      normalization == "Unnormalized")

pdf(paste0(plot_dir,"/allele_controls_OD_vs_elementconcentration_raw_to_pgpercellBC_unnormalized.pdf"), width=15,height=10)

ggplot(AE_Unnorm_OD, aes(x = SpecOD, 
                         y= value,
                         colour = plate))+
        geom_point(size=3,alpha=0.8)+
        geom_smooth(method="lm")+
        facet_wrap("element_measured",scales="free")+
        labs(title = "Raw PPB values")+
        theme_metallica()

ggplot(AE_Unnorm_OD, aes(x = SpecOD, 
                         y= ng_perwell,
                         colour = plate))+
        geom_point(size=3,alpha=0.8)+  
        geom_smooth(method="lm")+
        facet_wrap("element_measured",scales="free")+
        labs(title = "ng per well")+
        theme_metallica()

ggplot(AE_Unnorm_OD, aes(x = SpecOD, 
                         y= ng_perwell_BC,
                         colour = plate))+
        geom_point(size=3,alpha=0.8)+  
        geom_smooth(method="lm")+
        facet_wrap("element_measured",scales="free")+
        labs(title = "ng per well - Batch corrected")+
        theme_metallica()

ggplot(AE_Unnorm_OD, aes(x = SpecOD, 
                         y= pg_percell,
                         colour = plate))+
        geom_point(size=3,alpha=0.8)+
        geom_smooth(method="lm")+
        facet_wrap("element_measured",scales="free")+
        labs(title = "Picogram per cell")+
        theme_metallica()

ggplot(AE_Unnorm_OD, aes(x = SpecOD, 
                         y= pg_percell,
                         colour = plate))+
        geom_point(size=3,alpha=0.8)+
        geom_smooth(method="lm")+
        facet_wrap("element_measured",scales="free")+
        labs(title = "Picogram per cell")+
        theme_metallica()

ggplot(AE_Unnorm_OD, aes(x = SpecOD, 
                         y= pg_percell_BC,
                         colour = plate))+
        geom_point(size=3,alpha=0.8)+
        geom_smooth(method="lm")+
        facet_wrap("element_measured",scales="free")+
        labs(title = "Picogram per cell - Batch corrected")+
        theme_metallica()
dev.off()



#########################################################
## AE qc after batch correction ### P normalized Data ###
#########################################################

AE_Pnorm_OD<- filter(Sample_data_n,
                     element_perturbed == "AllEle" &
                     normalization == "Phosphorus")

pdf(paste0(plot_dir,"/allele_controls_OD_vs_elementconcentration_raw_to_pgpercellBC_phosphorus_normalized.pdf"), width=15,height=10)

ggplot(AE_Pnorm_OD, aes(x = SpecOD, 
                        y= value,
                        colour = plate))+
        geom_point(size=3,alpha=0.8)+
        geom_smooth(method="lm")+
        facet_wrap("element_measured",scales="free")+
        labs(title = "Raw PPB values P norm")+
        theme_metallica()

ggplot(AE_Pnorm_OD, aes(x = SpecOD, 
                        y= ng_perwell,
                        colour = plate))+
        geom_point(size=3,alpha=0.8)+  
        geom_smooth(method="lm")+
        facet_wrap("element_measured",scales="free")+
        labs(title = "ng per well P norm")+
        theme_metallica()

ggplot(AE_Pnorm_OD, aes(x = SpecOD, 
                        y= ng_perwell_BC,
                        colour = plate))+
        geom_point(size=3,alpha=0.8)+  
        geom_smooth(method="lm")+
        facet_wrap("element_measured",scales="free")+
        labs(title = "ng per well - Batch corrected P norm")+
        theme_metallica()

ggplot(AE_Pnorm_OD, aes(x = SpecOD, 
                        y= pg_percell,
                        colour = plate))+
        geom_point(size=3,alpha=0.8)+
        geom_smooth(method="lm")+
        facet_wrap("element_measured",scales="free")+
        labs(title = "Picogram per cell P norm")+
        theme_metallica()


ggplot(AE_Pnorm_OD, aes(x = SpecOD, 
                        y= pg_percell_BC,
                        colour = plate))+
        geom_point(size=3,alpha=0.8)+
        geom_smooth(method="lm")+
        facet_wrap("element_measured",scales="free")+
        labs(title = "Picogram per cell - Batch corrected P norm")+
        theme_metallica()
dev.off()

write.csv(Sample_data_n,paste0(outputtables_dir,"/metallomics_alldata_batchcorrected.csv"),row.names=F)

Sample_data_finalised <- Sample_data_n%>%
                      filter(normalization == "Phosphorus" & element_perturbed != "B" & !element_measured %in% c("Cu","Mo"))%>%
                      dplyr::select(element_measured,plate_position,element_perturbed,rel_element_concentration_th,
                                    SpecOD,ng_perwell_BC)


write.csv(Sample_data_finalised,paste0(outputtables_dir,"/metallomics_ngperwell_Pnorm_finalized.csv"),row.names = F)

#########################################################################
### Convert pg/cell to Atoms / cell and compare to published datasets ###
#########################################################################

AtMass <- read.csv(paste0(acc_file_dir,"/atomicmasses_metals.csv"),stringsAsFactors = F)

Sample_data_n <- merge(Sample_data_n,AtMass,by.x="element_measured",by.y="Element")%>%
                 mutate(Atoms_percell = (6.022*10^23 /AtomicMass)* (pg_percell_BC*10^-12) )

Atoms_percell_AllEle <- Sample_data_n %>%
                        filter(element_perturbed =="AllEle" )%>%
                        dplyr::select(Atoms_percell,element_measured,normalization,element_perturbed)%>%
                        group_by(element_measured,normalization)%>%
                        summarise(Mean_atoms_percell = mean(Atoms_percell,na.rm=T),
                                  SD_atoms_percell = sd(Atoms_percell,na.rm=T))%>%
                        ungroup()


colnames(Atoms_percell_AllEle)[which(colnames(Atoms_percell_AllEle)=="normalization")] <- "Study"

## Read in supplementary Data File 3 from CoFactor Yeast - has atoms/cell from various studies

published_atoms_percell <- vector()

for(i in 2:9){
  
  ds <- read_xlsx(paste0(published_dataset_dir,"/metallomics/CofactorYeast_SupplementaryDataset3.xlsx"),sheet = i)
  sn <- colnames(ds)[1]  
  ds <- ds[,c(1,which(colnames(ds)=="atom/cell"))]
  colnames(ds) <- c("element_measured","Mean_atoms_percell")
  ds$Study <- sn
  
  published_atoms_percell <- rbind(published_atoms_percell,ds)
}

published_atoms_percell$SD_atoms_percell <- NA
published_atoms_percell <- filter(published_atoms_percell, element_measured %in% unique(Atoms_percell_AllEle$element_measured))

Atoms_percell_AllEle <- rbind(Atoms_percell_AllEle,published_atoms_percell[,colnames(Atoms_percell_AllEle)])%>%
                        mutate( Study = ifelse(Study =="Phosphorus","This Study - P normalized",
                                               ifelse(Study =="Unnormalized","This Study - Unnormalized",
                                                      Study)),
                                Pub_or_EPP = ifelse(Study =="This Study - P normalized","This Study",
                                                 ifelse(Study == "This Study - Unnormalized","This Study",
                                                        "Other Study")))

pdf(paste0(plot_dir,"/atomspercell_comparison_with_published_data.pdf"),width=10,height=8)
ggplot(Atoms_percell_AllEle,
       aes(x=element_measured,
           y = log10(Mean_atoms_percell),
           colour = Study,
           group = Pub_or_EPP))+
  geom_point(aes(shape=Pub_or_EPP),
             size=3,alpha=0.8,
             position = position_jitterdodge(jitter.width = 0.3))+
  scale_colour_manual(values = c(brewer.pal(8,"Dark2"), "black","darkred"))+
  theme_metallica()
dev.off()

pdf(paste0(plot_dir,"/atomspercell_comparison_with_published_data_facet.pdf"),width=10,height=12)
ggplot(Atoms_percell_AllEle,
       aes(x = element_measured,
           y = Mean_atoms_percell,
           colour = Study,
           group = Pub_or_EPP))+
  facet_wrap("element_measured",scales="free")+
  geom_point(aes(shape = Pub_or_EPP),
             size = 3,
             alpha = 0.8,
             position = position_jitterdodge(jitter.width = 0.3))+
  scale_colour_manual(values = c(brewer.pal(8,"Dark2"), "black","darkred"))+
  theme_metallica()+
  theme(strip.text = element_text(size=13,colour="darkblue"),
        strip.background = element_rect(fill="white"))
dev.off()

############################
### Final Visualisations ###
############################

# Read in correct metal concentration axis 

meas_conc_axis <- read.csv(paste0(outputtables_dir,"/measured_element_concentration_axis.csv"),
                           stringsAsFactors = F)

colnames(meas_conc_axis) <- c("element_perturbed","BioSpecID","rel_element_concentration_th","rel_element_concentration_actual")

metallomics <- merge(Sample_data_finalised, meas_conc_axis,by.x=c("element_perturbed",
                                                               "rel_element_concentration_th"),
                                   by.y = c("element_perturbed",
                                            "rel_element_concentration_th"), all.x = T)
               
AE_ngperwell_vals <- Sample_data_finalised%>%
                     filter(element_perturbed == "AllEle")%>%
                     group_by(element_measured)%>%
                     dplyr::summarise(mean_ng_perwell_BC_AE = mean(ng_perwell_BC,na.rm=T))

metallomics_AEnorm <- merge(metallomics,AE_ngperwell_vals[,c("element_measured","mean_ng_perwell_BC_AE")],by="element_measured")%>%
                      mutate(Percentage_Change_ngpwell = (ng_perwell_BC-mean_ng_perwell_BC_AE)/mean_ng_perwell_BC_AE,
                             Ratio_to_AEngperwell = ng_perwell_BC/mean_ng_perwell_BC_AE,
                             BioSpecID = gsub("_"," ",BioSpecID))%>%
                      na.omit()

buffering_df_allpoints <- metallomics_AEnorm %>%
                          filter(element_perturbed == element_measured)%>%
                          group_by(BioSpecID, element_measured, rel_element_concentration_th,rel_element_concentration_actual)%>%
                          summarise(mean_Ratio_to_AEngperwell = mean(Ratio_to_AEngperwell,na.rm=T))%>%
                          ungroup()%>%
                          mutate(
                                 BioSpecID = factor(BioSpecID,
                                                            levels =
                                                              c(  "K 0.1","K 0.2","K 0.5","K 2","K 5","K 10",
                                                                  "Mg 0.05","Mg 0.1","Mg 0.2","Mg 0.5","Mg 2","Mg 5","Mg 10",
                                                                  "Fe 0","Fe 0.01","Fe 0.02","Fe 0.05","Fe 0.1","Fe 0.2","Fe 0.5","Fe 2","Fe 5","Fe 10","Fe 20","Fe 50",
                                                                  "Zn 0.01","Zn 0.02","Zn 0.05","Zn 0.1","Zn 0.2","Zn 0.5","Zn 2","Zn 5","Zn 10","Zn 20",
                                                                  "Cu 0","Cu 0.01","Cu 0.02","Cu 0.05","Cu 0.1","Cu 0.2","Cu 0.5","Cu 2",
                                                                  "Ca 0","Ca 0.01","Ca 0.02","Ca 0.05","Ca 0.1","Ca 0.2","Ca 0.5","Ca 2","Ca 5",
                                                                  "Mn 0","Mn 0.01","Mn 0.02","Mn 0.05","Mn 0.1","Mn 0.2","Mn 0.5","Mn 2","Mn 5","Mn 10","Mn 20","Mn 50","Mn 100",
                                                                  "Mo 0.2","Mo 0.5","Mo 2"
                                                              )))


pdf(paste0(plot_dir,"/metal_buffering.pdf"),width=13,height=6)

ggplot(buffering_df_allpoints,
       aes(x=BioSpecID,
           colour = BioSpecID,
           group = BioSpecID))+
  geom_segment( aes(x=BioSpecID, xend=BioSpecID, y=log2(rel_element_concentration_actual), 
                    yend=log2(mean_Ratio_to_AEngperwell)),
                alpha=0.8)+ 
  geom_point( aes(x=BioSpecID, y=log2(rel_element_concentration_actual)), size=3.5,alpha=0.8) +
  geom_point( aes(x=BioSpecID, y=log2(mean_Ratio_to_AEngperwell)),size=3,shape=8,colour="black") +
  scale_color_manual(values = colkey_BioSpecID)+
  theme_metallica()+
  labs(y="log2(Concentration relative to AE) \n circle - environmental \n star - intracellular)",
       x= "")+
  theme( panel.grid.major = element_line(size=0.15,colour="lightgray"), 
         panel.grid.minor = element_line(size=0.15,colour="lightgray"),
         legend.position = "none",
         axis.text.x = element_text(angle =90,hjust = 0)   )

dev.off()

####################################################
### write data with replicates for metallica app ###
####################################################

for_metallica_app <- metallomics_AEnorm%>%
                     dplyr::select(BioSpecID, element_perturbed,rel_element_concentration_actual,
                                   element_measured,ng_perwell_BC,mean_ng_perwell_BC_AE,Ratio_to_AEngperwell
                                   )
write.csv(for_metallica_app,paste0(metallica_app_dir,"/metallomica_app_metallomicsdata.csv"),row.names = F)

