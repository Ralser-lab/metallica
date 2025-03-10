#`---
#`  Title: "metpertWTgrowth_1_analyze_WTgrowth_WTgrowthdata.R"
#`  @author: Johannes Hartl and Simran Kaur Aulakh
#`  Date created: 8 April 2021 
#`  Description: Script to analyze all ICP-MS WTgrowthdata from intracellular metallomics measurements
#`---

#############################################
### source paths functions and libraries  ###
#############################################

# general
source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/layout_conversion.R"))

# experiment specific
source(paste0(code_dir,"/metpert_WTgrowth/metpertWTgrowth_0_libraries_functions.R"))

#################
### set paths ###
#################

# inputs
WTgrowthdata_dir <- paste0(proj_dir,"/experiment_data/metpert_WTgrowth/data_from_TECAN")

# outputs
plot_dir <- paste0(proj_dir,"/experiment_data/metpert_WTgrowth/output/plots")
outputtables_dir <- paste0(proj_dir,"/experiment_data/metpert_WTgrowth/output/tables")

dir.create(plot_dir, recursive = T)
dir.create(outputtables_dir, recursive = T)

#######################################
### read in WTgrowthdata from TECAN ###
#######################################

filename <- "metpertWTgrowthcurves_TECAN_data_BerlinLabMethod_20200819.xlsx"

### sheet names should follow the rules:
# note: tecan can only remove top plate from analysis; bottom plate (also empty) is acquired, call "empty" to remove

sheetnames <- c("empty",
                "SM0_P1-2", "SM0_P2-2", "SM0_P3-2", "SM0_P4-2","SM0_P5-5","SM0_P6-2", "SM0_P7-4", "SM0_P8-4",
                "AE_P1-2", "AE_P2-2", "AE_P3-2", "AE_P4-2","AE_P5-2","AE_P6-2", "AE_P7-4", "AE_P8-4"
               )
### any sheet called "empty" will be removed


num_plates <- length(sheetnames)

WTgrowthdata <- lapply(1:num_plates, function(x){ #1:n for n plates (sheets)
  
  x <<-x
  temp <- readxl::read_xlsx(paste0(WTgrowthdata_dir,"/",filename),
                            sheet = x,col_types = "text")
  
  starts <- which(grepl("^Cycles",as.character(unlist(temp[,1])))) # get all row numbers containing "Cycles"
  starts <- setNames(starts, nm = as.character(unlist(temp))[starts+1]) # get all plate row idx (found in column+1); with row idx  
  
  temp <- temp[starts+4,] %>% # get mean values (in Cycles idx+4)
          `colnames<-`(temp[starts[1]+2,]) %>%
          dplyr::rename("Pos" = "Time [s]") %>%
          dplyr::mutate(Pos = as.character(unlist(temp[starts+1,1]))) %>%
          #  tidyr::separate(Pos, c("Row", "Column"),sep = 1) %>%
          tidyr::gather("Time","OD600",-c(1)) %>%
          tidyr::spread(Pos,OD600) # %>%
        #    dplyr::mutate_all(as.numeric)
  
  return(temp)
  
})


WTgrowthdataset.name <- setNames(1:num_plates, nm = sheetnames) # set names of different sheets; use to identify; dimension 1:num_plates


WTgrowthdata <- lapply(names(WTgrowthdataset.name), function(x){
                WTgrowthdata[[WTgrowthdataset.name[x]]] %>%
                dplyr::mutate(WTgrowthdataset = x) %>%
                dplyr::select(c("WTgrowthdataset","Time",2:(ncol(.)-1)))
              })

WTgrowthdata <- do.call("rbind",WTgrowthdata)

WTgrowthdata <- WTgrowthdata %>%
                dplyr::mutate(Time = as.numeric(Time)) %>%
                dplyr::arrange(WTgrowthdataset, Time) %>%
                tidyr::gather("position","OD600",-c(1:2)) %>%
                dplyr::mutate_at(vars("OD600", "Time"), as.numeric) %>% 
                dplyr::mutate(Time = Time/3600.0)

### remove empty plate
WTgrowthdata <- WTgrowthdata %>% dplyr::filter(WTgrowthdataset != "empty")

WTgrowthdata$WTgrowthdataset <- gsub("-[0-9]","",WTgrowthdata$WTgrowthdataset)
WTgrowthdata <- WTgrowthdata%>%
                separate(WTgrowthdataset,into=c("preculture","plate"))%>%
                mutate(plate = gsub("P","",plate))%>%
                mutate(plate_position = paste(plate,position,sep="_"))

# Input layouts
layout_allplates <- read.csv(paste0(acc_file_dir,"/all_layouts_combined.csv"), stringsAsFactors = F)

WTgrowthdata$SampleName <- layout_allplates[match(WTgrowthdata$plate_position,layout_allplates$plate_position),"condition"]

## filter out "empty wells as many of these had old media and were just replaced in the LO files to avoid measuring ##

WTgrowthdata <- WTgrowthdata%>%
                filter(SampleName != "Empty")%>%
                filter(!grepl(" o",SampleName))%>% # old Fe media -- contaminated
                filter(SampleName != "Ca 10")%>%
                filter(!grepl("EGTA",SampleName))%>%
                filter(!(plate == 8 & position %in% c("B4","B8","D2","D5","G6")))%>%
                filter(!(position %in% c("C1","D1","E1","F1")))%>%
                filter(!(grepl("Na",SampleName) & plate ==5))%>%
                filter(SampleName != "AllEleprepQC")%>%
                filter(SampleName != "Mo 20")%>%
                filter(Time <= 36)%>%
                filter(!plate_position %in% c("6_C8","7_D6","3_B10"))%>%
                mutate(SampleName = gsub("Zn 0.010", "Zn 0.01", SampleName),
                       SampleName = gsub("Zn 0.0200","Zn 0.02",SampleName))

####################################################
### Plot all growth curves individual media wise ###
####################################################

pdf(paste0(plot_dir,"/growthcurves_allmediasamples.pdf"),width=40,height=40)

ggplot(WTgrowthdata,aes(x = Time, y = OD600, colour = preculture))+
  geom_point(alpha=0.4)+
  facet_wrap("SampleName")+
  geom_vline(xintercept = 18, colour = "red")+
  geom_vline(xintercept = 24, colour = "red")+
  geom_vline(xintercept = 30, colour = "red")+
  scale_color_viridis_d(option="D",begin =0.1,end=0.45)+
  theme_metallica()

dev.off()


#########################################
### Plot all growth curves metal wise ###
#########################################

WTgrowthdata_noChel <- WTgrowthdata%>%
                       data.frame()%>%
                       filter(!grepl("[+]", SampleName) & SampleName != "AllEle New" )%>%
                       separate(SampleName, into = c("metal","concentration"), sep =" ",remove = F)%>%
                       mutate(concentration = as.numeric(concentration),
                              direction = ifelse(concentration > 1,"excess",
                                                 ifelse(concentration < 1 , "depletion",
                                                        "control")))%>%
                       filter(metal != "B")

pdf(paste0(plot_dir,"/growth_curves_metalfacet.pdf"),width = 20,height = 12)
ggplot(unique(filter(WTgrowthdata_noChel, preculture == "SM0" &
                metal != "AllEle")[,c("Time",
                                      "OD600",
                                      "SampleName",
                                      "plate_position",
                                      "metal","direction")]),
       aes(x = Time,
           y = OD600,
           colour = SampleName))+
  geom_line(aes(group = plate_position), alpha = 0.5, linewidth = 0.5)+
  geom_point(size = 0.5, alpha = 0.5)+
  facet_wrap(c("metal", "direction"),ncol = 6)+
  scale_color_manual(values = colkey_BioSpecID)+
  theme_metallica()+
  theme(legend.position = "none")
dev.off()

###################################
### Calculate average OD change ###
###################################

# add T0 OD600 to WTgrowthdataframe

for(i in 1:nrow(WTgrowthdata)){
  
  T0od = filter(WTgrowthdata, Time == 0.0 )%>%
         filter(preculture == as.character(WTgrowthdata[i,"preculture"]) &
                                           plate_position == as.character(WTgrowthdata[i,"plate_position"]))
  T0od = as.numeric(T0od[,"OD600"])
  WTgrowthdata[i,"T0OD"] <- T0od
}

WTgrowthdata$ODchange = WTgrowthdata$OD600 - WTgrowthdata$T0OD

WTgrowthdata <- WTgrowthdata%>%
                separate(SampleName, into = c("element","chelator"), sep = "[+]")%>%
                separate(element, into = c("element","element_concentration"), sep = " ")%>%
                separate(chelator, into = c("chelator","chelator_concentration"),sep=" ")%>%
                mutate(element_concentration = as.numeric(element_concentration),
                       chelator = ifelse(is.na(chelator),"None",chelator) ,
                       chelator_concentration = gsub("uM","",chelator_concentration))%>%
                mutate(element_concentration = ifelse(grepl("AllEle",element),1,element_concentration))

## element Concentrations
MaxODdf <- WTgrowthdata%>%
           filter(Time > 35)

elep <- ggplot(filter(MaxODdf, chelator =="None" & !grepl("AllEle",element)),
             aes(x=factor(element_concentration,levels = sort(unique(element_concentration))),
                 y=ODchange,colour=preculture,group=preculture))+
        geom_point(size = 2, alpha = 0.7)+
        geom_smooth(se=F)+
        facet_wrap("element",ncol=4)+
        scale_colour_viridis_d(end=0.85)+
        labs(x="element Concentration")+
        theme_metallica(base_size = 25)+
        theme(axis.text.x = element_text(size = 9))


chelPl <- ggplot(filter(MaxODdf,chelator != "None"),
               aes(x=as.numeric(chelator_concentration),
                   y=ODchange, colour=preculture,group=preculture))+
          geom_point(size=2, alpha = 0.7)+
          stat_smooth(se=F)+
          facet_wrap(c("element","chelator"),scales = "free_x",ncol=4)+
          scale_colour_viridis_d(end=0.9)+
          labs(x="chelator Concentration")+
          theme_metallica(base_size = 25)

pl <- plot_layout(nrow = 2, heights = c(2,1)) 

pdf(paste0(plot_dir,"/concentration_vs_ODchange_allsamples.pdf"), width=22, height = 20)
(elep/chelPl)+pl
dev.off()

## Write file with all background subtracted growth curve WTgrowthdata

write.csv(WTgrowthdata,paste0(outputtables_dir,"/metpertWTgrowth_curves_allsamples.csv"),row.names = F)

##############################################################################
#### Plotting and extracting growth parameters with "growthcurver" package ###
##############################################################################

## Reshape WTgrowthdataframe with all OD WTgrowthdata for growthcurver


df_full_forGC <- WTgrowthdata[,c("preculture","Time","plate_position","ODchange")]%>%
                 filter(preculture == "SM0")%>%
                 mutate(time = round(Time,0))%>%
                 group_by(time,plate_position)%>%
                 summarise(ODchange = median(ODchange))%>%
                 reshape2::dcast(time~plate_position,value.var="ODchange")%>%
                 na.omit()


## Get the fits resulting from growth curver
GC_fit_results <- growthcurver::SummarizeGrowthByPlate(df_full_forGC) 

## Merge fit results with sample annotation from WTgrowthdata

GC_fit_results <- merge(GC_fit_results,
                        unique(WTgrowthdata[,c("plate_position","element",
                                       "element_concentration","chelator",
                                       "chelator_concentration")]
                        ),
                        by.x = "sample",
                        by.y = "plate_position")%>%
                   filter(!note %in% c("questionable fit","cannot fit data"))

## results df :
# k - carrying capacity / maximum possitble OD
# r - intrinsic growth rate / max growth rate in mid log phase
# t_mid - time taken to reach 1/2 of max od
# t_gen - fastest possible generation time
# auc_l - is the area under the logistic curve obtained by taking the integral of the logistic equation,
# auc_e - is the empirical area under the curve which is obtained by summing up the area under the experimental curve from the measurements in the input WTgrowthdata. 
# If using auc_l or auc_e, make sure that you specify the parameter t_trim so that these metrics are comparable across samples or plates that were grown for different lengths of time.

# sigma - goodness of fit -- smaller the better

## Plot distributions of sigma to see how well things were fit 

ggplot(GC_fit_results,
       aes(x=sigma,
           fill= paste(element,chelator)))+
  geom_histogram(alpha=0.6,
                 bins =50)+
  facet_wrap(c("element","chelator"))+
  theme_metallica()

## check for outliers of fit with PCA 
rownames(GC_fit_results) <- GC_fit_results$sample
pca.res <- prcomp(GC_fit_results %>% select(k:sigma), center=TRUE, scale=TRUE)


# Plot the results
as_data_frame(list(PC1=pca.res$x[,1],
                   PC2=pca.res$x[,2],
                   samples = GC_fit_results$sample)) %>% 
  ggplot(aes(x=PC1,y=PC2, label=samples)) + 
  geom_text(size = 3)

## Plot growth rates (r)

# elements
pdf(paste0(plot_dir,"/growthrate_aucl_elements_nochelators.pdf"),width=25,height = 25)
print(
  ggplot(filter(GC_fit_results,chelator == "None" & element_concentration != 0),
         aes(x=log10(element_concentration),
             y=r,
             colour=log10(element_concentration)))+
    geom_point(size=4,alpha=0.7)+
    labs(title="Growth Rate", colour = "element\nconcentration")+
    facet_wrap("element",ncol=4,scales = "free_x")+
    scale_colour_viridis_c()+
    stat_smooth(colour="black")+
    theme_metallica()
)

print(
  ggplot(filter(GC_fit_results,chelator=="None" & 
                  element !="AllEleNew" & element_concentration != 0),
         aes(x = log10(element_concentration),
             y = auc_l,
             colour = log10(element_concentration)))+
    geom_point(size = 4, alpha = 0.7)+
    labs(title="auc_l", colour = "element\nconcentration")+
    facet_wrap("element",ncol = 4,scales = "free_x")+
    scale_colour_viridis_c()+
    stat_smooth(colour = "black")+
    theme_metallica()
)
dev.off()

pdf(paste0(plot_dir,"/maxOD_aucs_tmids_elements_nochelators.pdf"), width = 25, height = 25)

ggplot(filter(GC_fit_results,chelator=="None" & 
                element !="AllEleNew" & element_concentration != 0),
       aes(x=log10(element_concentration),
           y=k,
           colour=log10(element_concentration)))+
  geom_point(size=4,alpha=0.7)+
  labs(title="Maximum OD", colour = "element\nConcentration")+
  facet_wrap("element",ncol=4,scales = "free_x")+
  scale_colour_viridis_c()+
  stat_smooth(colour="black")+
  theme_metallica()

ggplot(filter(GC_fit_results,chelator=="None" & 
                element !="AllEleNew" & element_concentration != 0),
       aes(x=log10(element_concentration),
           y=t_mid,
           colour=log10(element_concentration)))+
  geom_point(size=4,alpha=0.7)+
  labs(title="tmid", colour = "element\nConcentration")+
  facet_wrap("element",ncol=4,scales = "free_x")+
  scale_colour_viridis_c()+
  stat_smooth(colour="black")+
  theme_metallica()

ggplot(filter(GC_fit_results,chelator=="None" & 
                element !="AllEleNew" & element_concentration != 0),
       aes(x=log10(element_concentration),
           y=t_gen,
           colour=log10(element_concentration)))+
  geom_point(size=4,alpha=0.7)+
  labs(title="tgen", colour = "element\nConcentration")+
  facet_wrap("element",ncol=4,scales = "free_x")+
  scale_colour_viridis_c()+
  stat_smooth(colour="black")+
  theme_metallica()

dev.off()

# chelators

chels = na.omit(unique(GC_fit_results$chelator))
chels = chels[-grep("None",chels)]

for(c in 1:length(chels)){
  
  p1 <- ggplot(filter(GC_fit_results,chelator == chels[c] ),
             aes(x=as.numeric(chelator_concentration),
                 y=r,
                 colour=as.numeric(chelator_concentration)))+
    geom_point(size=3)+
    labs(title="Growth Rate", 
         colour = "",
         x="chelator Concentration")+
    facet_wrap(c("element","chelator"),ncol=4)+
    scale_colour_viridis_c(option="C",end=0.8)+
    stat_smooth(colour="black")+
    ylim(0,0.5)+
    theme_metallica()
  
  p2 <- ggplot(filter(GC_fit_results,chelator == chels[c] ),
             aes(x=as.numeric(chelator_concentration),
                 y=auc_l,
                 colour=as.numeric(chelator_concentration)))+
          geom_point(size=3)+
          labs(title="auc_l", 
               colour = "",
               x="chelator Concentration")+
          facet_wrap(c("element","chelator"),ncol=4)+
          scale_colour_viridis_c(option="C",end=0.8)+
          theme_metallica()
  
  
  if(chels[c] =="BCS"){
    
    pdf(paste0(plot_dir,"/BCS growth rate auc_l.pdf"),width=9,height=5)
    print(p1)
    print(p2)
    dev.off()
  }else if(chels[c]=="DfP"){
    pdf(paste0(plot_dir,"DfP growth rate auc_l.pdf"),width=14,height=5)
    print(p1)
    print(p2)
    dev.off()
  }else{
    pdf(paste0(plot_dir,"/",chels[c]," growth rate auc_l.pdf"),width=5,height=5)
    print(p1)
    print(p2)
    dev.off()
  }
}

## Median summary stat across all conditions

GC_fit_results_sum <- GC_fit_results%>%
                      group_by(element,element_concentration,chelator,chelator_concentration)%>%
                      summarize(Med.Growth.Rate = median(r),
                                Med.Max.OD = median(k,na.rm = T),
                                Med.Tmid = median(t_mid,na.rm = T),
                                Med.Tgen = median(t_gen,na.rm = T),
                                Med.aucl = median(auc_l,na.rm = T),
                                Med.auce = median(auc_e,na.rm = T))%>%
                      ungroup()

## plot distributions of all growth characteristics together

GC_fit_results_sum %>%
              filter(chelator=="None" & element != "AllEleNew")%>%
              reshape2::melt(id.vars=c("element","element_concentration",
                                       "chelator","chelator_concentration"))%>%
              ggplot(aes(x=value,
                         colour=element))+
              geom_density()+
              facet_wrap("variable",scales="free")+
              theme_metallica()


GC_fit_results_sum %>%
              filter(chelator!="None" & element != "AllEleNew")%>%
              reshape2::melt(id.vars=c("element","element_concentration",
                                       "chelator","chelator_concentration"))%>%
              ggplot(aes(x=value,
                         colour=chelator))+
              geom_density()+
              facet_wrap("variable",scales="free")+
              theme_metallica()


##########################
### Heatmap of summary ###
##########################

# elements

htmp_mat <- GC_fit_results_sum %>%
            filter(chelator=="None" & !element %in% c("AllEleNew","AllEle"))%>%
            reshape2::dcast(element~element_concentration, value.var = "Med.Growth.Rate")

ae_val <- filter(GC_fit_results_sum,chelator=="None" & element == "AllEle")$Med.Growth.Rate
htmp_mat$`1` <- ae_val

rownames(htmp_mat)<-htmp_mat$element
htmp_mat<-htmp_mat[,-1]

htmp_mat <- htmp_mat[order(as.numeric(as.character(colnames(htmp_mat)))) ]

pdf(paste0(plot_dir,"/heatmaps_growthrate_elesonly.pdf"),width=20,height=10)
pheatmap(htmp_mat,
         color = viridis(20,
                         option = "C"),
         cluster_rows = F,
         cluster_cols = F,
         na_col = NA,
         main = "Growth Rate")

plot.new()
pheatmap(htmp_mat,
         color = viridis(20,
                         option = "C"),
         cluster_rows = F,
         cluster_cols = F,
         na_col = NA,
         scale = "row",
         main = "Growth Rate Scaled by element")
dev.off()


htmp_mat <- GC_fit_results_sum %>%
            filter(chelator=="None" &
                     !element %in% c("AllEleNew","AllEle"))%>%
            reshape2::dcast(element~element_concentration,
                            value.var = "Med.aucl")

ae_val <- filter(GC_fit_results_sum,chelator=="None" & element == "AllEle")$Med.aucl
htmp_mat$`1` <- ae_val

rownames(htmp_mat)<-htmp_mat$element
htmp_mat<-htmp_mat[,-1]

htmp_mat <- htmp_mat[order(as.numeric(as.character(colnames(htmp_mat)))) ]


pdf(paste0(plot_dir,"/heatmaps_aucl_elesonly.pdf"),width=20,height=10)
print(
  pheatmap(htmp_mat,
           color = viridis(20,
                           option = "D"),
           cluster_rows = F,
           cluster_cols = F,
           na_col = NA,
           main = "auc_l")
)

plot.new()
print(
  pheatmap(htmp_mat,
           color = viridis(20,
                           option = "D"),
           cluster_rows = F,
           cluster_cols = F,
           na_col = NA,
           scale = "row",
           main = "auc_l Scaled by element")
)
dev.off()

#################################################################
## Plot all growth rate data in one plot for essential metals ###
#################################################################
# read in media metallomics file 
measured_conc_axis <- read.csv(paste0(metpert_WTmetallomics_dir,"/output/tables/measured_element_concentration_axis.csv"),
                               stringsAsFactors = F)

GC_fit_results_sum$BioSpecID <- paste(GC_fit_results_sum$element,
                                      GC_fit_results_sum$element_concentration)

GC_fit_results_sum_test <- merge(GC_fit_results_sum,measured_conc_axis, by = "BioSpecID")%>%
                           filter(chelator == "None")

pdf(paste0(plot_dir,"/log2_metalconc_vs_medgrowthrate.pdf"),width = 6,height = 3.5)
ggplot(GC_fit_results_sum_test,
       aes(x = log2(actual_AEnorm_conc),
           y = Med.Growth.Rate,
           colour = element,
           fill = element))+
  geom_smooth(se = T, linewidth = 0.5, alpha = 0.2)+
  geom_point(size = 2.5,alpha = 0.4)+
  geom_vline(xintercept = log2(5), colour = "red")+
  geom_vline(xintercept = -log2(5),colour ="red")+
  scale_colour_manual(values = colkey_Ele)+
  scale_fill_manual(values = colkey_Ele)+
  theme_metallica()+
  labs()
dev.off()
############################################
### Write summarized growth data to disk ###
############################################

write.csv(GC_fit_results_sum_test,paste0(outputtables_dir,"fitted_growthparams_summary.csv"),
          row.names = F)

write.csv(GC_fit_results,paste0(outputtables_dir,"fitted_growthparams_perreplicate.csv"),
          row.names = F)

##########################################################
### Write summarised data - small df for metallica app ###
##########################################################

for_metallica_app <- GC_fit_results %>%
               filter(chelator == "None", !element %in% c("B","AllEle"))%>%
               mutate(BioSpecID = paste(element,element_concentration))%>%
               dplyr::select(BioSpecID,element, element_concentration, r)%>%
               unique()

colnames(for_metallica_app) <- c("BioSpecID","metal","env_metal_concentration","growth_rate")

write.csv(for_metallica_app, paste0(metallica_app_dir,"/metallica_app_metpertWTgrowthrate.csv"),row.names = F)
