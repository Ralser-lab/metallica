#############################################################
####   Comparison with published metal depletion datasets ###
#############################################################


# general
source("/Users/aulakhs/Documents/Ralser Lab/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/layout_conversion.R"))


plot_dir <- paste0(metpert_WTproteomics_dir,"/output/plots/comparison_with_published_Fedep_Zndep_omcs_datasets")
dir.create(plot_dir, recursive = T)

outputtables_dir <- paste0(metpert_WTproteomics_dir,"/output/tables/comparison_with_published_Fedep_Zndep_omcs_datasets")
dir.create(outputtables_dir, recursive = T)

# Read the data
lmfitres_extracell <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/extracellularconc_vs_proteinabundance/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),stringsAsFactors = F) %>%
  dplyr::select(Genes,ORF, BioSpecID, Element, Element.Concentration,Significant, PValAdj_BH,Log2FC_vs_AE, LeastComplexModel)%>%
  unique()%>%
  mutate(gene_element = paste(Genes, Element))


lmfitres_intracell <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/intracellularconc_vs_proteinabundance/lmfit_intracellular_metalconc_vs_log2PQ_with_SigNotSig_AdjPVthresh_0.05_fcthresh_1.5.csv"),stringsAsFactors = F) %>%
  dplyr::select(Genes,ORF, Element,Significant, median_relative_intracellular_concentration,PValAdj_BH,median_log2_foldchangevsAE, LeastComplexModel)%>%
  mutate(gene_element = paste(Genes, Element))

######################
### Iron depletion ###
######################

Fe_depletion_results <- lmfitres_extracell%>%
                        dplyr::filter(Element == "Fe")%>%
                        dplyr::select(Genes,Element, Element.Concentration, Significant, Log2FC_vs_AE)%>%
                        dplyr::filter(Element.Concentration == min(Element.Concentration,na.rm = T))%>%
                        dplyr::group_by(Genes)%>%
                        dplyr::summarise(Log2FC_depletion = mean(Log2FC_vs_AE, na.rm = T))%>%
                        ungroup()%>%
                        dplyr::mutate(Gene = Genes)%>%
                        dplyr::select(Gene,Log2FC_depletion)

Fe_depletion_results$Dataset = "metal perturbation"


## Shakoury-Elizeh and Puig 

Fe_depletion_ShakPuig <- read.csv(paste0(published_dataset_dir,"/metal_perturbation_omics/Published data iron depletion omics.csv"), stringsAsFactors = F)%>%
                   dplyr::mutate(Log2FC_depletion = log2(FoldChange_vs_Control))%>%
                   dplyr::select(Gene,Log2FC_depletion,Dataset)
                   
## Navarrette-Perea

Fe_depletion_NP <- read.csv(paste0(published_dataset_dir,"/metal_perturbation_omics/BPStreatedtimecourse_foldchanges_NavarrettePerea2021.csv"),stringsAsFactors = F)%>%
                   dplyr::select(Gene.Symbol,FC_12h)%>%
                   dplyr::mutate(Log2FC_depletion = log2(FC_12h),
                                 Gene = Gene.Symbol)%>%
                   dplyr::select(Gene, Log2FC_depletion,)%>%
                   mutate(Dataset = "Navarrete-Perea")


Fe_depletion_all_ds_res <- rbind(Fe_depletion_ShakPuig,Fe_depletion_results,Fe_depletion_NP)%>%
                           as.data.frame()%>%
                           reshape2::dcast(Gene ~Dataset,value.var = "Log2FC_depletion",fun.aggregate = mean)

## Define function for calculating correlations between datasets

getCorrelations <- function(df) {
  # Select only the numeric columns from the dataframe
  df <- df[sapply(df, is.numeric)]
  
  # Get the names of the columns in the dataframe
  columnNames <- colnames(df)
  
  # Create an empty matrix to store the correlations
  correlationMatrix <- matrix(nrow = length(columnNames), ncol = length(columnNames))
  colnames(correlationMatrix) <- columnNames
  rownames(correlationMatrix) <- columnNames
  
  # Loop over the columns and calculate the correlations
  for (i in seq_along(columnNames)) {
    for (j in seq_along(columnNames)) {
      correlationMatrix[i, j] <- cor(df[[columnNames[i]]], df[[columnNames[j]]], use = "pairwise.complete.obs")
    }
  }
  
  return(correlationMatrix)
}

Fe_depletion_correlation <- getCorrelations(Fe_depletion_all_ds_res)

getCorrelationTests <- function(df) {
  # Select only the numeric columns from the dataframe
  df <- df[sapply(df, is.numeric)]
  
  columnNames <- colnames(df)
  
  # Create an empty data frame to store the correlation tests
  correlationTests <- data.frame(
    Column1 = character(),
    Column2 = character(),
    Correlation = numeric(),
    P.value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop over the columns and calculate the correlation tests
  for (i in seq_along(columnNames)) {
    for (j in seq_along(columnNames)) {
      if (i < j) {  # Avoid duplicates and self-correlations
        test <- cor.test(df[[columnNames[i]]], df[[columnNames[j]]], use = "pairwise.complete.obs")
        correlationTests <- rbind(correlationTests, data.frame(
          Column1 = columnNames[i],
          Column2 = columnNames[j],
          Correlation = test$estimate,
          P.value = test$p.value
        ))
      }
    }
  }
  
  return(correlationTests)
}

Fe_depletion_correlation_tests <- getCorrelationTests(Fe_depletion_all_ds_res)

library(corrplot)

col <- rev(corrplot::COL2("RdBu", n = 200))

library(corrplot)
col <- rev(corrplot::colorRampPalette(colors = c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(200))

pdf(paste0(plot_dir,"/Fe_depletion_correlation.pdf"), width = 6.5, height = 5.5)

corrplot(Fe_depletion_correlation, 
         method = "ellipse", 
         type = "lower", 
         title = "Fe Depletion Correlation",
         tl.col = "black", 
         tl.srt = 0, 
         col = col, 
         cl.pos = "r",
         diag = FALSE,   # switch off diagonal
         na.label = "-")

dev.off()


Fe_depletion_all_ds_res <- Fe_depletion_all_ds_res %>%
                           mutate( Colour = ifelse(`metal perturbation` > log2(1.5), "up",
                                                   ifelse(`metal perturbation` < log2(1/1.5),"down","nochange")),
                                   Colour_NP = ifelse(`Navarrete-Perea` > log2(1.5),"up",
                                                      ifelse(`Navarrete-Perea` < log2(1/1.5), "down","nochange")),
                                   Colour_ShakEl = ifelse(`Shakoury-Elizeh_2010` > log2(1.5),"up",
                                                          ifelse(`Shakoury-Elizeh_2010` < log2(1/1.5),"down","nochange")))


r2_cortest_Fe_NP_mp <- cor(Fe_depletion_all_ds_res$`metal perturbation`, Fe_depletion_all_ds_res$`Navarrete-Perea`, 
                     use = "pairwise.complete.obs")

pv_cortest_Fe_NP_mp <- cor.test(Fe_depletion_all_ds_res$`metal perturbation`, Fe_depletion_all_ds_res$`Navarrete-Perea`, 
                          use = "pairwise.complete.obs")[[3]]


pv_cortest_Fe_SE_mp <- cor.test(Fe_depletion_all_ds_res$`metal perturbation`, Fe_depletion_all_ds_res$`Shakoury-Elizeh_2010`, 
                                use = "pairwise.complete.obs")[[3]]
pv_cortest_Fe_Puig_mp <- cor.test(Fe_depletion_all_ds_res$`metal perturbation`, Fe_depletion_all_ds_res$`Puig_2005`, 
                                use = "pairwise.complete.obs")[[3]]



pdf(paste0(plot_dir,"/Fe depletion dataset comparison.pdf"),width = 6.3,height = 5)
ggplot(Fe_depletion_all_ds_res,
       aes(x = `metal perturbation`,
           y = Puig_2005))+
  geom_point(aes( colour = Colour), size = 4)+
  geom_hline(yintercept = log2(1.5),size =0.2, colour = "darkred")+
  geom_hline(yintercept = log2(1/1.5), size = 0.2, colour = "darkblue")+
  geom_vline(xintercept = log2(1.5), size = 0.2, colour = "darkred")+
  geom_vline(xintercept = log2(1/1.5), size = 0.2, colour = "darkblue")+
  geom_smooth(method = "lm", colour = "black")+
  scale_colour_manual(values = c("darkblue","gray","darkred"))+
  geom_text_repel(aes(label = Gene),size = 4)+
  theme_metallica()


ggplot(Fe_depletion_all_ds_res,
       aes(x = `metal perturbation`,
           y = `Shakoury-Elizeh_2010`))+
  geom_point(aes( colour = Colour), size = 4)+
  geom_hline(yintercept = log2(1.5),size =0.2, colour = "darkred")+
  geom_hline(yintercept = log2(1/1.5), size = 0.2, colour = "darkblue")+
  geom_vline(xintercept = log2(1.5), size = 0.2, colour = "darkred")+
  geom_vline(xintercept = log2(1/1.5), size = 0.2, colour = "darkblue")+
  geom_smooth(method = "lm", colour = "black")+
  scale_colour_manual(values = c("darkblue","gray","darkred"))+
  geom_text_repel(aes(label = Gene),size = 4)+
  theme_metallica()



ggplot(Fe_depletion_all_ds_res,
       aes(x = Puig_2005,
           y = `Shakoury-Elizeh_2010`
           ))+
  geom_point(aes(colour = Colour_ShakEl), size = 4)+
  geom_hline(yintercept = log2(1.5),size =0.2, colour = "darkred")+
  geom_hline(yintercept = log2(1/1.5), size = 0.2, colour = "darkblue")+
  geom_vline(xintercept = log2(1.5), size = 0.2, colour = "darkred")+
  geom_vline(xintercept = log2(1/1.5), size = 0.2, colour = "darkblue")+
  scale_colour_manual(values = c("darkblue","gray","darkred"))+
  geom_text_repel(aes(label = Gene),size = 4)+
  theme_metallica()


ggplot(Fe_depletion_all_ds_res,
       aes(x = `metal perturbation`,
           y = `Navarrete-Perea`))+
  geom_point(aes(colour = Colour), size = 4)+
  geom_smooth(method = "lm", colour = "black")+
  geom_hline(yintercept = log2(1.5),size =0.2, colour = "darkred")+
  geom_hline(yintercept = log2(1/1.5), size = 0.2, colour = "darkblue")+
  geom_vline(xintercept = log2(1.5), size = 0.2, colour = "darkred")+
  geom_vline(xintercept = log2(1/1.5), size = 0.2, colour = "darkblue")+
  scale_colour_manual(values = c("darkblue","gray","darkred"))+
 # labs(title = paste("p-value:",pv_cortest_Fe_NP_mp,"\nr2:",round(r2_cortest_Fe_NP_mp,2)))+
  theme_metallica()


ggplot(Fe_depletion_all_ds_res,
       aes(x = `Shakoury-Elizeh_2010`,
           y = `Navarrete-Perea`))+
  geom_point(aes(colour = Colour_NP),size = 4)+
  geom_hline(yintercept = log2(1.5),size =0.2, colour = "darkred")+
  geom_hline(yintercept = log2(1/1.5), size = 0.2, colour = "darkblue")+
  geom_vline(xintercept = log2(1.5), size = 0.2, colour = "darkred")+
  geom_vline(xintercept = log2(1/1.5), size = 0.2, colour = "darkblue")+
  scale_colour_manual(values = c("darkblue","gray","darkred"))+
  geom_text_repel(aes(label = Gene,colour = Colour_NP),size = 4)+
  theme_metallica()


ggplot(Fe_depletion_all_ds_res,
       aes(x = `Puig_2005`,
           y = `Navarrete-Perea`))+
  geom_point(aes(colour = Colour_NP),size = 4)+
  geom_hline(yintercept = log2(1.5),size =0.2, colour = "darkred")+
  geom_hline(yintercept = log2(1/1.5), size = 0.2, colour = "darkblue")+
  geom_vline(xintercept = log2(1.5), size = 0.2, colour = "darkred")+
  geom_vline(xintercept = log2(1/1.5), size = 0.2, colour = "darkblue")+
  scale_colour_manual(values = c("darkblue","gray","darkred"))+
  geom_text_repel(aes(label = Gene,colour = Colour_NP),size = 4)+
  theme_metallica()


dev.off()

mismatch_NP_metpert <- filter(Fe_depletion_all_ds_res,
                              `metal perturbation` < log2(1/1.5) & `Navarrete-Perea` > log2(1.5) |
                              `metal perturbation` > log2(1.5) & `Navarrete-Perea` < log2(1/1.5)
                             )

write.csv(mismatch_NP_metpert,paste0(outputtables_dir,"/mismatches navarette perea and Fe depletion metal perturbation.csv"))

######################
### Zinc depletion ###
######################


Zn_depletion_Wang2018 <- read.csv(paste0(published_dataset_dir,"/metal_perturbation_omics/Zn_depletion_foldchanges_WangEide2018.csv"), stringsAsFactors = F)%>%
                          dplyr::select(Gene.names,Fold.Change.4.hrs,Fold.Change.8.hrs,
                                        Fold.Change.12.hrs,Fold.Change.16.hrs)%>%
                          reshape2::melt(id.vars = "Gene.names", variable.name = "Dataset",
                                         value.name = "Log2FC_depletion")%>%
                          filter(Dataset == "Fold.Change.12.hrs")%>%
                          mutate(Gene = Gene.names)%>%
                          dplyr::select(Gene,Log2FC_depletion,Dataset)%>%
                          mutate(Dataset = "Wang 2018 FC 12h")
            

Zn_depletion_results <- lmfitres_extracell%>%
                        dplyr::filter(Element == "Zn")%>%
                        dplyr::select(Genes,Element, Element.Concentration, Significant, Log2FC_vs_AE)%>%
                        dplyr::filter(Element.Concentration == min(Element.Concentration,na.rm = T))%>%
                        dplyr::group_by(Genes)%>%
                        dplyr::summarise(Log2FC_depletion = mean(Log2FC_vs_AE, na.rm = T))%>%
                        ungroup()%>%
                        dplyr::mutate(Gene = Genes)%>%
                        dplyr::select(Gene,Log2FC_depletion)

Zn_depletion_results$Dataset = "metal perturbation"


Zn_depletion_all_ds_res <- rbind(Zn_depletion_Wang2018,Zn_depletion_results)%>%
                           as.data.frame()%>%
                           reshape2::dcast(Gene ~Dataset,value.var = "Log2FC_depletion",fun.aggregate = mean)%>%
                           na.omit()%>%
                           mutate( Colour = ifelse(`metal perturbation` > log2(1.5), "up",
                                                     ifelse(`metal perturbation` < log2(1/1.5),
                                                            "down","nochange")))

r2_cortest_zn <- cor(Zn_depletion_all_ds_res$`metal perturbation`, Zn_depletion_all_ds_res$`Wang 2018 FC 12h`, 
                     method = "pearson", 
                         use = "pairwise.complete.obs")

pv_cortest_zn <- cor.test(Zn_depletion_all_ds_res$`metal perturbation`, Zn_depletion_all_ds_res$`Wang 2018 FC 12h`,
                          method = "pearson", 
                          use = "pairwise.complete.obs")[[3]]


pdf(paste0(plot_dir,"/Zn depletion dataset comparison.pdf"),width = 6.3,height = 5)
ggplot(Zn_depletion_all_ds_res,
       aes(x = `metal perturbation`,
           y =`Wang 2018 FC 12h`))+
  geom_point(aes(colour = Colour))+
  geom_smooth(method = "lm", colour = "black")+
  geom_hline(yintercept = log2(1.5),size =0.2, colour = "darkred")+
  geom_hline(yintercept = log2(1/1.5), size = 0.2, colour = "darkblue")+
  geom_vline(xintercept = log2(1.5), size = 0.2, colour = "darkred")+
  geom_vline(xintercept = log2(1/1.5), size = 0.2, colour = "darkblue")+
  scale_colour_manual(values = c("darkblue","gray","darkred"))+
  labs(title = paste("p-value:",pv_cortest_zn,"\nr2:",round(r2_cortest_zn,2)))+
  theme_metallica()
dev.off()

mismatch_Wang_metpert_Zn<- filter(Zn_depletion_all_ds_res,
                              `metal perturbation` < log2(1/1.5) & `Wang 2018 FC 12h` > log2(1.5) |
                                `metal perturbation` > log2(1.5) & `Wang 2018 FC 12h` < log2(1/1.5)
)

write.csv(mismatch_NP_metpert,paste0(outputtables_dir,"/mismatches wang 2018 perea and Zn depletion metal perturbation.csv"))

