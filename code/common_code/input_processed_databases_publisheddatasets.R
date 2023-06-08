
#`---
#`  Title: "input_processed_databases_publisheddatasets.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 7 April 2023
#`  Description: inputs all databases and published genome wide datasets used for functional annotation
#`---


#############################################
### source paths functions and libraries  ###
#############################################

setwd(db_dir)

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
MetPDB_biosig <- filter(MetPDB,term %in% c("Ca","Cu","Fe","K","Mg","Mn","Mo","Na","Zn"))

gsets<-list(GOslim_BP,GOslim_CC,GOslim_MF,
            GO_gset_BP,GO_gset_CC,GO_gset_MF,
           # Reactome,
            KEGG,
            #SGDPhen,
            EC
            #Sc_GEM,
            )

gsetnames<-c("GOslim_BP","GOslim_CC","GOslim_MF",
             "GO_BP","GO_CC","GO_MF",
             #"Reactome",
             "KEGG",
             #"SGDPhen",
             "EC"
             #"Sc_GEM",
             )

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



# Relevant published datasets

StressResp_Brauer2007 <- read.csv("Brauer2007_SupplementaryTable1.csv",stringsAsFactors = F)

# Read in Common Metal Responsive Genes from Jin et Al 2008
CMRgenes_Jin2008 <- read.csv("Jin2008_CMRgenes.csv",stringsAsFactors = F)
