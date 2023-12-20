
#`---
#`  Title: "input_processed_databases_publisheddatasets.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 7 April 2023
#`  Description: inputs all databases and published genome wide datasets used for functional annotation
#`---


#############################################
### source paths functions and libraries  ###
#############################################

# general
source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")

setwd(db_dir)

GO_gset=read.csv("ORF2GO_all_gsets.csv",stringsAsFactors = F)
GO_gset_MF=unique(filter(GO_gset,ontology=="MF")[,c("ORF","term")])
GO_gset_BP=unique(filter(GO_gset,ontology=="BP")[,c("ORF","term")])
GO_gset_CC=unique(filter(GO_gset,ontology=="CC")[,c("ORF","term")])

GOslim_BP=read.delim("go_slim_mapping_BP.tab")[,c(1,5)]
colnames(GOslim_BP)<-c("ORF","term")
GOslim_BP <- filter(GOslim_BP, term != "biological_process")


GOslim_CC=read.delim("go_slim_mapping_CC.tab")[,c(1,5)]
colnames(GOslim_CC)<-c("ORF","term")
GOslim_CC <- filter(GOslim_CC, term != "cellular_component")

GOslim_MF=read.delim("go_slim_mapping_MF.tab")[,c(1,5)]
colnames(GOslim_MF)<-c("ORF","term")
GOslim_MF <- filter(GOslim_MF, term != "molecular_function")
rm(GO_gset)

Reactome<-read.csv("ORF2Reactome.csv",stringsAsFactors = F)
colnames(Reactome)<-c("ORF","term")

KEGG<-read.csv("ORF2KEGG.csv",stringsAsFactors = F)[,c("ORF","KEGG.Term.Name")]
colnames(KEGG)<-c("ORF","term")

KEGG$term<-gsub("[ ][-][ ]Saccharomyces cerevisiae[ ][(]budding yeast[)]","", KEGG$term)

KEGG <- filter(KEGG,!term %in%  c("Not Found in SCe KEGG","Metabolic pathways"))

SGDPhen<-read.csv("ORF2SGDPhenotypes.csv",stringsAsFactors = F)[,c("ORF","phenotype")]
colnames(SGDPhen)<-c("ORF","term")

EC <- read.delim("uniprot2EC_uniprotkb_download_2023_09_07.tsv",stringsAsFactors = F, sep = "\t")%>%
      mutate(ORF = as.character(lapply(X = Entry, FUN = convert_Uniprot2singleORF)))%>%
  tidyr::separate_rows(EC.number, sep = "; ") %>%
  mutate(EC.number = trimws(EC.number))%>%
  dplyr::select(ORF,EC.number)
colnames(EC)<-c("ORF","term")

EC_names <- list("Oxidoreductases", # 1
                 "Transferases", # 2
                 "Hydrolases", # 3
                 "Lyases", # 4
                 "Isomerases", # 5
                 "Ligases", # 6
                 "Translocases" # 7 
)
# Create the mapping vector
EC_map <- setNames(EC_names, 1:length(EC_names))

EC_df <- EC %>%
  tidyr::separate(term, into = c("EC_level_1","EC_level_2",NA,NA), sep = "[.]")%>%
  mutate(EC_level1_name = recode(EC_level_1, !!!EC_map))%>%
  filter(EC_level_1 != "")


Sc_GEM<-read.delim("Scerevisiae_GEM.txt",stringsAsFactors = F)
colnames(Sc_GEM)<-c("ORF","term")

MetPDB <- read.csv("ORF2MetPDB.csv",stringsAsFactors = F)
MetPDB_biosig <- filter(MetPDB,term %in% c("Ca","Cu","Fe","K","Mg","Mn","Mo","Na","Zn"))

GO_gset_MF_metalbinding <- GO_gset_MF %>%
  filter(grepl("binding",term)|
           grepl("protoheme IX farnesyltransferase activity",term))%>%
  filter(grepl("calcium",term) |
           grepl("copper",term) |
           grepl("cupric",term) |
           grepl("iron",term) |
           grepl("magnesium",term) |
           grepl("manganese",term) |
           grepl("molybdenum", term) |
           grepl("potassium", term) |
           grepl("sodium", term) |
           grepl("zinc", term) |
           grepl("heme",term)|
           grepl("metal",term))%>%
  filter(!grepl("calcium-dependent",term))


GO_gset_MF_binding_metalwise <- GO_gset_MF_metalbinding%>%
  mutate(term = ifelse(grepl("calcium",term),"Ca",
                       ifelse(grepl("copper", term),"Cu",
                              ifelse(grepl("iron",term),"Fe",
                                     ifelse(grepl("magnesium",term),"Mg",
                                            ifelse(grepl("manganese",term),"Mn",
                                                   ifelse(grepl("molybdenum",term),"Mo",
                                                          ifelse(grepl("potassium",term),"K",
                                                                 ifelse(grepl("sodium",term),"Na",
                                                                        ifelse(grepl("zinc",term),"Zn",
                                                                               ifelse(grepl("heme",term),"Fe",
                                                                                      "unspecific")))))))))))%>%
  dplyr::select(ORF,term)%>%
  unique()

metalwise_metalbinding_ORFs <- unique(filter(GO_gset_MF_binding_metalwise, term != "unspecific")$ORF)

GO_gset_MF_metalbinding_unspecific <- filter(GO_gset_MF_binding_metalwise, term == "unspecific" & 
                                               !(ORF %in% metalwise_metalbinding_ORFs))%>%
                                       unique()

GO_gset_MF_binding_metalwise <- filter(GO_gset_MF_binding_metalwise, term != "unspecific")%>%
                                unique()


rm(GO_gset_MF_metalbinding)

 
##########################################
### Metal transporters in S.cerevisiae ###
##########################################

philpott_metal_transporter_genes <- list(
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
# Convert the list into a dataframe
philpott_metal_transporter_df <- stack(philpott_metal_transporter_genes)

# Rename the columns to "Gene" and "Term"
colnames(philpott_metal_transporter_df) <- c("Gene", "term")

philpott_metal_transporter_df <- philpott_metal_transporter_df%>%
                                    mutate(ORF = as.character(lapply(Gene, convert_GeneName2ORF)))%>%
                                    dplyr::select(ORF,term)

metal_transporter_philpott_ORFs <- unique(as.character(lapply(unlist(philpott_metal_transporter_genes), convert_GeneName2ORF) ))

GO_gset_MF_metaltransporter <- unique(filter( GO_gset_MF,
                                              # unspecific annotations
                                                        grepl("metal ion transmembrane transporter activity",term) |
                                              # specific - Calcium 
                                                        grepl("calcium ion transmembrane transporter activity",term)|
                                                        grepl("calcium:proton antiporter activity",term)|
                                                        grepl("calcium:sodium antiporter",term) | 
                                                        grepl("P-type calcium transporter activity",term)|
                                                        grepl("calcium channel",term)|
                                                        
                                                      # Copper
                                                        grepl("copper chaperone activity",term) |
                                                        grepl("P-type divalent copper transporter activity",term)|
                                                        grepl("copper ion transmembrane transporter activity",term)|
                                                        
                                                      # Iron 
                                                        grepl("iron ion transmembrane transporter activity",term)|                      
                                                        grepl("ferrous iron transmembrane transporter activity",term)|                  
                                                        grepl("siderophore-iron transmembrane transporter activity",term)|
                                                        grepl("siderophore uptake transmembrane transporter activity",term)|
                                                        grepl("ferric-enterobactin transmembrane transporter activity",term)|
                                                        grepl("iron chaperone activity",term)|   
                                                      
                                                      # Potassium
                                                        grepl("high-affinity potassium ion transmembrane transporter activity",term)|
                                                        grepl("P-type potassium transmembrane transporter activity",term)|
                                                        grepl("potassium channel activity",term)|
                                                        grepl("potassium ion leak channel activity",term)|
                                                        grepl("potassium ion transmembrane transporter activity",term)|
                                                        grepl("potassium:chloride symporter activity",term)|
                                                        grepl("potassium:proton antiporter activity",term)|
                                                        grepl("voltage-gated potassium channel activity",term)|
                                                        
                                                      # Magnesium
                                                        grepl("magnesium ion transmembrane transporter activity",term)|
                                                        
                                                      # Manganese 
                                                        grepl("ABC-type manganese transporter activity",term)|
                                                        grepl("manganese ion transmembrane transporter activity",term)|
                                                        
                                                      # Molybdenum -- NO ANNOTATIONS
                                                        
                                                      # Sodium
                                                        grepl("calcium:sodium antiporter activity involved in regulation of cardiac muscle cell membrane potential",term)|
                                                        grepl("P-type sodium transporter activity",term)|
                                                        grepl("sodium channel activity",term)|
                                                        grepl("sodium:inorganic phosphate symporter activity",term)|
                                                        grepl("sodium:proton antiporter activity",term)|
                                                      
                                                      # Zinc
                                                       grepl("low-affinity zinc ion transmembrane transporter activity",term)|
                                                       grepl("high-affinity zinc transmembrane transporter activity",term)|
                                                       grepl("zinc ion transmembrane transporter activity",term)
))
                                                     
                                         
  
GO_gset_MF_metaltransporter_ORFs <- unique(GO_gset_MF_metaltransporter$ORF)


GO_gset_MF_metaltransporter <- GO_gset_MF_metaltransporter%>%
                                    mutate(metal = ifelse(grepl("calcium",term),"Ca",
                                                         ifelse(grepl("copper", term),"Cu",
                                                                ifelse(grepl("iron",term),"Fe",
                                                                       ifelse(grepl("ferrous",term),"Fe",
                                                                              ifelse(grepl("ferric",term),"Fe",
                                                                                     ifelse(grepl("siderophore",term),"Fe",
                                                                                           ifelse(grepl("magnesium",term),"Mg",
                                                                                                  ifelse(grepl("manganese",term),"Mn",
                                                                                                         ifelse(grepl("molybdenum",term),"Mo",
                                                                                                                ifelse(grepl("potassium",term),"K",
                                                                                                                       ifelse(grepl("sodium",term),"Na",
                                                                                                                              ifelse(grepl("zinc",term),"Zn","unspecific")))))))))))))

## add calcium sodium as sodium tooo

ca_sod_Na_annotation <- filter(GO_gset_MF_metaltransporter,   grepl("calcium:sodium antiporter",term))%>%
                        mutate(metal = "Na")

GO_gset_MF_metaltransporter <- rbind(GO_gset_MF_metaltransporter, ca_sod_Na_annotation)

GO_gset_MF_transporter_metalwise <- filter(GO_gset_MF_metaltransporter, metal != "unspecific")%>%
                                    dplyr::select(ORF,metal)

colnames(GO_gset_MF_transporter_metalwise)<- c("ORF","term")

GO_gset_MF_metaltransporter_unspecific <- filter(GO_gset_MF_metaltransporter, metal == "unspecific" &
                                                   !(ORF %in% GO_gset_MF_transporter_metalwise$ORF))
# NO UNSPECIFIC transporters

#############################################################
### Remove metal transporters from metal binding proteins ###
#############################################################

GO_gset_MF_binding_metalwise <- GO_gset_MF_binding_metalwise%>%
                                filter( !ORF %in% unique(c(GO_gset_MF_transporter_metalwise$ORF,
                                                           philpott_metal_transporter_df$ORF)))

GO_gset_MF_metalbinding_unspecific <- GO_gset_MF_metalbinding_unspecific%>%
                                      filter( !ORF %in% unique(c(GO_gset_MF_transporter_metalwise$ORF,
                                                                 philpott_metal_transporter_df$ORF)))


################################
### Add other metal related terms that are not binding or transporter

other_metal_terms <- c( # general - unspecific
                        "metalloaminopeptidase activity","metallocarboxypeptidase activity",
                        "metallodipeptidase activity","metalloendopeptidase activity",
                        "metallopeptidase activity","oxidoreductase activity, acting on metal ions",
                        "oxidoreductase activity, acting on metal ions, flavin as acceptor",
                        "protein tyrosine phosphatase activity, metal-dependent",
                        
                        # Calcium
                       "calcium-dependent cysteine-type endopeptidase activity",
                       "calcium-dependent phospholipid binding","calcium-dependent protein binding",
                       "calcium-dependent protein kinase C activity",
                       "calcium-dependent protein serine/threonine phosphatase activity",
                       "calcium-dependent protein serine/threonine phosphatase regulator activity",
                       "calcium-independent phospholipase A2 activity",
                       "calcium-dependent protein serine/threonine kinase activity",
                       
                       # Cu - none left after binding and transport
                       
                       # Iron
                       "acireductone dioxygenase [iron(II)-requiring] activity",
                       "oxidoreductase activity, acting on iron-sulfur proteins as donors",
                       "ferredoxin hydrogenase activity","ferric-chelate reductase (NADPH) activity",
                       "ferric-chelate reductase activity","ferrochelatase activity",
                       "ferroxidase activity",
                       "sirohydrochlorin ferrochelatase activity",
                       "sulfite reductase (ferredoxin) activity",
                       
                      # Potassium -- none left after binding and transport
                      # Magnesium -- none left after binding and transport
                      # Manganese -- none left after binding and transport
                      # Sodium -- none left after binding and transport
                      
                      # Zinc
                      "alcohol dehydrogenase activity, zinc-dependent")

GO_gset_MG_other_metalrelated <- filter(GO_gset_MF,
                                        term %in% other_metal_terms)%>%
                                 mutate(term = ifelse( grepl("calcium",term),"Ca",
                                                       ifelse(grepl("iron",term),"Fe",
                                                              ifelse(grepl("ferredoxin",term),"Fe",
                                                                     ifelse(grepl("ferric",term),"Fe",
                                                                            ifelse(grepl("ferroxidase",term),"Fe",
                                                                            ifelse(grepl("ferrochelatase",term),"Fe",
                                                                                   ifelse(grepl("zinc",term),"Zn",
                                                                                          "unspecific")))))))
                                        )
                                 

## write "metal related gene" annotation for python scripts

GO_MF_all_metal_related_anno <- rbind(GO_gset_MF_binding_metalwise,
                                      GO_gset_MF_metalbinding_unspecific,
                                      GO_gset_MF_transporter_metalwise,
                                      philpott_metal_transporter_df,
                                      GO_gset_MG_other_metalrelated)%>%
  unique()

write.csv(GO_MF_all_metal_related_anno, paste0(db_dir,"/metal_binding_transport_anno.csv"),
          row.names = F)


####################################
### Make a list of all gene sets ###
####################################


gsets<-list(GOslim_BP,GOslim_CC,GOslim_MF,
            GO_gset_BP,GO_gset_CC,GO_gset_MF,
            # Reactome,
            KEGG,
            GO_gset_MF_binding_metalwise,
            #SGDPhen,
            EC
            #Sc_GEM,
)

gsetnames<-c("GOslim_BP","GOslim_CC","GOslim_MF",
             "GO_BP","GO_CC","GO_MF",
             #"Reactome",
             "KEGG",
             "GO_MF_metalspecific_metalbinding",
             #"SGDPhen",
             "EC"
             #"Sc_GEM",
)






#######################################################################
### true positive for metal specific metal binding GOMF and MetPDB ####
#######################################################################

# Relevant published datasets

StressResp_Brauer2007 <- read.csv("Brauer2007_SupplementaryTable1.csv",stringsAsFactors = F)

# Read in Common Metal Responsive Genes from Jin et Al 2008
CMRgenes_Jin2008 <- read.csv("Jin2008_CMRgenes.csv",stringsAsFactors = F)
