
#`---
#`  Title: "database_identifier_conversion_functions.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 5 April 2023
#`  Description: reads in uniprot exported converter and contains functions for all ORF <> UniprotID <> Gene/Protein Name conversions
#`---

#####################################
### Set Paths used in all scripts ###
#####################################

# general
source("/Users/aulakhs/Documents/Ralser Lab/metallica/code/common_code/initialise_common_paths.R")
library(dplyr)
#######################################
### Read in uniprot conversion file ###
#######################################


GenProt_SGD <- read.csv(paste0(db_dir,"/uniprot_allScerevisiae_20230208.tsv"), sep="\t",stringsAsFactors = F)[c("Entry","Gene.Names..ordered.locus.","Gene.Names..primary.","Annotation")]%>%
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

##########################################################################
### Functions to convert between ORF <> UniProtID <> Gene/Protein Name ###
##########################################################################

convert_Uniprot2GeneName <- function(unp_id){
  
  uids <- unlist(strsplit(unp_id,";"))
  
  allgns <- ""
  
  for( i in 1:length(uids)){
    
    gn <- GenProt_SGD%>%
      filter(Uniprot.ID == uids[i])
    gn <- gn[["Gene.Name"]]
    
    if(length(gn) > 0 ){  
      gns <- ""
      
      for(j in 1:length(gn)){ 
        gns <- paste(gns,gn[j],sep=" ") }
      gns <- trimws(gns)
      allgns<-paste(allgns,gns,sep=" ")
      allgns <- trimws(allgns)}
    
    else{
      allgns <- unp_id
    }
    return(allgns)
  }
}

## function to convert protein group to Gene Name -- return uniprot for proteins without gene names -- and multiple gene names for protein groups with more than 1 uniprot id

convert_Uniprot2SingleGeneName <- function(unp_id){
  
  uids <- unlist(strsplit(unp_id,";"))
  
  allgns <- ""
  
  for( i in 1:length(uids)){
    
    gn <- GenProt_SGD%>%
      filter(Uniprot.ID == uids[i])
    gn <- gn[["Gene.Name"]]
    
    if(length(gn) > 0 ){  
      gnr <- gn[[1]]
    } else{
      gnr <- unp_id
    }
    return(gnr)
  }
}

convert_ORF2GeneName <- function(ORFname){
  
  gn <- GenProt_SGD%>%
    filter(ORF == ORFname)
  gn <- gn[["Gene.Name"]]
  
  if(length(gn) < 0 ){
    
    gn <- ORFname
  }
  return(gn)
}

convert_ORF2SingleGeneName<-function(ORFname){
  
  gn <- GenProt_SGD%>%
    filter(ORF == ORFname)
  gn <- unlist(strsplit(gn[["Gene.Name"]],"[;]"))
  
  if(length(gn) < 0 ){
    
    gn <- ORFname
  }else if(length(gn) > 1){
    
    gn <- trimws(gn[[1]])
  }
  return(gn)
}

convert_GeneName2ORF <- function(genename){
  
  orf <- GenProt_SGD%>%
    filter(Gene.Name == genename)%>%
    dplyr::select(Gene.Name,ORF)%>%
    unique()
  orf <- orf[["ORF"]]
  
  if(length(orf) < 0 ){
    
    orf <- genename
  }
  return(orf)
}

convert_GeneName2UNIPROT <- function(genename){
  
  unp <- GenProt_SGD%>%
    filter(Gene.Name == genename)%>%
    dplyr::select(Gene.Name,Uniprot.ID)%>%
    unique()
  unp <- unp[["Uniprot.ID"]]
  
  if(length(unp) < 0 ){
    
    unp <- genename
  }
  return(unp)
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


# Function to convert UniprotIDs to ORFs

convert_Uniprot2singleORF <- function(x){
  
  uid = unlist(strsplit(x,";"))[[1]]
  
  singleORF = unlist(filter(GenProt_SGD,Uniprot.ID==uid)$ORF[[1]])
  
  return(singleORF)
  
}

