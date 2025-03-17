
#`---
#`  Title: "create_database_mappings.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 7 April 2023
#`  Description: create ORF to database mappings or note down where they were downloaded from 
#`---


#############################################
### source paths functions and libraries  ###
#############################################

# general
source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")


#############
### KEGG  ###
#############

# Function to fetch pathway name for a given pathway ID
get_pathway_name <- function(pathway_id) {
  Sys.sleep(0.05) # Delay for 1 second between requests
  tryCatch({
    pathway_info <- keggGet(pathway_id)
    pathway_name <- pathway_info[[1]]$NAME
    return(pathway_name)
  }, error = function(e) {
    return(NA) # Return NA in case of an error
  })
}

library(KEGGREST)

# Assuming you already have the keggLink data
kegg_pathways <- keggLink("pathway", "sce")

# Convert the kegg_data into a dataframe
kegg_df <- data.frame(ORF = names(kegg_pathways), KEGG_Term_ID = kegg_pathways, stringsAsFactors = FALSE)%>%
           mutate(KEGG_Term_Name = as.character(lapply(KEGG_Term_ID, get_pathway_name)),
                  ORF = gsub("sce:","",ORF),
                  KEGG_Term_ID = gsub("path:","",KEGG_Term_ID),
                  KEGG_Term_Name = gsub(" - Saccharomyces cerevisiae \\(budding yeast\\)$", "", KEGG_Term_Name))
                  
write.csv(kegg_df, paste0(db_dir,"ORF2KEGG_new.csv"),row.names = F)



##################################
### Gene Ontology - full GO db ###
##################################
library(AnnotationDbi)
library(org.Sc.sgd.db)
library(GO.db)

# Get the mapping of ORFs to GO terms
orf_to_go <- as.data.frame(AnnotationDbi::select(org.Sc.sgd.db,
                                                 keys = keys(org.Sc.sgd.db),
                                                 columns = c("GO", "ORF"),
                                                 keytype = "ORF"))

# Function to get GO term names
get_go_term_name <- function(go_id) {
  term <- tryCatch({
    Term(GOTERM[[go_id]])
  }, error = function(e) NA)
  return(term)
}

# Apply the function to get GO term names
orf_to_go$GOTermName <- sapply(orf_to_go$GO, get_go_term_name)

# View the DataFrame with GO term names
orf_to_go <- na.omit(orf_to_go)%>%
             dplyr::select(ORF,ONTOLOGY,GOTermName)%>%
             unique()

colnames(orf_to_go)<- c("ORF","ontology","term")
write.csv(orf_to_go, paste0(db_dir,"/ORF2GO_all_gsets.csv"), row.names = F)

