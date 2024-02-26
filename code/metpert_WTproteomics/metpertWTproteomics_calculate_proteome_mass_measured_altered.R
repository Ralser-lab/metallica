#`---
#`  Title: "metpertWTproteomics_calculate_proteome_mass_measured_altered.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 4 Feb 2024 
#`  Description: Script to calculate what % of the proteome by mass we can quantify and what % is differentially abundant

#############################################
### source paths functions and libraries  ###
#############################################
# general
source("/Users/aulakhs/Documents/RalserLab/metallica/code/common_code/initialise_common_paths.R")
source(paste0(code_dir,"/common_code/graphics_parameters.R"))


plot_dir <- paste0(metpert_WTproteomics_dir,"/output/plots/proteome_mass_quant_da")
dir.create(plot_dir, recursive = T)

outputtables_dir <- paste0(metpert_WTproteomics_dir,"/output/tables/proteome_mass_quant_da")
dir.create(outputtables_dir, recursive = T)

################################################################
### Read in protein expression data () and protein mass data ###
################################################################

## copy number data

prot_copynum <- read.csv(paste0(db_dir,"/protein_copy_number_scerevisiae_ghaemmaghami_2003.csv"),stringsAsFactors = F)%>%
                filter(Observed.Expression != "None" & !grepl("#",Protein.Molecules.Cell) & !grepl("%", Protein.Molecules.Cell))%>%
                mutate(Protein.Molecules.Cell = as.numeric(Protein.Molecules.Cell))%>%
                dplyr::select(ORF, Protein.Molecules.Cell)%>%
                unique()%>%
                na.omit()

# mass
prot_mass <- read.table(paste0(db_dir, "/uniprotkb_scerevisiae_proteinmass_2024_02_04.tsv"), 
                        sep = "\t", 
                        header = TRUE, 
                        fill = TRUE)

# merge

prot_mass_cpn <- merge(na.omit(unique(prot_mass[,c("Gene.Names..ordered.locus.","Mass")])),
                       prot_copynum[,c("ORF","Protein.Molecules.Cell")],
                       by.x = ("Gene.Names..ordered.locus."),
                       by.y = ("ORF"))%>%
                 mutate(cpn_weighted_mass = Protein.Molecules.Cell*Mass,
                        frac_of_total_mass = cpn_weighted_mass/sum(cpn_weighted_mass))

##############################################################
### Read in proteomics measurements and stat analysis data ###
##############################################################

proteome_quant_85pcnt <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/completematrix/allmetals/allsamples/Protein Quantities wImpValues FullMatrix.csv"))%>%
                         dplyr::select(Protein.Ids)%>%
                         unique()%>%
                         na.omit()%>%
                         mutate(ORF = as.character(lapply(Protein.Ids, convert_Uniprot2singleORF)))%>%
                         pull(ORF)

# Read the data
extra_sig_ORFs <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/extracellularconc_vs_proteinabundance/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),stringsAsFactors = F) %>%
  dplyr::select(Genes,ORF ,Significant )%>%
  filter(Significant==1)%>%
  pull(ORF)%>%
  unique()


intra_sig_ORFs <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/intracellularconc_vs_proteinabundance/lmfit_intracellular_metalconc_vs_log2PQ_with_SigNotSig_AdjPVthresh_0.05_fcthresh_1.5.csv"),stringsAsFactors = F) %>%
  dplyr::select(Genes,ORF,Significant )%>%
  filter(Significant==1)%>%
  pull(ORF)%>%
  unique()

all_sig_ORFs <- unique(c(extra_sig_ORFs,intra_sig_ORFs))
length(all_sig_ORFs)

all_meas_ORFs <- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/extracellularconc_vs_proteinabundance/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),stringsAsFactors = F) %>%
                 pull(ORF)%>%
                 unique()

################################
### Proteome mass quantified ###
################################

prot_mass_quant_any <- prot_mass_cpn%>%
                       filter(Gene.Names..ordered.locus. %in% all_meas_ORFs)

print(paste("proteome mass quantified:",sum(prot_mass_quant_any$frac_of_total_mass)))

prot_mass_quant_85pcnt <- prot_mass_cpn%>%
  filter(Gene.Names..ordered.locus. %in% proteome_quant_85pcnt)

print(paste("proteome mass quantified 85% completeness :",sum(prot_mass_quant_85pcnt$frac_of_total_mass)))


#################################
### Proteome mass significant ###
#################################


prot_mass_sig <- prot_mass_cpn%>%
  filter(Gene.Names..ordered.locus. %in% all_sig_ORFs)

print(paste("proteome mass significant:",sum(prot_mass_sig$frac_of_total_mass)))

## adjust for what was measured

prot_mass_sig_adj <- prot_mass_cpn%>%
  filter(Gene.Names..ordered.locus. %in% all_meas_ORFs)%>%
  mutate(frac_of_total_mass = cpn_weighted_mass/sum(cpn_weighted_mass))%>%
  filter(Gene.Names..ordered.locus. %in% all_sig_ORFs)

print(paste("proteome mass significant, adj for measured only:",sum(prot_mass_sig_adj$frac_of_total_mass)))


############################################
### Classify proteins as highly abundant ###
############################################

ggplot(prot_mass,
       aes(x = Mass))+
  geom_histogram(bins = 300)+
  theme_metallica()
