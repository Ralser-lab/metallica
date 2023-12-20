#`---
#`  Title: "compare_ecYeast8simulation_with_measuredproteinabundance.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 6 June 2021 
#`  Description: Script to compare results of simulation of metal depletion and excess with experimentally obtained protein quantity results
#`---


#############################################
### source paths functions and libraries  ###
#############################################

source("/Users/aulakhs/Documents/Ralser Lab/metallica/code/common_code/initialise_common_paths.R")

# general

source(paste0(code_dir,"/common_code/graphics_parameters.R"))
source(paste0(code_dir,"/common_code/input_processed_databases_publisheddatasets.R"))
source(paste0(code_dir,"/common_code/database_identifier_conversion_functions.R"))

# specific

source(paste0(code_dir,"/metpert_ecYeast8simulation/metpert_ecYeast8simulation_0_libraries_functions.R"))


plot_dir <- paste0(metpert_ecYeast8simulation_dir,"/output/plots/")
dir.create(plots_dir,recursive = T)

output_tables_dir <- paste0(metpert_ecYeast8simulation_dir,"/output/tables/")
dir.create(output_tables_dir,recursive = T)

################################################################################
### Read in experimentally measured protein abundance and simulation results ###
################################################################################

protabun_data<- read.csv(paste0(metpert_WTproteomics_dir,"/output/tables/lmfit DE res with PQ and SigNotSig AdjPVthresh 0.05 FCthresh 1.5.csv"),
                        stringsAsFactors = F)


GEMsim_0pt9_grwth <- data.frame(read_excel(paste0(metpert_ecYeast8simulation_dir,"/matlab_simulations/cofactorYeast8_GEM_simulations_from_YuChen_08062023.xlsx"),sheet = 2))

colnames(GEMsim_0pt9_grwth) <- GEMsim_0pt9_grwth[1,]

GEMsim_0pt9_grwth <- GEMsim_0pt9_grwth[-1,]

GEMsim_0pt9_grwth <- GEMsim_0pt9_grwth%>%
  ## filter out metal transport fluxes -- these were the ones that were minimized or maximized to achieve simulation 
  filter( !Reaction.Name %in% c("Ca(2+) exchange","Ca(2+) transport", ## Ca  
                                "Uptake of Cu(+)","Cu2(+) transport", "Cu2(+) exchange",
                                "iron (II) transport","	iron(3+) exchange","iron(2+) exchange",
                                "Mg(2+) exchange",
                                "Mn(2+) exchange",
                                "potassium exchange","potassium transport",
                                "sodium exchange",
                                "Zn(2+) exchange"))

                     reshape2::melt(id.vars=c("ORF","Gene"),variable.name = "ElePert",value.name = "simulated_protein_fc")%>%
                     mutate(ElePert = gsub("_"," ",ElePert),
                            ElePert = gsub("depletion","Depletion",ElePert),
                            ElePert = gsub("excess","Excess",ElePert),
                            ORF_ElePert = paste(ORF,ElePert),
                            simulated_protein_fc = as.numeric(simulated_protein_fc))
    
protabun_data_forcomparison <- protabun_data%>%
                                dplyr::select(Genes,ORF,Significant,Element,Element.Concentration,Log2FC_vs_AE)%>%
                                mutate(Perturbation = ifelse(Element.Concentration > 1.2,"Excess",
                                                             ifelse(Element.Concentration < 0.8,"Depletion",
                                                                    "Control")))%>%
                                na.omit()%>%
                                mutate(ORF_ElePert = paste(ORF,Element,Perturbation))%>%
                                group_by(Genes, Perturbation,Element)%>%
                                mutate(MaxMag_log2_experimental_protein_fc = ifelse(abs(max(Log2FC_vs_AE,na.rm=T)) > abs(min(Log2FC_vs_AE,na.rm=T)),
                                                                  max(Log2FC_vs_AE,na.rm = T),
                                                                  min(Log2FC_vs_AE,na.rm = T)),
                                       MaxMag_experimental_protein_fc = 2^(MaxMag_log2_experimental_protein_fc))%>%
                                ungroup()%>%
                                dplyr::select(ORF_ElePert,Element,MaxMag_experimental_protein_fc)%>%
                                unique()


simulation_vs_experiment <- merge(GEMsim_0pt9_grwth[,c("ORF_ElePert","ElePert","ORF","Gene","simulated_protein_fc")],
                                  protabun_data_forcomparison,
                                  by = "ORF_ElePert")%>%
                            ## calculate deviation between exp and simulation and % absolute deviation
                            mutate(deviation = MaxMag_experimental_protein_fc - simulated_protein_fc,
                                   abs_deviation = abs(deviation),
                                   perc_abs_deviation = round(abs_deviation/abs(simulated_protein_fc)*100,2) ,
                                   perc_abs_deviation = ifelse(!is.infinite(perc_abs_deviation),perc_abs_deviation,NA))%>%
                            na.omit()%>%
                            group_by(Element)%>%
                            mutate(median_perc_abs_deviation_metal = median(perc_abs_deviation, na.rm = T))%>%
                            ungroup()%>%
                            group_by(ElePert)%>%
                            mutate(median_perc_abs_deviation_metalpert = median(perc_abs_deviation, na.rm = T))%>%
                            ungroup()%>%
                            group_by(ORF)%>%
                            mutate(median_perc_abs_deviation_ORF = median(perc_abs_deviation, na.rm = T))%>%
                            ungroup()

############################
### Visualise deviations ###
############################

# median per metal - perturbation 

pdf(paste0(plot_dir,"//median_percentage_absolute_deviation_per_metalpert.pdf"), width = 5,height = 6)
ggplot(unique(simulation_vs_experiment[,c("ElePert","median_perc_abs_deviation_metalpert")]),
       aes(x = ElePert,
           y = median_perc_abs_deviation_metalpert,
           colour = ElePert,
           fill = ElePert))+
  geom_bar(stat = "identity", width = 0.4)+
  scale_fill_manual(values = colkey_EleDir)+
  scale_colour_manual(values = colkey_EleDir)+
  theme_metallica()+
  theme(axis.text.x = element_text(angle = 90), legend.position = "none")+
  labs(x = "",y = "median (% abs(deviation_exp_vs_sim))",
       fill = "")
dev.off()

# per ORF across all metals -- prepare data for iPATH

ggplot(unique(simulation_vs_experiment[,c("ORF","median_perc_abs_deviation_ORF")]),
       aes(x = median_perc_abs_deviation_ORF))+
  geom_histogram(bins = 200)+
  xlim(0,200)

iPath_devperORF <- unique(simulation_vs_experiment[,c("ORF","median_perc_abs_deviation_ORF")]) %>%
                   mutate(Uniprot_ID = as.character(lapply(ORF,convert_ORF2Uniprot)))%>%
                   mutate(colour = as.character(lapply(median_perc_abs_deviation_ORF,convert_percabsdeviation_to_colour)),
                           width = "W12",
                           opacity = 1,
                           UNIPROT.ID = paste0("UNIPROT:",Uniprot_ID))%>%
                   dplyr::select(UNIPROT.ID,colour,width,opacity)

write_delim(iPath_devperORF, paste0(output_tables_dir,"/foriPATHmedian_perc_abs_deviation_perORF.tsv"),col_names = F)

## Save colour-key used in iPath

# Create a dataframe with pad values and corresponding colours

df <- data.frame(pad = 1:101, 
                 color = sapply(1:101, convert_percabsdeviation_to_colour))

# Plot color key with ggplot2
pdf(paste0(plot_dir,"/colour_key_iPath_medpercabsdeviation.pdf"),width = 3,height = 1)
ggplot(df, aes(x = pad, y = 1, fill = color)) +
  geom_tile(color = NA) +
  scale_fill_identity() +
  labs(x = NULL, y = NULL) +
  theme_metallica()+
  theme(panel.grid = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "none") +
  labs(title = "colour key iPATH - median percentage absolute deviation")
dev.off()



################################################################
### compare median deviation of metal binders vs non binders ###
################################################################

sim_vs_exp_metBanno <- merge(simulation_vs_experiment,GO_gset_MF_binding_metalwise, by = "ORF")%>%
                       group_by(term,ElePert)%>%
                       summarise(median_perc_abs_deviation_metBinder = median(perc_abs_deviation, na.rm = T))%>%
                       ungroup()%>%
                       mutate(colour = as.character(lapply(median_perc_abs_deviation_metBinder,convert_percabsdeviation_to_colour)))

pdf(paste0(plot_dir,"//median_percentage_absolute_deviation_metalbinders_only_metalwise.pdf"), width = 4,height = 7)
ggplot(sim_vs_exp_metBanno, aes(x = term, y = ElePert, fill = colour)) +
  geom_tile() +
  scale_fill_identity() +
  labs(x = "Term", y = "Element Perturbation", fill = "Median % Abs Deviation") +
  theme_metallica()+
  labs(x = "metal binding annotation in GO-MF", y = "")
dev.off()


sim_vs_exp_metBanno <- merge(simulation_vs_experiment,GO_gset_MF_metalbinding, by = "ORF")%>%
  group_by(term,ElePert)%>%
  summarise(median_perc_abs_deviation_metBinder = median(perc_abs_deviation, na.rm = T))%>%
  ungroup()%>%
  mutate(colour = as.character(lapply(median_perc_abs_deviation_metBinder,convert_percabsdeviation_to_colour)))

pdf(paste0(plot_dir,"/median_percentage_absolute_deviation_metalbinders_only_GOMFtermwise.pdf"), width = 7,height = 5)
ggplot(sim_vs_exp_metBanno, 
       aes(x = ElePert, 
           y = term, fill = colour)) +
  geom_tile() +
  scale_fill_identity() +
  labs(x = "Term", y = "Element Perturbation", fill = "Median % Abs Deviation") +
  theme_metallica()+
  theme(axis.text.x = element_text(angle = 90))+
  labs(x = "", y = "")
dev.off()
