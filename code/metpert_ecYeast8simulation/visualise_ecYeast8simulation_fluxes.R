#`---
#`  Title: "visualise_ecYeast8simulation_fluxes.R"
#`  @author: Simran Kaur Aulakh
#`  Date created: 6 June 2021 
#`  Description: Script to visualise flux distributions of metal depletin and excess simulations by Yu Chen
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

GEM_sim_fluxes <- data.frame(read_excel(paste0(metpert_ecYeast8simulation_dir,"/matlab_simulations/Simulated results 220602.xlsx"),sheet = 3))

# Separate control from the rest
control_fluxes <- GEM_sim_fluxes %>%
  select(Reaction.ID, Control)

# separate excess and depletion conditions
conditions <- c("Excess", "Depletion")

GEM_sim_fluxes_long <- map_dfr(conditions, ~ {
  condition_df <- GEM_sim_fluxes %>%
    select(Reaction.ID, Reaction.Name, Formula, grRules, Control, matches(paste0("\\.", .x))) %>%
    rename_with(~ stringr::str_replace(., paste0("\\.", .x), ""), -Control) %>%
    pivot_longer(cols = -c(Reaction.ID, Reaction.Name, Formula, grRules, Control),
                 names_to = "metal_perturbation", 
                 values_to = "flux_condition") %>%
    mutate(perturbation = .x,
           flux_change = abs(flux_condition - Control),
           perc_flux_change = 100 * flux_change / abs(Control))
})


# filter out infinite and NA values
GEM_sim_fluxes_long <- GEM_sim_fluxes_long %>%
  mutate(metal_perturbation = gsub("\\.", " ", metal_perturbation))%>%
  filter(is.finite(perc_flux_change), !is.na(perc_flux_change))


## sumarise per metal and per ORF
flux_change_summary <- GEM_sim_fluxes_long %>%
                        # Average flux deviation per metal-perturbation
                        group_by(metal_perturbation) %>%
                        mutate(median_perc_flux_change_permetal = ifelse(is.finite(median(perc_flux_change, na.rm = TRUE)),
                                                                         median(perc_flux_change, na.rm = TRUE), NA))%>%
                        ungroup()%>%
                        # Average flux deviation per reaction - across all metal-perturbations
                        group_by(Reaction.ID,Reaction.Name,Formula,grRules)%>%
                        mutate(median_perc_flux_change_perRxn = ifelse(is.finite(median(perc_flux_change, na.rm = TRUE)), 
                                                                       median(perc_flux_change, na.rm = TRUE), NA))%>%
                        na.omit()


pdf(paste0(plot_dir,"/flux_change_density_distributions.pdf"), width = 7,height = 4)

ggplot(filter(unique(GEM_sim_fluxes_long[,c("Reaction.ID","metal_perturbation","perc_flux_change")]),
              perc_flux_change < 100 & perc_flux_change > 1),
       aes(x = perc_flux_change,
           colour = metal_perturbation,
           fill = metal_perturbation))+
  geom_density(alpha = 0.55)+
  scale_fill_manual(values = colkey_EleDir)+
  scale_colour_manual(values = colkey_EleDir)+
  theme_metallica()+
  labs(x = "",
       fill = "")+
  scale_x_log10()
dev.off()

pdf(paste0(plot_dir,"/median_flux_change_per_metalpert.pdf"), width = 5,height = 6)
ggplot(unique(flux_change_summary[,c("metal_perturbation","median_perc_flux_change_permetal")]),
       aes(x = metal_perturbation,
           y = median_perc_flux_change_permetal,
           colour = metal_perturbation,
           fill = metal_perturbation))+
  geom_bar(stat = "identity", width = 0.4)+
  scale_fill_manual(values = colkey_EleDir)+
  scale_colour_manual(values = colkey_EleDir)+
  theme_metallica()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")+
  labs(x = "",y = "median(% flux change metal-pert vs cntrl)",
       fill = "")
dev.off()

























# Gather and reshape
GEM_sim_fluxes_long <- GEM_sim_fluxes %>%
                        select(-Control) %>%
                        gather(key = "condition", value = "value", 
                               c(Ca.Excess, Cu.Excess, Fe.Excess, K.Excess, 
                                 Mg.Excess, Mn.Excess, Na.Excess, Zn.Excess, 
                                 Ca.Depletion, Cu.Depletion, Fe.Depletion, K.Depletion, 
                                 Mg.Depletion, Mn.Depletion, Na.Depletion, Zn.Depletion)) %>%
                        separate(condition, into = c("metal", "state"), sep = "\\.") %>%
                        group_by(Reaction.ID, metal) %>%
                        mutate(flux_change = abs(diff(value)))

# Join back the control_fluxes
GEM_sim_fluxes_long <- GEM_sim_fluxes_long %>%
                        left_join(control_fluxes, by = "Reaction.ID")%>%
# Calculate perc_flux_change
                        mutate(perc_flux_change= flux_change / Control * 100)


# load tidyverse
library(tidyverse)

# read data


# filter out infinite and NA values
GEM_sim_fluxes_long <- GEM_sim_fluxes_long %>%
  filter(is.finite(perc_flux_change), !is.na(perc_flux_change))

# summarise
flux_change_summary <- GEM_sim_fluxes_long %>%
  group_by(metal) %>%
  summarise(median_perc_flux_change_permetal = median(perc_flux_change, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(metal = factor(metal, levels = unique(metal)))

# plot
ggplot(flux_change_summary, 
       aes(x = metal, y = median_perc_flux_change_permetal, fill = metal)) +
  geom_col(width = 0.6, colour = "white") +
  scale_fill_manual(values = colkey_Ele) +
  theme_metallica() +
  labs(x = "", y = "median (% flux change metal excess vs depletion)",
       fill = "") +
  theme(axis.text.x = element_text(angle = 90))
