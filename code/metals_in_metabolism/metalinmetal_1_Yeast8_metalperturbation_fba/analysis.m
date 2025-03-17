%% Analysis of the simulations by CofactorYeast

load('perturb_res.mat');
simdata = perturb_res;

load('CofactorYeast.mat');
% printRxnFormula(model,'rxnAbbrList',model.rxns,'metNameFlag',true);


% process the reference condition
load('fluxes_ref.mat');
protein_conc_ref = calculateProteinConc(model,model.genes,fluxes_ref);
% remove low absolute protein level in reference
cutoff_low_abs = 0.05;
low_abs_value = quantile(protein_conc_ref(protein_conc_ref>0),cutoff_low_abs);
gene_list = model.genes(protein_conc_ref>low_abs_value);
protein_abund = protein_conc_ref(protein_conc_ref>low_abs_value);

[~,txt,~] = xlsread('Yeast8_Modification.xlsx','SGDgeneNames');
gname_1 = txt(2:end,1);
gname_2 = txt(2:end,2);
clear txt;
gname_2(ismember(gname_2,'')) = gname_1(ismember(gname_2,''));
[~,b] = ismember(gene_list,gname_1);
protein_list = gname_2(b);
clear b;

ref_res = struct();
ref_res.gene_list = gene_list;
ref_res.protein_list = protein_list;
ref_res.protein_abund = protein_abund;


% process the perturbation conditions
simdata.fluxes = [simdata.fluxes_depletion simdata.fluxes_excess];
labels1 = cellfun(@(x) [x '_depletion'],simdata.metals,'UniformOutput',false);
labels2 = cellfun(@(x) [x '_excess'],simdata.metals,'UniformOutput',false);
simdata.labels = [labels1 labels2]; clear labels1 labels2;

% estimate concentrations only for the proteins with non-zero
% concentrations under the control condition
protein_conc = calculateProteinConc(model,ref_res.gene_list,simdata.fluxes);

sim_tmp = struct();
sim_tmp.gene_list = ref_res.gene_list;
sim_tmp.protein_list = ref_res.protein_list;
sim_tmp.labels = simdata.labels;
sim_tmp.FC = protein_conc./ref_res.protein_abund;

sim_tmp.regulation = sim_tmp.FC;
sim_tmp.regulation(sim_tmp.FC > 1) = 1;
sim_tmp.regulation(sim_tmp.FC < 1) = -1;
sim_tmp.regulation(sim_tmp.FC == 1) = 0;

sim = sim_tmp;
save('sim.mat','sim');


