%% simulation for the reference and perturbation conditions
% Timing: ~ 350 s
load('CofactorYeast.mat');
load('enzymedata.mat');
tic;

soplexpath = '/Users/cheyu/build/bin/soplex'; % change this to the soplex path on your PC

metals = {'Ca' 'Cu' 'Fe' 'K' 'Mg' 'Mn' 'Na' 'Zn'};
metalrxns = {'Ca(2+) exchange' 'Cu2(+) exchange' 'iron(2+) exchange' 'potassium exchange' ...
  'Mg(2+) exchange' 'Mn(2+) exchange' 'sodium exchange' 'Zn(2+) exchange'};

max_mu = 0.379;


%% Set model
% set medium
model = setMedia(model,1);% minimal media (Delft media) (default)
model = changeRxnBounds(model,'r_1714',-1000,'l');% glucose
% set oxygen
model = changeRxnBounds(model,'r_1992',-1000,'l');
% block reactions
model = blockRxns(model);


%% Set optimization

tot_protein = 0.46; % g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'r_4041','s_3717[c]'); % g/gProtein
% r_4041 is pseudo_biomass_rxn_id in the GEM, and s_3717[c] is protein id

f = tot_protein * f_modeled_protein;
f_mito = 0.1;
clear tot_protein f_modeled_protein;

factor_k_withoutcofator = 0; %%%%%%%%% this could also be changed to 0.5

%% Reference
disp('Simulating reference condition...');
model_tmp = changeRxnBounds(model,'r_2111',max_mu,'b');
rxnID = 'dilute_dummy';
fileName = writeLP(model_tmp,max_mu,f,f_mito,'Maximize',rxnID,enzymedata,factor_k_withoutcofator);
command = sprintf([soplexpath,' -s0 -g5 -t300 -f1e-20 -o1e-20 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-9 --real:fpopttol=1e-9 %s > %s.out %s'],fileName,fileName);
system(command,'-echo');
[~,~,fluxes_ref] = readSoplexResult('Simulation.lp.out',model_tmp);



%% Metal perturbation

fluxes_excess = zeros(length(model.rxns),length(metals));
fluxes_depletion = zeros(length(model.rxns),length(metals));

file = fopen('Growth Paramters for GEM.csv');
c = textscan(file,[repmat('%s',1,3),'%f',repmat('%s',1,8)],'Delimiter',',','HeaderLines',1);
listmu = c{4};
listcond = c{end};
measured_max_mu = listmu(contains(listcond,'AllEle'));

for i = 1:length(metals)
    rxnID = model.rxns(ismember(model.rxnNames,metalrxns(i)));
    
    % simulate metal excess by maximizing metal uptake
    disp(['Simulating ' metals{i} ' excess...']);
    relative_mu = listmu(contains(listcond,metals{i}) & contains(listcond,'excess'))/measured_max_mu;
    if relative_mu > 1
        relative_mu = 1;
    end
    mu_tmp = max_mu * relative_mu;
    model_tmp = changeRxnBounds(model,'r_2111',mu_tmp,'b');
    fileName = writeLP(model_tmp,mu_tmp,f,f_mito,'Minimize',rxnID,enzymedata,factor_k_withoutcofator);
    command = sprintf([soplexpath,' -s0 -g5 -t300 -f1e-20 -o1e-20 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-9 --real:fpopttol=1e-9 %s > %s.out %s'],fileName,fileName);
    system(command,'-echo');
    [~,sol_status,sol_full] = readSoplexResult('Simulation.lp.out',model_tmp);
    disp(['solution status: ' sol_status]);
    if strcmp(sol_status,'optimal')
        fluxes_excess(:,i) = sol_full;
    end
    
    % simulate metal depletion by minimizing metal uptake
    disp(['Simulating ' metals{i} ' depletion...']);
    relative_mu = listmu(contains(listcond,metals{i}) & contains(listcond,'depletion'))/measured_max_mu;
    if relative_mu > 1
        relative_mu = 1;
    end
    mu_tmp = max_mu * relative_mu;
    model_tmp = changeRxnBounds(model,'r_2111',mu_tmp,'b');
    fileName = writeLP(model_tmp,mu_tmp,f,f_mito,'Maximize',rxnID,enzymedata,factor_k_withoutcofator);
    command = sprintf([soplexpath,' -s0 -g5 -t300 -f1e-20 -o1e-20 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-9 --real:fpopttol=1e-9 %s > %s.out %s'],fileName,fileName);
    system(command,'-echo');
    [~,sol_status,sol_full] = readSoplexResult('Simulation.lp.out',model_tmp);
    disp(['solution status: ' sol_status]);
    if strcmp(sol_status,'optimal')
        fluxes_depletion(:,i) = sol_full;
    end
    
end


perturb_res = struct();
perturb_res.metals = metals;
perturb_res.fluxes_excess = fluxes_excess;
perturb_res.fluxes_depletion = fluxes_depletion;

save('fluxes_ref.mat','fluxes_ref');
save('perturb_res.mat','perturb_res');

clear;

toc;







