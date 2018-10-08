clear

%initCobraToolbox
%changeCobraSolver('gurobi')

model = readCbModel('iJO1366.mat');
Biomass_fxn = 'BIOMASS_Ec_iJO1366_core_53p95M';
[tmp, Biomass_num] = ismember(Biomass_fxn, model.rxns);
mod1=changeRxnBounds(model, {'EX_mg2_e', 'EX_o2_e'},[0 -20], 'l');

% Orth et al., 2011, Heshiki 2013
mod1 = changeRxnBounds(mod1,'CAT',0,'b');
mod1 = changeRxnBounds(mod1,'DHPTDNR',0,'b');
mod1 = changeRxnBounds(mod1,'DHPTDNRN',0,'b');
mod1 = changeRxnBounds(mod1,'FHL',0,'b');
mod1 = changeRxnBounds(mod1,'SPODM',0,'b');
mod1 = changeRxnBounds(mod1,'SPODMpp',0,'b');
mod1 = changeRxnBounds(mod1,'SUCASPtpp',0,'b');
mod1 = changeRxnBounds(mod1,'SUCFUMtpp',0,'b');
mod1 = changeRxnBounds(mod1,'SUCMALtpp',0,'b');
%mod1 = changeRxnBounds(mod1,'SUCOAS',0,'b');
%mod1 = changeRxnBounds(mod1, 'AKGDH',0,'b');
%mod1 = changeRxnBounds(mod1,'SUCDi',0,'b');

% mod1 = alphabetizeModel(mod1);

substrates = {'EX_glc__D_e', 'EX_fru_e', 'EX_glyc_e', 'EX_xyl__D_e','EX_ac_e', 'EX_lac__L_e', 'EX_glcn_e'};

Growth = zeros(1,length(substrates));

% Optimization per single carbon source/substrate
for i = 1:length(substrates)
    mod1 = changeRxnBounds(mod1, 'EX_glc__D_e', 0, 'l');
    mod1 = changeRxnBounds(mod1, substrates(i), -20, 'l');
    
    mod_B = changeObjective(mod1, Biomass_fxn);
    FBAsolveB = optimizeCbModel(mod_B, 'max', 0, false);
    Growth(i) = FBAsolveB.x(Biomass_num);

    mod1 = changeRxnBounds(mod1, substrates(i), 0, 'l');
end


potale = {'EX_ac_e', 'EX_ppa_e', 'EX_lac__L_e'};

mod1 = changeRxnBounds(mod1, 'EX_glc__D_e', 0, 'l');
mod1 = changeRxnBounds(mod1, potale, -20, 'l');
    
mod_B = changeObjective(mod1, Biomass_fxn);
FBAsolveB = optimizeCbModel(mod_B, 'max', 0, false);
Growth = [Growth,FBAsolveB.x(Biomass_num)];
    
mod1 = changeRxnBounds(mod1, potale, 0, 'l');

substrates{1,end+1} = 'pot ale';
dat = cell(length(substrates), 10);

for i = 1:(length(substrates))
    dat{i,1} = substrates(i);
    dat{i,2} = Growth(i);
end    

