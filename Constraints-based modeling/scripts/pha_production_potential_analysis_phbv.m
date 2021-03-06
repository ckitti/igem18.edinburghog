clear

%initCobraToolbox
%changeCobraSolver('gurobi')

model = readCbModel('iJO1366.mat');
Biomass_fxn = 'BIOMASS_Ec_iJO1366_core_53p95M';
[tmp, Biomass_num] = ismember(Biomass_fxn, model.rxns);
model=changeRxnBounds(model, {'EX_pnto__R_e', 'EX_o2_e'},[0 -20], 'l');

%% Modifications

% Add Metabolites

mod1 = addMetabolite(model, '3kvcoa_c', 'metFormula', 'C26H42N7O18P3S','Charge', -4);
mod1 = addMetabolite(mod1, '3hvcoa_c', 'metFormula', 'C26H44N7O18P3S','Charge', -4);
mod1 = addMetabolite(mod1, 'PHB_c','Charge', -4);
mod1 = addMetabolite(mod1, 'PHV_c','Charge', -4);

% Add Reactions

% Add PhaA (acetyl-CoA)
mod1 = addReaction(mod1, 'PhaA_ac', {'accoa_c', 'aacoa_c', 'coa_c'},...
    [-2 1 1], 1, 0, 1000);
% Add PhaA (propionyl-CoA)
mod1 = addReaction(mod1, 'PhaA_pr', {'accoa_c', 'ppcoa_c', '3kvcoa_c', 'coa_c'},...
    [-1 -1 1 1], 1, 0, 1000);

% Add PhaB (acetoacetyl-CoA)
mod1 = addReaction(mod1, 'PhaB_ac', {'h_c','aacoa_c', 'nadph_c', '3hbcoa_c', 'nadp_c'},...
    [-1 -1 -1 1 1], 1, 0, 1000);
% Add PhaB (3-ketovaleryl-CoA)
mod1 = addReaction(mod1, 'PhaB_pr', {'h_c','3kvcoa_c', 'nadph_c', '3hvcoa_c', 'nadp_c'},...
    [-1 -1 -1 1 1], 1, 0, 1000);

% Add PhaC (PHB/V)
mod1 = addReaction(mod1, 'PhaC_ac', {'3hbcoa_c', 'PHB_c', 'coa_c'}, [-1 1 1], 1, 0, 1000);
mod1 = addReaction(mod1, 'PhaC_pr', {'3hvcoa_c', 'PHV_c', 'coa_c'}, [-1 1 1], 1, 0, 1000);
mod1 = addReaction(mod1, 'PhaC_phbv', {'PHB_c', 'PHV_c', 'PHBV_c'}, [-1 -1 1], 1, 0, 1000);

% Add PHBV extraction

mod1 = addReaction(mod1,'PHBV_p', {'PHBV_c', 'PHBV_p'}, [-1 1], 1, 0, 1000);
mod1 = addReaction(mod1,'PHBV_e', {'PHBV_p', 'PHBV_e'}, [-1 1], 1, 0, 1000);
mod1 = addReaction(mod1,'EX_PHBV_e', 'PHBV_e <=>');

mod1 = changeRxnBounds(mod1, Biomass_fxn, 0.1, 'l');


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
mod1 = changeRxnBounds(mod1,'SUCTARTtpp',0,'b');
mod1 = changeRxnBounds(mod1,'SUCOAS',0,'b');

% mod1 = alphabetizeModel(mod1);

substrates = {'EX_glc__D_e', 'EX_fru_e', 'EX_glyc_e', 'EX_xyl__D_e', 'EX_succ_e','EX_ac_e', 'EX_ppa_e', 'EX_lac__L_e', 'EX_glcn_e'};

PHBVrate = zeros(1,length(substrates));
GrowthBV = zeros(1,length(substrates));

% Optimization per single carbon source/substrate
for i = 1:length(substrates)
    mod1 = changeRxnBounds(mod1, 'EX_glc__D_e', 0, 'l');
    mod1 = changeRxnBounds(mod1, substrates(i), -20, 'l');
    
    mod_BV = changeObjective(mod1, {'EX_PHBV_e'});
    FBAsolveBV = optimizeCbModel(mod_BV, 'max', 0, false);
    PHBVrate(i) = FBAsolveBV.f;
    GrowthBV(i) = FBAsolveBV.x(Biomass_num);

    mod1 = changeRxnBounds(mod1, substrates(i), 0, 'l');
end

% Sinking sort
for i = 1:length(substrates)
    for j = 1:length(substrates) - 1
        if PHBVrate(j) < PHBVrate(j+1) || PHBVrate(j) > PHBVrate(j+1)
            if PHBVrate(j) < PHBVrate(j+1)
                tPHBV = PHBVrate(j);
                PHBVrate(j) = PHBVrate(j+1);
                PHBVrate(j+1) = tPHBV;
                
                tGR = GrowthBV(j);
                GrowthBV(j) = GrowthBV(j+1);
                GrowthBV(j+1) = tGR;
                
                tSub = substrates(j);
                substrates(j) = substrates(j+1);
                substrates(j+1) = tSub;

            end 
        else
            
            if GrowthBV(j) < GrowthRate(j+1)
                tGR = GrowthBV(j);
                GrowthBV(j) = GrowthBV(j+1);
                GrowthBV(j+1) = tGR;
                
                tSub = substrates(j);
                substrates(j) = substrates(j+1);
                substrates(j+1) = tSub;
                
                tPHBV = PHBVrate(j);
                PHBVrate(j) = PHBVrate(j+1);
                PHBVrate(j+1) = tPHBV;
            end
        end
    end
end

potale = {'EX_ac_e', 'EX_ppa_e', 'EX_lac__L_e'};

mod1 = changeRxnBounds(mod1, 'EX_glc__D_e', 0, 'l');
mod1 = changeRxnBounds(mod1, potale, -20, 'l');
    
mod_BV = changeObjective(mod1, 'EX_PHBV_e');
FBAsolveBV = optimizeCbModel(mod_BV, 'max', 0, false);
PHBVrate = [PHBVrate, FBAsolveBV.f];
GrowthBV = [GrowthBV, FBAsolveBV.x(Biomass_num)];

mod1 = changeRxnBounds(mod1, potale, 0, 'l');

% Calculate product yield
productYieldBV = PHBVrate/20;

% Calculate substrate specific productivity
sspBV = productYieldBV .* GrowthBV;

substrates{1,end+1} = 'pot ale';
dat = cell(length(substrates), 4);

for i = 1:(length(substrates))
    dat{i,1} = substrates(i);
    dat{i,2} = PHBVrate(i);
    dat{i,3} = productYieldBV(i);
    dat{i,4} = sspBV(i);
end    

