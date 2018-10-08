clear

%initCobraToolbox
%changeCobraSolver('gurobi')

model = readCbModel('iJO1366.mat');
Biomass_fxn = 'BIOMASS_Ec_iJO1366_core_53p95M';
[tmp, Biomass_num] = ismember(Biomass_fxn, model.rxns);
model=changeRxnBounds(model, {'EX_mg2_e', 'EX_o2_e'},[-20 -20], 'l');

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

% Add PHBV extraction

mod1 = addReaction(mod1,'PHB_p', {'PHB_c', 'PHB_p'}, [-1 1], 1, 0, 1000);
mod1 = addReaction(mod1,'PHB_e', {'PHB_p', 'PHB_e'}, [-1 1], 1, 0, 1000);
mod1 = addReaction(mod1,'EX_PHB_e', 'PHB_e <=>');

mod1 = addReaction(mod1,'PHV_p', {'PHV_c', 'PHV_p'}, [-1 1], 1, 0, 1000);
mod1 = addReaction(mod1,'PHV_e', {'PHV_p', 'PHV_e'}, [-1 1], 1, 0, 1000);
mod1 = addReaction(mod1,'EX_PHV_e', 'PHV_e <=>');

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

% mod1 = alphabetizeModel(mod1);

substrates = {'EX_glc__D_e', 'EX_fru_e', 'EX_glyc_e', 'EX_xyl__D_e', 'EX_succ_e','EX_ac_e', 'EX_ppa_e', 'EX_lac__L_e', 'EX_glcn_e'};

PHBrate = zeros(1,length(substrates));
GrowthB = zeros(1,length(substrates));
PHVrate = zeros(1,length(substrates));
GrowthV = zeros(1,length(substrates));
PHBVrate= zeros(1,length(substrates));
GrowthBV= zeros(1,length(substrates));

% Optimization per single carbon source/substrate
for i = 1:length(substrates)
    mod1 = changeRxnBounds(mod1, 'EX_glc__D_e', 0, 'l');
    mod1 = changeRxnBounds(mod1, substrates(i), -20, 'l');
    
    mod_B = changeObjective(mod1, 'EX_PHB_e');
    FBAsolveB = optimizeCbModel(mod_B, 'max', 0, false);
    PHBrate(i) = FBAsolveB.f;
    GrowthB(i) = FBAsolveB.x(Biomass_num);
    
    mod_V = changeObjective(mod1, 'EX_PHV_e');
    FBAsolveV = optimizeCbModel(mod_V, 'max', 0, false);
    PHVrate(i) = FBAsolveV.f;
    GrowthV(i) = FBAsolveV.x(Biomass_num);
    
    mod_BV = changeObjective(mod1, {'EX_PHB_e', 'EX_PHV_e'})  ;
    FBAsolveBV = optimizeCbModel(mod_BV, 'max', 0, false);
    PHBVrate(i) = FBAsolveBV.f;
    GrowthBV(i) = FBAsolveBV.x(Biomass_num);

    mod1 = changeRxnBounds(mod1, substrates(i), 0, 'l');
end

% Sinking sort
for i = 1:length(substrates)
    for j = 1:length(substrates) - 1
        if PHBrate(j) < PHBrate(j+1) || PHBrate(j) > PHBrate(j+1)
            if PHBrate(j) < PHBrate(j+1)
                tPHB = PHBrate(j);
                PHBrate(j) = PHBrate(j+1);
                PHBrate(j+1) = tPHB;
                
                tGR = GrowthB(j);
                GrowthB(j) = GrowthB(j+1);
                GrowthB(j+1) = tGR;
                
                tSub = substrates(j);
                substrates(j) = substrates(j+1);
                substrates(j+1) = tSub;

            end 
        else
            
            if GrowthB(j) < GrowthRate(j+1)
                tGR = GrowthB(j);
                GrowthB(j) = GrowthB(j+1);
                GrowthB(j+1) = tGR;
                
                tSub = substrates(j);
                substrates(j) = substrates(j+1);
                substrates(j+1) = tSub;
                
                tPHB = PHBrate(j);
                PHBrate(j) = PHBrate(j+1);
                PHBrate(j+1) = tPHB;
            end
        end
    end
end

for i = 1:length(substrates)
    for j = 1:length(substrates) - 1
        if PHVrate(j) < PHVrate(j+1) || PHVrate(j) > PHVrate(j+1)
            if PHVrate(j) < PHVrate(j+1)
                tPHV = PHVrate(j);
                PHVrate(j) = PHVrate(j+1);
                PHVrate(j+1) = tPHV;
                
                tGR = GrowthV(j);
                GrowthV(j) = GrowthV(j+1);
                GrowthV(j+1) = tGR;
                
                tSub = substrates(j);
                substrates(j) = substrates(j+1);
                substrates(j+1) = tSub;

            end 
        else
            
            if GrowthV(j) < GrowthRate(j+1)
                tGR = GrowthV(j);
                GrowthV(j) = GrowthV(j+1);
                GrowthV(j+1) = tGR;
                
                tSub = substrates(j);
                substrates(j) = substrates(j+1);
                substrates(j+1) = tSub;
                
                tPHV = PHVrate(j);
                PHVrate(j) = PHVrate(j+1);
                PHVrate(j+1) = tPHV;
            end
        end
    end
end

potale = {'EX_ac_e', 'EX_ppa_e', 'EX_lac__L_e'};

mod1 = changeRxnBounds(mod1, 'EX_glc__D_e', 0, 'l');
mod1 = changeRxnBounds(mod1, potale, -20, 'l');
    
mod_B = changeObjective(mod1, 'EX_PHB_e');
FBAsolveB = optimizeCbModel(mod_B, 'max', 0, false);
PHBrate = [PHBrate,FBAsolveB.f];
GrowthB = [GrowthB,FBAsolveB.x(Biomass_num)];
    
mod_V = changeObjective(mod1, 'EX_PHV_e');
FBAsolveV = optimizeCbModel(mod_V, 'max', 0, false);
PHVrate = [PHVrate, FBAsolveV.f];
GrowthV = [GrowthV, FBAsolveV.x(Biomass_num)];
    
mod_BV = changeObjective(mod1, {'EX_PHB_e', 'EX_PHV_e'})  ;
FBAsolveBV = optimizeCbModel(mod_BV, 'max', 0, false);
PHBVrate = [PHBVrate, FBAsolveBV.f];
GrowthBV = [GrowthBV, FBAsolveBV.x(Biomass_num)];

mod1 = changeRxnBounds(mod1, potale, 0, 'l');

% Calculate product yield
productYieldB = PHBrate/20;
productYieldV = PHVrate/20;
productYieldBV = PHBVrate/20;

% Calculate substrate specific productivity
sspB = productYieldB .* GrowthB;
sspV = productYieldV .* GrowthV;
sspBV = productYieldBV .* GrowthBV;

substrates{1,end+1} = 'pot ale';
dat = cell(length(substrates), 10);

for i = 1:(length(substrates))
    dat{i,1} = substrates(i);
    dat{i,2} = PHBrate(i);
    dat{i,3} = PHVrate(i);
    dat{i,4} = PHBVrate(i);
    dat{i,5} = productYieldB(i);
    dat{i,6} = productYieldV(i);
    dat{i,7} = productYieldBV(i);
    dat{i,8} = sspB(i);
    dat{i,9} = sspV(i);
    dat{i,10} = sspBV(i);
end    


