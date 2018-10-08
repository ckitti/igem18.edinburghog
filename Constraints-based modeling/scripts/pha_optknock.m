clear

initCobraToolbox
changeCobraSolver('gurobi')

model = readCbModel('iJO1366.mat');
Biomass_fxn = 'BIOMASS_Ec_iJO1366_core_53p95M';
[tmp, Biomass_num] = ismember(Biomass_fxn, model.rxns);
model=changeRxnBounds(model, {'EX_glyc_e','EX_ac_e', 'EX_ppa_e', 'EX_lac__L_e', 'EX_o2_e'},[-20 -20 -20 -20 -20], 'l');
% substrates = {'EX_glc__D_e', 'EX_fru_e', 'EX_glyc_e', 'EX_xyl__D_e', 'EX_succ_e','EX_ac_e', 'EX_ppa_e', 'EX_lac__L_e', 'EX_glcn_e'};

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
%mod1 = addReaction(mod1, 'PhaC_phbv', {'PHB_c', 'PHV_c', 'PHBV_c'}, [-1 -1 1], 1, 0, 1000);

% Add PHBV extraction

mod1 = addReaction(mod1,'PHB_p', {'PHB_c', 'PHB_p'}, [-1 1], 1, 0, 1000);
mod1 = addReaction(mod1,'PHB_e', {'PHB_p', 'PHB_e'}, [-1 1], 1, 0, 1000);
mod1 = addReaction(mod1,'EX_PHB_e', 'PHB_e <=>');

mod1 = addReaction(mod1,'PHV_p', {'PHV_c', 'PHV_p'}, [-1 1], 1, 0, 1000);
mod1 = addReaction(mod1,'PHV_e', {'PHV_p', 'PHV_e'}, [-1 1], 1, 0, 1000);
mod1 = addReaction(mod1,'EX_PHV_e', 'PHV_e <=>');

%mod1 = addReaction(mod1,'PHBV_p', {'PHBV_c', 'PHBV_p'}, [-1 1], 1, 0, 1000);
%mod1 = addReaction(mod1,'PHBV_e', {'PHBV_p', 'PHBV_e'}, [-1 1], 1, 0, 1000);
%mod1 = addReaction(mod1,'EX_PHBV_e', 'PHBV_e <=>');

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

reducedRxns = {'3HAD120', '3HAD140', '3HAD160', '3OAR120', '3OAR160',...
    '3OAS120', '3OAS160', 'ACKr', 'ADK1', 'ADNK1', 'ASNS2', 'ASPO3', 'CBMKr',...
    'CPPPGO', 'CYTK1', 'DHAPT', 'DHORD2', 'DMPPS', 'EAR1OOx', 'EAR120x', 'EAR121x'...
    'EAR140x', 'EAR141x', 'EAR160x', 'EAR161x', 'EAR40x', 'EAR60x', 'EAR80x', 'ENO',...
    'F6PA', 'FBA', 'FLDR', 'FUM', 'G5SD', 'G6PDH2r', 'GAPD', 'GARFT', 'GART',...
    'GHMT2r', 'GLU5K', 'GLUDy', 'GND', 'IMPD', 'IPDPS', 'KAS15', 'MDH', 'MTHFC', ...
    'MTHFD', 'NPDK1', 'NPDK3', 'NPDK5', 'NPDK7', 'OPHHX', 'PAPSR', 'PDH', 'PDX5PO2', ...
    'PFK', 'PFL', 'PGCD', 'PGI', 'PGK', 'PGL', 'PGM', 'PPC', 'RNDR3', 'RNDR4', 'RNTR1c',...
    'RNTR2c', 'RNTR3c', 'RPE', 'RPI', 'SUCDi', 'SUCOAS', 'TALA', 'TKT1', 'TPI', 'TRDR',...
    'TRPAS2', 'TRPS3', 'URIDK2r', 'VALTA', 'ABUTt2pp', 'ACCOAL', 'ACOATA', 'ACS', 'ADK3', 'ADPT',...
    'ASPO4', 'ASPO5', 'DHORD5', 'FRD2', 'FRD3', 'GLBRAN2', 'GLCP', 'GLCP2', 'GLCS1', 'GLDBRAN2',...
    'GLGC', 'GRXR', 'GTHOr', 'HPYRI', 'HYPYRRx', 'ICHORS', 'ICHORSi', 'IPDDI', 'KAS14', 'KAT1', ...
    'NPDK6', 'NPDK8', 'PAPSR2', 'PPAKr', 'PPCSCT', 'PTA2', 'PYK', 'R15BPK', 'R1PK', 'RNDR1', 'RNDR1b',...
    'RNDR2b', 'RNDR3b', 'RNDR4b', 'RNTR4c', 'TKT2', 'TRSARr', 'VPAMT'};

%OptKnock analysis

options.targetRxn = {'PhaC_pr'};
options.vMax = 1000;
options.numDel = 5;
options.numDelSense = 'L';
constrOpt.rxnList = {Biomass_fxn, 'ATPM'};
constrOpt.values = [0.01, 3.15];
constrOpt.sense = 'GE';

%[optKnockSol, bilevelMILPproblem] = OptKnock(mod1, reducedRxns, options,
%constrOpt, {}, false);
optKnockSol = OptKnock(mod1, reducedRxns, options, constrOpt);

%[growthRate, minProd, maxProd] = testOptKnockSol(mod1, 'EX_PHB_e','optKnockSol.rxnList);

for i = 1:length(optKnockSol.rxnList)
    modelKO = changeRxnBounds(mod1,optKnockSol.rxnList(i), 0, 'b');
end

FBAsolution = optimizeCbModel(modelKO, 'max', 0, false);
disp('rxns to KO');
disp([optKnockSol.rxnList]);
disp('Biomass');
disp(FBAsolution.f);
disp('Production');
disp(optKnockSol.obj)
