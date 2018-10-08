%% Flux Balance Analysis (FBA)
% Cheewin Kittikunapong
% MSc Synthetic Biology and Biotechnology
% The University of Edinburgh iGEM Overgraduate 2018
%% Initialise The Cobra Toolbox and set the solver.
% If necessary, initialise the cobra toolbox:
clear 
% initCobraToolbox
%% 
% Setup the appropriate solver for the machine you are using by removing 
% the "%" (comment) sign for only the desired solver.

% changeCobraSolver('glpk','all');
% changeCobraSolver('tomlab_cplex','all');
% changeCobraASolver('ibm_cplex','all');
 changeCobraSolver ('gurobi', 'all');
%% Model Setup
% 
% Before proceeding with the simulations, load the model into the workspace:

global CBTDIR
modelFileName = 'iJO1366.mat';
modelDirectory = getDistributedModelFolder(modelFileName); %Look up the folder for the distributed Models.
modelFileName= [modelDirectory filesep modelFileName]; % Get the full path. Necessary to be sure, that the right model is loaded
model = readCbModel(modelFileName);
%% 
% Assumption: the cellular objectives include energy production or
% optimisation of uptake rates and by-product secretion for various 
% physiological functions of the human body.
%% PROCEDURE
%% 1. Standard FBA
% Standard FBA predicts an optimal solution for a cellular objective within 
% a given set of constraints on a metabolic network (see Equation 1). Constraints 
% on the network are set by assigning limits on the uptake, consumption or production 
% of metabolites in reactions.
% 
% *Timing:*
% 
% The time to determine a FBA solution depends on the size of the genome-scale 
% model and is commonly less than a second for a medium sized model.
% 
% For each new simulation, the original model will be copied to a new variable. 
% This preserves the constraints of the original model to perform further simulations 
% with new constraints. Additionally, this method of renaming the model avoids 
% confusion while performing multiple simulations at the same time.

modelaerobic = model;
%% 
% biomass objective function

 printRxnFormula(model, 'BIOMASS_Ec_iJO1366_WT_53p95M');
%% 
% We will set this reaction as our objective with the |'changeObjective'| 
% command. Maximising the flux through the biomass function will result in maximal growth conditions
% (up to the limit of the reaction).

modelaerobic = changeObjective (modelaerobic, 'BIOMASS_Ec_iJO1366_WT_53p95M');
%% 
% The glucose and oxygen, in this case, are provided in high amounts for 
% calculating the flux through ATP demand. 
% 
% The |'changeRxnBounds' |function changes the flux constraints of the lower 
% ('|l|'), upper ('|u|'), or both the bounds ('|b|'), of the specified reaction. 
% Here, we will change the maximal uptake of glucose to 20 $<math xmlns="http://www.w3.org/1998/Math/MathML" 
% display="inline"><mrow><mi>&mu;</mi></mrow></math>$mol/min/gDW and of oxygen 
% to 1000 $<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mrow><mi>&mu;</mi></mrow></math>$mol/min/gDW. 
% The uptake of oxygen is effectively unconstrainted (i.e. infinity). 

%modelaerobic = changeRxnBounds (modelaerobic, 'EX_gal_e', -20, 'l');  
modelaerobic = changeRxnBounds (modelaerobic, 'EX_o2_e', -1000, 'l'); 
%modelaerobic = changeRxnBounds (modelaerobic, 'SUCOAS', 0, 'l'); 
%modelaerobic = changeRxnBounds (modelaerobic, 'AKGDH', 0, 'l'); 
%modelaerobic = changeRxnBounds (modelaerobic, 'SUCDi', 0, 'b'); 

%% 
% The function |optimizeCbModel| calculates one of the optimal solutions 
% for a (maximum or minimum) objective reaction within the defined solution space. 
% In the above example, the maximal flux through the |<http://vmh.life/#human/all/DM_atp_c_ 
% DM_atp_c_>| is desired. 

FBAaerobic = optimizeCbModel (modelaerobic, 'max');
%% 
% * *Anticipated results*
% 
% When oxygen and all carbon sources (internal and external) are provided 
% the flux through ATP demand reaction can reach its maximum rate of 1000 $<math 
% xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mrow><mi>&mu;</mi></mrow></math>$mol/min/gDW. 
% 
% * *Troubleshooting*
% 
% IMPORTANT ADVICE!!!!!!!
% If there are multiple carbon sources available in the model, it may be 
% necessary to specify more constraints in order to examine the effect of a single 
% carbon source on biosynthesis.  
% 
% To avoid this issue, all external carbon sources need to be closed with 
% the exception of the single carbon source of interest.

%Closing the uptake of all energy and oxygen sources
[exchBool,uptBool] = findExcRxns(model);
uptakes = model.rxns(uptBool);

% IMPORTANT ADVICE!!!!
% The alternative way to do that, in case you were using another large model, 
% that does not contain defined Subsystem is
% to find uptake exchange reactions with following codes:
% [selExc, selUpt] = findExcRxns(model);
% uptakes1 = model.rxns(selUpt);

% Selecting from the exchange uptake reactions those 
% which contain at least 1 carbon in the metabolites included in the reaction:
 subuptakeModel = extractSubNetwork(model, uptakes);
 hiCarbonRxns = findCarbonRxns(subuptakeModel,1);
 
% Closing the uptake of all the carbon sources
 modelalter = model;
 modelalter = changeRxnBounds(modelalter, hiCarbonRxns, 0, 'b');
 
% Closing other oxygen and energy sources. 
% changes = {'R1',   'R2',    'R3'};
% modelalter = changeRxnBounds (modelalter, changes, 0, 'l');
%% 
% *Calculating maximum ATP energy production under anaerobic and glucose 
% only conditions:*

% modelanaerobic = modelalter;
% modelanaerobic = changeRxnBounds(modelanaerobic, 'EX_glc__D_e',-20,'l');
% modelanaerobic = changeRxnBounds (modelanaerobic, 'EX_o2_e', 0, 'l');
% modelanaerobic = changeObjective(modelanaerobic,'BIOMASS_Ec_iJO1366_WT_53p95M');
% FBAanaerob = optimizeCbModel(modelanaerobic,'max')


%% testing cell growth of E coli grown in acetate and propionate media
modelaerobic = changeRxnBounds (modelaerobic, 'EX_glc__D_e', -20, 'l');  
modelaerobic = changeRxnBounds (modelaerobic, 'EX_ac_e', -20, 'l');  
modelaerobic = changeRxnBounds (modelaerobic, 'EX_ppa_e', -20, 'l');
modelaerobic = changeRxnBounds (modelaerobic, 'EX_o2_e', -1000, 'l'); 
FBA1 = optimizeCbModel (modelaerobic, 'max')

modelaerobic = changeRxnBounds (modelaerobic, 'EX_ac_e', 0, 'l');  
modelaerobic = changeRxnBounds (modelaerobic, 'EX_ppa_e', -20, 'l');
FBA2 = optimizeCbModel (modelaerobic, 'max')

modelaerobic = changeRxnBounds (modelaerobic, 'EX_ac_e', -20, 'l');  
modelaerobic = changeRxnBounds (modelaerobic, 'EX_ppa_e', 0, 'l');
FBA3 = optimizeCbModel (modelaerobic, 'max')

modelaerobic = changeRxnBounds (modelaerobic, 'EX_ac_e', 0, 'l');  
modelaerobic = changeRxnBounds (modelaerobic, 'EX_ppa_e', 0, 'l');
FBA4 = optimizeCbModel (modelaerobic, 'max')

modelaerobic = changeRxnBounds (modelaerobic, 'EX_glc__D_e', 0, 'l');  

modelaerobic = changeRxnBounds (modelaerobic, 'EX_ac_e', -20, 'l');  
modelaerobic = changeRxnBounds (modelaerobic, 'EX_ppa_e', -20, 'l');
modelaerobic = changeRxnBounds (modelaerobic, 'EX_o2_e', -1000, 'l'); 
FBA5 = optimizeCbModel (modelaerobic, 'max')

modelaerobic = changeRxnBounds (modelaerobic, 'EX_ac_e', 0, 'l');  
modelaerobic = changeRxnBounds (modelaerobic, 'EX_ppa_e', -20, 'l');
FBA6 = optimizeCbModel (modelaerobic, 'max')

modelaerobic = changeRxnBounds (modelaerobic, 'EX_ac_e', -20, 'l');  
modelaerobic = changeRxnBounds (modelaerobic, 'EX_ppa_e', 0, 'l');
FBA7 = optimizeCbModel (modelaerobic, 'max')

modelaerobic = changeRxnBounds (modelaerobic, 'EX_ac_e', 0, 'l');  
modelaerobic = changeRxnBounds (modelaerobic, 'EX_ppa_e', 0, 'l');
FBA8 = optimizeCbModel (modelaerobic, 'max')

trial = 1:1:8
output = [FBA1.obj FBA2.obj FBA3.obj FBA4.obj FBA5.obj FBA6.obj FBA7.obj FBA8.obj]
%% 
% * *Anticipated results*
% 
% Compared to the aerobic condition, anaerobic condition with only glucose 
% as an energy source has reduced flux through ATP demand (82 $<math xmlns="http://www.w3.org/1998/Math/MathML" 
% display="inline"><mrow><mi>&mu;</mi></mrow></math>$mol/min/gDW), signifying 
% the need to oxygen to run the oxidative phosphorylation. The results are dependant 
% on the model you are using. For Recon 3.0, under anaerobic conditions with only 
% glucose as an energy source, the flux for ATP demand is 40 $<math xmlns="http://www.w3.org/1998/Math/MathML" 
% display="inline"><mrow><mi>&mu;</mi></mrow></math>$mol/min/gDW.