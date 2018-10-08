%%  Kinetic Model of C. necator PHA Pathway
%   Cheewin (Big) Kittikunapong
%   MSc Synthetic Biology and Biotechnology 2017-18
%   The University of Edinburgh (Overgraduate iGEM Team)

%% Kinetic parameters from literature

% acs (6.2.1.1)

% Km acetate = 200 uM
% Km CoA = 200 uM
% kcat = 2.63 s^-1 (Brown et al. 1977)

Km_ac = 200; Km_coa = 200; k0 = 2.63;

% phaA (EC 2.3.1.9)

% degradative (acetoacetyl-CoA)
% Km = 38.38 uM (Kim and Kim 2015)
% kcat = 3.03 x 10^5 s^-1 (Kim and Kim 2015) [k-1]
% k-2 = 6.894 x 10^-1 s^-1 (SB2015 kcat derivation)

% biosynthetic (acetyl-CoA)
% Km = 390 uM (Oeding and Schlegel 1972)
% kcat = 3 x 10^6 s^-1 (CK) / 26.46 s^-1 (SB2015) [k2]
% k1 = 6.785 x 10^-2 s^-1 (SB2015 kcat derivation)

Km1 = 390; k1 = 3 * 10^(-6);

% bktB (EC 2.3.1.16)

% degradative
% Km = 11.58 uM (Kim et al. 2014)
% Vmax = 1.5 umol/min (Kim et al. 2014)
% kcat = 102.18 s^-1 (Kim et al. 2014)

% Km = 394 uM (acetoacetyl-CoA, Haywood et al., 1988)
% Km = 93 uM (CoA, Haywood et al. 1988)

% biosynthetic
% Km = 1100 uM (acetyl-CoA, Haywood et al. 1988)
% kcat = 10 s^-1 (CK)

% Km1b = 1100; k1b = 10;

% phaB (EC 1.1.1.36)

% reduction
% kcat = 102 s^-1 (Matsumoto et al. 2013)
% Km (acetoacetyl-CoA) = 5.7  uM (Matsumoto et al. 2013)
% Km (NADPH) = 149 uM (Matsumoto et al. 2013), 19 uM (BRENDA)

% oxidation
% kcat = 29 s^-1 (Qin et al. 1999, Candida tropicalis)
% Km (3HB-CoA) = 0.055	(Qin et al. 1999, Candida tropicalis)
% Km (NADP+) = 31 (Belova et al. 1997, Cupriavidus necator)

Km2 = 5.7; k2 = 102;

% phaC (EC 2.3.1.B2)

% 3-hydroxybutyryl-CoA

% kcat = 18.67 s^-1 (Yuan et al. 2001)
% Km = 450 uM (Yuan et al. 2001)

% 3-hydroxyvaleryl-CoA

% kcat = 1.43 s^-1 (Yuan et al. 2001)
% Km = 320 uM (Yuan et al. 2001)

Km3 = 450; k3 = 18.67;

%% Species parameters for model

 E0 = 100;     % acs
 E1 = 100;  % phaA
 E2 = 100;   % phaB
 E3 = 100;   % phaC
 
 ES1 = 0;
 ES2 = 0;

Ac = 2000;      % acetate
CoA = 500;      % CoA
NADPH = 1000;
NADP = 0;

S1 = 0;     % acetyl-CoA
S2 = 0;     % propionyl-CoA
S3 = 0;     % acetoacetyl-CoA
S4 = 0;     % 3-ketovaleryl-CoA
S5 = 0;     % 3-hydroxybutyryl-CoA
S6 = 0;     % 3-hydroxyvaleryl-CoA
S7 = 0;     % 3HB
S8 = 0;     % 3HV

%% Reaction parameters for model

% acetyl-CoA synthase (acs)

% v0 = E0 * Ac * CoA * k0 / (

% beta-ketothiolase (phaA)

% v1a = E1a * S1 * k1a / (Km1a + S1);

% 3-oxoacyl-CoA reductase (phaB)

% v2 = E2 * S3 * k2 / (Km2 + S3);

% PHA polymerase (phaC)

% v3 = E3 * S5 * k3 / (Km3 + S5);

%% Functions

time = 259200; % runtime 72 hours, 3 days

% PHB biosynthesis with only PhaA, irreversible, SB2015
[t1, x1] = ode15s(@fPHBdt, [0 time], [Ac CoA S1 S3 S5 S7]);
% [t, x2] = ode15s(@fPHBdt, [0 1000], [1000 200 0 0 0 0]);
x1(end,6)

% PHB biosynthesis with only PhaA, reversible PhaA
[t2, x2] = ode15s(@fPHBdt_PhaA_rr, [0 time], [Ac CoA S1 S3 S5 S7]);
x2(end,6)

% PHB biosynthesis with only PhaA, reversible PhaA, PhaB
[t3, x3] = ode15s(@fPHBdt_rr, 0:1800:time, [Ac CoA NADPH NADP S1 S3 S5 S7]);
x3(end,8)
%plot(t3, x3(:,5:8))

% PHB biosynthesis ", variable E/ES concentrations
[t4, x4] = ode15s(@fPHBdt_e, 0:1800:time,[Ac CoA NADPH NADP S1 S3 S5 S7 E0 E1 E2 E3 ES1 ES2]);
x4(end,8)
%plot(t4, x4(:,5:8))