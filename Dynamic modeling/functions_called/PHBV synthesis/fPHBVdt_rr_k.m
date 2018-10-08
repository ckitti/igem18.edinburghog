function [ f ] = fPHBVdt_rr_k( t, x, k1, k2 )
% PHBV biosynthesis with only PHA, reversible

%% Species parameters for model

c = 1;

e0 = c;    %acs
e1 = c;    %phaA
e1b = c;   %bktB
e2 = c;    %phaB
e3 = c;    %phaC

ac = x(1);    %acetate
pr = x(2);    %propionate
coa = x(3);   %CoA
s1 = x(4);    %acetyl-CoA
s2 = x(5);    %propionyl-CoA
s3 = x(6);    %acetoacetyl-CoA
s4 = x(7);    %3KVCoA
s5 = x(8);   %3HBCoA
s6 = x(9);   %3HVCoA
s7 = x(10);   %PHB
s8 = x(11);   %PHV    

%% Kinetic parameters from literature

% Michaelis constants

% acs (acetate)
Km_ac = 200; %(similar to Sc = 208)
Km_coa = 200;
% acs (propionate)
Kmb_pr = 1000; % uM (BRENDA ref 550, Sc)

% phaA (acetate)
Km1_accoa = 390; %390
Km1_a2coa = 38.38;

% bktB (propionate)
% Km1b_a2coa = 11.58;
Km1b_kvcoa = 60; % (Candida BRENDA)
Km1b_accoa = k1;
Km1b_prcoa = k2; % assumption
Km1b_coa = 0.093; % (C. necator, Haywood)

Km2 = 5.7;
% Km2_nadph = 149;
% Km2_nadp = 0.031;
% Km2_3hbcoa = 0.055;

Km3 = 450;
Kmb3 = 320; % Yuan 2011

% Turnover number and rate constants

% acs
k0 = 2.63;
k0b = 6; %s^-1 (BRENDA methanothermobacter)

% phaA reversible reaction
k1_f2 = 26.46;        % kcat forward (biosynthetic) 26.46
k1_b1 = 3.30E5;       % kcat reverse (degradative)
% kf1 = (kb1 + kf2)/Km1_accoa;
% kb2 = (kb1 + kf2)/Km1_a2coa;

k1_s = k1_f2/Km1_accoa;
k1_p = k1_b1/Km1_a2coa;

% bktB reversible reaction
k1b_f2 = 10; % assumption 10
k1b_b1 = 102.18;

k1b_s = k1b_f2/(2*Km1b_prcoa);
k1b_p = k1b_b1/Km1b_kvcoa;

% phaB irreversible reaction
k2 = 102;
% k2_b1 = 29;
% k2_s = k2_f2/(Km2_a2coa * Km2_nadph);
% k2_p = k2_b1/(Km2_nadp * Km2_3hbcoa);

k3 = 18.67;
k3b = 1.433; % (Yuan 2011, kcat 3HV)

%% Reaction parameters for model

% acs

% acetate
v0 = e0 * ac * coa* k0 / (1 + (Km_ac * coa) + (Km_coa * ac) + (ac * coa));
% propionate
v0b= e0 * pr * coa* k0b / (1 + (Kmb_pr * coa) + (Km_coa * pr) + (pr * coa));

% phaA, multi-substrate kinetics
v1 = e1 * (s1 * s1 * k1_s - coa * s3 * k1_p) / ...
    (1 + (s1 * s1)/Km1_accoa + (coa * s3)/Km1_a2coa);

% bktB, multi-substrate kinetics
v1b = e1b * (s2 * s1 * k1b_s - coa * s4 * k1b_p) /...
    (1 + (s2 * s1)/(Km1b_prcoa + Km1b_accoa) + ...
    (coa * s4)/(Km1b_kvcoa + Km1b_coa));

% phaB
% 3HB production
v2 = e2 * s3 * k2 / (Km2 + s3);
% 3HV production
v2b = e2 * s4 * k2 / (Km2 + s4);

% phaC
v3 = e3 * s5 * k3 / (Km3 + s5);
v3b = e3 * s6 * k3b /(Kmb3 + s6);


%% Differential equations

dAcdt = -v0;
dPrdt = -v0b;
dCoAdt = v3 + v1 - (v0 + v0b);
dS1dt = v0 - v1;
dS2dt = v0b - v1b;
dS3dt = v1 - v2;
dS4dt = v1b- v2b;
dS5dt = v2 - v3;
dS6dt = v2b - v3b;
dS7dt = v3;
dS8dt = v3b;

f = [dAcdt; dPrdt; dCoAdt; dS1dt; dS2dt; dS3dt; dS4dt; dS5dt; dS6dt; dS7dt; dS8dt];