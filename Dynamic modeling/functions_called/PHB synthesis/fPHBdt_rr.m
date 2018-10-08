function [ f ] = fPHBdt_rr( t, x )
% PHB biosynthesis with only PHA, reversible

%% Species parameters for model

c = 100;

e0 = c;    %acs
e1 = c;      %phaA
e2 = c;    %phaB
e3 = c;    %phaC

ac = x(1);    %acetate
coa = x(2);   %CoA
nadph = x(3);
nadp = x(4);
s1 = x(5);    %acetyl-CoA
s3 = x(6);    %acetoacetyl-CoA
s5 = x(7);    %3HBCoA
s7 = x(8);    %PHB

%% Kinetic parameters from literature

% Michaelis constants

Km_ac = 200;
Km_coa = 200;

Km1_accoa = 390;
Km1_a2coa = 38.38;

Km2_a2coa = 5.7;
Km2_nadph = 149;
Km2_nadp = 0.031;
Km2_3hbcoa = 0.055;

Km3 = 450;

% Turnover number and rate constants

% acs
kcat0 = 2.63;

% phaA reversible reaction
k1_f2 = 26.46;        % kcat forward (biosynthetic)
k1_b1 = 3.03E5;       % kcat reverse (degradative)
% kf1 = (kb1 + kf2)/Km1_accoa;
% kb2 = (kb1 + kf2)/Km1_a2coa;

k1_s = k1_f2/Km1_accoa;
k1_p = k1_b1/Km1_a2coa;

% phaB reversible reaction
k2_f2 = 102;
k2_b1 = 29;

k2_s = k2_f2/(Km2_a2coa * Km2_nadph);
k2_p = k2_b1/(Km2_nadp * Km2_3hbcoa);

k3 = 18.67;

%% Reaction parameters for model

% acs
v0 = e0 * ac * coa* kcat0 / (1 + (Km_ac * coa) + (Km_coa * ac) + (ac * coa));
% phaA, multi-substrate kinetics
v1 = e1 * (s1 * s1 * k1_s - coa * s3 * k1_p) / ...
    (1 + (s1 * s1)/Km1_accoa + (coa + s3)/Km1_a2coa);
% phaB
v2 = e2 * (nadph * s3 * k2_s - nadp * s5 * k2_p) /...
    (1 + (s3 * nadph)/(Km2_a2coa*Km2_nadph) + (nadp*s5)/...
    (Km2_nadp*Km2_3hbcoa));
% phaC
v3 = e3 * s5 * k3 / (Km3 + s5);

%% Differential equations

dAcdt = -v0;
dCoAdt = v3 + v1 - v0;
dNADPHdt = -v2;
dNADPdt = v2;
dS1dt = v0 - v1;
dS3dt = v1 - v2;
dS5dt = v2 - v3;
dS7dt = v3;

f = [dAcdt; dCoAdt; dNADPHdt; dNADPdt; dS1dt; dS3dt; dS5dt; dS7dt];