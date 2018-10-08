function [ f ] = fPHBdt_e( t, x )
% PHB biosynthesis with only PHA, reversible

%% Species parameters for model

ac = x(1);    %acetate
coa = x(2);   %CoA
nadph = x(3);
nadp = x(4);
s1 = x(5);    %acetyl-CoA
s3 = x(6);    %acetoacetyl-CoA
s5 = x(7);    %3HBCoA
s7 = x(8);    %PHB
e0 = x(9);    %acs
e1 = x(10);    %phaA
e2 = x(11);    %phaB
e3 = x(12);    %phaC
es1 = x(13);
es2 = x(14);

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
k1_f1 = (k1_b1 + k1_f2)/Km1_accoa;
k1_b2 = (k1_b1 + k1_f2)/Km1_a2coa;

k1_s = k1_f2/Km1_accoa;
k1_p = k1_b1/Km1_a2coa;

% phaB reversible reaction
k2_f2 = 102;
k2_b1 = 29;
k2_f1 = (k2_b1 + k2_f2)/(Km2_a2coa*Km2_nadph);
k2_b2 = (k2_b1 + k2_f2)/(Km2_nadp*Km2_3hbcoa);

k2_s = k2_f2/(Km2_a2coa * Km2_nadph);
k2_p = k2_b1/(Km2_nadp * Km2_3hbcoa);

k3 = 18.67;

%% Reaction parameters for model

% acs
v0 = e0 * ac * coa* kcat0 / (1 + (Km_ac * coa) + (Km_coa * ac) + (ac * coa));
% phaA, multi-substrate kinetics
v1 = es1 * (s1 * s1 * k1_s - coa * s3 * k1_p) / ...
    (1 + (s1 * s1)/Km1_accoa + (coa + s3)/Km1_a2coa);
% phaB
v2 = es2 * (nadph * s3 * k2_s - nadp * s5 * k2_p) /...
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

de0dt = 0;
de1dt = (k1_f2 + k1_b1) * es1 ...
    - k1_f1 * s1 * s1 * e1 - k1_b2 * e1 * s3 * coa;
des1dt = e1 * (k1_f1 * s1 * s1 + k1_b2 * s3 * coa) ...
    - (k1_f2 + k1_b1) * es1;
de2dt = (k2_f2 + k2_b1) * es2 ...
    - k2_f1 * s3 * nadph * e1 - k2_b2 * e2 * s5 * nadp;
des2dt = k2_f1 * s3 * nadph * e1 + k2_b2 * e2 * s5 * nadp ...
    - (k2_f2 + k2_b1) * es2;
de3dt = 0;

f = [dAcdt; dCoAdt; dNADPHdt; dNADPdt; dS1dt; dS3dt; dS5dt; dS7dt; ...
    de0dt; de1dt; de2dt; de3dt; des1dt; des2dt];