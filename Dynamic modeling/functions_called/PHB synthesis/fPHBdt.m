function [ f ] = fPHBdt( t, x )
% PHB biosynthesis with only PHA, irreversible, SB2015

%% Species parameters for model

c = 100;

e0 = c;    %acs
e1 = c;    %phaA
e2 = c;    %phaB
e3 = c;    %phaC

ac = x(1);    %acetate
coa = x(2);   %CoA
s1 = x(3);    %acetyl-CoA
s3 = x(4);    %acetoacetyl-CoA
s5 = x(5);    %3HBCoA
s7 = x(6);    %PHB

%% Kinetic parameters from literature

k0 = 2.63;
k1 = 26.46;    %SB2015
k2 = 102;
k3 = 18.67;

Km_ac = 200;
Km_coa = 200;
Km1 = 390;
Km2 = 5.7;
Km3 = 450;

%% Reaction parameters for model

v0 = e0 * ac * coa* k0 / (1+ (Km_ac * coa) + (Km_coa * ac) + (ac * coa));
v1 = e1 * s1 * k1 / (Km1 + s1);
v2 = e2 * s3 * k2 / (Km2 + s3);
v3 = e3 * s5 * k3 / (Km3 + s5);

%% Differential equations

dAcdt = -v0;
dCoAdt = v3 + v1 - v0;
dS1dt = v0 + (-2 * v1);
dS3dt = v1 - v2;
dS5dt = v2 - v3;
dS7dt = v3;

f = [dAcdt; dCoAdt; dS1dt; dS3dt; dS5dt; dS7dt];