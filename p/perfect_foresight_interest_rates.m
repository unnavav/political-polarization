%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MIT SHOCK AND THEN CHECKING THE TRANSITION
% november 2024 
%
% vaasavi
% 
% this code takes the decision and voting rules of the household from pop-
% ulism to liberalism, then saves that as the expectation

restoredefaultpath;
clear all; clc;
addpath(genpath(pwd));

%% importing liberalist household

cd ..\d\liberalism_populism

load results_t0.1500_eta1.0100.mat

VP = V;
GP = G;
adistrP = adistr;
rP = r;

load terms_struct_t0.1500_eta1.0000.mat
termsP = terms;

load results_t0.0500_eta1.0500.mat

VL = V;
GL = G;
adistrL = adistr;
rL = r;

load terms_struct_t0.0500_eta1.1000.mat
termsL = terms;

load results_transition_t0.0500_eta1.1000.mat
VLP = V;
GLP = G;
adistrLP = adistr;
rLP = r;

load results_transition_t0.1500_eta1.0000.mat
VPL = V;
GPL = G;
adistrPL = adistr;
rPL = r;

cd ../../p

%% now to do perfect foresight transitions

% give a series of interest rates that switch partway through

%lagg will change; first period and then the growth rate changes, 
% which means new prices for everything

lP=1.2145; lL = 1.3359; % pulling these manually bc i forgot to save them

lt = [repelem(lP, 25), repelem(lL, 200)];

rstart = rP;
rend = rL;

kst = ((rstart+delta)/alpha)^(1/(alpha-1))*lP;
k0 = ((rend+delta)/alpha)^(1/(alpha-1))*lL;

