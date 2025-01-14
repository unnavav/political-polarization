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

load results_t0.0000_eta1.0000.mat

VP = V;
GP = G;
adistrP = adistr;
rP = r;
lP = growth_lagg;

mkdir("results_t0.0000_eta1.0000")
cd results_t0.0000_eta1.0000\
writematrix(V,"V.csv");
writematrix(G,"G.csv");
writematrix(adistr, "adistr.csv");
writematrix(r, "r.csv");
writematrix(growth_lagg, "l.csv");
cd ..

load terms_struct_t0.0000_eta1.0000.mat
termsP = terms;

load results_t0.0000_eta1.1000.mat
VL = V;
GL = G;
adistrL = adistr;
rL = r;
lL = growth_lagg;

mkdir("results_t0.0000_eta1.1000")
cd results_t0.0000_eta1.1000\
writematrix(V,"V.csv");
writematrix(G,"G.csv");
writematrix(adistr, "adistr.csv");
writematrix(r, "r.csv");
writematrix(growth_lagg, "l.csv");
cd ..

load terms_struct_t0.0000_eta1.1000.mat
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
lt = repelem(lL, 500);

r0 = rP;
r1 = rL;

alpha = 0.36; delta = 0.06;

terms.alpha = alpha;
terms.delta = delta;

lambda = 0.5;

transition(r0, r1, lt, terms, dTol, lambda)

lt = repelem(lP, 500);

r0 = rL;
r1 = rP;

lambda = 0.5;

transition(r0, r1, lt, terms, dTol, lambda)

save transition_L_to_P
