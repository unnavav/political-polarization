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

load results_t0.0500_eta1.0500.mat

VL = V;
GL = G;
adistrL = adistr;
rL = r;

cd ../../p

%% now to do perfect foresight transitions

% give a series of interest rates that switch partway through

R = [repelem(rP, 50), repelem(rL, 400)];