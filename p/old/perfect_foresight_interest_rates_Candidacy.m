%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MIT SHOCK AND THEN CHECKING THE TRANSITION
% november 2024 
%
% vaasavi
% 
% this code takes the decision and voting rules of the household from pop-
% ulism to liberalism, then saves that as the expectation

cd C:\VAASAVI\Dropbox\Education\OSU\Ongoing_Research\Populism\political-polarization\p
restoredefaultpath;
clear all; clc;
addpath(genpath(pwd));

%% importing liberalist household

cd ..\d\rsteadystates
load results_rho85sig2_t_35_eta13.mat
cd ../../p

%% now to do perfect foresight transitions

% we're doing perfect foresight over regime transitions, so it has to be 
% consistent with that. In particular 

terms.etagrid = [etagrid(1) etagrid(3)];
terms.taugrid = [taugrid(5) taugrid(3)];

T = 300;                 % number of periods
Kss = kval;               % steady-state capital level
rho = 0.85;              % persistence
shock_size = 0.02;        % 30% deviation from Kss

K_neg = zeros(T,1);
K_neg(1) = Kss * (1 - shock_size);  % start with 30% drop

for t = 2:T
    K_neg(t) = Kss + rho * (K_neg(t-1) - Kss);  % mean-revert
end

kt = [Kss K_neg'];

% % give a series of capital bumpsthat switch partway through
% T = 100;
% shock_range = 0.01;
% rng(7644870);             
% Kss = Karray{3, 1};
% X = (2*shock_range).*rand(T-1,1) - shock_range;  % U(-0.1, +0.1)
% K = Kss .* (1 + X);
% 
% kt = [Kss K' repelem(Kss, 100)];

terms.alpha = alpha;
terms.delta = delta;
terms.beta = beta;

%increasing the 
agrid = compute.logspace(agrid(1), 100, na);
amu = compute.logspace(agrid(1), agrid(na), nmu);
terms.agrid = agrid;
terms.amu = amu';

p0 = p;

dTol = 1e-4;
[Rguess, pt, EVarray, Garray, Warray] = predict.perfectForesight(kt, p0, terms, dTol);

cd ../d/
save kwiggles_guess_eta13_tau35_april2025