%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementing Krussel Smith
% vaasavi
% oct 2025
% 
% this code aims to implement krussel smith by several aggregate states:
% capital depreciation (delta), regime (R), aggregate capital (K).
% algorithm detailed in paper, etc. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

restoredefaultpath;
clear all; clc;
addpath(genpath(pwd));


%% pulling in steady states of interest and relevant policies

cd ../d/steadystates/
load resultsna100ah50rho90sig3.mat

clearvars -except Warray Karray etagrid taugrid captax

g1 = Warray{1,5}; g2 = Warray{3,3};
K_ss_pop = Karray{1,5}; K_ss_lib = Karray{3,3};

etap = etagrid(1); etal = etagrid(3);
taup = taugrid(5); taul = taugrid(3);

taugrid = [taup taul]; etagrid = [etap etal];

cd ../../p

%% set up params to load into the KS algo

vTol = 1e-6;

% model params
alpha = 0.36; delta = 0.06; beta = 0.96; sigma = 3; phi = 0;

% grid sizes
nl = 7; na = 100; nmu = na*10; nr = 2; nk = 13;

al = 0; ah = 50;
% get labor distribution and aggregate values
%step 1: make labor grid and labor transition matrix

mu = 0;
rho = .9;
sig_l = 0.2;
range = 2.575;

% need to back out sigma^2_e given sigma^2_l
sigx = sig_l*sqrt(1 - rho^2);

[pil, lgrid] = compute.getTauchen(nl,  mu, sigx, rho, range);

agrid = compute.logspace(al, ah, na)';
nmu = na*10;
amu = linspace(agrid(1), agrid(na), nmu);

verbose = true;

% set size of aggregate moment grid
kl = .85;
kh = 1.15;
Kgrid = linspace(K_ss_pop*kl, K_ss_lib*kh, nk); 

% we also need the scaling grid for lambda to go in:
stationary_pil = asymptotics(dtmc(pil));      % nl x 1
Eeps     = stationary_pil*lgrid';              % sum ε Ω
Eeps_sum = sum(Eeps);

Eeps_pow = sum((lgrid'.^(1 - taugrid)) .* stationary_pil'); % sum ε^(1-τ) Ω  (τ regime-specific)

lambda_ratio = Eeps_sum ./ Eeps_pow;   

% start with a forecast that's basically steady state:

Kfore = [log(K_ss_pop)*0.01, 0.99;   % Small persistence, close to SS
         log(K_ss_lib)*0.01, 0.99];

Rfore = [.95 0.05; ...
    0.05 .95]; %REMEMBER TO FIX THIS WHEN YOU BEGIN ACTUAL RFORE
%% prep VFI

terms = struct('alpha', alpha, ...
    'beta', beta, ...
    'delta', delta, ...
    'sigma', sigma, ...
    'phi', phi, ...
    'agrid', agrid, ...
    'lgrid', lgrid, ...
    'pil', pil, ...
    'Kfore', Kfore, ...
    'Rfore', Rfore, ...
    'captax', captax, ...
    'taugrid', taugrid, ...
    'G', 0, ...
    'Kgrid', Kgrid, ...
    'etagrid', etagrid, ...
    'lamval', lambda_ratio, ...
    'rnseed', 1234567);

rng("default")
T = 1000;
% Rt = predict.sim(T,2,"default",Rguess);

verbose = true;
forearray = cell(50,1);

newfore = Kfore;
forearray{1} = Kfore;
array_ind = 1;
b1 = Kfore(1,:); b2 = Kfore(2,:);

foredist = 10;
kforedist = 10;
%% begin iteration

% closing in on this bitch
iter_ct = 1;
Varray = cell(1);
EVarray = cell(1);
Garray = cell(1);
Kforearray = cell(1);
Rforearray = cell(1);

terms.starter_distr = g1;

while foredist > vTol

    if foredist > 1e-1
        vTol = 1e-4;
    else
        vTol = 1e-6;
    end

    fprintf("\n ===== Generating K forecast Rule ===== \n")


    fprintf("\nGetting Regression Data\n")
    [Kprdat Rdata Prdata distr_array V G EV] = ...
        ks.getRegData(T, terms, vTol, verbose);

    fprintf('\n\nRegime 1: %d periods, Regime 2: %d periods\n', ...
        sum(Rdata==1), sum(Rdata==2));

    Varray{iter_ct} = V;
    EVarray{iter_ct} = EV;
    Garray{iter_ct} = G;
        
    K_next = log(Kprdat(2:end));
    K_curr = log(Kprdat(1:end-1));
    R_curr = Rdata(1:end-1);
        
    % Capital law: log K_{t+1} = a(R_t) + b(R_t) log K_t
    if sum(R_curr == 1) > 10
        kmdl1 = fitlm(K_curr(R_curr==1), K_next(R_curr==1));
        b1 = kmdl1.Coefficients.Estimate;
    end
    
    if sum(R_curr == 2) > 10
        kmdl2 = fitlm(K_curr(R_curr==2), K_next(R_curr==2));
        b2 = kmdl2.Coefficients.Estimate;
    end
    Kfore_new = [b1'; b2'];
    foredist = norm(Kfore_new-Kfore, 'inf');

    % Update forecast
    Kfore = 0.8*Kfore + 0.2*Kfore_new;
    terms.Kfore = Kfore;

    fprintf('K(1) = %0.4f\n', Kprdat(1));
    fprintf('K range: [%0.4f, %0.4f]\n', min(Kprdat), max(Kprdat));

    fprintf('Regime 1: log K'' = %0.4f + %0.4f log K\n', b1(1), b1(2));
    fprintf('Regime 2: log K'' = %0.4f + %0.4f log K\n', b2(1), b2(2));

    Kforearray{iter_ct} = Kfore;

    fprintf("\nCapital: R2 for pop = %0.6f\nR2 for lib = %0.6f", ...
        kmdl1.Rsquared.Adjusted, ...
        kmdl2.Rsquared.Adjusted)
    fprintf("\nForedist = %0.6f\n\n", foredist)

    iter_ct = iter_ct + 1;
end

save ../d/ks_rmseed_in_pr.mat