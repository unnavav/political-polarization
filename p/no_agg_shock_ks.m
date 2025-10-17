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

%% set up params to load into the KS algo

vTol = 1e-6;

% model params
alpha = 0.36; delta = 0.06; beta = 0.96; sigma = 3; phi = 0;

% grid sizes
nl = 7; na = 250; nmu = na*10; nr = 2; nk = round(na/10);

alpha = 0.36; delta = 0.06; beta = 0.96; sigma = 1; phi = 0;

al = 0; ah = 100;
% get labor distribution and aggregate values
%step 1: make labor grid and labor transition matrix

mu = 0;
rho = .9;
sig_l = 0.5327;
range = 2.575;

% need to back out sigma^2_e given sigma^2_l
sigx = sig_l*sqrt(1 - rho^2);

[pil, lgrid] = compute.getTauchen(nl,  mu, sigx, rho, range);

% also need to create a grid for capital depreciation shock
mu_d = 0.025; rho_d = .9; sig_mu = .01;
sigmx = sig_mu*sqrt(1 - rho_d^2);
range = 3.5;

% [pid, dgrid] = compute.getTauchen(2, mu_d, sigmx, rho_d, range);

pid = [.95 0.05; 0.05 .95];
dgrid = [0.02 -0.02];

agrid = compute.logspace(al, ah, na)';

verbose = true;

% set size of aggregate moment grid
Kgrid = linspace(8, 11, nk);

% temp placement, two coefficients each from model
Kfore = [-0.6568, 0.8283;
    -0.6812, 0.8213];
Rfore = [.7 .3; ...
          .3 .7]; % deterministic

% definitions for the probability of changing regime
probscale1 = 10;
probscale2 = 10;

phis = [probscale1 probscale2];

% finally set policies
taul = 0.08; taup = 0.18; etal = .45; etap = 0;

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
    'dgrid', dgrid, ...
    'captax', zeros(1,nl), ...
    'taugrid', [taup taul], ...
    'G', 0, ...
    'Kgrid', Kgrid, ...
    'etagrid', [etap etal], ...
    'lamval', 0, ...
    'phis', phis, ...
    'rnseed', 1234567);

rng("default")
T = 1000;
% Rt = predict.sim(T,2,"default",Rguess);

verbose = true;
forearray = cell(50,1);

newfore = Kfore;
forearray{1} = Kfore;
array_ind = 1;

foredist = 10;
modelfun = @(b, x)(1+exp(b(1)+b(2)*x))^-1;
beta0 = [.5 5];
%% begin iteration

Rt = rand(T,1);

% closing in on this bitch
iter_ct = 1;
Varray = cell(1);
EVarray = cell(1);
Garray = cell(1);
Kforearray = cell(1);
Rforearray = cell(1);

while foredist > vTol

    Rforearray{iter_ct} = Rfore;
    Kforearray{iter_ct} = Kfore;

    if foredist > 1e-1
        vTol = 1e-4;
    else
        vTol = 1e-6;
    end
    fprintf("\nGetting Regression Data\n")
    [Kprdat Rdata Prdata V G EV] = ks.getRegData(Rt, terms, vTol, verbose);

    Varray{iter_ct} = V;
    EVarray{iter_ct} = EV;
    Garray{iter_ct} = G;

    K1 = Kprdat(Rdata == 1);
    K2 = Kprdat(Rdata == 2);
    Pr1 = Prdata(Rdata == 1);
    Pr2 = Prdata(Rdata == 2);

    logit1 = fitnlm(K1, Pr1, beta0, modelfun);

end