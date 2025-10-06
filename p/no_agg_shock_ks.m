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

alpha = 0.36; delta = 0.06; beta = 0.96; sigma = 3; phi = 0;

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
Kgrid = linspace(0, 15, nk);

Kfore = zeros(nr,nk,2);
Rfore = zeros(nr,nk,2);

%% prep VFI

terms = struct('alpha', alpha, ...
    'beta', beta, ...
    'delta', delta, ...
    'sigma', sigma, ...
    'phi', phi, ...
    'agrid', agrid, ...
    'lgrid', lgrid, ...
    'pil', pil, ...
    'RDmat', RDmat, ...
    'dgrid', dgrid, ...
    'captax', zeros(1,nl), ...
    'taugrid', [taul taup], ...
    'G', 0, ...
    'Kfore', Kfore, ...
    'Kgrid', Kgrid, ...
    'etagrid', [etal etap], ...
    'lamval', 0);

rng("default")
T = 1000;
Rt = predict.sim(T,2,"default",Rguess);

verbose = true;
forearray = cell(50,1);

newfore = Kfore;
forearray{1} = Kfore;
array_ind = 1;

foredist = 10;

%% begin iteration