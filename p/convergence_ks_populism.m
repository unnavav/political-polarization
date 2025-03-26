%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setting up original household problem and solving it
% vaasavi
% may 2024
% 
% aiyagari + progressive tax policy
%
% testing different scenarios to see when republicans
% win vs when democrats win
% INPUTS: 
%
% OUTPUTS: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

restoredefaultpath;
clear all; clc;
addpath(genpath(pwd));

% first set up some grids, pulling a lot from aiyagari
% project

vTol = 1e-4;

%% Import Initialized Values
% cd ../d
% terms = load("terms.m");
% cd ../p

%% params
alpha = 0.36; delta = 0.06; beta = 0.96; sigma = 1; phi = 0;

etap = 0;
etal = .1;
taup = 0;
taul = 0;
 
nl = 7;
na = 250;
nmu = na*10;
np = 2;

al = 0+phi; ah = 50+phi;

% initial r guess: 
r = 0.04;

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
dgrid = [0.00 -0.02];

agrid = compute.logspace(al, ah, na)';

verbose = true;

% set size of aggregate moment grid
nm = 25;
Kgrid = linspace(5, 20, nm);

p = .9;

% forecast variables (LL, LH, PL, PH)
Kfore = [0.0, 0.85, -0.1;  % LL
         0.0, 0.85, -0.1;  % LH
         0.0, 0.85, -0.1;  % PL
         0.0, 0.85, -0.1]; % PH


Rguess = [.9 .1; ...
          .1 .9];

RDmat = kron(Rguess, pid);

%prepare for VFI
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
    'lambda', 0);

rng("default")
T = 400;
jt = rand(T,1);
jt = jt > .5;
jt(1) = 1;
foredist = 10;

verbose = true;
[Varray, Garray, EVarray] = ks.parsolve(terms, vTol, verbose);

while foredist > vTol

    [reg_data VP VL GP GL V0] = ks.getRegData(jt, terms, vTol, verbose);
    
    treg = 101:T;
    rd = reg_data(treg);
    % regression one--populist regime
    P_ind = find(jt(treg)) ;
    Ppr_ind = P_ind+1;
    P_ind = P_ind(Ppr_ind < length(treg));
    Ppr_ind = Ppr_ind(Ppr_ind < length(treg));
    L_ind = find(~jt(treg));
    Lpr_ind = L_ind+1;
    L_ind = L_ind(Lpr_ind < length(treg));
    Lpr_ind = Lpr_ind(Lpr_ind < length(treg));
    
    
    dat = [ones(length(P_ind),1) log(rd(P_ind))];
    pfore = regress(log(rd(Ppr_ind)), dat);
    
    dat = [ones(length(L_ind),1) log(rd(L_ind))];
    lfore = regress(log(rd(Lpr_ind)), dat);

    newfore = [pfore'; lfore'];
    foredist = compute.dist(newfore, terms.Kfore,2);
    terms.Kfore = .5*newfore + .5*terms.Kfore;
    fprintf("Foredist = %1.4f\n", foredist);
    disp(newfore)
end