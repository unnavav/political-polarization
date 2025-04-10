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

% local = parcluster('local');
% fprintf("Numbe of workers available: %i\n", local.NumWorkers)
% pool = local.parpool(23);

%% first set up some grids, pulling a lot from aiyagari
% project

vTol = 1e-4;

% Import Initialized Values
% cd ../d
% terms = load("terms.m");
% cd ../p

% params
alpha = 0.36; delta = 0.06; beta = 0.96; sigma = 3; phi = 0;

etap = 0;
etal = 0.25;
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
dgrid = [0.02 -0.02];

agrid = compute.logspace(al, ah, na)';

verbose = true;

% set size of aggregate moment grid
nm = 25;
Kgrid = linspace(0, 15, nm);

p = .9;

% forecast variables (LL, LH, PL, PH)
% guess from previous runs of the code
Kfore = [-3.3110, 0.1402;  % LL
         -0.6568, 0.8283;  % LH
         -0.8760, 0.7702;  % PL
         -0.6812, 0.8213]; % PH
       
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
    'lamval', 0);

rng("default")
T = 1000;
Rt = predict.sim(T,2,"default",Rguess);
dt = predict.sim(T,2,"default",pid);
jt = [Rt dt];

verbose = true;
forearray = cell(50,1);

newfore = Kfore;
forearray{1} = Kfore;
array_ind = 1;

foredist = 10;

cd ../d/

while foredist > vTol

    [Kprdata Varray Garray EVarray] = ks.getRegData(jt, terms, vTol, verbose);
    
    treg = 50:T;

    for Rstate = 1:2
        for dstate = 1:2

            state_index = Rstate*2-1 + dstate-1;

            data_inds = find(jt(treg,1) == Rstate & jt(treg,2) == dstate);
            Kpr_inds = data_inds+1;
            data_inds = data_inds(Kpr_inds < length(treg));
            Kpr_inds = Kpr_inds(Kpr_inds < length(treg));
            
            regdata = [ones(length(Kpr_inds),1) log(Kprdata(data_inds))];

            state_fore = regress(log(Kprdata(Kpr_inds)), regdata);
            newfore(state_index,:) = state_fore';
        end
    end

    
    array_ind = array_ind+1;
    forearray{array_ind} = newfore;
    foredist = compute.dist(newfore, terms.Kfore, 2);

    fprintf("New forecast:\n")
    disp(newfore)

    fprintf("REGRESSION DIST = %1.4f     ---------------------", foredist)
    terms.Kfore = .5*newfore + .5*terms.Kfore;
end
filename = sprintf("KS_migration_results_eta_%0.2f_.mat", etal);
save(filename)
