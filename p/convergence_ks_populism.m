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

vTol = 1e-6; dTol = 1e-4;

%% Import Initialized Values
% cd ../d
% terms = load("terms.m");
% cd ../p

%% params
alpha = 0.36; delta = 0.06; beta = 0.96; sigma = 1; phi = 0;

etap = 0;
etal = .1;
 
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
rho = .5;
sig_l = 0.5327;
 
% need to back out sigma^2_e given sigma^2_l
sigx = sig_l*sqrt(1 - rho^2);

[pil, lgrid] = compute.getTauchen(nl,  mu, sigx, rho);

agrid = compute.logspace(al, ah, na)';

verbose = true;

% set size of aggregate moment grid
nm = 25;
Kgrid = linspace(5, 20, nm);

p = .9;

% forecast variables
Kfore = [ .05 .95; ...
           .03 .9];

%prepare for VFI
terms = struct('alpha', alpha, ...
    'beta', beta, ...
    'delta', delta, ...
    'sigma', sigma, ...
    'phi', phi, ...
    'agrid', agrid, ...
    'lgrid', lgrid, ...
    'pil', pil, ...
    'captax', zeros(1,nl), ...
    'tau', 0, ...
    'G', 0, ...
    'p', p, ...
    'Kfore', Kfore, ...
    'Kgrid', Kgrid, ...
    'etap', etap, ...
    'etal', etal);

reg_data = ks.getRegData(jt, terms, rule_guess, nl, na, nK, vTol, verbose);