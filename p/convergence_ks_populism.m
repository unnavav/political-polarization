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

terms = import terms;

verbose = true;

% set size of aggregate moment grid
nK = 25;

% because I'm too lazy to change my household mod, I'll instead do the
% thing where I do an cell array of value functions

Varray = cell(nK,1);

terms.Varray = Varray;

reg_data = ks.getRegData(jt, terms, rule_guess, nl, na, nK, vTol, verbose);