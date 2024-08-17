%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perfect Foresight Paths in the Aiyagari Model
% vaasavi
% nov 2023
% 
% this code solves the perfect foresight paths in the aiyagari model, then
% creates figures for TFP, capital, output, and consumption
% 
% dependencies: 
%     compute.m
% 
% outputs:
%     TBD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

restoredefaultpath;
clear all; clc;
addpath(genpath(pwd));

%% parameters
%giving five states for z
nz = 25;
%using PS2 vals, assuming they've been calibrated correctly
rho = 0.859;
sig = 0.014;

alpha = 0.36;
delta = 0.06;
beta = 0.96;
sigma = 2.0;

T = 500;
seed = 8999305; %seed for random number generator

lambda = 0.90;
dTol = 1e-7;

% tax wedge : TODO: import tau
tau = .15;

%calculate kst
kst = ((1-beta*(1-tau))/(alpha*beta*(1-tau))-delta/alpha)^(1/(alpha-1));
k0 = 1*kst;

ist = delta*kst;
cst = kst^alpha - ist;
%% first, generate sequence of z's 

[piz, zgrid] = compute.getTauchen(nz, 0, sig, rho);

% z_deviations = [13:20 flip(10:2:18) 8:5:13];
% 
% zi = [repelem(13,15) z_deviations repelem(13,T-31)];

zt = repelem(zgrid(13), T-1);
zt(30) = (nz-median(1:nz))*1.5;

%% begin the iterating

kguess = repelem(kst, T-1);
kguess(1) = k0;

kprguess = repelem(0, T-1);

dist = max(abs(kguess - kprguess));

iter_ct = 1;
while dist > dTol

    yt = zt.*(kguess).^alpha;
    rt = alpha.*zt.*(kguess.^(alpha-1)) - delta;

    ct = repelem(cst,T-1);
    for i = flip(1:T-2)
        ct(i) = ct(i+1)*((beta*(1+rt(i+1)))^(-1/sigma));
    end
    ct(T-1) = cst;

    kt = repelem(kst, T);
    kt(1) = k0;
    for i = 1:T-1
        kt(i+1) = yt(i) - ct(i) + (1-delta)*kt(i);
    end

    kprguess = lambda*kguess + (1-lambda)*kt(1:T-1);

    dist = max(abs(kprguess - kguess)); 

    fprintf("Iteration %i: ||K' - K|| = %4.8f\n", iter_ct, dist);

    iter_ct = iter_ct + 1;
    kguess = kprguess;
end

%% make the charts

col1 = [.159 .072 .130]*(1/(.255));
col2 = [.255 .217 .120]*(1/(.255));
col3 = [.110 .182 .095]*(1/(.255));
col4 = [.086 .145 .128]*(1/(.255));

%plotting
tiledlayout(2,2)
nexttile
plot(1:T-1, zt, 'Color', col1)
title("TFP")
nexttile
plot(1:T-1, yt, 'Color', col2)
title("Output")
nexttile
plot(1:T, kt, 'Color', col3)
title("Capital")
nexttile
plot(1:T-1, ct, 'Color', col4)
title("Consumption")

exportgraphics(gcf, "./graphs/Q1_charts.png", 'Resolution', 500)
