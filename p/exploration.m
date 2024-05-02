%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setting up original household problem and solving it
% vaasavi
% march 2024
% 
% solving via grid search,because there's noise and 
% ergo i donut trust splines (not sure if this instinct)
% is correct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

restoredefaultpath;
clear all; clc;
addpath(genpath(pwd));

% first set up some grids, pulling a lot from aiyagari
% project
% delete(gcp('nocreate'));
% parpool('local',4);

vTol = 1e-4; gTol = 1e-6; kTol = 1e-4;
%% params
alpha = 0.36; delta = 0.08; beta = 0.96173; sigma = 1; phi = 2;

nl = 7;
na = 250;
nmu = 2500;
np = 2;

al = 1; ah = 101;

wage = 1.1; r = 0.04;

%political params
identity = 0.01; pctDem = .5;
% get labor distribution and aggregate values
%step 1: make labor grid and labor transition matrix

mu = 0;
rho = .9554;
sig_l = 0.5327;

% need to back out sigma^2_e given sigma^2_l
sigx = sig_l*sqrt(1 - rho^2);

[pil, lgrid] = compute.getTauchen(nl,  mu, sigx, rho);

% get steady state labor distribution

ldist = asymptotics(dtmc(pil));

lagg = ldist*lgrid';

% get bounds for capital guesses
rst = 1.0/beta - 1.0;

klwrbnd = (rst + delta)/(alpha);
klwrbnd = klwrbnd/(lagg^(1-alpha));
klwrbnd = klwrbnd^(1/(alpha-1));
klmult = 1.025;
khmult = 20;

kl = klwrbnd*klmult;
kh = klwrbnd*khmult;

% make asset grid

agrid = linspace(al, ah, na);

% make distribution grid

amu = linspace(al, ah, nmu);

% make party policy grid

tgrid = [.1 .5];
lamgrid = [.07 .21];

p = .5; % initial 50/50 chance of getting any party

%% solving for wages and r by getting aggregate capital

kDist = 10;
while kDist > kTol

    kval = .5*(kl + kh);

    fprintf("A guess: %4.4f. Begin iteration for solution...\n", kval)
    fprintf("\t Solving value function:\n")

    r = alpha*(kval^(alpha - 1)*(lagg^(1-alpha))) - delta;
    wage = (1-alpha)*((kval^(alpha))*(lagg^(-alpha)));

    %reassign whole thing just bc i'm debugging and idk what's wrong
    terms = struct('beta', beta, ...
        'sigma', sigma, ...
        'phi', phi, ...
        'identity', identity, ...
        'r', r, ...
        'wage', wage, ...
        'agrid', agrid, ...
        'lgrid', lgrid, ...
        'tgrid', tgrid, ...
        'lamgrid', lamgrid, ...
        'pil', pil, ...
        'p', p);

    [V, g, VOTES] = HH.solve(nl, na, np, terms, vTol);

    fprintf("\t Solving asset distribution:\n")
    [adistr, kagg] = HH.getDist(g, amu, agrid, pil, pctDem);

    f = kagg - kval;

    flatdistr=squeeze(sum(adistr,2));

    if f > 0
        fprintf("Aggregate capital is too low.\n\n")
        pause(1)
        figure
        mesh(flatdistr) 
        kl = .5*(kval+kl);
    else
        fprintf("Aggregate capital is too high.\n\n")
        pause(1)
        figure
        mesh(flatdistr) 
        kh = .5*(kval+kh);
    end

    kdist = abs(kagg - kval);
end

tiledlayout(4,1);
nexttile
mesh(V(:,:,1,1))
nexttile
mesh(V(:,:,2,1))
nexttile
mesh(VOTES(:,:,2))
nexttile
mesh(EV(:,:,1,1) - EV(:,:,1,2))


%% (B) getting perfect foresight impulse responses

T = 150;
f = @(rho, t, e0) rho^t*e0;
rho = .9;
e0 = 0.01;
k0 = 1;
lambda = .9;

At = ones(T,1);
for t = 1:T
    At(t) = f(rho, t, e0);
end

At = exp(At);
At = [ones(10,1); At];
At = At';
At = At(1:T);

params = [alpha beta delta sigma lagg];

[At, yt, ct, kt] = predict.perfectForesight(At,kagg,k0,lambda,params,vTol);

col1 = [.159 .072 .130]*(1/(.255));
col2 = [.255 .217 .120]*(1/(.255));
col3 = [.110 .182 .095]*(1/(.255));
col4 = [.086 .145 .128]*(1/(.255));

%plotting
tiledlayout(2,2)
nexttile
plot(At, 'Color', col1)
title("TFP")
nexttile
plot(yt, 'Color', col2)
title("Output")
nexttile
plot(kt, 'Color', col3)
title("Capital")
nexttile
plot(ct, 'Color', col4)
title("Consumption")

exportgraphics(gcf, "../v/Q2b_charts.png", 'Resolution', 500)

% need to figure out impulse responses 