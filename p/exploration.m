%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setting up original household problem and solving it
% vaasavi
% may 2024
% 
% aiyagari + progressive tax policy
%
% The progressive tax policy is given and calibrated by
% Heathcote et al. 2017. I use the values of 
% progressivity already calibrated from the US tax
% schedule
% 
% The code is structured as follows: 
% 1. Make guess of interest rate and lambda
% 2. Back out aggregate capital and government
% 3. If not correct, adjust.
%
% In particular, the choice of lambda comes from 
% determining whether Y - C - I = G = 0 in equilibrium.
% I add up tax revenue with a given lambda and then 
% determine if the lambda is the correct value. If
% revenue is too high, lower lambda, and vice versa.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

restoredefaultpath;
clear all; clc;
addpath(genpath(pwd));

% first set up some grids, pulling a lot from aiyagari
% project
% delete(gcp('nocreate'));
% parpool('local',4);

vTol = 1e-6; gTol = 1e-8; dTol = 1e-3;
%% params
alpha = 0.36; delta = 0.08; beta = 0.96173; sigma = 1; phi = 0;

nl = 7;
na = 250;
nmu = 2500;

al = 0+phi; ah = 100+phi;

wage = 1.1; r = 0.04;

%political params
identity = 0.01; pctDem = .5;
goal = .1;
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
klmult = -0.3;
khmult = 2;

kl = klwrbnd*klmult;
kh = klwrbnd*khmult;

% make asset grid

agrid = linspace(al, ah, na);

% make distribution grid

amu = linspace(al, ah, nmu);

% initial tau guess
tau = 0.181; %heathcote et al 2017

% lambda upper lower bounds
ll = 0; lh = 1;
lamval = .8; %based on some guess from a 3D plot I made; this gives transfers
adj = .05;

%% solving for wages and r by getting aggregate capital

kDist = 10;
gDist = 10;
f=10;

kval = (kl+ kh)/2;

DIST = max(kDist,gDist);

while DIST > dTol

    fprintf("\n\n\n--------------------------------- > Lambda = %4.4f", lamval)

    while kDist > dTol
    
        fprintf("\n*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*")
        fprintf("\nA guess: %4.4f. Begin iteration for solution...\n", kval)
        fprintf("\t Solving value function:\n")
    
        r = alpha*(kval^(alpha - 1)*(lagg^(1-alpha))) - delta;
        wage = (1-alpha)*((kval^(alpha))*(lagg^(-alpha)));
    
        % making tax schedule
        
        wage_inc = repmat(wage*lgrid,na,1)';
        cap_inc = repmat(r*(agrid-1),nl,1);
        Y = wage_inc+cap_inc;
        T = gov.tax(Y,lamval,tau)-2;
    
        wage_inc_mu = repmat(wage*lgrid, nmu,1)';
        cap_inc_mu = repmat(r*amu, nl, 1);
        Ymu = wage_inc_mu + cap_inc_mu;
        Tmu = gov.tax(Ymu, lamval, tau);
    
        %prepare for VFI
        terms = struct('beta', beta, ...
            'sigma', sigma, ...
            'phi', phi, ...
            'r', r, ...
            'wage', wage, ...
            'agrid', agrid, ...
            'lgrid', lgrid, ...
            'pil', pil, ...
            'T', T);
    
        [V, g] = HH.solve(nl, na, terms, vTol);
    
        % asset distr
        fprintf("\tSolving asset distribution:\n")
        [adistr, kagg] = HH.getDist(g, amu, agrid, pil);
 
%         % CHECK:
%         % before the tax schedule is right, this is not going to converge to
%         % the correct equilibrium. So what I do instead is to see the
%         % amount the capital value has changed, and when it's no longer
%         % changing, exit the loop and update the tax schedule. 
%         pct_chg_diff = (kagg-kval)/f;
        f = kagg - kval;
    
        if f > 0
            fprintf("\n||Kguess - Kagg|| = %4.4f. \tAggregate capital is too low.", abs(f))
            kl = .5*(kval+kl);
        else
            fprintf("\n||Kguess - Kagg|| = %4.4f. \tAggregate capital is too high.", abs(f))
            kh = .5*(kval+kh);
        end

        kval = .5*(kl + kh);
    
        kDist = abs(f);  % check whether the capital diff
                                        % is changing at all
    end
    
    t = adistr.*Tmu; %getting all taxes collected
    t = sum(sum(t));

    if t>goal
       fprintf("\nGov't rev collected = %4.4f. Lam = %4.4f. " + ...
           "\tTax rate is too high.\n\n", t, lamval)
       lh = (lamval*adj+(1-adj)*lh);
    else
       fprintf("\nGov't rev collected = %4.4f. Lam = %4.4f. " + ...
           "\tTax rate is too low.\n\n", t, lamval);
       ll = (lamval*adj+(1-adj)*ll);
    end

    gDist = abs(t);
    lamval = .5*(ll + lh);

    DIST = max(kDist, gDist);
    fprintf("||DIST|| = %4.4f\n", DIST)
    kl = klwrbnd*klmult;
    kh = klwrbnd*khmult;
    kDist = 10;
end

tiledlayout(3,1);
nexttile
mesh(V)
nexttile
mesh(EV)
nexttile
mesh(g)
% nexttile
% mesh(EV(:,:,1,1) - EV(:,:,1,2))

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