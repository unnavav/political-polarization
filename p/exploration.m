%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setting up original household problem and solving it
% vaasavi
% may 2024
% 
% aiyagari + progressive tax policy
%
% Gonna use endogenous grid method to figure out the 
% HH problem . then from there, try to converge on a 
% lambda - r pair solution. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

restoredefaultpath;
clear all; clc;
addpath(genpath(pwd));

% first set up some grids, pulling a lot from aiyagari
% project
delete(gcp('nocreate'));
parpool('local',4);
addAttachedFiles(gcp, "dependencies")

load handel

vTol = 1e-5; gTol = 1e-6; dTol = 1e-3;
%% params
alpha = 0.36; delta = 0.08; beta = 0.96173; sigma = 7; phi = 0; goal = .3626;

nl = 7;
na = 100;
nmu = na*10;

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
klmult = 0;
khmult = 1.2;

kl = klwrbnd*klmult;
kh = klwrbnd*khmult;

% make asset grid

agrid = compute.logspace(al, ah, na)';

% make distribution grid

amu = linspace(al, ah, nmu);

% initial tau guess
tau = 0.181; %heathcote et al 2017

% lambda upper lower bounds
ll = 0; lh = 1;
lamval = .01; %based on some guess from a 3D plot I made; this gives transfers
adj = .5;

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
        cap_inc = repmat(r*(agrid),nl,1);
        Y = wage_inc+cap_inc;
        T = gov.tax(Y,lamval,tau);
    
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
    
        [V, g] = HH.solve(nl, na, terms, vTol, gTol);
    
        tiledlayout(2,1);
        nexttile
        mesh(V)
        nexttile
        mesh(g)

        % asset distr
        fprintf("\tSolving asset distribution:\n")
        [adistr, kagg] = HH.getDist(g, amu, agrid, pil);
 
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

    Y = kagg^(alpha)*lagg^(1-alpha);
    G_Y = t/Y;
    fprintf("\n Government/Output = %4.4f", G_Y);

    if G_Y<goal
       fprintf("\nGov't rev collected = %4.4f. Lam = %4.4f. " + ...
           "\tTax rate is too high.\n\n", t, lamval)
       lh = (lamval*adj+(1-adj)*lh);
    else
       fprintf("\nGov't rev collected = %4.4f. Lam = %4.4f. " + ...
           "\tTax rate is too low.\n\n", t, lamval);
       ll = (lamval*adj+(1-adj)*ll);
    end

    gDist = abs(G_Y - goal);
    lamval = .5*(ll + lh);

    DIST = max(kDist, gDist);
    fprintf("||DIST|| = %4.4f\n", DIST)
    kl = kval*klmult;
    kh = kval*2;
    kDist = 10;
end

sound(y, Fs)

%%
date = string(datetime("today"));
filename = strcat("results_",date);
cd ..\d
save(filename)