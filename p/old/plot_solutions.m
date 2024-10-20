% This code looks at the various solutions to the problem given a grid of
% capital and tax values. It's to help me determine if the problem has
% concavity or whether it's lumpy.

restoredefaultpath;
clear all; clc;
addpath(genpath(pwd));

load handel.mat

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
klmult = .5;
khmult = 2;

kl = klwrbnd*klmult;
kh = klwrbnd*khmult;

% make asset grid

agrid = linspace(al, ah, na);

% make distribution grid

amu = linspace(al, ah, nmu);

% initial tau guess
tau = 0.08; %heathcote et al 2017, actual tax rate

% lambda upper lower bounds
ll = 0; lh = .9;

%% set grids and run

np = 10;

lamgrid = linspace(ll, lh, np);
kapgrid = linspace(kl, kh, np);

Gvals = zeros(np, np);

Kvals = zeros(np);

kDist = 10;

for ip = 1:np

    lamval = lamgrid(ip);
    fprintf("\n\n\n--------------------------------- > Lambda = %4.4f", lamval)

    for ip2 = 1:np
        kval = kapgrid(ip2);
    
    
        fprintf("\n*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*")
        fprintf("\nA guess: %4.4f. Begin iteration for solution...\n", kval)
        fprintf("\t Solving value function:\n")
    
        r = alpha*(kval^(alpha - 1)*(lagg^(1-alpha))) - delta;
        wage = (1-alpha)*((kval^(alpha))*(lagg^(-alpha)));
    
        % making tax schedule
        
        wage_inc = repmat(wage*lgrid,na,1)';
        cap_inc = repmat(r*(agrid-1),nl,1);
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
    
        [V, g] = HH.solve(nl, na, terms, vTol);
    
        % asset distr
        fprintf("\tSolving asset distribution:\n")
        [adistr, kagg] = HH.getDist(g, amu, agrid, pil);
 
        t = adistr.*Tmu; %getting all taxes collected
        t = sum(sum(t));


        Y = kagg^(alpha)*lagg^(1-alpha);
        G_Y = t/Y;

        Gvals(ip, ip2) = G_Y;
    end
end

save ..\d\GYvals_pt08 Gvals 

sound(y, Fs);

%%

for ip = 1:np

    lamval = lamgrid(ip);
    fprintf("\n\n\n--------------------------------- > Lambda = %4.4f", lamval)

    while kDist > dTol
        kval = kapgrid(ip2);
    
    
        fprintf("\n*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*")
        fprintf("\nA guess: %4.4f. Begin iteration for solution...\n", kval)
        fprintf("\t Solving value function:\n")
    
        r = alpha*(kval^(alpha - 1)*(lagg^(1-alpha))) - delta;
        wage = (1-alpha)*((kval^(alpha))*(lagg^(-alpha)));
    
        % making tax schedule
        
        wage_inc = repmat(wage*lgrid,na,1)';
        cap_inc = repmat(r*(agrid-1),nl,1);
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
    
        [V, g] = HH.solve(nl, na, terms, vTol);
    
        % asset distr
        fprintf("\tSolving asset distribution:\n")
        [adistr, kagg] = HH.getDist(g, amu, agrid, pil);
    
        t = adistr.*Tmu; %getting all taxes collected
        t = sum(sum(t));
    
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

    Y = kagg^(alpha)*lagg^(1-alpha);
    G_Y = t/Y;

    Gvals(ip, ip2) = G_Y;
end