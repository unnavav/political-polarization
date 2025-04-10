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

% setting up the parallelized pool
% parpool('local', 16)

[ETA, DELTA] = ndgrid(linspace(0,.5, 10), linspace(0,.1, 10));

EGRID = ETA+DELTA;

[TAU, DELTA] = ndgrid(linspace(0, .3, 10), linspace(0.1,.6, 10));

TGRID = TAU + DELTA;

% first set up some grids, pulling a lot from aiyagari
% project

vTol = 1e-5; dTol = 1e-3;
%% params
alpha = 0.36; delta = 0.06; beta = 0.96; sigma = 2; phi = 0;
 
nl = 7;
na = 250;
nmu = na*10;
np = 2;

al = 0+phi; ah = 100+phi;

% initial r guess: 
r = 0.04;

% get labor distribution and aggregate values
%step 1: make labor grid and labor transition matrix

mu = 0;
rho = .5;
sig_l = 0.65; %0.5327;
 
% need to back out sigma^2_e given sigma^2_l
sigx = sig_l*sqrt(1 - rho^2);

[pil, lgrid] = compute.getTauchen(nl,  mu, sigx, rho);

% get steady state labor distribution

ldist = asymptotics(dtmc(pil));

lagg = ldist*lgrid';

kagg = ((r+delta)/(alpha*lagg^(1-alpha)))^(1/(alpha-1));

% get bounds for capital guesses
rst = 1.0/beta - 1.0;

klwrbnd = (rst + delta)/(alpha);
klwrbnd = klwrbnd/(lagg^(1-alpha));
klwrbnd = klwrbnd^(1/(alpha-1));
klmult = .3;
khmult = 2;

kl = klwrbnd*klmult;
kh = klwrbnd*khmult;

% make asset grid

agrid = compute.logspace(al, ah, na)';

% make distribution grid

amu = linspace(al, ah, nmu);

% lambda upper lower bounds
ll = 0; lh = 1;
lamval = (lh + ll)/2; % people are really reactive to taxes
adj = 1/3;

%party tax regimes
tgrid = [0 .2]; 
% tgrid = [.15 .05]; 
eta = [1.00 1.1];

captax = repelem(0, nl);

G = [0 0];

[ni, nj] = size(EGRID);

K = cell(nj);
G_array = cell(nj);
adistr_array = cell(nj);
wage_array = cell(nj);
r_array = cell(nj);
EV_array = cell(nj);
%% eta 
folname = "eta_test";
mkdir("../d/",folname)
cd ../d/eta_test/

% doing something that will literally only work for 0-9
[ni, nj] = size(EGRID); 

dataset = zeros(nj, 12);

colnames = ["eta1" "eta2" "p" "K" "K_over_N" "r" "w" "Gini" "avg_a_l" "avg_a_p" ...
    "avg_inc_l" "avg_inc_p"];

for i = 1:nj

    fprintf("\n*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*\n")

    % first solve for individual steady states 

    adistr = zeros(nl, nmu);

    for im = 1:nmu
        for il = 1:nl
            adistr(il,im) = 1.0/(nmu*nl);
        end 
    end
    
    growth_lagg = lagg*(1+EGRID(1,i));

    kDist = 10;
    gDist = 10; 
    f=10;  
    
    %init from eqm i've already solved for
    kval = 8.75;
    kl = klwrbnd*klmult;
    kh = klwrbnd*khmult;
    
    DIST = max(kDist,gDist);

    while kDist > 0.05
    
        fprintf("\nA guess: %4.8f:", kval)
        
        % first getting equilibrium objects like prices, govt spending, etc
        r = alpha*(kval/growth_lagg)^(alpha - 1) - delta;
        wage = (1-alpha)*(kval/growth_lagg)^(alpha);
    
        fprintf("\t implied r = %1.4f", r);

        terms = struct('beta', beta, ...
            'sigma', sigma, ...
            'phi', phi, ...
            'agrid', agrid, ...
            'lgrid', lgrid, ...
            'pil', pil, ...
            'lamval', lamval, ...
            'captax', zeros(1,nl), ...
            'tau', 0, ...
            'G', 0, ...
            'K', kagg, ...
            'L', growth_lagg, ...
            'r', r, ...
            'w', wage);
       
        [V, G, ~, EV] = HH.solve(nl, na, terms, vTol, false);
     
        [adistr, kagg] = HH.getDist(G, amu, agrid, pil, false);
    
%         % CONDENSE DISTR
        acond = compute.condense(adistr, amu, agrid);
            
        kdist = kagg - kval;
    
        if kdist > 0
            fprintf(" ||Kguess - Kagg|| = %4.5f. \n\tAggregate capital is too low.\n", abs(kdist))
            kl = (kval+kl)/2;
        else
            fprintf(" ||Kguess - Kagg|| = %4.5f. \n\tAggregate capital is too high.\n", abs(kdist))
            kh = (kval+kh)/2;
        end
    
        kDist = abs(kdist);  
        kval = .5*(kl + kh);

    end

    K{i} = kagg;
    EV_array{i} = EV;
    G_array{i} = G;
    adistr_array{i}  = adistr;
    wage_array{i} = wage;
    r_array{i} = r;

    % use these steady states to back out decision rules
    VOTESp = EV_array{1} > EV_array{i};
    VOTESl = EV_array{i} > EV_array{1};

    acond_distr = compute.condense(adistr_array{1},amu,agrid);
    ap = acond_distr.*VOTESp; %weighted votes
    pp = sum(sum(ap));

    if pp >= .5 
        index = 1;
    else 
        index = i;
    end

    al = acond_distr.*VOTESl;
    algrid = sum(al,1).*agrid;
    welgrid = sum(ap,2).*(wage_array{1}.*lgrid)';
    albar = sum(algrid);
    welbar = sum(welgrid);

    apgrid = sum(ap,1).*agrid;
    wepgrid = sum(ap,2).*(wage_array{1}.*lgrid)';
    apbar = sum(apgrid);
    wepbar = sum(wepgrid);

    % gini coeff - stolen from aubhik bc i am dissociating rn
    muall = sum(adistr_array{index});
    a1 = amu.*muall;
    b1 = cumsum(a1);
    b1 = b1./b1(nmu);
    c1 = cumsum(muall);
    
    dc = c1(2:nmu) - c1(1:nmu-1);
    db = (b1(2:nmu) + b1(1:nmu-1))/2;
    B = sum(dc.*db);
    A = 0.5 - B;
    GiniC1 = 2.0*A;

    dataset(i,:) = [EGRID(1,1) EGRID(1,i) pp K{i} K{i}/growth_lagg ...
        r_array{i} wage_array{i} GiniC1 albar apbar welbar wepbar];

    dataset_table = array2table(dataset, 'VariableNames', colnames);
    writetable(dataset_table,"results_sigma_2_sigx_0pt65.csv")
end

cd ../../p

%% now changing growth paths

folname = "tau_test";
mkdir("../d/",folname)

% doing something that will literally only work for 0-9
[ni, nj] = size(EGRID); 

dataset = zeros(nj, 11);

colnames = ["tau1" "tau2" "p" "K" "r" "w" "Gini" "avg_a_l" "avg_a_p" ...
    "avg_inc_l" "avg_inc_p"];

for i = 1:nj

    fprintf("\n*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*\n")

    % first solve for individual steady states 

    adistr = zeros(nl, nmu);

    for im = 1:nmu
        for il = 1:nl
            adistr(il,im) = 1.0/(nmu*nl);
        end 
    end
    
    growth_lagg = lagg*(1+EGRID(1,i));

    kDist = 10;
    gDist = 10; 
    
    %init from eqm i've already solved for
    kval = 8.75;
    kl = klwrbnd*klmult;
    kh = klwrbnd*khmult;
    
    tau = TGRID(1, i);

    DIST = max(kDist,gDist);

    while kDist > 0.05
    
        fprintf("\nA guess: %4.8f:", kval)
        
        % first getting equilibrium objects like prices, govt spending, etc
        r = alpha*(kval^(alpha - 1)*(growth_lagg^(1-alpha))) - delta;
        wage = (1-alpha)*((kval^(alpha))*(growth_lagg^(-alpha)));
    
        wage_inc_mu = repmat(wage*lgrid, nmu,1)';
        cap_inc_mu = zeros(size(wage_inc_mu));
        for il = 1:nl
            cap_inc_mu(il, :) = repmat(r*amu, 1, 1)*(1-captax(il));
        end
        Tmu = zeros(nl, nmu);
        Tmu(:,:) = wage_inc_mu - gov.tax(wage_inc_mu,lamval,tau) + ...
            cap_inc_mu;

        g = sum(sum(Tmu(:,:).*adistr(:,:)));

        fprintf("\t implied r = %1.4f", r);

        terms = struct('beta', beta, ...
            'sigma', sigma, ...
            'phi', phi, ...
            'agrid', agrid, ...
            'lgrid', lgrid, ...
            'pil', pil, ...
            'lamval', lamval, ...
            'captax', zeros(1,nl), ...
            'tau', tau, ...
            'G', 0, ...
            'K', kagg, ...
            'L', growth_lagg, ...
            'r', r, ...
            'w', wage);
       
        [V, G, ~, EV] = HH.solve(nl, na, terms, vTol, false);
     
        [adistr, kagg] = HH.getDist(G, amu, agrid, pil, false);
    
%         % CONDENSE DISTR
        acond = compute.condense(adistr, amu, agrid);
            
        kdist = kagg - kval;
    
        if kdist > 0
            fprintf(" ||Kguess - Kagg|| = %4.5f. \n\tAggregate capital is too low.\n", abs(kdist))
            kl = (kval+kl)/2;
        else
            fprintf(" ||Kguess - Kagg|| = %4.5f. \n\tAggregate capital is too high.\n", abs(kdist))
            kh = (kval+kh)/2;
        end
    
        kDist = abs(kdist);  
        kval = .5*(kl + kh);

    end

    K{i} = kagg;
    EV_array{i} = EV;
    G_array{i} = G;
    adistr_array{i}  = adistr;
    wage_array{i} = wage;
    r_array{i} = r;

    % use these steady states to back out decision rules
    VOTESp = EV_array{1} > EV_array{i};
    VOTESl = EV_array{i} > EV_array{1};

    acond_distr = compute.condense(adistr_array{1},amu,agrid);
    ap = acond_distr.*VOTESp; %weighted votes
    pp = sum(sum(ap));

    if pp >= .5 
        index = 1;
    else 
        index = i;
    end

    al = acond_distr.*VOTESl;
    algrid = sum(al,1).*agrid;
    welgrid = sum(ap,2).*(wage_array{1}.*lgrid)';
    albar = sum(algrid);
    welbar = sum(welgrid);

    apgrid = sum(ap,1).*agrid;
    wepgrid = sum(ap,2).*(wage_array{1}.*lgrid)';
    apbar = sum(apgrid);
    wepbar = sum(wepgrid);

    % gini coeff - stolen from aubhik bc i am dissociating rn
    muall = sum(adistr_array{index});
    a1 = amu.*muall;
    b1 = cumsum(a1);
    b1 = b1./b1(nmu);
    c1 = cumsum(muall);
    
    dc = c1(2:nmu) - c1(1:nmu-1);
    db = (b1(2:nmu) + b1(1:nmu-1))/2;
    B = sum(dc.*db);
    A = 0.5 - B;
    GiniC1 = 2.0*A;

    dataset(i,:) = [EGRID(1,1) EGRID(1,i) pp K{i} ...
        r_array{i} wage_array{i} GiniC1 albar apbar welbar wepbar];

    dataset_table = array2table(dataset, 'VariableNames', colnames);
    writetable(dataset_table,"../d/tau_test/results_no_remittance.csv")
end
