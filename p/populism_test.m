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

vTol = 1e-5; dTol = 1e-3;
%% params
alpha = 0.36; delta = 0.06; beta = 0.96; sigma = 1; phi = 0;

eta = 1.0;
 
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
sig_l = 0.5327;
 
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
khmult = 6;

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
etagrid = [1.00 1.1];

p = .7;

% captax = [repelem(0, ceil(nl/2)) repelem(.15, floor(nl/2))]; %IRS
% captax = compute.logspace(0, 20, nl)/100;
% captax = linspace(0, .20, nl);
% captax = repelem(.15, nl);
% captax = repelem(0, nl);
captax = repelem(0, nl);
goal = .3652; % IMF for US govt spending (G/Y)

G = [0 0];

folname = "detrended";
mkdir("../d/",folname)

for i = 1:2

    adistr = zeros(nl, nmu);
    
    for im = 1:nmu
        for il = 1:nl
            adistr(il,im) = 1.0/(nmu*nl);
        end 
    end
    
    eta = etagrid(i);

    kDist = 10;
    gDist = 10; 
    f=10;  
    
    %init from eqm i've already solved for
    kval = 8.75;
    lamval = 0;
    ll = .25; lh = 1;
    
    DIST = max(kDist,gDist);
    
    while DIST > dTol*10
    
        fprintf("\n\n\n--------------------------------- > Lambda = %4.4f", lamval)
    
        while kDist > dTol
        
            fprintf("\n*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*")
            fprintf("\nA guess: %4.8f. Begin iteration for solution...\n", kval)
            fprintf("\t Solving value function:\n")
            
            % first getting equilibrium objects like prices, govt spending, etc
            r = alpha*(kval/(1+eta))^(alpha - 1) - delta;
            wage = (kval/(1+eta))^(alpha)*(1-alpha*(1+eta));
    
            wage_inc_mu = repmat(wage*lgrid, nmu,1)';
            cap_inc_mu = zeros(size(wage_inc_mu));
            for il = 1:nl
                cap_inc_mu(il, :) = repmat(r*amu, 1, 1)*(1-captax(il));
            end
            Tmu = zeros(nl, nmu, np);
            for ip = 1:np
                Tmu(:,:,ip) = wage_inc_mu - gov.tax(wage_inc_mu,lamval,tgrid(ip)) + ...
                    cap_inc_mu;
            end
    
%            G = sum(sum(Tmu(:,:,i).*adistr(:,:)));
    
            %prepare for VFI
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
                'p', p, ...
                'K', kagg, ...
                'L', 1, ...
                'r', r, ...
                'w', wage);
        
            fprintf("Solving value function:\n")
    
            [V, G, ~, EV] = HH.solve(nl, na, terms, vTol, false);
         
    
            fprintf("\tSolving asset distribution:\n")
            [adistr, kagg] = HH.getDist(G, amu, agrid, pil, true);
        
    %         % CONDENSE DISTR
            acond = compute.condense(adistr, amu, agrid);
                
            kdist = kagg - kval;
        
            if kdist > 0
                fprintf("\n||Kguess - Kagg|| = %4.5f. \tAggregate capital is too low.\n", abs(kdist))
                kl = (kval+kl)/2;
            else
                fprintf("\n||Kguess - Kagg|| = %4.5f. \tAggregate capital is too high.\n", abs(kdist))
                kh = (kval+kh)/2;
            end
        
            kDist = abs(kdist);  
            kval = .5*(kl + kh);
        end
        
        t = adistr.*Tmu(:,:,i); %getting all taxes collected in A govt
        t = sum(sum(t));
        Y = kval^(alpha)*growth_lagg^(1-alpha);
        gy = t/Y;
    
        if gy<goal
           fprintf("\nGov't rev collected = %4.4f. Lam = %4.4f. " + ...
               "\tTax rate is too high.\n\n", gy, lamval)
           lh = (lamval*adj+(1-adj)*lh);
        else
           fprintf("\nGov't rev collected = %4.4f. Lam = %4.4f. " + ...
               "\tTax rate is too low.\n\n", gy, lamval);
           ll = (lamval*adj+(1-adj)*ll);
        end
    
        gDist = abs(gy-goal);
        lamval = .5*(ll + lh);
    
        DIST = max(kDist, gDist);
        fprintf("||DIST|| = %4.4f\n", DIST)
        kl = klwrbnd*klmult;
        kh = klwrbnd*khmult;
        kDist = 10;
    end
    
    newdir = strcat("../d/", folname);
    cd(newdir)

    filename = strcat("results_t",sprintf('%0.4f', 0), ...
        "_eta", sprintf('%0.4f', eta(i)), ".mat");
    save(filename, 'adistr', "G","V", "wage", "r", "kagg", ...
        "lamval", "eta", 'tgrid', 'growth_lagg', 'delta', 'alpha')

    filename = strcat("terms_struct_t", sprintf('%0.4f', 0), ...
        "_eta", sprintf('%0.4f', eta(i)), ".mat");
    save(filename, "terms")

    cd ../../p/

end

%% now changing growth paths

for i = 1:2

    adistr = zeros(nl, nmu);
    
    for im = 1:nmu
        for il = 1:nl
            adistr(il,im) = 1.0/(nmu*nl);
        end
    end
    
    growth_lagg = lagg*eta(1)*eta(2);


    kDist = 10;
    gDist = 10; 
    f=10;  
    
    %init from eqm i've already solved for
    kval = 8.75;
    lamval = .6221;
    ll = .25; lh = 1;
    
    DIST = max(kDist,gDist);
    
    while DIST > dTol*10
    
        fprintf("\n\n\n--------------------------------- > Lambda = %4.4f", lamval)
    
        while kDist > dTol
        
            fprintf("\n*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*")
            fprintf("\nA guess: %4.8f. Begin iteration for solution...\n", kval)
            fprintf("\t Solving value function:\n")
            
            % first getting equilibrium objects like prices, govt spending, etc
            r = alpha*(kval^(alpha - 1)*(growth_lagg^(1-alpha))) - delta;
            wage = (1-alpha)*((kval^(alpha))*(growth_lagg^(-alpha)));
    
            wage_inc_mu = repmat(wage*lgrid, nmu,1)';
            cap_inc_mu = zeros(size(wage_inc_mu));
            for il = 1:nl
                cap_inc_mu(il, :) = repmat(r*amu, 1, 1)*(1-captax(il));
            end
            Tmu = zeros(nl, nmu, np);
            for ip = 1:np
                Tmu(:,:,ip) = wage_inc_mu - gov.tax(wage_inc_mu,lamval,tgrid(ip)) + ...
                    cap_inc_mu;
            end
    
            G = sum(sum(Tmu(:,:,i).*adistr(:,:)));
    
            %prepare for VFI
            terms = struct('alpha', alpha, ...
                'beta', beta, ...
                'sigma', sigma, ...
                'phi', phi, ...
                'agrid', agrid, ...
                'lgrid', lgrid, ...
                'pil', pil, ...
                'lamval', lamval, ...
                'captax', captax, ...
                'tau', tgrid(i), ...
                'G', G, ...
                'p', p, ...
                'K', kagg, ...
                'L', growth_lagg, ...
                'r', r, ...
                'w', wage);
        
            fprintf("Solving value function:\n")
    
            [V, G, EV] = HH.solve(nl, na, terms, vTol);
         
    
            fprintf("\tSolving asset distribution:\n")
            [adistr, kagg] = HH.getDist(G, amu, agrid, pil);
        
    %         % CONDENSE DISTR
            acond = compute.condense(adistr, amu, agrid);
                
            kdist = kagg - kval;
        
            if kdist > 0
                fprintf("\n||Kguess - Kagg|| = %4.5f. \tAggregate capital is too low.\n", abs(kdist))
                kl = (kval+5*kl)/6;
            else
                fprintf("\n||Kguess - Kagg|| = %4.5f. \tAggregate capital is too high.\n", abs(kdist))
                kh = (kval+5*kh)/6;
            end
        
            kDist = abs(kdist);  
            kval = .5*(kl + kh);
        end
        
        t = adistr.*Tmu(:,:,i); %getting all taxes collected in A govt
        t = sum(sum(t));
        Y = kval^(alpha)*growth_lagg^(1-alpha);
        gy = t/Y;
    
        if gy<goal
           fprintf("\nGov't rev collected = %4.4f. Lam = %4.4f. " + ...
               "\tTax rate is too high.\n\n", gy, lamval)
           lh = (lamval*adj+(1-adj)*lh);
        else
           fprintf("\nGov't rev collected = %4.4f. Lam = %4.4f. " + ...
               "\tTax rate is too low.\n\n", gy, lamval);
           ll = (lamval*adj+(1-adj)*ll);
        end
    
        gDist = abs(gy-goal);
        lamval = .5*(ll + lh);
    
        DIST = max(kDist, gDist);
        fprintf("||DIST|| = %4.4f\n", DIST)
        kl = klwrbnd*klmult;
        kh = klwrbnd*khmult;
        kDist = 10;
    end
    
    newdir = strcat("../d/", folname);
    cd(newdir)

    filename = strcat("results_transition_t",sprintf('%0.4f', tgrid(i)), ...
        "_eta", sprintf('%0.4f', eta(i)), ".mat");
    save(filename, 'adistr', "G","V", "wage", "r", "kagg", ...
        "lamval", "eta", 'tgrid')

    save("params", "terms")

    cd ../../p/

end

