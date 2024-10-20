%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setting up original household problem and solving it
% vaasavi
% may 2024
% 
% aiyagari + progressive tax policy
%
% testing different scenarios to see when republicans
% win vs when democrats win
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

restoredefaultpath;
clear all; clc;
addpath(genpath(pwd));

% first set up some grids, pulling a lot from aiyagari
% project

load handel

vTol = 1e-5; dTol = 1e-3;
%% params
alpha = 0.36; delta = 0.06; beta = 0.96; sigma = 1; phi = 0;
 
nl = 7;
na = 250;
nmu = na*10;
np = 2;

al = 0+phi; ah = 100+phi;

wage = 1.1; r = 0.04; % original guesses

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
klmult = .3;
khmult = 1.5;

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
% tgrid = [1-.181 1-.086]; 
tgrid = [.086 .181]; 

p = [1 0];

num_sens = 5;

pctA = .5;

% save for later
cd ../d/
save agrid.mat agrid
save lgrid.mat lgrid
save amu.mat amu
cd ../p

captax = [repelem(0, ceil(nl/2)) repelem(.15, floor(nl/2))]; %IRS
% captax = compute.logspace(0, 20, nl)/100;
% captax = linspace(0, .20, nl);
% captax = repelem(.15, nl);
goal = .3652; % IMF for US govt spending (G/Y)

G = [0 0];

ul = .15;
uh = 0.25;
uvals = linspace(ul, uh, num_sens);

%% scenario 1: varying tax rates

adistrA = zeros(nl, nmu);
adistrB = zeros(nl, nmu);
adistr = zeros(nl, nmu, np);

for i = 1:(nmu/2)
    for il = 1:nl
            adistrA(il,i) = 1.0/(nmu*nl/2);
            adistrB(il,i) = 1.0/(nmu*nl/2);
    end
end

adistr(:,:,1) = adistrA;
adistr(:,:,2) = adistrB;

shares = zeros(num_sens, np+3);

tl = .95;
th = .975;
tvals = linspace(tl, th, num_sens);

ubonus = 0.1;
ubonusA = [ubonus 0];
ubonusB = [0 ubonus];

for is = 1:num_sens

    tgrid_i = [0 .95];

    kDist = 10;
    gDist = 10; 
    f=10;  

    %init from eqm i've already solved for
    kval = 8.75;
    kh = ah;
    kl = al;
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
            r = alpha*(kval^(alpha - 1)*(lagg^(1-alpha))) - delta;
            wage = (1-alpha)*((kval^(alpha))*(lagg^(-alpha)));
    
            wage_inc_mu = repmat(wage*lgrid, nmu,1)';
            cap_inc_mu = zeros(size(wage_inc_mu));
            for il = 1:nl
                cap_inc_mu(il, :) = repmat(r*amu, 1, 1)*(1-captax(il));
            end
            Tmu = zeros(nl, nmu, np);
            for ip = 1:np
                Tmu(:,:,ip) = gov.tax(wage_inc_mu,lamval,tgrid_i(ip)) + ...
                    cap_inc_mu;
            end
    
            adistr(:,:,1) = adistrA;
            adistr(:,:,2) = adistrB;
            for ip = 1:np
                G(ip) = sum(sum(Tmu(:,:,ip).*adistr(:,:,ip)));
            end
    
    
            %prepare for VFI
            terms = struct('beta', beta, ...
                'sigma', sigma, ...
                'phi', phi, ...
                'r', r, ...
                'wage', wage, ...
                'agrid', agrid, ...
                'lgrid', lgrid, ...
                'pil', pil, ...
                'lamval', lamval, ...
                'captax', captax, ...
                'tau', tgrid_i, ...
                'G', G, ...
                'p', p);
        
            fprintf("Solving value functions:\n")
    
            fprintf("\tHH Type A:\n")
            terms.ubonus = ubonusA;
            [Va, Ga, EVa] = HH.solve(nl, na, np, terms, vTol);
        
            fprintf("\tHH Type B:\n")
            terms.ubonus = ubonusB;
            [Vb, Gb, EVb] = HH.solve(nl, na, np, terms, vTol);
            
            mesh(EVa(:,:,1) - EVa(:,:,2))
            VOTESa = (EVa(:,:,1) > EVa(:,:,2));
            VOTESb = (EVb(:,:,1) > EVb(:,:,2));
        
            test = EVa(:,:,1) - EVb(:,:,1);
        
            %TO DO: 
            % SOLVE ASSET DISTRIBUTION 
            % SOLVE VOTE COUNTS
    
            fprintf("\tSolving asset distribution:\n")
            [adistrA, kaggA, adistrB, kaggB] = HH.getDist(Ga, Gb, amu, agrid, pil, pctA);
        
    %         % CONDENSE DISTR
            acondA = compute.condense(adistrA, amu, agrid);
            acondB = compute.condense(adistrB, amu, agrid);
        
            share2A = pctA*sum(sum(VOTESa.*acondA)) + ...
                (1-pctA)*sum(sum(VOTESb.*acondA));
            share2B = pctA*sum(sum(VOTESa.*acondB)) + ...
                (1-pctA)*sum(sum(VOTESb.*acondB));
    %         p(1) = share2A > share2B;
    %         p(2) = share2A > share2B;
    
            
%             f = kaggA - kval;
            f = kaggB - kval;
        
            if f > 0
                fprintf("\n||Kguess - Kagg|| = %4.5f. \tAggregate capital is too low.\n", abs(f))
                kl = (kval+2*kl)/3;
            else
                fprintf("\n||Kguess - Kagg|| = %4.5f. \tAggregate capital is too high.\n", abs(f))
                kh = (kval+2*kh)/3;
            end
        
%             kDist = abs(kaggA - kval);  
            kDist = abs(kaggB - kval);  
            kval = .5*(kl + kh);
        end
        
        tA = adistrA.*Tmu(:,:,1); %getting all taxes collected in A govt
        tA = sum(sum(tA));
        tB = adistrB.*Tmu(:,:,2); %getting all taxes collected in B govt
        tB = sum(sum(tB));
        Y = kval^(alpha)*lagg^(1-alpha);
        gy = tB/Y;
    
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

    sound(y, Fs)

    folname = strcat("tax_scheme_Beqm_",sprintf('%0.4f', tgrid_i(1)), ...
        sprintf('%0.4f', tgrid_i(2)));
    mkdir("../d/",folname)

    newdir = strcat("../d/", folname);
    cd(newdir)

    filename = 'results_pval10';
    save(filename, 'adistrA', 'adistrB', 'share2A', 'share2B', ...
        "Ga","Gb", "Va", "Vb", "VOTESa", "VOTESb", "wage", "r", "kaggA", ...
        "shares", "lamval", "ubonus", 'acondA', 'acondB')

    cd ../../p/

end

%% scenario 2: varying personal gains

adistrA = zeros(nl, nmu);
adistrB = zeros(nl, nmu);
adistr = zeros(nl, nmu, np);

for i = 1:(nmu/2)
    for il = 1:nl
            adistrA(il,i) = 1.0/(nmu*nl/2);
            adistrB(il,i) = 1.0/(nmu*nl/2);
    end
end

adistr(:,:,1) = adistrA;
adistr(:,:,2) = adistrB;

shares = zeros(num_sens, np+3);

for is = 1:num_sens


    ubonus = uvals(is);
    ubonusA = [ubonus 0];
    ubonusB = [0 ubonus];
    
    kDist = 10;
    gDist = 10; 
    f=10;  
    
    kval = 8.75;
    kh = kval*1.5;
    kl = kval*.5;
    lamval = .7023;
    ll = lamval*.7; lh = lamval*1.3;

    DIST = max(kDist,gDist);
    
    while DIST > dTol*10
    
        fprintf("\n\n\n--------------------------------- > Lambda = %4.4f", lamval)
    
        while kDist > dTol
        
            fprintf("\n*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*")
            fprintf("\nA guess: %4.8f. Begin iteration for solution...\n", kval)
            fprintf("\t Solving value function:\n")
            
            % first getting equilibrium objects like prices, govt spending, etc
            r = alpha*(kval^(alpha - 1)*(lagg^(1-alpha))) - delta;
            wage = (1-alpha)*((kval^(alpha))*(lagg^(-alpha)));
    
            wage_inc_mu = repmat(wage*lgrid, nmu,1)';
            cap_inc_mu = zeros(size(wage_inc_mu));
            for il = 1:nl
                cap_inc_mu(il, :) = repmat(r*amu, 1, 1)*(1-captax(il));
            end
            Tmu = zeros(nl, nmu, np);
            for ip = 1:np
                Tmu(:,:,ip) = gov.tax(wage_inc_mu,lamval,tgrid(ip)) + ...
                    cap_inc_mu;
            end
    
            adistr(:,:,1) = adistrA;
            adistr(:,:,2) = adistrB;
            for ip = 1:np
                G(ip) = sum(sum(Tmu(:,:,ip).*adistr(:,:,ip)));
            end
    
    
            %prepare for VFI
            terms = struct('beta', beta, ...
                'sigma', sigma, ...
                'phi', phi, ...
                'r', r, ...
                'wage', wage, ...
                'agrid', agrid, ...
                'lgrid', lgrid, ...
                'pil', pil, ...
                'lamval', lamval, ...
                'captax', captax, ...
                'tau', tgrid, ...
                'G', G, ...
                'p', p);
        
            fprintf("Solving value functions:\n")
    
            fprintf("\tHH Type A:\n")
            terms.ubonus = ubonusA;
            [Va, Ga, EVa] = HH.solve(nl, na, np, terms, vTol);
        
            fprintf("\tHH Type B:\n")
            terms.ubonus = ubonusB;
            [Vb, Gb, EVb] = HH.solve(nl, na, np, terms, vTol);
            
            mesh(EVa(:,:,1) - EVa(:,:,2))
            VOTESa = (EVa(:,:,1) > EVa(:,:,2));
            VOTESb = (EVb(:,:,1) > EVb(:,:,2));
        
            test = EVa(:,:,1) - EVb(:,:,1);
        
            %TO DO: 
            % SOLVE ASSET DISTRIBUTION 
            % SOLVE VOTE COUNTS
    
            fprintf("\tSolving asset distribution:\n")
            [adistrA, kaggA, adistrB, kaggB] = HH.getDist(Ga, Gb, amu, agrid, pil, pctA);
        
    %         % CONDENSE DISTR
            acondA = compute.condense(adistrA, amu, agrid);
            acondB = compute.condense(adistrB, amu, agrid);
        
            share2A = pctA*sum(sum(VOTESa.*acondA)) + ...
                (1-pctA)*sum(sum(VOTESb.*acondA));
            share2B = pctA*sum(sum(VOTESa.*acondB)) + ...
                (1-pctA)*sum(sum(VOTESb.*acondB));
    %         p(1) = share2A > share2B;
    %         p(2) = share2A > share2B;
    
            % 
            f = kaggA - kval;
%             f = kaggB - kval;
        
            if f > 0
                fprintf("\n||Kguess - Kagg|| = %4.5f. \tAggregate capital is too low.\n", abs(f))
                kl = (kval+2*kl)/3;
            else
                fprintf("\n||Kguess - Kagg|| = %4.5f. \tAggregate capital is too high.\n", abs(f))
                kh = (kval+2*kh)/3;
            end
        
            kDist = abs(kaggA - kval);  
            kval = .5*(kl + kh);
        end
        
        tA = adistrA.*Tmu(:,:,1); %getting all taxes collected in A govt
        tA = sum(sum(tA));
        tB = adistrB.*Tmu(:,:,2); %getting all taxes collected in B govt
        tB = sum(sum(tB));
        Y = kval^(alpha)*lagg^(1-alpha);
        gy = tA/Y;
    
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

    sound(y, Fs)

    folname = strcat("personal_return",sprintf('%0.2f', uvals(is)));
    mkdir('../d/',folname)

    newdir = strcat("../d/", folname);
    cd(newdir)

    filename ='results';
    save(filename, 'adistrA', 'adistrB', 'share2A', 'share2B', ...
        "Ga","Gb", "Va", "Vb", "VOTESa", "VOTESb", "wage", "r", "kaggA", ...
        "kaggB","lamval")

    cd ../../p/
end




