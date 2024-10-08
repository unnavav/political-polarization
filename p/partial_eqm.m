%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARTIAL EQIULIBRIA 
% vaasavi
% sept 2024
% 
% this just looks at partial equilibrium, open economy outcomes for voting 
% rules of A&B houesholds to see if we need to throw the baby out with the
% bathwater
%
% capital market doesn't clear, just labor market (which is 'degenerate', 
% idk the right word) keeping a fixed r
% 
% outputs:
%     - graphs of voting rules without any social benefits
%     - graphs of voting rules with social benefits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

restoredefaultpath;
clear all; clc;
addpath(genpath(pwd));

%% params
vTol = 1e-5; dTol = 1e-3;

alpha = 0.36; delta = 0.06; beta = 0.96; sigma = 1; phi = 0;
 
nl = 7;
na = 250;
nmu = na*10;
np = 2;

al = 0+phi; ah = 100+phi;


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

r = rst; % original guesses
kagg = (r/alpha)^(1/(alpha-1));
wage = (1-alpha)*kagg^(alpha);

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
% captax = [repelem(0, ceil(nl/2)) repelem(.15, floor(nl/2))]; %IRS
captax = [repelem(0, nl)]; %IRS
goal = .3652; % IMF for US govt spending (G/Y

p = [1 0];

num_sens = 25;

tl = .0;
th = .8;
tvals = linspace(tl, th, num_sens);

pctA = .5;

%%% now do the partial eqm analysis


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

ul = 0;
uh = .15;
uvals = linspace(ul, uh, num_sens);

folname = "partial_eqm_test";
mkdir("../d/",folname)

newdir = strcat("../d/", folname);
cd(newdir)

for is = 1:num_sens

    tgrid_i = [.0860 .181];
    
    ubonus = uvals(is);
    ubonusA = [ubonus 0];
    ubonusB = [0 ubonus];

    gDist = 10; 
    f=10;  
    
    kval = kagg;
    
    DIST = max(0,gDist);
    
    % lambda upper lower bounds
    lamval = .7026; % OG guess
    ll = lamval*.6; lh = lamval*1.4;

    while DIST > dTol*10
    
        fprintf("\n\n\n--------------------------------- > Lambda = %4.4f\n", lamval)
       
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

%         G = [0 0];

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
        
        tA = adistrA.*Tmu(:,:,1); %getting all taxes collected in A govt
        tA = sum(sum(tA));
        Y = kagg^(alpha)*lagg^(1-alpha);
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
    
        DIST = gDist;
        fprintf("||DIST|| = %4.4f\n", DIST)
    end

    shares(is,:) = [uvals(is) share2A share2B wage r];

    filename = strcat('results_',sprintf("%0.2f",uvals(is)),'_noG.mat');
    save(filename, 'adistrA', 'adistrB', 'share2A', 'share2B', ...
        "Ga","Gb", "Va", "Vb", "VOTESa", "VOTESb", "wage", "r", "kaggA", ...
        "shares")

end