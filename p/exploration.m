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

load handel

vTol = 1e-5; gTol = 1e-8; dTol = 1e-3;
%% params
alpha = 0.36; delta = 0.08; beta = 0.96173; sigma = 1; phi = 0; gamma = .5;
 
nl = 7;
na = 250;
nmu = na*10;
np = 2;

al = 0+phi; ah = 100+phi;

wage = 1.1; r = 0.04;

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
adj = .5;

%party tax regimes
% tgrid = [1-.181 1-.086]; 
tgrid = [.086 .181]; 

p = [1 0];

ubonus = .005;

ubonusA = [ubonus 0];
ubonusB = [0 ubonus];

pctA = .5;

captax = .15; %from US tax code
goal = .3652;

G = [0 0];

%% solving for wages and r by getting aggregate capital

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

kDist = 10;
gDist = 10; 
f=10; 

kval = (kl+ kh)/2;

DIST = max(kDist,gDist);

VOTES = zeros(na, nl);
    
while DIST > dTol

    fprintf("\n\n\n--------------------------------- > Lambda = %4.4f", lamval)

    while kDist > dTol
    
        fprintf("\n*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*")
        fprintf("\nA guess: %4.4f. Begin iteration for solution...\n", kval)
        fprintf("\t Solving value function:\n")
        
        % first getting equilibrium objects like prices, govt spending, etc
        r = alpha*(kval^(alpha - 1)*(lagg^(1-alpha))) - delta;
        wage = (1-alpha)*((kval^(alpha))*(lagg^(-alpha)));

        wage_inc_mu = repmat(wage*lgrid, nmu,1)';
        cap_inc_mu = repmat(r*amu, nl, 1);
        Tmu = zeros(nl, nmu, np);
        for ip = 1:np
            Tmu(:,:,ip) = gov.tax(wage_inc_mu,lamval,tgrid(ip)) + ...
                cap_inc_mu.*(1-captax);
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
            'tau', tgrid, ...
            'G', G, ...
            'p', p, ...
            'ubonus',ubonus, ...
            'captax', captax);
    
        fprintf("\t Solving value function for HH a:\n")
        Va = zeros(nl, na, np);
        Ga = Va;
        EVa = Va;
        for ip = 1:np
            fprintf("\t\t Party %i:\n", ip)

            terms.tau = tgrid(ip);
            terms.G = G(ip);
            terms.ubonus = ubonusA(ip);

            [Vaip, Gaip, EVaip] = HH.solve(nl, na, np, terms, vTol, gTol);
            Va(:,:,ip) = Vaip;
            Ga(:,:,ip) = Gaip;
            EVa(:,:,ip) = EVaip;
        end

        fprintf("\t Solving value function for HH b:\n")
        Vb = zeros(nl, na, np);
        Gb = Vb;
        EVb = Vb;
        for ip = 1:np
            fprintf("\t\t Party %i:\n", ip)

            terms.tau = tgrid(ip);
            terms.G = G(ip);
            terms.ubonus = ubonusB(ip);

            [Vbip, Gbip, EVbip] = HH.solve(nl, na, np, terms, vTol, gTol);
            Vb(:,:,ip) = Vbip;
            Gb(:,:,ip) = Gbip;
            EVb(:,:,ip) = EVbip;
        end
        
        VOTESa = (EVa(:,:,1) > EVa(:,:,2));
        VOTESb = (EVb(:,:,1) > EVb(:,:,2));
    
        test = EVa(:,:,1) - EVb(:,:,1);
    
        %TO DO: 
        % SOLVE ASSET DISTRIBUTION 
        % SOLVE VOTE COUNTS

        fprintf("\tSolving asset distribution:\n")
        [adistrA, kaggA, adistrB, kaggB] = HH.getDist(Ga, Gb, amu, agrid, pil, pctA);
    
%         % CONDENSE DISTR
%         acondA = compute.condense(adistrA, amu, agrid);
%         acondB = compute.condense(adistrB, amu, agrid);
%     
%         share2A = pctA*sum(sum(VOTESa.*acondA)) + ...
%             (1-pctA)*sum(sum(VOTESb.*acondA));
%         share2B = pctA*sum(sum(VOTESa.*acondB)) + ...
%             (1-pctA)*sum(sum(VOTESb.*acondB));
%         p(1) = share2A;
%         p(2) = share2B;
%     
        % 
        f = kaggA - kval;
    %     f = kaggB - kval;
    
        if f > 0
            fprintf("\n||Kguess - Kagg|| = %4.4f. \tAggregate capital is too low.\n", abs(f))
            kl = .5*(kval+kl);
        else
            fprintf("\n||Kguess - Kagg|| = %4.4f. \tAggregate capital is too high.\n", abs(f))
            kh = .5*(kval+kh);
        end
    
        kDist = abs(kaggA - kval);  % check whether the capital diff
                                    % is changing at all
    
        kval = .5*(kl + kh);
    end
    
    tA = adistrA.*Tmu(:,:,1); %getting all taxes collected
    tA = sum(tA);
    tB = adistrB.*Tmu(:,:,2); %getting all taxes collected
    tB = sum(tB);
    Y = kval^(alpha)*lagg^(1-alpha);
    gy = tA/Y;

    taxdist = sum(tA-goal);

    if taxdist>goal
       fprintf("\nGov't rev collected = %4.4f. Lam = %4.4f. " + ...
           "\tTax rate is too high.\n\n", taxdist, lamval)
       lh = (lamval*adj+(1-adj)*lh);
    else
       fprintf("\nGov't rev collected = %4.4f. Lam = %4.4f. " + ...
           "\tTax rate is too low.\n\n", taxdist, lamval);
       ll = (lamval*adj+(1-adj)*ll);
    end

    gDist = abs(taxdist);
    lamval = .5*(ll + lh);

    DIST = max(kDist, gDist);
    fprintf("||DIST|| = %4.4f\n", DIST)
    kl = klwrbnd*klmult;
    kh = klwrbnd*khmult;
    kDist = 10;
end

sound(y, Fs)

%%
date = string(datetime("today"));
filename = strcat("..\d\results_wage_tax_EGM",date);
save(filename)

%% graphing
[xgrid, ygrid] = meshgrid(lgrid,agrid);

lrep = repmat(lgrid,1,na)';
arep = repmat(agrid,nl,1);
arep = arep(:);
dat = ga(:,:,1);
dat = dat(:);
G1dat = [lrep arep dat];

tiledlayout(1,1)
graphdata = griddata(G1dat(:,1),G1dat(:,2), G1dat(:,3), xgrid, ygrid);
ax1 = nexttile;
mesh(graphdata);
title('Household A Savings Choices', 'FontSize', 16)
yticks(1:15:75)
yticklabels([agrid(1:15:75)])
ylabel('Assets', 'FontSize', 14)
xticks(1:2:7)
xticklabels(lgrid(1:2:7))
xlabel('Labor Productivity','FontSize', 14)
colormap(ax1,winter)
% exportgraphics(gcf,"..\v\hha_decision.png", 'Resolution',300)

dat = gb(:,:,1);
dat = dat(:);
G1dat = [lrep arep dat];
tiledlayout(1,1)
graphdata = griddata(G1dat(:,1),G1dat(:,2), G1dat(:,3), xgrid, ygrid);
ax1 = nexttile;
mesh(graphdata);
title('Household B Savings Choices', 'FontSize', 16)
yticks(1:15:75)
yticklabels([agrid(1:15:75)])
ylabel('Assets', 'FontSize', 14)
xticks(1:2:7)
xticklabels(lgrid(1:2:7))
xlabel('Labor Productivity','FontSize', 14)
colormap(ax1,winter)

dat = VOTESa;
dat = dat(:);
G1dat = [lrep arep dat];
tiledlayout(1,1)
graphdata = griddata(G1dat(:,1),G1dat(:,2), G1dat(:,3), xgrid, ygrid);
ax1 = nexttile;
contourf(graphdata);
title('Household A Voting Choices', 'FontSize', 16)
yticks(1:15:75)
yticklabels([agrid(1:15:75)])
ylabel('Assets', 'FontSize', 14)
xticks(1:2:7)
xticklabels(lgrid(1:2:7))
xlabel('Labor Productivity','FontSize', 14)
colormap(ax1,winter)

dat = VOTESb;
dat(1,1) = 1e-2;
dat = dat(:);
G1dat = [lrep arep dat];
tiledlayout(1,1)
graphdata = griddata(G1dat(:,1),G1dat(:,2), G1dat(:,3), xgrid, ygrid);
ax1 = nexttile;
contourf(graphdata);
title('Household B Voting Choices', 'FontSize', 16)
yticks(1:15:75)
yticklabels([agrid(1:15:75)])
ylabel('Assets', 'FontSize', 14)
xticks(1:2:7)
xticklabels(lgrid(1:2:7))
xlabel('Labor Productivity','FontSize', 14)
colormap(ax1,winter)