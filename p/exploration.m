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

vTol = 1e-5; gTol = 1e-8; dTol = 1e-3;
%% params
alpha = 0.36; delta = 0.08; beta = 0.96173; sigma = 1; phi = 0; gamma = .5;
 
nl = 7;
na = 75;
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
klmult = 0;
khmult = .;

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
lamval = .01; % people are really reactive to taxes
adj = .5;

%party tax regimes
tgrid = [.45 .5]; 

p = [.3 .15];

ubonus = .005;
pctA = .5;

%% solving for wages and r by getting aggregate capital

adistr = zeros(nl, nmu, np);

for i = 1:nmu
    for il = 1:nl
        for ip = 1:np
            adistr(il,i, ip) = 1.0/(nmu*nl*np);
        end
    end
end

kDist = 10;
gDist = 10; 
f=10; 

kval = (kl+ kh)/2;

DIST = max(kDist,gDist);

VOTES = zeros(na, nl);

while kDist > dTol

    fprintf("\n*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*")
    fprintf("\nA guess: %4.4f. Begin iteration for solution...\n", kval)
    fprintf("\t Solving value function:\n")

    p = [.3 .7];

    r = alpha*(kval^(alpha - 1)*(lagg^(1-alpha))) - delta;
    wage = (1-alpha)*((kval^(alpha))*(lagg^(-alpha)));

    % making tax schedule
    
    wage_inc = repmat(wage*lgrid,na,1)';
    cap_inc = repmat(r*(agrid),nl,1);
    Y = wage_inc+cap_inc;
    T = zeros(nl, na, np);
    for i = 1:np
        T(:,:,i) = gov.tax(Y,lamval,tgrid(i));
    end

    wage_inc_mu = repmat(wage*lgrid, nmu,1)';
    cap_inc_mu = repmat(r*amu, nl, 1);
    Ymu = wage_inc_mu + cap_inc_mu;
    Tmu = zeros(nl, nmu, np);
    for ip = 1:np
        Tmu(:,:,ip) = gov.tax(Ymu,lamval,tgrid(ip));
    end

    G = [0 0];
    for ip = 1:np
        G(ip) = sum(sum(Tmu(:,:,ip).*adistr(:,:,ip)));
        Y(:,:,ip) = T(:,:,ip)+G(ip);
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
        'Y', Y, ...
        'G', G, ...
        'p', p, ...
        'ubonus',ubonus);

    [Va, Vb, EVa, EVb, ga, gb] = HH.solve(nl, na, np, terms, vTol, gTol);

    VOTESa = (EVa(:,:,2) > EVa(:,:,1));
    VOTESb = (EVb(:,:,2) > EVb(:,:,1));

    test = EVa(:,:,1) - EVb(:,:,1);

    tiledlayout(5,1)
    nexttile
    mesh(ga(:,:,1)-ga(:,:,2));
    nexttile
    mesh(gb(:,:,1)-gb(:,:,2));
    nexttile
    mesh(ga(:,:,1) - gb(:,:,1));
    nexttile
    mesh(VOTESa - VOTESb);
    nexttile
    mesh(EVb(:,:,2) - EVb(:,:,1))

        % asset distr
    fprintf("\tSolving asset distribution:\n")
    [adistr, kagg] = HH.getDist(ga, gb, amu, agrid, pil, pctA);

    % CONDENSE DISTR
    acond = compute.condense(adistr, amu, agrid);

    for ip = 1:np
        share2 = pctA*sum(sum(VOTESa.*acond(:,:,ip))) + ...
            (1-pctA)*sum(sum(VOTESb.*acond(:,:,ip)));
        p(ip) = exp(share2/gamma)/(exp(share2/gamma) + exp((1-share2)/gamma));
    end

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
    
sound(y, Fs)

%%
date = string(datetime("today"));
filename = strcat("..\d\results_old_params",date);
save(filename)

%% graphing
[xgrid, ygrid] = meshgrid(lgrid,agrid);



set(gca, 'XTick', [agrid(1) agrid(30:15:75)], 'XTickLabels', xlab)

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
exportgraphics(gcf,"..\v\hha_decision.png", 'Resolution',300)

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

%% now we do the stochastic responses 

zt = [repelem(3, 100) 5 repelem(3, 1000)];
params = [alpha beta delta sigma lagg];
k0 = kval;
[zt yt ct kt] = predict.perfectForesight(zt, kval, k0, .95, params, dTol);