%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementing Krussel Smith
% vaasavi
% oct 2025
% 
% this code aims to implement krussel smith by several aggregate states:
% capital depreciation (delta), regime (R), aggregate capital (K).
% algorithm detailed in paper, etc. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

restoredefaultpath;
clear all; clc;
addpath(genpath(pwd));

%% set up params to load into the KS algo

vTol = 1e-6;

% model params
alpha = 0.36; delta = 0.06; beta = 0.96; sigma = 3; phi = 0;

% grid sizes
nl = 7; na = 100; nmu = na*10; nr = 2; nk = 13;

alpha = 0.36; delta = 0.06; beta = 0.96; sigma = 3; phi = 0;

al = 0; ah = 200;
% get labor distribution and aggregate values
%step 1: make labor grid and labor transition matrix

mu = 0;
rho = .9;
sig_l = 0.5327;
range = 2.575;

% need to back out sigma^2_e given sigma^2_l
sigx = sig_l*sqrt(1 - rho^2);

[pil, lgrid] = compute.getTauchen(nl,  mu, sigx, rho, range);

% also need to create a grid for capital depreciation shock
mu_d = 0.025; rho_d = .9; sig_mu = .01;
sigmx = sig_mu*sqrt(1 - rho_d^2);
range = 3.5;

% [pid, dgrid] = compute.getTauchen(2, mu_d, sigmx, rho_d, range);

pid = [.95 0.05; 0.05 .95];
dgrid = [0.02 -0.02];

agrid = compute.logspace(al, ah, na)';
nmu = na*10;
amu = linspace(agrid(1), agrid(na), nmu);

verbose = true;

% set size of aggregate moment grid
Kgrid = linspace(7, 45, nk);

% definitions for the probability of changing regime
probscale1 = 10;
probscale2 = 10;

phis = [probscale1 probscale2];

% finally set policies
taul = 0.086; taup = 0.181; etal = .5; etap = 0;

taugrid = [taup taul]; etagrid = [etap etal];

captax = [.15 .15 .15 .15 .20 .20 .20];

% we also need the scaling grid for lambda to go in:
stationary_pil = asymptotics(dtmc(pil));      % nl x 1
Eeps     = stationary_pil*lgrid';              % sum ε Ω
Eeps_sum = sum(Eeps);

Eeps_pow = sum((lgrid'.^(1 - taugrid)) .* stationary_pil'); % sum ε^(1-τ) Ω  (τ regime-specific)

lambda_ratio = Eeps_sum ./ Eeps_pow;   


% start with a forecast that's basically steady state:
Kfore = [.1 .99; ...
          .08 .99];

Rfore = [2.2 8; ...
          -2.2 8];% trying a different form 
%% prep VFI

terms = struct('alpha', alpha, ...
    'beta', beta, ...
    'delta', delta, ...
    'sigma', sigma, ...
    'phi', phi, ...
    'agrid', agrid, ...
    'lgrid', lgrid, ...
    'pil', pil, ...
    'Kfore', Kfore, ...
    'Rfore', Rfore, ...
    'dgrid', dgrid, ...
    'captax', captax, ...
    'taugrid', taugrid, ...
    'G', 0, ...
    'Kgrid', Kgrid, ...
    'etagrid', etagrid, ...
    'lamval', lambda_ratio, ...
    'phis', phis, ...
    'rnseed', 1234567);

rng("default")
T = 1000;
% Rt = predict.sim(T,2,"default",Rguess);

verbose = true;
forearray = cell(50,1);

newfore = Kfore;
forearray{1} = Kfore;
array_ind = 1;
b1 = Kfore(1,:); b2 = Kfore(2,:);

foredist = 10;
kforedist = 10;
modelfun = @(b, x)(1+exp(b(1)+b(2)*x))^-1;
beta0 = [.5 5];
%% begin iteration
Rt = rand(T,1);

% closing in on this bitch
iter_ct = 1;
Varray = cell(1);
EVarray = cell(1);
Garray = cell(1);
Kforearray = cell(1);
Rforearray = cell(1);

while foredist > vTol

    Rforearray{iter_ct} = Rfore;
    Kforearray{iter_ct} = Kfore;

    if foredist > 1e-1
        vTol = 1e-4;
    else
        vTol = 1e-6;
    end

    fprintf("\n ===== Generating K forecast Rule ===== \n")
    while kforedist > vTol

        fprintf("\nGetting Regression Data\n")
        [Kprdat Rdata Prdata distr_array V G EV] = ...
            ks.getRegData(T, terms, vTol, verbose);

        Varray{iter_ct} = V;
        EVarray{iter_ct} = EV;
        Garray{iter_ct} = G;
        
        K1 = Kprdat(Rdata == 1);
        K2 = Kprdat(Rdata == 2);
        Pr1 = Prdata(Rdata == 1);
        Pr2 = Prdata(Rdata == 2);
        
        K_next = log(Kprdat(2:end));
        K_curr = log(Kprdat(1:end-1));
        R_curr = Rdata(1:end-1);
        
        X = [ones(size(K_curr)), K_curr];
        
        % Capital law: log K_{t+1} = a(R_t) + b(R_t) log K_t
        if sum(R_curr == 1) > 0
            kmdl1 = fitlm(K_curr(R_curr==1,:), K_next(R_curr==1));
            b1 = kmdl1.Coefficients.Estimate;
        end
        
        if sum(R_curr == 2) > 10
            kmdl2 = fitlm(K_curr(R_curr==2), K_next(R_curr==2));
            b2 = kmdl1.Coefficients.Estimate;
        end
        Kfore_new = [b1'; b2'];

        % adjusting forecast slowly
        Kfore = 0.5*Kfore + 0.5*Kfore_new;
        
        fprintf("\nCapital: R2 for pop = %0.6f\nR2 for lib = %0.6f", ...
            kmdl1.Rsquared.Adjusted, ...
            kmdl2.Rsquared.Adjusted)

        kforedist = norm(Kfore_new-Kfore, 'inf');
        fprintf("\nForedist = %0.6f\n\n", kforedist)

    end

    Kpr = ks.forecastK(Kfore, Kgrid); % nk x r
    [EV1, EV2, Votes_EV] = gov.getVotingExpectations(V, pil, Kpr, Kgrid);

    Rt = zeros(T,1);
    pr_logit = Rt;
    Rt(1) = Rdata(1);
    for t = 2:T
        Kt = Kprdata(t-1);
        R_t = Rt(t-1);
        [ix we] = compute.weight(Kgrid, Kt);

        g_today = distr_array{t};
        acond = compute.condense(g_today, amu, agrid);
    
        % weighting over whichever Kpr state
        Votes = we*Votes_EV(ix,R_t, :, :) + ...
            (1-we)*Votes_EV(ix+1,R_t, :, :);
        Votes = squeeze(Votes);
    
        pr = sum(sum(Votes.*acond));
        % now feed it through logit
        pr_logit(t) = 1./(1+exp(-1*(pr-.5)));
        Rt(t) = 2 - (pr_logit(t)>0.5);
    end
%     [w1 kagg1] = HH.getDist(squeeze(G(1,1,:,:)), amu, agrid, pil, verbose);
%     [w2 kagg2] = HH.getDist(squeeze(G(25,1,:,:)), amu, agrid, pil, verbose);
% 
%     [wl kaggl] = HH.getDist(squeeze(G(1,2,:,:)), amu, agrid, pil, verbose);
%     [wp kaggp] = HH.getDist(squeeze(G(1,1,:,:)), amu, agrid, pil, verbose);


    % Regime law: Pr(R_{t+1}=1 | R_t, K_t) = logit(d(R_t)+e(R_t) log K_t)
    % Do two logits conditioned on current regime:
    mdl1 = fitglm(K_curr(Rt(1:999)==1), pr_logit(Rt(1:999)==1), ...
                  'Distribution','binomial','Link','logit','Intercept',true);
    mdl2 = fitglm(K_curr(Rt(1:999)==2), pr_logit(Rt(1:999)==2), ...
                  'Distribution','binomial','Link','logit','Intercept',true);
    Rfore_new = [mdl1.Coefficients.Estimate'; mdl2.Coefficients.Estimate'];  % [d1 e1; d2 e2]
    

    % Damping
    Rfore = 0.5*Rfore + 0.5*Rfore_new;
    
    % Forecaster distance
    foredist = max( [norm(Kfore_new-Kfore, 'inf'), norm(Rfore_new-Rfore, 'inf')] );

    fprintf("\nCapital: R2 for pop = %0.6f\nR2 for lib = %0.6f", ...
        kmdl1.Rsquared.Adjusted, ...
        kmdl2.Rsquared.Adjusted)
    fprintf("\nRegime: R-sq for pop = %0.6f\nR2 for lib = %0.6f", mdl1.Rsquared.ordinary, ...
        mdl2.Rsquared.ordinary)
    fprintf("\nForedist = %0.6f\n\n", foredist)
%     % only once we start getting to a place where beliefs show up, I guess?
%     if ~(isempty(K1) || isempty(K2))
%         break
%         logit1 = fitnlm(K1, Pr1, beta0, modelfun);
%     else
%         fprintf("\n Kgrid not correct maybe?\n")
%         K_min = min(Kprdat);
%         K_max = max(Kprdat);
%         buffer = 0.15;  % 15% slack
%         Kgrid = linspace((1-buffer)*K_min, (1+buffer)*K_max, nk);
%         terms.Kgrid = Kgrid;
%         disp(Kgrid);
%     end

    iter_ct = iter_ct + 1;
end

save ../d/ks_rmseed_in_pr