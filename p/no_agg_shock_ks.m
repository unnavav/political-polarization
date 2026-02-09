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


%% pulling in steady states of interest and relevant policies

cd ../d/steadystates/
load resultsna100ah50rho90sig3.mat

clearvars -except Warray Karray etagrid taugrid captax

g1 = Warray{1,5}; g2 = Warray{3,3};
K_ss_pop = Karray{1,5}; K_ss_lib = Karray{3,3};

etap = etagrid(1); etal = etagrid(3);
taup = taugrid(5); taul = taugrid(3);

taugrid = [taup taul]; etagrid = [etap etal];

cd ../../p

%% set up params to load into the KS algo

vTol = 1e-6;

% model params
alpha = 0.36; delta = 0.06; beta = 0.96; sigma = 3; phi = 0;

% grid sizes
nl = 7; na = 100; nmu = na*10; nr = 2; nk = 13;

al = 0; ah = 50;
% get labor distribution and aggregate values
%step 1: make labor grid and labor transition matrix

mu = 0;
rho = .9;
sig_l = 0.2;
range = 2.575;

% need to back out sigma^2_e given sigma^2_l
sigx = sig_l*sqrt(1 - rho^2);

[pil, lgrid] = compute.getTauchen(nl,  mu, sigx, rho, range);

agrid = compute.logspace(al, ah, na)';
nmu = na*10;
amu = linspace(agrid(1), agrid(na), nmu);

verbose = true;

% set size of aggregate moment grid
kl = .85;
kh = 1.15;
Kgrid = linspace(K_ss_pop*kl, K_ss_lib*kh, nk); 

% we also need the scaling grid for lambda to go in:
stationary_pil = asymptotics(dtmc(pil));      % nl x 1
Eeps     = stationary_pil*lgrid';              % sum ε Ω
Eeps_sum = sum(Eeps);

Eeps_pow = sum((lgrid'.^(1 - taugrid)) .* stationary_pil'); % sum ε^(1-τ) Ω  (τ regime-specific)

lambda_ratio = Eeps_sum ./ Eeps_pow;   

% start with a forecast that's basically steady state:

Kfore = [log(K_ss_pop)*0.01, 0.99;   % Small persistence, close to SS
         log(K_ss_lib)*0.01, 0.99];

Rfore = [-.5 0; ...
    -.5 0]; %start them at some random Rfore

Rswitch = [.95 .05; ...
    .05 .95];
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
    'Rswitch', Rswitch, ... 
    'captax', captax, ...
    'taugrid', taugrid, ...
    'G', 0, ...
    'Kgrid', Kgrid, ...
    'etagrid', etagrid, ...
    'lamval', lambda_ratio, ...
    'rnseed', 1234567);

rng("default")
T = 1000;
% Rt = predict.sim(T,2,"default",Rguess);

verbose = true;
forearray = cell(50,1);

forearray{1} = Kfore;
array_ind = 1;
b1 = Kfore(1,:); b2 = Kfore(2,:);

foredist = 10;
kforedist = 10;
rforedist = 10;

% setting up predicitions for Regime change
% start from assuming capital has no effect, then start to build out
% forecast
mymodelfun = @(beta,x) 1./(1 + exp(-(beta(1) + beta(2).*x)));
beta0_1 = [2.9957; 0];
beta0_2 = [0.0513; 0];

%% begin iteration

% closing in on this bitch
iter_ct = 1;
Varray = cell(1);
EVarray = cell(1);
Garray = cell(1);
Kforearray = cell(1);
Rforearray = cell(1);

terms.starter_distr = g1;

while foredist > vTol

    Kforearray{iter_ct} = Kfore;
    Rforearray{iter_ct} = Rfore;

    if foredist > 1e-1
        vTol = 1e-4;
    else
        vTol = 1e-6;
    end

    fprintf("\n ===== Generating forecast rules ===== \n")


    fprintf("\nGetting Regression Data\n")
    [Kprdat Rdata Prdata distr_array V G EV] = ...
        ks.getRegData(T, terms, vTol, verbose);

    fprintf('\n\nRegime 1: %d periods, Regime 2: %d periods\n', ...
        sum(Rdata==1), sum(Rdata==2));

    Varray{iter_ct} = V;
    EVarray{iter_ct} = EV;
    Garray{iter_ct} = G;
        
    K_next = log(Kprdat(2:end));
    K_curr = log(Kprdat(1:end-1));
    R_curr = Rdata(1:end-2);
    R_next = Rdata(2:end-1);

    %           updating K coefficients
    % Capital law: log K_{t+1} = a(R_t) + b(R_t) log K_t
    if sum(R_curr == 1) > 10
        kmdl1 = fitlm(K_curr(R_curr==1), K_next(R_curr==1));
        b1 = kmdl1.Coefficients.Estimate;
    end
    
    if sum(R_curr == 2) > 10
        kmdl2 = fitlm(K_curr(R_curr==2), K_next(R_curr==2));
        b2 = kmdl2.Coefficients.Estimate;
    end
    Kfore_new = [b1'; b2'];
    kforedist = norm(Kfore_new-Kfore, 'inf');

    % Update K forecast
    Kfore = 0.8*Kfore + 0.2*Kfore_new;
    terms.Kfore = Kfore;

    %              Updating R Coefficients
    % R law = Pr(R' = 1 | R) = exp(1/(d(R) + e(R)lnK))

    % Binary indicator for next regime being 1
    Y = (R_next == 1);  
    if sum(R_curr == 1) > 10
        X1   = K_curr(R_curr == 1);
        X1_std = (X1 - mean(X1)) / std(X1);
        Y1   = Y(R_curr == 1);       % 0/1
        opts = statset('Display','iter');
        nlm1 = fitnlm(X1_std, Y1, mymodelfun, beta0_1);
        br1  = nlm1.Coefficients.Estimate;   % [beta0; beta1]
    end

%     scatter(X1 , Y1)
%     hold on
%     fplot(@(kk) mymodelfun(nlm1.Coefficients.Estimate, kk), [min(X1), max(X1)])
%     hold off

    if sum(R_curr == 2) > 10
        X2   = K_curr(R_curr == 2);
        X2_std = (X2 - mean(X2)) / std(X2);
        Y2   = Y(R_curr == 2);
        nlm2 = fitnlm(X2_std, Y2, mymodelfun, beta0_2);
        br2  = nlm2.Coefficients.Estimate;
    end
    Rfore_new = [br1'; br2'];
    % Update K forecast
    testK = linspace(min(log(K_curr)), max(log(K_curr)), nk);
    p_old = ks.forecastR(Rfore,testK);
    p_new = ks.forecastR(Rfore_new,testK);
    rforedist = norm(p_new - p_old, 'inf');

    Rfore = 0.8*Rfore + 0.2*Rfore_new;
    terms.Rfore = Rfore;

    fprintf('K(1) = %0.4f\n', Kprdat(1));
    fprintf('K range: [%0.4f, %0.4f]\n', min(Kprdat), max(Kprdat));

    figure;

    % 1) Two-panel plot
    subplot(2,1,1);
    plot(Prdata,'LineWidth',1.6,'Color',[0.80 0.20 0.20]);
    yline(0.5,'--','Color',[0.4 0.4 0.4],'LineWidth',1);  % majority threshold
    xlabel('Time');
    ylabel('P(populism)');
    title('Probability of Voting for Populism');
    ylim([0 1]);
    grid on;
    set(gca,'FontSize',12);

    subplot(2,1,2);
    plot(Kprdat,'LineWidth',1.6,'Color',[0.20 0.30 0.75]);
    xlabel('Time');
    ylabel('Aggregate Capital K');
    title('Aggregate Capital Path');
    grid on;
    set(gca,'FontSize',12);

    fprintf('Regime 1: log K'' = %0.4f + %0.4f log K\n', b1(1), b1(2));
    fprintf('Regime 1: Pr R pr = 1 = exp((%0.4f + %0.4f log K)^-1)\n', br1(1), br1(2));

    fprintf('Regime 2: log K'' = %0.4f + %0.4f log K\n', b2(1), b2(2));
    fprintf('Regime 2: Pr R pr = 1 = exp((%0.4f + %0.4f log K)^-1)\n', br2(1), br2(2));

    fprintf("\nCapital: R2 for pop = %0.6f\nR2 for lib = %0.6f", ...
        kmdl1.Rsquared.Adjusted, ...
        kmdl2.Rsquared.Adjusted)

    fprintf("\nRegime Guessing Distance: %1.4f", ...
        rforedist)
    foredist = max(rforedist, kforedist);
    fprintf("\nForedist = %0.6f\n\n", foredist)

    iter_ct = iter_ct + 1;

    % After running fitnlm and having:
    % br1, br2  % 2x1 coefficient vectors [beta0; beta1]
    % K_curr, R_curr, R_next
    
    % % Build function handle for logit
    % mymodelfun = @(beta,x) 1./(1 + exp(-(beta(1) + beta(2).*x)));
    % 
    % % Use log K if that’s what you estimated on
    % X1 = log(K_curr(R_curr == 1));
    % X2 = log(K_curr(R_curr == 2));
    % 
    % p1 = mymodelfun(br1, X1);   % P(R_{t+1}=1 | R_t=1, K_t)
    % p2 = mymodelfun(br2, X2);   % P(R_{t+1}=1 | R_t=2, K_t)
    % 
    % figure;
    % scatter(X1, p1, 10, 'b', 'filled'); hold on;
    % scatter(X2, p2, 10, 'r', 'filled');
    % 
    % ylim([0 1]);
    % xlabel('log K_t');
    % ylabel('P(R_{t+1} = 1)');
    % legend('Current R_t = 1','Current R_t = 2','Location','best');
    % grid on;

end

save ../d/ks_rfore_endo.mat