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
load 20260417_etatest_results_rho90sig2.mat

clearvars -except Warray Karray etagrid taugrid captax dgrid pid parray

g1 = Warray{1,3}; g2 = Warray{4,2};
K_ss_pop = Karray{1,3}; K_ss_lib = Karray{4,2};

etap = etagrid(1); etal = etagrid(4);
taup = taugrid(3); taul = taugrid(2);

taugrid = [taup taul]; etagrid = [etap etal];

cd ../../p

%% set up params to load into the KS algo

vTol = 1e-6;
T = 4000;

% model params
alpha = 0.36; delta = 0.06; beta = 0.96; sigma = 3; phi = 0;

% grid sizes
nl = 7; na = 100; nmu = na*10; nr = 2; nk = 13; nd = length(dgrid);

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

%you're about to see something horrible. Don't read it. I'm tired
Kfored = zeros(nd, nr, 2); Rfored = Kfored;
Kfored(1,:,:) = Kfore; Kfored(2,:,:) = Kfore;
Rfored(1,:,:) = Rfore; Rfored(2,:,:) = Rfore;
Kfore = Kfored; Rfore = Rfored;

Kmeans = zeros(nr, nd); Kstds = ones(nr, nd);

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
    'pid', pid, ...
    'captax', captax, ...
    'taugrid', taugrid, ...
    'G', 0, ...
    'Kgrid', Kgrid, ...
    'etagrid', etagrid, ...
    'lamval', lambda_ratio, ...
    'rnseed', 1234567, ...
    'Kmeans', Kmeans, ...
    'Kstds', Kstds);

rng("default")

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

br11_prev = beta0_1; br12_prev = beta0_1;
br21_prev = beta0_2; br22_prev = beta0_2;

%% begin iteration

% closing in on this bitch
iter_ct = 1;
Varray = cell(1);
EVarray = cell(1);
Garray = cell(1);
Kforearray = cell(1);
Kmeansarray = cell(1);
Kstdsarray = cell(1);
Rforearray = cell(1);
nlms = cell(1,2,2);



terms.starter_distr = g2;

% Get today's date as a datetime object
t = datetime('today');

% Convert the datetime object to a string with the specified format
todayDateStr = string(t, 'yyyyMMdd');

filename = strcat('ks_rfore_endo_all_',todayDateStr ,'2.mat');

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
    [Kprdat Rdata Prdata ddata distr_array V G EV] = ...
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
    d_curr = ddata(1:end-1);

    % indices for four (R_t, d_t) states
    ix11 = (R_curr == 1 & d_curr == 1);
    ix12 = (R_curr == 1 & d_curr == 2);
    ix21 = (R_curr == 2 & d_curr == 1);
    ix22 = (R_curr == 2 & d_curr == 2);

    %           updating K coefficients
    % Capital law: log K_{t+1} = a(R_t) + b(R_t) log K_t
    bK11 = NaN(2,1);  % [a11; b11]
    bK12 = NaN(2,1);
    bK21 = NaN(2,1);
    bK22 = NaN(2,1);
    
    if sum(ix11) > 10
        mdl11   = fitlm(K_curr(ix11), K_next(ix11));
        bK11    = mdl11.Coefficients.Estimate;   % [const; slope]
    end
    if sum(ix12) > 10
        mdl12   = fitlm(K_curr(ix12), K_next(ix12));
        bK12    = mdl12.Coefficients.Estimate;
    end
    if sum(ix21) > 10
        mdl21   = fitlm(K_curr(ix21), K_next(ix21));
        bK21    = mdl21.Coefficients.Estimate;
    end
    if sum(ix22) > 10
        mdl22   = fitlm(K_curr(ix22), K_next(ix22));
        bK22    = mdl22.Coefficients.Estimate;
    end
    
    % stack as 4×2: rows = (R,d) pairs in fixed order
    % e.g. row1:(1,1), row2:(1,2), row3:(2,1), row4:(2,2)
    Kfore1 = [bK11'; bK21'];
    Kfore2 = [bK12'; bK22'];
    Kfore_new(1,:,:) = Kfore1; Kfore_new(2,:,:) = Kfore2;

    kforedist = compute.dist(Kfore_new,Kfore, 3);
    Kfore     = 0.8*Kfore + 0.2*Kfore_new;
    terms.Kfore = Kfore;
    %              Updating R Coefficients
    % R law = Pr(R' = 1 | R) = exp(1/(d(R) + e(R)lnK))

    % Binary indicator for next regime being 1
    Y = (R_next == 1);  % 0/1
    
    br11 = NaN(2,1); br12 = NaN(2,1);
    br21 = NaN(2,1); br22 = NaN(2,1);
    
    opts = statset('Display','off');  % silence if you want
    
    % I have to scale the logK because otherwise the difference is
    % not big enough to get regression coefficients that work. 
    % (R_t,d_t) = (1,1)
    if sum(ix11) > 10
        Kmean = mean(K_curr(ix11));
        Kstd = std(K_curr(ix11));
        X1 = (K_curr(ix11) - Kmean) / Kstd;

        Kmeans(1,1) = Kmean; Kstds(1,1) = Kstd;

        Y1    = Y(ix11);
        nlm11 = fitnlm(X1, Y1, mymodelfun, br11_prev, 'Options', opts);
        br11  = nlm11.Coefficients.Estimate;   % [beta0; beta1]
        br11_prev = br11;
    end
    %preventing degenerate coeff when there's a low number of obs
    if any(abs(br11) > 20), br11 = br11_prev; end 

    
    % (1,2)
    if sum(ix12) > 10
        Kmean = mean(K_curr(ix12));
        Kstd = std(K_curr(ix12));
        X2 = (K_curr(ix12) - Kmean) / Kstd;

        Kmeans(1,2) = Kmean; Kstds(1,2) = Kstd;
        
        Y2    = Y(ix12);
        nlm12 = fitnlm(X2, Y2, mymodelfun, br12_prev, 'Options', opts);
        br12  = nlm12.Coefficients.Estimate;
        br12_prev = br12;
    end
    if any(abs(br12) > 20), br12 = br12_prev; end

    % (2,1)
    if sum(ix21) > 10
        Kmean = mean(K_curr(ix21));
        Kstd = std(K_curr(ix21));
        X3 = (K_curr(ix21) - Kmean) / Kstd;

        Kmeans(2,1) = Kmean; Kstds(2,1) = Kstd;
        Y3    = Y(ix21);
        nlm21 = fitnlm(X3, Y3, mymodelfun, br21_prev, 'Options', opts);
        br21  = nlm21.Coefficients.Estimate;
        br21_prev = br21;
    end
    if any(abs(br21) > 20), br21 = br21_prev; end

    % (2,2)
    if sum(ix22) > 10
        Kmean = mean(K_curr(ix22));
        Kstd = std(K_curr(ix22));
        X4 = (K_curr(ix22) - Kmean) / Kstd;

        Kmeans(2,2) = Kmean; Kstds(2,2) = Kstd;
        Y4    = Y(ix22);
        nlm22 = fitnlm(X4, Y4, mymodelfun, br22_prev, 'Options', opts);
        br22  = nlm22.Coefficients.Estimate;
        br22_prev = br22;
    end
    if any(abs(br22) > 20), br22 = br22_prev; end

    terms.Kmeans = Kmeans;
    terms.Kstds = Kstds;
    Kmeansarray{iter_ct} = Kmeans; Kstdsarray{iter_ct} = Kstds;


    Rfore1 = [br11'; br21'];
    Rfore2 = [br12'; br22'];
    p_old = ks.forecastR(Rfore,Kgrid,Kmeans, Kstds);
    Rfore_new(1,:,:) = Rfore1; Rfore_new(2,:,:) = Rfore2;
    Rfore = 0.8*Rfore + 0.2*Rfore_new;
    p_new = ks.forecastR(Rfore_new,Kgrid, Kmeans, Kstds);

    % Update R forecast
    terms.Rfore = Rfore;

    testK = linspace(min(K_curr), max(K_curr), nk);
    foredist = compute.dist(p_new, p_old, 3);
    nlms{iter_ct, 1, 1} = nlm11; nlms{iter_ct, 1, 2} = nlm11;
    nlms{iter_ct, 2, 1} = nlm21; nlms{iter_ct, 2, 2} = nlm22;
    
    fprintf('K(1) = %0.4f\n', Kprdat(1));
    fprintf('K range: [%0.4f, %0.4f]\n', min(Kprdat), max(Kprdat));

    figure;
        
    set(groot, 'defaultTextInterpreter', 'latex')
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
    set(groot, 'defaultLegendInterpreter', 'latex')
% 1) Three-panel plot
    subplot(3,1,1);
    plot(Prdata(1500:2000),'LineWidth',1.6,'Color',[0.80 0.20 0.20]);
    yline(0.5,'--','Color',[0.4 0.4 0.4],'LineWidth',1);
    xlabel('Time');
    ylabel('P(populism)');
    title('Probability of Voting for Populism');
    xlim([0 500]);
    grid on;
    set(gca,'FontSize',18);
    
    subplot(3,1,2);
    plot(Kprdat(1500:2000),'LineWidth',1.6,'Color',[0.20 0.30 0.75]);
    xlabel('Time');
    ylabel('Aggregate Capital $K$', 'Interpreter', 'latex');
    title('Aggregate Capital Path');
    xlim([0 500]);
    grid on;
    set(gca,'FontSize',18);
    
    subplot(3,1,3);
    plot(ddata(1500:2000),'LineWidth',1.6,'Color',[0.5 0.25 0.42]);
    xlabel('Time');
    ylabel('$\delta$', 'Interpreter', 'latex');
    title('Capital Depreciation Shock');
    xlim([0 500]);
    grid on;
    set(gca,'FontSize',18);
    
% Kfore_new and Rfore_new are 4×2: [a  b]
% rows: (R,d) = (1,1),(1,2),(2,1),(2,2)

    fprintf('\nCapital law (log K'' = a + b log K):\n');
    fprintf('  R=1, d=1: a = %0.4f, b = %0.4f\n', Kfore_new(1,1,1), Kfore_new(1,1,2));
    fprintf('  R=2, d=1: a = %0.4f, b = %0.4f\n', Kfore_new(1,2,1), Kfore_new(1,2,2));
    fprintf('  R=1, d=2: a = %0.4f, b = %0.4f\n', Kfore_new(2,1,1), Kfore_new(2,1,2));
    fprintf('  R=2, d=2: a = %0.4f, b = %0.4f\n', Kfore_new(2,2,1), Kfore_new(2,2,2));
    
    fprintf('\nRegime transition (Pr(R''=1 | R,d,K) = exp((beta0 + beta1 log K)^-1)):\n');
    fprintf('  R=1, d=1: beta0 = %3.4f, beta1 = %3.4f\t (%i Periods)\n', ...
        Rfore_new(1,1,1), Rfore_new(1,1,2), sum(ix11));
    fprintf('  R=2, d=2: beta0 = %3.4f, beta1 = %3.4f\t (%i Periods)\n', ...
        Rfore_new(1,2,1), Rfore_new(1,2,2), sum(ix12));
    fprintf('  R=2, d=1: beta0 = %3.4f, beta1 = %3.4f\t (%i Periods)\n', ...
        Rfore_new(2,1,1), Rfore_new(2,1,2), sum(ix21));
    fprintf('  R=2, d=2: beta0 = %3.4f, beta1 = %3.4f\t (%i Periods)\n',...
        Rfore_new(2,2,1), Rfore_new(2,2,2), sum(ix22));
    
    fprintf('\nCapital regressions (log K'' on log K):\n');
    
    if exist('mdl11','var')
        fprintf('  R=1, d=1: R2 = %0.6f\n', mdl11.Rsquared.Adjusted);
    else
        fprintf('  R=1, d=1: R2 =   n/a (too few obs)\n');
    end
    
    if exist('mdl12','var')
        fprintf('  R=1, d=2: R2 = %0.6f\n', mdl12.Rsquared.Adjusted);
    else
        fprintf('  R=1, d=2: R2 =   n/a (too few obs)\n');
    end
    
    if exist('mdl21','var')
        fprintf('  R=2, d=1: R2 = %0.6f\n', mdl21.Rsquared.Adjusted);
    else
        fprintf('  R=2, d=1: R2 =   n/a (too few obs)\n');
    end
    
    if exist('mdl22','var')
        fprintf('  R=2, d=2: R2 = %0.6f\n', mdl22.Rsquared.Adjusted);
    else
        fprintf('  R=2, d=2: R2 =   n/a (too few obs)\n');
    end
        
    fprintf('Regime guessing distance: %1.4f\n', foredist);

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

    cd ../d/
    save(filename)
end

