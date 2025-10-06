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

cd C:\VAASAVI\Dropbox\Education\OSU\Ongoing_Research\Populism\political-polarization\p
restoredefaultpath;
clear all; clc;
addpath(genpath(pwd));

% first set up some grids, pulling a lot from aiyagari
% project

%% params
vTol = 1e-5; dTol = 1e-2;
alpha = 0.36; delta = 0.06; beta = 0.96; sigma = 1; phi = 0;

neta = 20;
ntau = 20;
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
rho = .9;
sig_l = 0.2;
 
% need to back out sigma^2_e given sigma^2_l
sigx = sig_l*sqrt(1 - rho^2);
range = 2.575;

[pil, lgrid] = compute.getTauchen(nl,  mu, sigx, rho, range);

% get steady state labor distribution

ldist = asymptotics(dtmc(pil));

lagg = ldist*lgrid';

kagg = ((r+delta)/(alpha*lagg^(1-alpha)))^(1/(alpha-1));

% get bounds for capital guesses
rst = 1.0/beta - 1.0;

klwrbnd = (rst + delta)/(alpha);
klwrbnd = klwrbnd/(lagg^(1-alpha));
klwrbnd = klwrbnd^(1/(alpha-1));
klmult = .05;
khmult = 6;

kl = klwrbnd*klmult;
kh = klwrbnd*khmult;

% make asset grid

agrid = compute.logspace(al, ah, na)';
amu = linspace(al, ah, nmu);

adistr = zeros(nl, nmu);

for im = 1:nmu
    for il = 1:nl
        adistr(il,im) = 1.0/(nmu*nl);
    end 
end

etagrid = linspace(.0,.9,neta);
taugrid = linspace(.0,.9,ntau);
% captax = [repelem(0, ceil(nl/2)) repelem(.15, floor(nl/2))]; %IRS
% captax = compute.logspace(0, 20, nl)/100;
% captax = linspace(0, .20, nl);
% captax = repelem(.15, nl);
% captax = repelem(0, nl);
captax = repelem(0, nl);

G = [0 0];

kDist = 10;
gDist = 10; 

kval = .5;

DIST = max(kDist,gDist);

Varray = cell(1,2);
Garray = cell(1,2);
EVarray = cell(1,2);
Warray = cell(1,2);
Karray = cell(1,2);
allkp = cell(5,2);

TV = zeros(nl, na);
TG = TV; V= TV;

EV = zeros(nl, na);

%prepare for VFI
terms = struct('beta', beta, ...
    'sigma', sigma, ...
    'phi', phi, ...
    'agrid', agrid, ...
    'lgrid', lgrid, ...
    'pil', pil, ...
    'captax', zeros(1,nl));

folname = "rsteadystates";
mkdir("../d/",folname)
newdir = strcat("../d/", folname);
cd(newdir)

eta_pop = etagrid(1);
tau_pop = taugrid(6);

% Liberal policy grid to loop over
eta_L_vals = etagrid(2:7);       % Liberal migration rates
tau_L_vals = flip(taugrid(1:5));       % Liberal tax rates

% now set up for actual VFI:
r = 0.04; w = 1.3;
scale = 1;
for ia = 1:na
    k_val = agrid(ia);
    for il = 1:nl
        yval = scale*(1+r)*k_val + w*lgrid(il) - r*phi;
        ymin = max(1e-10, yval);
        V(il, ia) = log(ymin);
    end
end


EVarray{1,1} = predict.getExpectation(V, pil);
EVarray{1,2} = predict.getExpectation(V, pil);

kl = 0;
kh = 15;
kval = (kl + kh)/2;
kDist = 10; adj = .3;

p = .5;
pnew = p;

%% TAXATION

fprintf("TAXATION BITCHES -----------------------------------------------\n\n")

for i = 1:length(tau_L_vals)
    tau_lib = tau_L_vals(i);
    eta_lib = eta_pop;

    % Set policies
    tpolicygrid = [tau_pop tau_lib];
    epolicygrid = [eta_pop eta_lib];

    % Run your kp-solver here...
    % (insert solver block from earlier, using epolicygrid/tpolicygrid)

    kl = 0;
    kh = 15;
    kval = (kl + kh)/2;
    kDist = 10; adj = .5;
    kp = [0 0];
    while kDist > dTol

        K = kval;
        fprintf("\n*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*")
        fprintf("\nA guess: %4.8f. Begin iteration for solution...\n", kval)
        fprintf("Guessing Probabilities:\n")
        
        pdist = 10;
        while pdist > dTol
            EV = p*EVarray{1} + (1-p)*(EVarray{2});
            % solve for VT:
            terms.r = vaas.calcr(alpha, delta, K, epolicygrid(1));
            terms.w = vaas.calcw(alpha, K, epolicygrid(1));
    
            terms.tau = tpolicygrid(1);
            % we have to get the value of lambda such that taxation 
            % is redistributing everything. aka BB
            tot_inc = (terms.w*lgrid).*ldist;
            tot_inc = sum(tot_inc);
            denom = (terms.w*lgrid).^(1-terms.tau).*ldist;
            denom = sum(denom);
            terms.lamval = tot_inc/denom;
    
            [Varray{1}, Garray{1}] = solve.interpV(nl, na, EV, terms, dTol);
    
            terms.r = vaas.calcr(alpha, delta, K, epolicygrid(2));
            terms.w = vaas.calcw(alpha, K, epolicygrid(2));
    
            terms.tau = tpolicygrid(2);
            % we have to get the value of lambda such that taxation 
            % is redistributing everything. aka BB
            tot_inc = (terms.w*lgrid).*ldist;
            tot_inc = sum(tot_inc);
            denom = (terms.w*lgrid).^(1-terms.tau).*ldist;
            denom = sum(denom);
            terms.lamval = tot_inc/denom;
    
            [Varray{2}, Garray{2}] = solve.interpV(nl, na, EV, terms, dTol);
    
            EVarray{1} = predict.getExpectation(Varray{1}, pil);
            EVarray{2} = predict.getExpectation(Varray{2}, pil);
    
            [Warray{1}, Karray{1}] = HH.getDist(Garray{1}, amu, agrid, ...
                pil, false);  
    
            acond = compute.condense(Warray{1}, amu, agrid);
            VOTES = EVarray{1,1} > EVarray{1,2};
            VOTES = acond .* VOTES;
            pnew = sum(sum(VOTES));
            
            pdist = abs(pnew-p);
            fprintf("\tProbability diff: |%0.2f - %0.2f| = %0.2f\n", pnew, ...
                p, pdist);
            p = (p+pnew)/2;
        end
    
        kp = [kp; kval p];
    
        [Warray{2}, Karray{2}] = HH.getDist(Garray{2}, amu, agrid, ...
              pil, false);            
        kdist = Karray{1} - kval;
    
        if kdist > 0
            fprintf("\n||Kguess - Kagg|| = %4.5f. \tAggregate capital is too low.\n", abs(kdist))
            kl = (1-adj)*kval+adj*kl;
        else
            fprintf("\n||Kguess - Kagg|| = %4.5f. \tAggregate capital is too high.\n", abs(kdist))
            kh = (1-adj)*kval+(adj)*kh;
        end
    
        kDist = abs(kdist);  
        kval = .5*(kl + kh);
    end

    filename = strcat("results_rho85sig2_t_",sprintf("%i",i),"6_eta11.mat");
    save(filename)

    kp = kp(2:end,:);
    all_kp{i,1} = kp;  % store result
end

fprintf("Percentage Voting for Populists: %0.2f\n", p);

%% now migration
fprintf("MIGRATION BITCHES------------------------------------------------")

p = .5;
pnew = p;

for i = 1:length(eta_L_vals)
    tau_lib = tau_pop;
    eta_lib = eta_L_vals(i);

    % Set policies
    tpolicygrid = [tau_pop tau_lib];
    epolicygrid = [eta_pop eta_lib];

    % Run your kp-solver here...
    % (insert solver block from earlier, using epolicygrid/tpolicygrid)

    kl = 0;
    kh = 15;
    kval = (kl + kh)/2;
    kDist = 10; adj = .5;
    kp = [0 0];
    while kDist > dTol

        K = kval;
        fprintf("\n*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*")
        fprintf("\nA guess: %4.8f. Begin iteration for solution...\n", kval)
        fprintf("Guessing Probabilities:\n")
        
        pdist = 10;
        while pdist > dTol
            EV = p*EVarray{1} + (1-p)*(EVarray{2});
            % solve for VT:
            terms.r = vaas.calcr(alpha, delta, K, epolicygrid(1));
            terms.w = vaas.calcw(alpha, K, epolicygrid(1));
    
            terms.tau = tpolicygrid(1);
            % we have to get the value of lambda such that taxation 
            % is redistributing everything. aka BB
            tot_inc = (terms.w*lgrid).*ldist;
            tot_inc = sum(tot_inc);
            denom = (terms.w*lgrid).^(1-terms.tau).*ldist;
            denom = sum(denom);
            terms.lamval = tot_inc/denom;
    
            [Varray{1}, Garray{1}] = solve.interpV(nl, na, EV, terms, dTol);
    
            terms.r = vaas.calcr(alpha, delta, K, epolicygrid(2));
            terms.w = vaas.calcw(alpha, K, epolicygrid(2));
    
            terms.tau = tpolicygrid(2);
            % we have to get the value of lambda such that taxation 
            % is redistributing everything. aka BB
            tot_inc = (terms.w*lgrid).*ldist;
            tot_inc = sum(tot_inc);
            denom = (terms.w*lgrid).^(1-terms.tau).*ldist;
            denom = sum(denom);
            terms.lamval = tot_inc/denom;
    
            [Varray{2}, Garray{2}] = solve.interpV(nl, na, EV, terms, dTol);
    
            EVarray{1} = predict.getExpectation(Varray{1}, pil);
            EVarray{2} = predict.getExpectation(Varray{2}, pil);
    
            [Warray{1}, Karray{1}] = HH.getDist(Garray{1}, amu, agrid, ...
                pil, false);  
    
            acond = compute.condense(Warray{1}, amu, agrid);
            VOTES = EVarray{1,1} > EVarray{1,2};
            VOTES = acond .* VOTES;
            pnew = sum(sum(VOTES));
            
            pdist = abs(pnew-p);
            fprintf("\tProbability diff: |%0.2f - %0.2f| = %0.2f\n", pnew, ...
                p, pdist);
            p = (p+pnew)/2;
        end
    
        kp = [kp; kval p];
    
        [Warray{2}, Karray{2}] = HH.getDist(Garray{2}, amu, agrid, ...
              pil, false);            
        kdist = Karray{1} - kval;
    
        if kdist > 0
            fprintf("\n||Kguess - Kagg|| = %4.5f. \tAggregate capital is too low.\n", abs(kdist))
            kl = (1-adj)*kval+adj*kl;
        else
            fprintf("\n||Kguess - Kagg|| = %4.5f. \tAggregate capital is too high.\n", abs(kdist))
            kh = (1-adj)*kval+(adj)*kh;
        end
    
        kDist = abs(kdist);  
        kval = .5*(kl + kh);
    end

    filename = strcat("results_rho85sig2_t_66_eta1",sprintf("%i",i),".mat");
    save(filename)

    kp = kp(2:end,:);
    all_kp{i,2} = kp;  % store result
end

fprintf("Percentage Voting for Populists: %0.2f\n", p);


cd ../../p/

%% just graphing this prettily

cd ../v/Paper/

skp = sort(kp);
smooth= smoothdata(skp(:,2), 'gaussian', 5);

hexColors = ["#11302a", "#036264", "#8f5774", "#e0829d", "#dac4d0"];
rgbColors = sscanf(join(hexColors,''),'#%2x%2x%2x',[3,length(hexColors)]).'/255;

% Interpolate to 256 steps
n = 5;
darkAcademiaChic = interp1(linspace(0,1,size(rgbColors,1)), rgbColors, linspace(0,1,n));

fig = figure;
set(fig, 'Units', 'inches', 'Position', [1, 1, 8.53, 4.8]);

set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [8.53, 4.8], ...
         'PaperPosition', [0, 0, 8.53, 4.8]);
hold on;
plot(skp(:,1),smooth, 'Color','#dac4d0','LineWidth',2);
% Add original points as circles
plot(skp(:,1),skp(:,2), 'o', 'MarkerSize', 3, ...
    'MarkerEdgeColor', '[0.2 0.2 0.6]', 'MarkerFaceColor', '[0.2 0.2 0.6]');
yline(.5, '--', 'Regime Switch', 'LabelVerticalAlignment','bottom', ...
    'Alpha', 0.6, 'Interpreter', 'latex');
xlabel('Aggregate Capital ($K$)', ...
      'FontSize', 18, 'Interpreter', 'latex');
ylabel('Populist Votes ($p$)', ...
      'FontSize', 18, 'Interpreter', 'latex');
title('Populist Voting Curve', ...
      'FontSize', 22, 'Interpreter', 'latex');
hold off;
legend('Smoothed Vote Share', 'Model Data Points', 'Location', 'southeast');
print(fig, 'populist_voting_curve.pdf', '-dpdf', '-vector');
