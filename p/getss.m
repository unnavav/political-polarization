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

Varray = cell(neta,ntau);
Garray = cell(neta,ntau);
EVarray = cell(neta,ntau);
Warray = cell(neta,ntau);
Karray = cell(neta,ntau);
parray = cell(neta,ntau);


meanEV_diffs = zeros(neta, ntau);
stdEV_diffs = zeros(neta, ntau);
maxEV_diffs = zeros(neta, ntau);
minEV_diffs = zeros(neta, ntau);

%% migration ss's

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

folname = "steadystates";
mkdir("../d/",folname)
newdir = strcat("../d/", folname);
cd(newdir)

for i = 1:neta
    eta = etagrid(i);
    kl = 0;
    kh = 15;
    kval = (kl + kh)/2;
    kDist = 10;
    fprintf("Migration Rate: %0.4f\n", eta);

    for j = 1:ntau

        terms.tau = taugrid(j);
    
        while kDist > dTol
        
            fprintf("\n*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*")
            fprintf("\nA guess: %4.8f. Begin iteration for solution...\n", kval)
            fprintf("\t Solving value function:\n")
    
            terms.r = vaas.calcr(alpha, delta, kval, eta);
            terms.w = vaas.calcw(alpha, kval, eta);
    
            % we have to get the value of lambda such that taxation 
            % is redistributing everything. aka BB
            tot_inc = (terms.w*lgrid).*ldist;
            tot_inc = sum(tot_inc);
            denom = (terms.w*lgrid).^(1-terms.tau).*ldist;
            denom = sum(denom);
            terms.lamval = tot_inc/denom;
    
            iter_ct = 1;
            dist = 10;
            G = zeros(nl,na);
    
            % set up V so that it doesn't start empty
            scale = 1;
            for ia = 1:na
                k_val = agrid(ia);
                for il = 1:nl
                    yval = scale*(1+terms.r)*k_val + terms.w*lgrid(il) - r*phi;
                    ymin = max(1e-10, yval);
                    V(il, ia) = log(ymin);
                end
            end
            
            %init expected vals
            for ia = 1:na
                for il = 1:nl
                    EV(il, ia) = pil(il,:)*V(:,ia);
                end
            end
    
            while dist > vTol && iter_ct < 500
                
    %             if iter_ct < 400
    %                 [TV, TG, ~, EV] = egm.solve(nl, na, terms, EV, V);
    %             else 
                    if iter_ct < 200
                        sTol = 1e-4;
                    else 
                        sTol = 1e-6;
                    end
                    [TV, TG, EV] = compute.interpV(terms, V, EV, sTol);
    %             end
        
                vdist = compute.dist(TV,V,2);
                gdist = compute.dist(TG,G,2);
                dist = max(vdist,gdist);
        
                if mod(iter_ct, 50) == 0
                    fprintf("\n\tIteration %i: \n\t\t||TV - V|| = %4.6f" + ...
                        "\n\t\t||TG - G|| = %4.6f", iter_ct, vdist, gdist);
                end
                iter_ct = iter_ct+1;
    
                V = TV;
                G = TG;
            end
    
            cap_inc = repmat((1+terms.r*agrid), nl,1);
            wages = [terms.w*lgrid - gov.tax(terms.w*lgrid, terms.lamval, terms.tau)]';
            wage_inc = repmat(wages,1,na);
            y = cap_inc + wage_inc;
       
            Varray{i,j}= V;
            Garray{i,j} = G;
            EVarray{i,j} = EV; 
    
            [Warray{i,j}, Karray{i,j}] = HH.getDist(Garray{i,j}, amu, agrid, ...
                pil, false);      
        
            kdist = Karray{i,j} - kval;
        
            %if you're looking at this, ignore it. it is chaos incarnate.
            if kdist < 10
                adj = .5;
            else
                adj = .7;
            end
    
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
    
        acond = compute.condense(Warray{i,j}, amu, agrid);
        if i > 1
            Votes_Pop = EVarray{i,j} < EVarray{1,j};
            p = Votes_Pop.*acond;
            ev_diff = EVarray{i,j} - EVarray{1,j};
            meanEV_diffs(i,j) = mean(ev_diff(:));
            stdEV_diffs(i,j)  = std(ev_diff(:));
            maxEV_diffs(i,j)  = max(ev_diff(:));
            minEV_diffs(i,j)  = min(ev_diff(:));
    
        else
            p = zeros(size(acond));
        end

        parray{i,j} = sum(sum(p));
        fprintf("Percentage Voting for Populists: %0.2f\n", parray{i,j});

        filename = strcat("results_rho85sig2_t",sprintf('%0.4f', taugrid(j)),"_eta", sprintf('%0.4f', etagrid(i)), ".mat");
        save(filename)

    end
end 

cd ../../p/
