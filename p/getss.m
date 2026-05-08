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
alpha = 0.36; delta = 0.04; beta = 0.96; sigma = 2; 
neta = 5;
ntau = 5;
nl = 7;
na = 100;
nmu = na*10;
np = 2;
nd = 2;

al = 0; ah = 50;

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

% delta shock grid
dgrid = [0.04 0.08];
pid = [.94 0.06;
    .9 .1]; %totally made this up

etagrid = linspace(.0,.10,neta);
taugrid = linspace(.09,.19,ntau);
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

TV = zeros(nd, nl, na);
TG = TV; V= TV;

EV = zeros(nd, nl, na);

%prepare for VFI
terms = struct('beta', beta, ...
    'sigma', sigma, ...
    'agrid', agrid, ...
    'lgrid', lgrid, ...
    'pil', pil, ...
    'captax', zeros(1,nl));

folname = "steadystates";
mkdir("../d/",folname)
newdir = strcat("../d/", folname);
cd(newdir)


% Get today's date as a datetime object
t = datetime('today');

% Convert the datetime object to a string with the specified format
todayDateStr = string(t, 'yyyyMMdd');


for i = 1:neta
    eta = etagrid(i);

    for j = 1:ntau

        terms.tau = taugrid(j);
        fprintf("Migration Rate: %0.4f, Progressivity: %0.4f\n", eta, terms.tau);

        kl = 5;
        kh = 20;
        kval = (kl + kh)/2;
        kDist = 10;
        while kDist > dTol
        
            fprintf("\n*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*")
            fprintf("\nA guess: %4.8f. Begin iteration for solution...\n", kval)
            fprintf("\t Solving value function:\n")
    
            rgrid = vaas.calcr(alpha, dgrid, kval, eta);
            terms.r = rgrid(1); % just an initialization value
            terms.w = vaas.calcw(alpha, kval, eta);
            if terms.r > 0
                %10 pct of natural borrowing limit
                phi = 0.1*(terms.w*min(lgrid))/terms.r;
            else
                phi = 0;
            end

            terms.phi  = 0; %terms.phi = phi when there's a bc, right now not playing that game 
    
            % we have to get the value of lambda such that taxation 
            % is redistributing everything. aka BB
            tot_inc = (terms.w*lgrid).*ldist;
            tot_inc = sum(tot_inc);
            denom = (terms.w*lgrid).^(1-terms.tau).*ldist;
            denom = sum(denom);
            terms.lamval = tot_inc/denom;
    
            iter_ct = 1;
            dist = 10; kdist = 10;
            G = zeros(nd, nl,na);
    
            % set up V so that it doesn't start empty
            scale = 1;
            for id = 1:nd
                for ia = 1:na
                    k_val = agrid(ia);
                    for il = 1:nl
                        yval = scale*(1+terms.r)*k_val + ...
                            terms.w*lgrid(il) - r*terms.phi;
                        ymin = max(1e-10, yval);
                        V(id,il, ia) = log(ymin);
                    end
                end
            end
            
            %init expected vals
            for id = 1:nd
                for ia = 1:na
                    for il = 1:nl
                        EV(id,il, ia) = pid(id,1)*dot(pil(il,:),V(1,:,ia)) + ...
                            (1-pid(id,1))*dot(pil(il,:),V(2,:,ia));
                    end
                end
            end
    
            while kdist > vTol
                
                for id = 1:nd
                    for ia = 1:na
                        for il = 1:nl
                            EV(id,il, ia) = pid(id,1)*dot(pil(il,:),V(1,:,ia)) + ...
                                (1-pid(id,1))*dot(pil(il,:),V(2,:,ia));
                        end
                    end
                end
                % now converging on the value function and decision rule
                % for every capital-regime combo (EGM bc this is 50
                % convergences)
                for id = 1:nd
                    terms.r = rgrid(id);
                    dV = squeeze(V(id, :,  :));
                    dEV = squeeze(EV(id, :,  :));
                    [TVd, TGd]= egm.solve(terms, dEV, dV);
                    TV(id,:,:) = TVd; TG(id,:,:) = TGd;
                end

                % check distance
                dist = compute.dist(V, TV, 3);
                kdist = compute.dist(G, TG, 3);
            
                if mod(iter_ct, 250) == 0
                    fprintf("\n\tIteration %i: \n\t\t||TV - V|| = %4.6f" + ...
                        "\n\t\t||TG - G|| = %4.6f", iter_ct, dist, kdist);
% %                     fprintf("\nInitial Values:");
% %                     fprintf("\nMin Kgrid: %2.4f, Max Kgrid: %2.4f", min(Kgrid), max(Kgrid));
% %                     fprintf("\nInitial Capital Forecasts: Kp = %2.4f, Kl = %2.4f", Kp, Kl);
% %                     fprintf("\nInitial Policy Function: Min Gp = %2.4f, Max Gp = %2.4f", ...
% %                         min(TGP(:)), max(TGP(:)));
% %                      fprintf("\nInitial Policy Function: Min Gl = %2.4f, Max Gl = %2.4f", ...
% %                         min(TGL(:)), max(TGL(:)));
% %                    fprintf("\nInitial Value Function: Min VP = %2.4f, Max VP = %2.4f", ...
% %                        min(TVP(:)), max(TVP(:)));
% %                    fprintf("\nInitial Value Function: Min VL = %2.4f, Max VL = %2.4f", ...
% %                        min(TVL(:)), max(TVL(:)));
                end
            
                iter_ct = iter_ct + 1;
            
                G = 0.4 * TG + 0.6 * G;
                V = 0.4 * TV + 0.6 * V;

            end  
            fprintf("\n\tIteration %i: \n\t\t||TV - V|| = %4.6f" + ...
                "\n\t\t||TG - G|| = %4.6f", iter_ct, dist, kdist);
            diff = TV - V;
            fprintf('\n\t\tMin diff: %4.6f, Max diff: %4.6f', min(diff(:)), max(diff(:)));
            fprintf('\n\t\tMean diff: %4.6f, Std diff: %4.6f', mean(diff(:)), std(diff(:)));

            Varray{i,j}= V;
            Garray{i,j} = G;
            EVarray{i,j} = EV; 
    
            
            [Warray{i,j}, Karray{i,j}] = HH.getDist(Garray{i,j}, amu, agrid, ...
                pil, pid, phi, false);      
        
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
    
        filename = strcat(todayDateStr, "_etatest_results_rho90sig2.mat");
        save(filename)

    end
end 

%populist voting
parray = zeros(neta,ntau);
% choose the ss to compare against
popi = 1; popj = 5;  % eta=0, tau=0.19
libi = 3; libj = 1;  % eta=0.05, tau=0.09

for i = 1:neta
    for j = 1:ntau
    Votes_Pop = EVarray{popi,popj} > EVarray{libi,libj};
    acond = compute.condense(Warray{i,j}, amu, agrid); 
    p = Votes_Pop.*acond;
    % ev_diff = EVarray{i,j} - EVarray{1,j};
    % meanEV_diffs(i,j) = mean(ev_diff(:));
    % stdEV_diffs(i,j)  = std(ev_diff(:));
    % maxEV_diffs(i,j)  = max(ev_diff(:));
    % minEV_diffs(i,j)  = min(ev_diff(:));
    
    parray(i,j) = sum(sum(sum(p)));
    fprintf('Eta = %.2f; Tau = %.2f\n', etagrid(i), taugrid(j))
    fprintf('\tVote share for populism at populist SS: %.4f\n', parray(popi,popj))
    fprintf('\tVote share for populism at liberal SS: %.4f\n', parray(libi,libj))

    % fprintf("Percentage Voting for Populists: %0.2f\n", parray{i,j});
    end
end
parray


popi = 1; popj = 3;  % eta=0, tau=0.19 — populist candidate
libi = 5; libj = 2;  % eta=0.05, tau=0.09 — liberal candidate

% At each SS distribution, share preferring populist over liberal
parray_switch_to_pop = zeros(neta, ntau);
% At each SS distribution, share preferring liberal over populist  
parray_switch_to_lib = zeros(neta, ntau);

for i = 1:neta
    for j = 1:ntau
        acond = compute.condense(Warray{i,j}, amu, agrid);

        % Who prefers populist platform
        Votes_Pop = EVarray{popi,popj} > EVarray{libi,libj};
        parray_switch_to_pop(i,j) = sum(sum(sum(Votes_Pop.*acond)));

        % Who prefers liberal platform
        Votes_Lib = EVarray{libi,libj} > EVarray{popi,popj};
        parray_switch_to_lib(i,j) = sum(sum(sum(Votes_Lib.*acond)));
    end
end

fprintf('Populist SS — share preferring populist: %.4f\n', parray_switch_to_pop(popi,popj))
fprintf('Populist SS — share preferring liberal: %.4f\n', parray_switch_to_lib(popi,popj))
fprintf('Liberal SS — share preferring populist: %.4f\n', parray_switch_to_pop(libi,libj))
fprintf('Liberal SS — share preferring liberal: %.4f\n', parray_switch_to_lib(libi,libj))

filename = strcat(todayDateStr, "_etatest_results_rho90sig2.mat");
save(filename)

cd ../../p/
