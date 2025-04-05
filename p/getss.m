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

vTol = 1e-5; dTol = 1e-2;
%% params
alpha = 0.36; delta = 0.06; beta = 0.96; sigma = 3; phi = 0;

neta = 6;
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

% make distribution grid

amu = linspace(al, ah, nmu);

% lambda upper lower bounds
ll = 0; lh = 1;
lamval = (lh + ll)/2; % people are really reactive to taxes
adj = 1/3;

%party tax regimes
tgrid = [0 .2]; 
proggrid = [.70 .25 .05 0 0 0 0];
% tgrid = [.15 .05]; 
etagrid = linspace(.0,.01,neta);

% captax = [repelem(0, ceil(nl/2)) repelem(.15, floor(nl/2))]; %IRS
% captax = compute.logspace(0, 20, nl)/100;
% captax = linspace(0, .20, nl);
% captax = repelem(.15, nl);
% captax = repelem(0, nl);
captax = repelem(0, nl);

G = [0 0];

folname = "steadystates";
mkdir("../d/",folname)


adistr = zeros(nl, nmu);

for im = 1:nmu
    for il = 1:nl
        adistr(il,im) = 1.0/(nmu*nl);
    end 
end


% first init equilibrium objects like prices, govt spending, etc
kval = (kh+kl)/2;
eta = etagrid(1);
r = vaas.calcr(alpha, delta, kval/lagg, eta);
wage = vaas.calcw(alpha, kval/lagg, eta);

wage_inc_mu = repmat(wage*lgrid, nmu,1)';
cap_inc_mu = zeros(size(wage_inc_mu));

for il = 1:nl
    cap_inc_mu(il, :) = repmat(r*amu, 1, 1)*(1-captax(il));
end
Tmu = zeros(nl, nmu, np);
for ip = 1:np
Tmu(:,:,ip) = gov.tax(wage_inc_mu,0,tgrid(ip)) + ...
    0*r*cap_inc_mu;
end

Gvals = sum(sum(Tmu(:,:,:).*adistr(:,:)));
Gvals = squeeze(Gvals)';
lumpsum = .7.*Gvals;
progressive =.3.*Gvals;
g(1,:) = lumpsum(1) + progressive(1).*proggrid;
g(2,:) = lumpsum(2) + progressive(2).*proggrid;    

kDist = 10;
gDist = 10; 

kval = .5;
g = zeros(np, nl);
terms.lamval = 0; % only migration in this place

DIST = max(kDist,gDist);

Varray = cell(neta,1);
Garray = cell(neta,1);
EVarray = cell(neta,1);
Warray = cell(neta,1);
Karray = cell(neta,1);
parray = cell(neta,1);

%% migration ss's

TV = zeros(nl, na);
TG = TV; V= TV;

EV = zeros(nl, na);


for i = 5:neta
    eta = etagrid(i);
    kl = 6.83;
    kh = 6.9;
    kval = (kl + kh)/2;
    kDist = 10;
    fprintf("Migration Rate: %0.4f\n", eta);

    while kDist > dTol
    
        fprintf("\n*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*")
        fprintf("\nA guess: %4.8f. Begin iteration for solution...\n", kval)
        fprintf("\t Solving value function:\n")
        
       %prepare for VFI
        terms = struct('beta', beta, ...
            'sigma', sigma, ...
            'phi', phi, ...
            'agrid', agrid, ...
            'lgrid', lgrid, ...
            'pil', pil, ...
            'lamval', 0, ...
            'captax', zeros(1,nl), ...
            'tau', 0, ...
            'g', g);
    

        terms.r = vaas.calcr(alpha, delta, kval, eta);
        terms.w = vaas.calcw(alpha, kval, eta);
    
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
%                 [TV, TG, ~, EV] = egm.solve(nl, na, 1, terms, EV, V);
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
            if dTol < 1
                V = V + 0.1*(TV-V);
            else
                V = TV;
            end
             
            G = TG;
        end

        cap_inc = repmat((1+terms.r*agrid), nl,1);
        wages = [terms.w*lgrid - gov.tax(terms.w*lgrid, terms.lamval, terms.tau)]';
        wage_inc = repmat(wages,1,na);
        y = cap_inc + wage_inc;
    


        Varray{i}= V;
        Garray{i} = G;
        EVarray{i} = EV; 

        [Warray{i}, Karray{i}] = HH.getDist(Garray{i}, amu, agrid, ...
            pil, false);      
    %     adistr = Warray{1};
    %     Gvals = sum(sum(Tmu(:,:,:).*adistr(:,:)));
    %     Gvals = squeeze(Gvals)';
    %     lumpsum = .7.*Gvals;
    %     progressive =.3.*Gvals;
    %     g(1,:) = lumpsum(1) + progressive(1).*proggrid;
    %     g(2,:) = lumpsum(2) + progressive(2).*proggrid;    
    
        kdist = Karray{i} - kval;
    
        %if you're looking at this, ignore it. it is chaos incarnate.
        if kdist < 10
            adj = .5;
            vTol = 1e-4;
        else
            adj = .7;
            vTol = 1e-6;
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

    acond = compute.condense(Warray{i}, amu, agrid);
    Votes_Pop = EVarray{i} < EVarray{1};
    p = Votes_Pop.*acond;
    parray{i} = sum(sum(p));
    fprintf("Percentage Voting for Populists: %0.2f", parray{i});
end 
temp  = EVarray{i};
newdir = strcat("../d/", folname);
cd(newdir)

filename = strcat("results_t",sprintf('%0.4f', 0), ...
    "_eta", sprintf('%0.4f', eta(2)), ".mat");
save(filename, 'adistr', "G","V", "wage", "r", "kagg", ...
    "lamval", "eta", 'tgrid', 'growth_lagg', 'delta', 'alpha')

filename = strcat("terms_struct_t", sprintf('%0.4f', 0), ...
    "_eta", sprintf('%0.4f', eta(2)), ".mat");
save(filename, "terms")

cd ../../p/


%% now changing growth paths

for i = 1:2

    adistr = zeros(nl, nmu);
    
    for im = 1:nmu
        for il = 1:nl
            adistr(il,im) = 1.0/(nmu*nl);
        end
    end
    
    growth_lagg = lagg*eta(1)*eta(2);


    kDist = 10;
    gDist = 10; 
    f=10;  
    
    %init from eqm i've already solved for
    kval = .5;
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
            r = alpha*(kval^(alpha - 1)*(growth_lagg^(1-alpha))) - delta;
            wage = (1-alpha)*((kval^(alpha))*(growth_lagg^(-alpha)));
    
            wage_inc_mu = repmat(wage*lgrid, nmu,1)';
            cap_inc_mu = zeros(size(wage_inc_mu));
            for il = 1:nl
                cap_inc_mu(il, :) = repmat(r*amu, 1, 1)*(1-captax(il));
            end
            Tmu = zeros(nl, nmu, np);
            for ip = 1:np
                Tmu(:,:,ip) = wage_inc_mu - gov.tax(wage_inc_mu,lamval,tgrid(ip)) + ...
                    cap_inc_mu;
            end
    
            G = sum(sum(Tmu(:,:,i).*adistr(:,:)));
    
            %prepare for VFI
            terms = struct('alpha', alpha, ...
                'beta', beta, ...
                'sigma', sigma, ...
                'phi', phi, ...
                'agrid', agrid, ...
                'lgrid', lgrid, ...
                'pil', pil, ...
                'lamval', 0, ...
                'captax', captax, ...
                'tau',0, ...
                'G', G, ...
                'p', p, ...
                'K', kagg, ...
                'L', growth_lagg, ...
                'r', r, ...
                'w', wage);
        
            fprintf("Solving value function:\n")
    
            [V, G, EV] = HH.solve(nl, na, terms, vTol);
         
    
            fprintf("\tSolving asset distribution:\n")
            [adistr, kagg] = HH.getDist(G, amu, agrid, pil);
        
    %         % CONDENSE DISTR
            acond = compute.condense(adistr, amu, agrid);
                
            kdist = kagg - kval;
        
            if kdist > 0
                fprintf("\n||Kguess - Kagg|| = %4.5f. \tAggregate capital is too low.\n", abs(kdist))
                kl = (kval+5*kl)/6;
            else
                fprintf("\n||Kguess - Kagg|| = %4.5f. \tAggregate capital is too high.\n", abs(kdist))
                kh = (kval+5*kh)/6;
            end
        
            kDist = abs(kdist);  
            kval = .5*(kl + kh);
        end
        
        t = adistr.*Tmu(:,:,i); %getting all taxes collected in A govt
        t = sum(sum(t));
        Y = kval^(alpha)*growth_lagg^(1-alpha);
        gy = t/Y;
    
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
    
    newdir = strcat("../d/", folname);
    cd(newdir)

    filename = strcat("results_transition_t",sprintf('%0.4f', tgrid(i)), ...
        "_eta", sprintf('%0.4f', eta(i)), ".mat");
    save(filename, 'adistr', "G","V", "wage", "r", "kagg", ...
        "lamval", "eta", 'tgrid')

    save("params", "terms")

    cd ../../p/

end

