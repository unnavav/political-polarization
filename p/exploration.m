%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setting up original household problem and solving it
% vaasavi
% march 2024
% 
% solving via grid search,because there's noise and 
% ergo i donut trust splines (not sure if this instinct)
% is correct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

restoredefaultpath;
clear all; clc;
addpath(genpath(pwd));

% first set up some grids, pulling a lot from aiyagari
% project
% delete(gcp('nocreate'));
% parpool('local',4);

vTol = 1e-4; gTol = 1e-6;
%% params
alpha = 0.36; delta = 0.08; beta = 0.96173; sigma = 1; phi = 1; lambda = .07;

nl = 7;
na = 250;
nmu = 2500;
np = 2;

al = 0; ah = 100;

wage = 1.1; r = 0.04;
%% get labor distribution and aggregate value
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
klmult = 1.025;
khmult = 1.2;

kl = klwrbnd*klmult;
kh = klwrbnd*khmult;

%% make asset grid

agrid = compute.logspace(al, ah, na);

%% make distribution grid

amu = linspace(al, ah, nmu);

%% make party policy grid

tgrid = [.7 .3];

%% initialize value functions
V = zeros(nl, na);

for il = 1:nl
    l = lgrid(il);
    for ia = 1:na
        a = agrid(ia);

        max_cons = wage*l + (1+r)*a - r*phi;
            
        if sigma == 1
            V(il, ia) = log(max_cons)/(1-beta); 
        else
            V(il, ia) = (max_cons^(1-sigma))/((1-beta)*(1-sigma));
        end
        
    end
end

% p = .5; % initial 50/50 chance of getting any party

%% solve value function: grid lookup

dist = 10;
EV = zeros(nl, na);
g = zeros(nl, na);
V1 = V;
g1 = g;
DEV = EV;
iter_ct = 1;

% start initial r guess at a different 
kval = kl;
r = alpha*(kval^(alpha - 1)*(lagg^(1-alpha))) - delta;
wage = (1-alpha)*((kval^(alpha))*(lagg^(-alpha)));

while dist > vTol
   
    % get expected value function 

    %TODO ADD PARTY DIFFERENCES
    for il = 1:nl
        EV(il, :) = pil(il,:)*V(:,:);
    end

    %DEV approx:

    DEV = egm.numdev(EV, agrid);

    % now just make choices for assets

    E = zeros(size(DEV));

    for il = 1:nl
        l = lgrid(il);
        for ia = 1:na
            apr = agrid(ia);

            c = (beta*DEV(il, ia))^(-1/sigma);
            a = (c + apr - (wage*l - r*phi))/(1+r);
    
            E(il, ia) = a; 
        end
    end

    tiledlayout(3,1);
    nexttile
    mesh(V)
    nexttile
    mesh(DEV)
    nexttile
    mesh(E)

    for il = 1:nl
        achoices = E(il,:)';

        for ia = 1:na

            lb = achoices(1);
            
            ahat = agrid(ia);

            if ahat < lb
                g(il, ia) = 0;
            else
                [ix, we] = compute.weight(achoices, ahat);
                g(il, ia) = we*agrid(ix) + (1-we)*agrid(ix + 1);
            end

            c = (1+r)*ahat + wage*lgrid(il) - r*phi - g(il, ia);

            if c < 0
                disp([il ia])
            end

            if sigma == 1
                V(il, ia) = log(c) + beta*V1(il, ia);
            else
                V(il, ia) = (c^(1-sigma))/(1-sigma) + ...
                    beta*V1(il, ia);
            end 

        end
    end

    dist = compute.dist(V1, V,2);

    if mod(iter_ct, 10) == 0
        fprintf("Iteration %i: ||TV - V|| = %4.7f\n", iter_ct, dist);
    end

    iter_ct = iter_ct + 1;

    V1 = V;

    tiledlayout(3,1);
    nexttile
    mesh(V)
    nexttile
    mesh(E)
    nexttile
    mesh(g)

end

tiledlayout(4,1);
nexttile
mesh(V)
% nexttile
% mesh(V(:,:,2))
nexttile
mesh(g)
% nexttile
% mesh(g(:,:,2))