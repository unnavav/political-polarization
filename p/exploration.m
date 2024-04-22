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
alpha = 0.36; delta = 0.08; beta = 0.96173; sigma = 1; phi = 2;

identity = 0.01;

nl = 7;
na = 250;
nmu = 2500;
np = 2;

al = 1; ah = 101;

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

tgrid = [.1 .5];
lamgrid = [.07 .21];
%% initialize value functions
V = zeros(nl, na, np, np); %(labor, assets, my type, goverment type)

for il = 1:nl
    l = lgrid(il);
    for ia = 1:na
        a = agrid(ia);

        max_cons = wage*l + (1+r)*a - r*phi;
            
        for ip = 1:np

            max_cons = gov.tax(max_cons, lamgrid(ip), tgrid(ip));
            
            if sigma == 1
                V(il, ia, ip, ip) = log(max_cons)/(1-beta); 
            else
                V(il, ia, ip, ip) = (max_cons^(1-sigma))/((1-beta)*(1-sigma));
            end
        end
    end
end

p = .5; % initial 50/50 chance of getting any party

%% solve value function: grid lookup

dist = 10;
EV = zeros(nl, na, np, np);
g = zeros(nl, na, np, np);
V1 = V;
g1 = g;
iter_ct = 1;

VOTES = zeros(nl, na, np);

% start initial r guess at a different 
kval = kl;
r = alpha*(kval^(alpha - 1)*(lagg^(1-alpha))) - delta;
wage = (1-alpha)*((kval^(alpha))*(lagg^(-alpha)));

% get expected value function to start with

%TODO ADD PARTY DIFFERENCES
for il = 1:nl
    for ip = 1:np
        EV(il, :, 1, ip) = (p*pil(il,:)*V(:,:,1,1) + ...
            (1-p)*pil(il,:)*V(:,:,1,2))';
        EV(il, :, 2, ip) = (p*pil(il,:)*V(:,:,2,1) + ...
            (1-p)*pil(il,:)*V(:,:,2,2))';
    end
end

while dist > vTol

 
    %g choice: maximize expected value

    for il = 1:nl
        l = lgrid(il);
        for ia = 1:na

            a = agrid(ia);

            y = (1+r)*a + wage*l - r*phi;

            yvec = ones(250,1)*y;

            C = yvec - agrid;
            C = C(C>0);
    
            for ip = 1:np

                nchoices = max(size(C));
                for ic = 1:nchoices
                    C(ic) = gov.tax(agrid(ic),lamgrid(ip), tgrid(ip));
                end

                %C = (1-tgrid(ip))*C;

                if sigma == 1
                    Val = log(C)' + beta*EV(il, 1:size(C,1),ip);
                else
                    Val = (C.^(1-sigma))/(1-sigma)' + beta*EV(il, 1:size(C,1), ip);
                end
                
                [val, ai] = max(Val);
    
                if ip - 1 > 0
                    V(il, ia, 1, ip) = val;
                    g(il, ia, 1, ip) = agrid(ai);
                    V(il, ia, 2, ip) = val + identity;
                    g(il, ia, 2, ip) = agrid(ai);
                else 
                    V(il, ia, 1, ip) = val + identity;
                    g(il, ia, 1, ip) = agrid(ai);
                    V(il, ia, 2, ip) = val;
                    g(il, ia, 2, ip) = agrid(ai);
                end
            end
        end
    end

    % get expected value function after getting value function

    %TODO ADD PARTY DIFFERENCES
    for il = 1:nl
        for ip = 1:np
            EV(il, :, 1, ip) = pil(il,:)*V(:,:,1,ip);
            EV(il, :, 2, ip) = pil(il,:)*V(:,:,2,ip);
        end
    end

    % calculate voting decision

    VOTES(:,:,1) = EV(:,:,1,1) >= EV(:,:,1,2);
    VOTES(:,:,2) = EV(:,:,2,1) >= EV(:,:,2,2);

    dist = compute.dist(V, V1, 4);

    if mod(iter_ct, 10) == 0
        fprintf("Iteration %i: ||TV - V|| = %4.7f\n", iter_ct, dist);
    end

    iter_ct = iter_ct + 1;

    V1 = V;
    g1 = g;

end

tiledlayout(4,1);
nexttile
mesh(V(:,:,1,1))
nexttile
mesh(V(:,:,2,1))
nexttile
mesh(VOTES(:,:,2))
nexttile
mesh(EV(:,:,1,1) - EV(:,:,1,2))