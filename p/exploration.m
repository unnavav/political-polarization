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

 
    %g choice: maximize expected value

    for il = 1:nl
        l = lgrid(il);
        for ia = 1:na

            a = agrid(ia);

            y = (1+r)*a + wage*l - r*phi;

            yvec = ones(250,1)*y;

            C = yvec - agrid;
            C = C(C>0);

            if sigma == 1
                Val = log(C)' + beta*EV(il, 1:size(C,1));
            else
                Val = (C.^(1-sigma))/(1-sigma)' + beta*EV(il, 1:size(C,1));
            end
            
            [val, ai] = max(Val);

            V(il, ia) = val;
            g(il, ia) = agrid(ai);

        end
    end

    dist = compute.dist(V, V1, 2);

    if mod(iter_ct, 10) == 0
        fprintf("Iteration %i: ||TV - V|| = %4.7f\n", iter_ct, dist);
    end

    iter_ct = iter_ct + 1;

    V1 = V;
    g1 = g;

end

tiledlayout(2,1);
nexttile
mesh(V)
% nexttile
% mesh(V(:,:,2))
nexttile
mesh(g)
% nexttile
% mesh(g(:,:,2))