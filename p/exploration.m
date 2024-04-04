%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setting up original household problem and solving it
% vaasavi
% march 2024
% 
% solving via grid search,because there's noise and 
% ergo i donut trust splines (not sure if this instinct)
% is correct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% first set up some grids, pulling a lot from aiyagari
% project
addpath(genpath(pwd));

% delete(gcp('nocreate'));
% parpool('local',4);

vTol = 1e-6; gTol = 1e-8;
%% params
alpha = 0.36; delta = 0.08; beta = 0.96173; sigma = 1; phi = 1; lambda = .07;

nl = 7;
na = 250;
nmu = 2500;
np = 2;

al = 0; ah = 100;

wage = 1.2482; r = 0.02;
%% get labor distribution and aggregate value
%step 1: make labor grid and labor transition matrix

mu = 0;
rho = .9554;
sig_l = 0.5327;

% need to back out sigma^2_e given sigma^2_l
sigx = sig_l*sqrt(1 - rho^2);

[pil, lgrid] = compute.getTauchen(nl,  mu, sigx, rho);

% get steady state labor distribution

ldist = ones(1,nl)/nl;

dist = 10;
while dist > gTol
    ldist_1 = ldist*pil;
    dist = compute.dist(ldist_1, ldist, 1);
    ldist = ldist_1;
end

lagg = ldist*lgrid';

% get bounds for capital guesses
rst = 1.0/beta - 1.0;

klwrbnd = (((rst + delta)/(alpha))/(lagg^(1.0 - alpha)))^(1.0 - alpha);
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
V = zeros(nl, na, np);

for il = 1:nl
    l = lgrid(il);
    for ia = 1:na
        a = agrid(ia);

        max_cons = wage*l + (1+r)*a - r*phi;

        for ip = 1:np

            tau = tgrid(ip);
            
%             max_cons = gov.tax(max_cons, lambda, tau);

            if sigma == 1
                V(il, ia,ip) = log(max_cons)/(1-beta); 
            else
                V(il, ia,ip) = (max_cons^(1-sigma))/((1-beta)*(1-sigma));
            end
        end
    end
end

p = .5; % initial 50/50 chance of getting any party

%% solve value function: grid lookup

dist = 10;
EV = zeros(nl, na, np);
g = zeros(nl, na, np);
V1 = V;
g1 = g;
DEV = EV;
iter_ct = 1;

while dist > vTol
   
    % get expected value function 

    %TODO ADD PARTY DIFFERENCES
    for ia = 1:na
        for il = 1:nl
            for ip = 1:np
                pl = pil(il,:);
                val = p*V(:,  ia, 1) + (1-p)*V(:, ia, 2);

                EV(il, ia, ip) = pl*val;
            end
        end
    end

    %DEV approx:

    DEV = egm.numdev(EV, agrid);

    % now just make choices for assets

    for il = 1:nl
        l = lgrid(il);
        for ia = 1:na
            a = agrid(ia);
            for ip = 1:np
                tau = tgrid(ip);

                y = wage*l + (1+r)*a - r*phi;
%                 y = gov.tax(y, lambda, tau);

                params = [y beta sigma];

                [vguess, aguess, v] = compute.gss(EV(il, :, ip), params, ...
                    agrid, gTol);

                figure
                x = v(:,1);
%                 y = v(:,2);
%                 plot(x,y)
%                 hold on
%                 plot(x, v(:,3))
%                 hold off
%                 hold on
                plot(x, v(:,4))
%                 hold off
%                 legend('fd', 'vd', 'aiyagari')

                V(il, ia, ip) = vguess;
                g(il, ia, ip) = aguess;
            end
        end
    end

    dist = compute.dist(V1, V, 3);

    if mod(iter_ct, 10) == 0
        fprintf("Iteration %i: ||TV - V|| = %4.7f\n", iter_ct, dist);
    end

    iter_ct = iter_ct + 1;

    V1 = V;

    tiledlayout(3,1);
    nexttile
    mesh(EV(:,:,1))
    nexttile
    mesh(V(:,:,1))
    nexttile
    mesh(g(:,:,1))

end

tiledlayout(4,1);
nexttile
mesh(V(:,:,1))
nexttile
mesh(V(:,:,2))
nexttile
mesh(g(:,:,1))
nexttile
mesh(g(:,:,2))