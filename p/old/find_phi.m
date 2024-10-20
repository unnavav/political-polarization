%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phi triangulation
% june 2024
% vaasavi
%
% Given that the share of votes has to be equal to .54, we need to have a
% phi value that's not high enough to have people only vote for their own
% party. So this iterates to a bonus that hits that number, based on what
% we know about decision rules. In particular, we know that if we expand
% the lump sum remittance while decreasing the bonus of having one's own
% party in power, eventually the lowest B voters should switch over to A.
% I attempt to iterate to that answer here.

restoredefaultpath;
clear all; clc;
addpath(genpath(pwd));

load ..\d\results-20240703\results_changedy_flatcaptax_EGM_phi32604-Jul-2024.mat

dTol = 15e-4; %it's just being weird on the .001 level, so raising it slightly

uh = 10; ul = 0;

ubonus = (uh+ul)/2;

ubonusA = [ubonus 0];
ubonusB = [0 ubonus];

goal = 0.540321873; % average of post war elections results

voteDist = 10;
iter_ct = 1;

while voteDist > dTol

        fprintf("\n\nGuess %i\nSolving value functions for phi level %1.8f:\n", ...
            iter_ct, ubonus)

        fprintf("\tHH Type A:\n")
        terms.ubonus = ubonusA;
        [Va, Ga, EVa] = HH.solve(nl, na, np, terms, vTol);
    
        fprintf("\tHH Type B:\n")
        terms.ubonus = ubonusB;
        [Vb, Gb, EVb] = HH.solve(nl, na, np, terms, vTol);
        
        VOTESa = (EVa(:,:,1) > EVa(:,:,2));
        VOTESb = (EVb(:,:,1) > EVb(:,:,2));
    
        test = EVa(:,:,1) - EVb(:,:,1);
    
        %TO DO: 
        % SOLVE ASSET DISTRIBUTION 
        % SOLVE VOTE COUNTS

        fprintf("\tSolving asset distribution:\n")
        [adistrA, kaggA, adistrB, kaggB] = HH.getDist(Ga, Gb, amu, agrid, pil, pctA);
    
        % CONDENSE DISTR for asset 
        acondA = compute.condense(adistrA, amu, agrid);
        acondB = compute.condense(adistrB, amu, agrid);
    
        share2A = pctA*sum(sum(VOTESa.*acondA)) + ...
            (1-pctA)*sum(sum(VOTESb.*acondA));
        share2B = pctA*sum(sum(VOTESa.*acondB)) + ...
            (1-pctA)*sum(sum(VOTESb.*acondB));

        diff = share2A - goal;
        voteDist = abs(diff);

        if diff > 0
            fprintf("\nToo many votes for A: %1.5f percent off mark.\n", ...
                voteDist)
            ul = (ubonus+2*ul)/3;
            fprintf("Increase phi guess.")
        else
            fprintf("\nToo few votes for A: %1.5f percent off mark.\n", ...
                voteDist)
            uh = (ubonus+2*uh)/3;
            fprintf("Derease phi guess.")
        end

        ubonus = (uh+ul)/2;

        ubonusA = [ubonus 0];
        ubonusB = [0 ubonus];

        iter_ct = iter_ct + 1;
end

[xgrid, ygrid] = meshgrid(lgrid,agrid);
lrep = repmat(lgrid,1,na)';
arep = repmat(agrid,nl,1);
arep = arep(:);

dat = VOTESa;
dat = dat(:);
G1dat = [lrep arep dat];
tiledlayout(1,1)
graphdata = griddata(G1dat(:,1),G1dat(:,2), G1dat(:,3), xgrid, ygrid);
ax1 = nexttile;
contourf(graphdata);
yticks(1:50:250)
yticklabels([agrid(1:50:250)])
ylabel('Assets', 'FontSize', 14)
xticks(1:2:7)
xticklabels(lgrid(1:2:7))
xlabel('Labor Productivity','FontSize', 14)
colormap(ax1,winter)
