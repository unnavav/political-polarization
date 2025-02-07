%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting the Voting Rules
% jan 2025
%
% vaasavi
% 
% this code takes the value functions from each perfect foresight
% transition to back out our voting rules

restoredefaultpath;
clear all; clc;
addpath(genpath(pwd));

cd ../d

load transition_L_to_P_jan2025.mat

Vpswitch = Varray{1};

load transition_P_to_L_jan2025.mat

Vlswitch = Varray{1};

cd liberalism_populism\

load results_t0.0000_eta1.0000.mat

Vp = V;

load results_t0.0000_eta1.1000.mat

Vl = V;

EVPswitch = V;
EVLswitch = V; 
EVp = V;
EVl = V;
for ia = 1:na
    for il = 1:nl
        EVPswitch(il, ia) = pil(il,:)*Vpswitch(:,ia);
        EVLswitch(il, ia) = pil(il,:)*Vlswitch(:,ia);
        EVp(il, ia) = pil(il,:)*Vp(:,ia);
        EVl(il, ia) = pil(il,:)*Vl(:,ia);
    end
end

mesh(Vpswitch > Vp)
mesh(Vlswitch > Vl)

%% plot voting rules true to size

tiledlayout(1,2)
%rescale to linear instead of log
nmu = length(amu);
rule2500 = zeros(nl, nmu);

votes = (EVPswitch>EVp);

for im = 1:nmu
    kval = amu(im);
    for il = 1:nl  
        [ix, we] = compute.weight(agrid, kval);
            
        %split between rep and dem capital choices
        kdval = votes(il,ix)*we + votes(il,ix+1)*(1.0 - we);

        rule2500(il, im) = kdval;
    end
end

rule = rule2500(:,1:25:2500)';


yticks = amu(250:250:2500);
nexttile
hold on
contourf(rule)
title("Switiching to a Liberalist Regime")
ylabel("Assets")
yticklabels(10:10:100)
xlabel("Productivity")
xticklabels(round(lgrid,2))
hold off


votes = (EVl>EVLswitch);

for im = 1:nmu
    kval = amu(im);
    for il = 1:nl  
        [ix, we] = compute.weight(agrid, kval);
            
        %split between rep and dem capital choices
        kdval = votes(il,ix)*we + votes(il,ix+1)*(1.0 - we);

        rule2500(il, im) = kdval;
    end
end

rule = rule2500(:,1:25:2500)';

nexttile
hold on
contourf(rule)
title("Switiching to a Populist Regime")
ylabel("Assets")
yticklabels(10:10:100)
xlabel("Productivity")
xticklabels(round(lgrid,2))
hold off

hL = legend('Switch', 'Stay');
hL.Layout.Tile = 'East';

%% plot voting rules

tiledlayout(1,2)


yticks = agrid(50:50:250);
nexttile
hold on
contourf((EVPswitch>EVp)')
colormap(gray)
title("Vote for Liberalism Under Populism")
ylabel("Assets")
yticklabels(yticks)
xlabel("Productivity")
xticklabels(round(lgrid,2))
hold off


votes = (EVLswitch>EVl)';
nexttile
hold on
contourf(votes)
colormap(gray)
title("Vote for Populism Under Liberalism")
ylabel("Assets")
yticklabels(yticks)
xlabel("Productivity")
xticklabels(round(lgrid,2))
hold off

% Create legend with reference to actual plots
ax = nexttile(2); % Attach legend to an existing tile
h1 = plot(NaN, NaN, 'ow');
h2 = plot(NaN, NaN, 'ob');
legend(ax, [h1, h2], {'Against','For'}, 'Location', 'eastoutside');

%% chatgpt sol
t = tiledlayout(1,2);

yticks = agrid(50:50:250);

poplib = (EVPswitch > EVp)';
% First subplot
ax1 = nexttile;
hold on
contourf(ax1, poplib)
colormap(gray)
title("Vote for Liberalism Under Populism")
ylabel("Assets")
yticklabels(yticks)
xlabel("Productivity")
xticklabels(round(lgrid,2))
hold off

% Second subplot
ax2 = nexttile;
hold on
votes = (EVLswitch > EVl)';
contourf(ax2, votes)
colormap(gray)
title("Vote for Populism Under Liberalism")
ylabel("Assets")
yticklabels(yticks)
xlabel("Productivity")
xticklabels(round(lgrid,2))
hold off

% Create a dummy plot for the legend (ensuring it belongs to one of the axes)
hold on;
h1 = scatter(NaN, NaN, 50, 'w', 'filled', 'MarkerEdgeColor', 'k'); % White marker for "Against"
h2 = scatter(NaN, NaN, 50, 'ok', 'filled'); % Black marker for "For"

% Attach legend to one of the axes
legend(ax2, [h1, h2], {'For', 'Against'}, 'Location', 'eastoutside');

%% some other checks

libpop = votes;
liblib = EVl > EVLswitch;
liblib = liblib';
mesh(liblib - poplib)