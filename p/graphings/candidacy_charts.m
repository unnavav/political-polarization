restoredefaultpath;
clear all; clc;
addpath(genpath(pwd));

cd ../d

load transition_L_to_P_jan2025.mat

lvals = linspace(lgrid(1), lgrid(length(lgrid)), na);
avals = linspace(agrid(1), agrid(50), na);
K = .1;
eta = 0;

wages5 = vaas.calcw(alpha, K, eta);
intr5 = vaas.calcr(alpha, delta, K, eta);

wealth5 = avals'*intr5;
income5 = lvals*wages5;

wealth5 = repmat(wealth5, 1, na);
income5 = repmat(income5, na, 1);

total = wealth5+income5;

pct_income = income5./total*100;

pct_inc_binary = pct_income < 50;
hold on
contourf(pct_income)

% chatgpt cutoff stuff

% Define a-grid
a_vals = linspace(min(agrid), max(agrid), 250);

% Functional cutoff: l = (r / w) * a
cutoff_lvals = (intr5 / wages5) * a_vals;

% Plot the heatmap underneath
imagesc(pct_income);
axis xy;
colormap(parula);
colorbar;
title('Wage Earnings as a Percentage of Total Income');
xlabel('Wage Earnings');
ylabel('Capital Income');

% Add the functional threshold
hold on;
plot(cutoff_lvals, a_vals, 'w-', 'LineWidth', 2);  % white line, bold

% Axis beautification (optional)
xticks(linspace(1, na, 5));
yticks(linspace(1, na, 5));
xticklabels(round(linspace(min(lvals), max(lvals), 5), 2));
yticklabels(round(linspace(min(agrid), max(agrid), 5), 2));


w = vaas.calcw(alpha, K, eta);
r = vaas.calcr(alpha, delta, K, eta);

% Define e (labor supply) grid
e_vals = linspace(min(lvals), max(lvals), 250);

% Corresponding a values that satisfy 50/50 income split:
a_vals = (w / r) * e_vals;

% Plot cutoff line
plot(e_vals, a_vals, 'LineWidth', 2);
xlabel('Wage Earnings');
ylabel('Capital Income');
title('Cutoff: Wage Income = Capital Income (50/50 Split)'); 

cd ../p

%% testing ground

e_vals = lvals;                      % Wage earnings axis (x)

% Overlay the 50% cutoff line
K = 1;
eta = 0;
w = vaas.calcw(alpha, K, eta);
r = vaas.calcr(alpha, delta, K, eta);

K = 3;
eta = 0;
w2 = vaas.calcw(alpha, K, eta);
r2 = vaas.calcr(alpha, delta, K, eta);


K = 5;
eta = 0;
w3 = vaas.calcw(alpha, K, eta);
r3 = vaas.calcr(alpha, delta, K, eta);


a_cutoff = (w / r) * e_vals;         
a_cutoff2 = (w2/r2)* e_vals; 
a_cutoff3 = (w3/r3)* e_vals;

atab = [a_cutoff; a_cutoff3; a_cutoff2];
colors = {'#034C53', '#007074', '#F38C79'};

hold on;
plot(e_vals, a_cutoff, 'LineWidth', 2, 'Color', '#034C53');  % overlay smooth cutoff curve
plot(e_vals, a_cutoff2, 'LineWidth', 2, 'Color', '#007074');
plot(e_vals, a_cutoff3, 'LineWidth', 2, 'Color', '#F38C79');
xlabel('Idiosyncratic Productivity $\varepsilon$', ...
      'FontSize', 18, 'Interpreter', 'latex');
ylabel('Capital Savings ($a$)', ...
      'FontSize', 18, 'Interpreter', 'latex');
title('Wage Earnings = Capital Income', ...
      'FontSize', 22, 'Interpreter', 'latex');

% text(e_vals(end-10), e_vals(end-10), sprintf('K = %i', K), 'FontSize', 14, ...
%      'FontName', 'CMU Serif', 'Interpreter', 'latex', 'Color', colors{1});

% yticks(linspace(1, length(agrid(agrid<10)), 5));  % Pick 5 nicely spaced ticks
% % yticklabels(round(linspace(min(agrid), max(agrid(agrid<10)), 5), 2));
% xticks(linspace(1, length(e_vals(e_vals<2)), 5));
% xticklabels(round(linspace(min(e_vals), max(e_vals(e_vals<2)), 5), 2));

% text(e_vals(end-40), e_vals(end-40), 'K = 1', ...
%      'Rotation', angle, ...
%      'HorizontalAlignment', 'center', ...
%      'VerticalAlignment', 'bottom', ...
%      'FontSize', 12, ...
%      'FontWeight', 'bold', ...
%      'Color', '#034C53');

legendEntries = cell(3,1); % Preallocate legend text
legendEntries{1} = "$K = 1$";
legendEntries{2} = "$K = 3$";
legendEntries{3} = "$K = 5$";
legend(legendEntries, 'Location', 'bestoutside', 'Interpreter', 'latex');
hold off

%% increase eta

e_vals = lvals;                      % Wage earnings axis (x)

% Overlay the 50% cutoff line
K = 5;
eta = .1;
w = vaas.calcw(alpha, K, eta);
r = vaas.calcr(alpha, delta, K, eta);

eta = .3;
w2 = vaas.calcw(alpha, K, eta);
r2 = vaas.calcr(alpha, delta, K, eta);


eta = .5;
w3 = vaas.calcw(alpha, K, eta);
r3 = vaas.calcr(alpha, delta, K, eta);


a_cutoff = (w/r) * e_vals;         
a_cutoff2 = (w2/r2)* e_vals; 
a_cutoff3 = (w3/r3)* e_vals;

atab = [a_cutoff; a_cutoff3; a_cutoff2];
colors = {'#034C53', '#007074', '#F38C79'};

hold on;
plot(e_vals, a_cutoff, 'LineWidth', 2, 'Color', '#034C53');  % overlay smooth cutoff curve
plot(e_vals, a_cutoff2, 'LineWidth', 2, 'Color', '#007074');
plot(e_vals, a_cutoff3, 'LineWidth', 2, 'Color', '#F38C79');
xlabel('Idiosyncratic Productivity $\varepsilon$', ...
      'FontSize', 18, 'Interpreter', 'latex');
ylabel('Capital Savings ($a$)', ...
      'FontSize', 18, 'Interpreter', 'latex');
title('Wage Earnings = Capital Income', ...
      'FontSize', 22, 'Interpreter', 'latex');

% yticks(linspace(1, length(agrid(agrid<10)), 5));  % Pick 5 nicely spaced ticks
% yticklabels(round(linspace(min(agrid), max(agrid(agrid<10)), 5), 2));
xticks(linspace(1, length(e_vals(e_vals<5)), 5));
xticklabels(round(linspace(min(e_vals), max(e_vals(e_vals<5)), 5), 2));

legendEntries = cell(3,1); % Preallocate legend text
legendEntries{1} = "$\eta = .1$";
legendEntries{2} = "$\eta = .3$";
legendEntries{3} = "$\eta = .5$";
legend(legendEntries, 'Location', 'bestoutside', 'Interpreter', 'latex');
hold off

%% tau

e_vals = lvals;                      % Wage earnings axis (x)
tau = .5;
lamval = .7;

% Overlay the 50% cutoff line
eta = 0;
w = vaas.calcw(alpha, K, eta);
r = vaas.calcr(alpha, delta, K, eta)*(1-.15);

w2 = vaas.calcw(alpha, K, eta);
r2 = vaas.calcr(alpha, delta, K, eta)*(1-.15);


w3 = vaas.calcw(alpha, K, eta);
r3 = vaas.calcr(alpha, delta, K, eta)*(1-.15);

tau = .1;
a_cutoff = ((w*e_vals - gov.tax(w*e_vals, lamval, tau)) / r);   

tau = .3;
a_cutoff2 = ((w2*e_vals - gov.tax(w2*e_vals, lamval, tau)) / r2); 

tau = .5;
a_cutoff3 = ((w3*e_vals - gov.tax(w3*e_vals, lamval, tau)) / r3);

colors = {'#034C53', '#007074', '#F38C79'};

hold on;
plot(a_cutoff, e_vals, 'LineWidth', 2, 'Color', '#034C53');  % overlay smooth cutoff curve
plot(a_cutoff2, e_vals, 'LineWidth', 2, 'Color', '#007074');
plot(a_cutoff3, e_vals, 'LineWidth', 2, 'Color', '#F38C79');
xlabel('Idiosyncratic Productivity $\varepsilon$', ...
      'FontSize', 18, 'Interpreter', 'latex');
ylabel('Capital Savings ($a$)', ...
      'FontSize', 18, 'Interpreter', 'latex');
title('Wage Earnings = Capital Income', ...
      'FontSize', 22, 'Interpreter', 'latex');

% yticks(linspace(1, length(agrid(agrid<10)), 5));  % Pick 5 nicely spaced ticks
% yticklabels(round(linspace(min(agrid), max(agrid(agrid<10)), 5), 2));
xticks(linspace(1, length(e_vals(e_vals<2)), 5));
xticklabels(round(linspace(min(e_vals), max(e_vals(e_vals<2)), 5), 2));

legendEntries = cell(3,1); % Preallocate legend text
legendEntries{1} = "$\tau = .1$";
legendEntries{2} = "$\tau = .3$";
legendEntries{3} = "$\tau = .5$";
legend(legendEntries, 'Location', 'bestoutside', 'Interpreter', 'latex');
hold off

%% eta + tau

e_vals = lvals;                      % Wage earnings axis (x)
lamval = .7;

% Overlay the 50% cutoff line
eta = 0;
w = vaas.calcw(alpha, K, eta);
r = vaas.calcr(alpha, delta, K, eta)*(1-.15);

eta = .5;
w2 = vaas.calcw(alpha, K, eta);
r2 = vaas.calcr(alpha, delta, K, eta)*(1-.15);


w3 = vaas.calcw(alpha, K, eta);
r3 = vaas.calcr(alpha, delta, K, eta)*(1-.15);

tau = .5;
a_cutoff = ((w*e_vals - gov.tax(w*e_vals, lamval, tau)) / r);   

tau = .1;
a_cutoff2 = ((w2*e_vals - gov.tax(w2*e_vals, lamval, tau)) / r2); 

colors = {'#034C53', '#007074', '#F38C79'};

hold on;
plot(a_cutoff, e_vals, 'LineWidth', 2, 'Color', '#007074');  % overlay smooth cutoff curve
plot(a_cutoff2, e_vals, 'LineWidth', 2, 'Color', '#F38C79');
% plot(a_cutoff3, e_vals, 'LineWidth', 2, 'Color', '#F38C79');
xlabel('Idiosyncratic Productivity $\varepsilon$', ...
      'FontSize', 18, 'Interpreter', 'latex');
ylabel('Capital Savings ($a$)', ...
      'FontSize', 18, 'Interpreter', 'latex');
title('Wage Earnings = Capital Income', ...
      'FontSize', 22, 'Interpreter', 'latex');

% yticks(linspace(1, length(agrid(agrid<10)), 5));  % Pick 5 nicely spaced ticks
% yticklabels(round(linspace(min(agrid), max(agrid(agrid<10)), 5), 2));
xticks(linspace(1, length(e_vals(e_vals<2)), 5));
xticklabels(round(linspace(min(e_vals), max(e_vals(e_vals<2)), 5), 2));

legendEntries = cell(2,1); % Preallocate legend text
legendEntries{1} = "$\tau = .5, \eta = 0$";
legendEntries{2} = "$\tau = .1, \eta = .5$";
% legendEntries{3} = "$\tau = .5$";
legend(legendEntries, 'Location', 'bestoutside', 'Interpreter', 'latex');
hold off

%% dems and republicans

e_vals = lvals;                      % Wage earnings axis (x)
lamval = .7;

% Overlay the 50% cutoff line
eta = .5;
w = vaas.calcw(alpha, K, eta);
r = vaas.calcr(alpha, delta, K, eta)*(1-.15);

eta = .0;
w2 = vaas.calcw(alpha, K, eta);
r2 = vaas.calcr(alpha, delta, K, eta)*(1-.15);


w3 = vaas.calcw(alpha, K, eta);
r3 = vaas.calcr(alpha, delta, K, eta)*(1-.15);

tau = .5;
a_cutoff = ((w*e_vals - gov.tax(w*e_vals, lamval, tau)) / r);   

tau = .1;
a_cutoff2 = ((w2*e_vals - gov.tax(w2*e_vals, lamval, tau)) / r2); 

colors = {'#034C53', '#007074', '#F38C79'};

hold on;
plot(a_cutoff, e_vals, 'LineWidth', 2, 'Color', '#007074');  % overlay smooth cutoff curve
plot(a_cutoff2, e_vals, 'LineWidth', 2, 'Color', '#F38C79');
% plot(a_cutoff3, e_vals, 'LineWidth', 2, 'Color', '#F38C79');
xlabel('Idiosyncratic Productivity $\varepsilon$', ...
      'FontSize', 18, 'Interpreter', 'latex');
ylabel('Capital Savings ($a$)', ...
      'FontSize', 18, 'Interpreter', 'latex');
title('Wage Earnings = Capital Income', ...
      'FontSize', 22, 'Interpreter', 'latex');

% yticks(linspace(1, length(agrid(agrid<10)), 5));  % Pick 5 nicely spaced ticks
% yticklabels(round(linspace(min(agrid), max(agrid(agrid<10)), 5), 2));
xticks(linspace(1, length(e_vals(e_vals<2)), 5));
xticklabels(round(linspace(min(e_vals), max(e_vals(e_vals<2)), 5), 2));

legendEntries = cell(2,1); % Preallocate legend text
legendEntries{1} = "$\tau = .5, \eta = .5$";
legendEntries{2} = "$\tau = .1, \eta = .0$";
% legendEntries{3} = "$\tau = .5$";
legend(legendEntries, 'Location', 'bestoutside', 'Interpreter', 'latex');
hold off

%% 
load KS_migration_results.mat
VLL = Varray{1};
VLH = Varray{2};
VPL = Varray{3};
VPH = Varray{4};

EVLL = Varray{1};
EVLH = Varray{2};
EVPL = Varray{3};
EVPH = Varray{4};


deltaV = squeeze(EVPH(25,:,:) - EVLH(25,:,:));
deltaV = deltaV>0;
plotmap = deltaV(:,1:100);
hold on;
contourf(plotmap')
xlabel('Labor Productivity (ε)');
ylabel('Assets (a)');
title('Votes Cast');
yticklabels([0 round(agrid(10:10:100),2)])
xticklabels(round(lgrid,2))
hold off;

figure;
hold on;
for ia = [1, round(na/2), na]
    plot(lgrid, deltaV(:,ia), 'DisplayName', sprintf('a = %.2f', agrid(ia)));
end
legend show;
xlabel('Labor Productivity (ε)');
ylabel('ΔV = V_{pop} - V_{lib}');
title('Voting Differential by ε');


alpha = 0.36;
k = 5:0.5:10;  % Sample K values

etal = 0.25;
etap = 0;

w_lib = vaas.calcw(alpha, k, etal);  % Higher migration
w_pop = vaas.calcw(alpha, k, etap);  % Lower migration

plot(k, w_pop, 'r-', 'DisplayName', 'Populist Regime');
hold on;
plot(k, w_lib, 'b--', 'DisplayName', 'Liberalist Regime');
legend;
xlabel('Capital (K)');
ylabel('Wage (w)');
title('Wage Comparison Across Regimes');
