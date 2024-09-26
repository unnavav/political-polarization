
restoredefaultpath;
clear all; clc;

%cd "C:\Users\unnav\Dropbox\Education\OSU\Ongoing_Research\2nd Year Paper\political-polarization\p\graphings"

addpath(genpath(pwd));

set(gca,'fontname','century')
cd ../d

load agrid
load lgrid


legend_dat = round(lgrid([1 4 7]),2);
leg = string(legend_dat);

cd ../d/tax_scheme_0.00000.9500/

load results_pval10.mat

cd ..\
cd ../v/

amu = linspace(min(agrid), max(agrid),max(size(adistrA)));



hold on
set(gca,'fontname','century')
loglog(Ga(1,:,1),agrid,'LineWidth',2.0)
loglog(Ga(4,:,1), agrid,'LineWidth',2.0)
loglog(Ga(7,:,1), agrid,'LineWidth',2.0)
xlabel('Assets')
ylabel('Asset Choice')
title({'Decision Rule for Household A', 'in Steady State'})
subtitle(' ')
ylim([0 100])
xlim([0 100])
ylabel('Asset Choice', 'FontSize', 12)
xlabel('Current Assets', 'FontSize', 12)
hold off
legend(leg)
exportgraphics(gcf,'A_t95.png','Resolution',300)

hold on
set(gca,'fontname','century')
loglog(Gb(1,:,1),agrid,'LineWidth',2.0)
loglog(Gb(4,:,1), agrid,'LineWidth',2.0)
loglog(Gb(7,:,1), agrid,'LineWidth',2.0)
xlabel('Assets')
ylabel('Asset Choice')
title({'Decision Rule for Household B', 'in Steady State'})
subtitle(' ')
ylim([0 100])
xlim([0 100])
ylabel('Asset Choice', 'FontSize', 12)
xlabel('Current Assets', 'FontSize', 12)
hold off
legend(leg)
exportgraphics(gcf,'B_t95.png','Resolution',300)

diffsA = Ga(:,:,2) - Gb(:,:,2);

% t = tiledlayout(1,1,'Padding','tight');
% t.Units = 'inches';
% t.OuterPosition = [0 0 4 4];
% t.InnerPosition = [1 1 3 3];
% nexttile;
hold on
set(gca,'fontname','century')
semilogx(agrid, diffsA(1,:),'LineWidth',2.0)
semilogx(agrid, diffsA(4,:),'LineWidth',2.0)
semilogx(agrid, diffsA(7,:),'LineWidth',2.0)
xlabel('Assets')
ylabel('Asset Choice')
title({'Difference in Choices Between Party A', 'and Party B Households'})
subtitle(' ')
xlim([0 100])
ylabel('Difference in Asset Choice', 'FontSize', 12)
xlabel('Assets', 'FontSize', 12)
hold off
legend(leg)
exportgraphics(gcf,'AminusB_t95.png','Resolution',300)

hold on
set(gca,'fontname','century')
semilogx(amu, adistrB(1,:),'LineWidth',2.0)
semilogx(amu, adistrB(4,:),'LineWidth',2.0)
semilogx(amu, adistrB(7,:),'LineWidth',2.0)
xlim([0 100])
xlabel('Assets')
ylabel('Asset Choice')
title({'Distribution of Asset Choices in Steady State'})
subtitle(' ')
% yticks(0:.005:.08)
% yticklabels(0:.005:.08)
ylabel('Probability', 'FontSize', 12)
xlabel('Assets', 'FontSize', 12)
hold off
legend(leg)
exportgraphics(gcf,'Adistr_t95.png','Resolution',300)


cd ../d/tax_scheme_0.00000.9000/

load results_pval10.mat

cd ..
cd ../v/

hold on
set(gca,'fontname','century')
loglog(Ga(1,:,1),agrid,'LineWidth',2.0)
loglog(Ga(4,:,1), agrid,'LineWidth',2.0)
loglog(Ga(7,:,1), agrid,'LineWidth',2.0)
xlabel('Assets')
ylabel('Asset Choice')
title({'Decision Rule for Household A', 'in Steady State'})
subtitle(' ')
ylim([0 100])
xlim([0 100])
ylabel('Asset Choice', 'FontSize', 12)
xlabel('Current Assets', 'FontSize', 12)
hold off
legend(leg)
exportgraphics(gcf,'A_t90.png','Resolution',300)

diffsA = Ga(:,:,1) - Gb(:,:,1);

hold on
set(gca,'fontname','century')
plot(agrid, diffsA(1,:),'LineWidth',2.0)
plot(agrid, diffsA(4,:),'LineWidth',2.0)
plot(agrid, diffsA(7,:),'LineWidth',2.0)
xlabel('Assets')
ylabel('Asset Choice')
title({'Difference in Choices Between Party A', 'and Party B Households'})
subtitle(' ')
yticks(0:2e-7:1.2e-6)
yticklabels(0:2e-7:1.2e-6)
ylabel('Difference in Asset Choice', 'FontSize', 12)
xticks(1:15:250)
xticklabels([agrid(1:15:250)])
xlabel('Assets', 'FontSize', 12)
hold off
legend(leg)
exportgraphics(gcf,'AminusB_t90.png','Resolution',300)

hold on
set(gca,'fontname','century')
semilogx(amu, adistrA(1,:),'LineWidth',2.0)
semilogx(amu, adistrA(4,:),'LineWidth',2.0)
semilogx(amu, adistrA(7,:),'LineWidth',2.0)
ylim([0 4e-3])
xlim([0 100])
xlabel('Probability')
ylabel('Asset Choice')
title({'Distribution of Asset Choices in Steady State'})
subtitle(' ')
ylabel('Probability', 'FontSize', 12)
xlabel('Assets', 'FontSize', 12)
hold off
legend(leg)
exportgraphics(gcf,'Adistr_t90.png','Resolution',300)

% voting results
custom = [0.780000000000000	0.365000000000000	0.671000000000000; ...
    0.731241830065360	0.402875816993464	0.701891067538127; ...
    0.682483660130719	0.440751633986928	0.732782135076253; ...
    0.633725490196079	0.478627450980392	0.763673202614379; ...
    0.584967320261438	0.516503267973856	0.794564270152505; ...
    0.536209150326797	0.554379084967320	0.825455337690632; ...
    0.487450980392157	0.592254901960784	0.856346405228758; ...
    0.438692810457516	0.630130718954248	0.887237472766885; ...
    0.389934640522876	0.668006535947712	0.918128540305011; ...
    0.341176470588235	0.705882352941177	0.949019607843137];

tiledlayout(1,2)

t1 = nexttile;
hold on
set(gca,'fontname','century')
contourf(transpose(VOTESa),'edgecolor','none')
title("Household A's Voting Rule")
subtitle(' ')
yticks([50 100 150 200 250])
yticklabels(round(agrid([50 100 150 200 250]), 2))
ylabel('Household Wealth')
xticklabels(round(lgrid, 2))
xlabel('Idiosyncratic Productivity (ε)')
% set(gca, 'YScale', 'log')
colormap(t1, custom)

t2 = nexttile;
hold on
set(gca,'fontname','century')
contourf(transpose(VOTESb),'edgecolor','none')
title("Household B's Voting Rule")
subtitle(' ')
yticks([50 100 150 200 250])
yticklabels(round(agrid([50 100 150 200 250]), 2))
ylabel('Household Wealth')
xticklabels(round(lgrid, 2))
xlabel('Idiosyncratic Productivity (ε)')
% set(gca, 'YScale', 'log')
colormap(t2, custom)

cb = colorbar;
cb.Layout.Tile = 'east';


cd ../d/personal_return0.25/

load results.mat

cd ..

cd ../v/

hold on
set(gca,'fontname','century')
loglog(Ga(1,:,1),agrid,'LineWidth',2.0)
loglog(Ga(4,:,1), agrid,'LineWidth',2.0)
loglog(Ga(7,:,1), agrid,'LineWidth',2.0)
xlabel('Assets')
ylabel('Asset Choice')
title({'Decision Rule for Household A', 'in Steady State'})
subtitle(' ')
ylim([0 100])
xlim([0 100])
ylabel('Asset Choice', 'FontSize', 12)
xlabel('Current Assets', 'FontSize', 12)
hold off
legend(leg)
exportgraphics(gcf,'A_u25.png','Resolution',300)

diffsA = Ga(:,:,1) - Gb(:,:,1);

hold on
set(gca,'fontname','century')
plot(agrid, diffsA(1,:),'LineWidth',2.0)
plot(agrid, diffsA(4,:),'LineWidth',2.0)
plot(agrid, diffsA(7,:),'LineWidth',2.0)
xlabel('Assets')
ylabel('Asset Choice')
title({'Difference in Choices Between Party A', 'and Party B Households'})
subtitle(' ')
ylabel('Difference in Asset Choice', 'FontSize', 12)
xlabel('Assets', 'FontSize', 12)
hold off
legend(leg)
exportgraphics(gcf,'AminusB_u25.png','Resolution',300)

hold on
set(gca,'fontname','century')
semilogx(amu, adistrA(1,:),'LineWidth',2.0)
semilogx(amu, adistrA(4,:),'LineWidth',2.0)
semilogx(amu, adistrA(7,:),'LineWidth',2.0)
ylim([0 4e-3])
xlim([0 100])
xlabel('Probability')
ylabel('Asset Choice')
title({'Distribution of Asset Choices in Steady State'})
subtitle(' ')
ylabel('Probability', 'FontSize', 12)
xlabel('Assets', 'FontSize', 12)
hold off
legend(leg)
exportgraphics(gcf,'Adistr_u25.png','Resolution',300)


cd ../d/personal_return0.05/

load results.mat

cd ..
cd ../v/

hold on
set(gca,'fontname','century')
loglog(Ga(1,:,1),agrid,'LineWidth',2.0)
loglog(Ga(4,:,1), agrid,'LineWidth',2.0)
loglog(Ga(7,:,1), agrid,'LineWidth',2.0)
xlabel('Assets')
ylabel('Asset Choice')
title({'Decision Rule for Household A', 'in Steady State'})
subtitle(' ')
ylim([0 100])
xlim([0 100])
ylabel('Asset Choice', 'FontSize', 12)
xlabel('Current Assets', 'FontSize', 12)
hold off
legend(leg)
exportgraphics(gcf,'A_u05.png','Resolution',300)


diffsA = Ga(:,:,1) - Gb(:,:,1);

hold on
set(gca,'fontname','century')
plot(agrid, diffsA(1,:),'LineWidth',2.0)
plot(agrid, diffsA(4,:),'LineWidth',2.0)
plot(agrid, diffsA(7,:),'LineWidth',2.0)
xlabel('Assets')
ylabel('Asset Choice')
title({'Difference in Choices Between Party A', 'and Party B Households'})
subtitle(' ')
yticks(0:5e-7:3e-6)
yticklabels(0:5e-7:3e-6)
ylabel('Difference in Asset Choice', 'FontSize', 12)
xticks(1:15:250)
xticklabels([agrid(1:15:250)])
xlabel('Assets', 'FontSize', 12)
hold off
legend(leg)
exportgraphics(gcf,'AminusB_u05.png','Resolution',300)

hold on
set(gca,'fontname','century')
set(gca,'fontname','century')
semilogx(amu, adistrA(1,:),'LineWidth',2.0)
semilogx(amu, adistrA(4,:),'LineWidth',2.0)
semilogx(amu, adistrA(7,:),'LineWidth',2.0)
ylim([0 4e-3])
xlim([0 100])
xlabel('Probability')
ylabel('Asset Choice')
title({'Distribution of Asset Choices in Steady State'})
subtitle(' ')
ylabel('Probability', 'FontSize', 12)
xlabel('Assets', 'FontSize', 12)
hold off
legend(leg)
exportgraphics(gcf,'Adistr_u05.png','Resolution',300)


cd ../d/personal_return0.15/

load results.mat

tiledlayout(1,2)

t1 = nexttile;
hold on
set(gca,'fontname','century')
contourf(transpose(VOTESa),'edgecolor','none')
title("Household A's Voting Rule")
subtitle(' ')
yticks([50 100 150 200 250])
yticklabels(round(agrid([50 100 150 200 250]), 2))
ylabel('Household Wealth')
xticklabels(round(lgrid, 2))
xlabel('Idiosyncratic Productivity (ε)')
% set(gca, 'YScale', 'log')
colormap(t1, custom)

tiledlayout(1,1)
t2 = nexttile;
hold on
set(gca,'fontname','century')
contourf(transpose(VOTESb),'edgecolor','none')
title("Household B's Voting Rule")
subtitle(' ')
yticks([50 100 150 200 250])
yticklabels(round(agrid([50 100 150 200 250]), 2))
ylabel('Household Wealth')
xticklabels(round(lgrid, 2))
xlabel('Idiosyncratic Productivity (ε)')
% set(gca, 'YScale', 'log')
colormap(t2, custom)

cb = colorbar;
cb.Layout.Tile = 'east';


cd ../d/partial_eqm_test/

load results_0.00

tiledlayout(3,2)

t1 = nexttile;
hold on
set(gca,'fontname','century')
contourf(transpose(VOTESa),'edgecolor','none')
title("Household Dem's Voting Rule")
subtitle(' ')
yticks([50 100 150 200 250])
yticklabels(round(agrid([50 100 150 200 250]), 2))
ylabel('Household Wealth')
xticklabels(round(lgrid, 2))
xlabel('Idiosyncratic Productivity (ε)')
% set(gca, 'YScale', 'log')
colormap(t1, custom)

t2 = nexttile;
hold on
set(gca,'fontname','century')
contourf(transpose(VOTESb),'edgecolor','none')
title("Household Rep's Voting Rule")
subtitle(' ')
yticks([50 100 150 200 250])
yticklabels(round(agrid([50 100 150 200 250]), 2))
ylabel('Household Wealth')
xticklabels(round(lgrid, 2))
xlabel('Idiosyncratic Productivity (ε)')
% set(gca, 'YScale', 'log')
colormap(t2, custom)

load results_0.01.mat

t3 = nexttile;
hold on
set(gca,'fontname','century')
contourf(transpose(VOTESa),'edgecolor','none')
title("Household Dem's Voting Rule")
subtitle(' ')
yticks([50 100 150 200 250])
yticklabels(round(agrid([50 100 150 200 250]), 2))
ylabel('Household Wealth')
xticklabels(round(lgrid, 2))
xlabel('Idiosyncratic Productivity (ε)')
% set(gca, 'YScale', 'log')
colormap(t3, custom)

t4 = nexttile;
hold on
set(gca,'fontname','century')
contourf(transpose(VOTESb),'edgecolor','none')
title("Household Rep's Voting Rule")
subtitle(' ')
yticks([50 100 150 200 250])
yticklabels(round(agrid([50 100 150 200 250]), 2))
ylabel('Household Wealth')
xticklabels(round(lgrid, 2))
xlabel('Idiosyncratic Productivity (ε)')
% set(gca, 'YScale', 'log')
colormap(t4, custom)

load results_0.02.mat

t5 = nexttile;
hold on
set(gca,'fontname','century')
contourf(transpose(VOTESa),'edgecolor','none')
title("Household Dem's Voting Rule")
subtitle(' ')
yticks([50 100 150 200 250])
yticklabels(round(agrid([50 100 150 200 250]), 2))
ylabel('Household Wealth')
xticklabels(round(lgrid, 2))
xlabel('Idiosyncratic Productivity (ε)')
% set(gca, 'YScale', 'log')
colormap(t5, custom)

t6 = nexttile;
hold on
set(gca,'fontname','century')
contourf(transpose(VOTESb),'edgecolor','none')
title("Household Rep's Voting Rule")
subtitle(' ')
yticks([50 100 150 200 250])
yticklabels(round(agrid([50 100 150 200 250]), 2))
ylabel('Household Wealth')
xticklabels(round(lgrid, 2))
xlabel('Idiosyncratic Productivity (ε)')
% set(gca, 'YScale', 'log')
colormap(t6, custom)

cb = colorbar;
cb.Layout.Tile = 'east';


load results_0.04_noG

tiledlayout(3,2)

t1 = nexttile;
hold on
set(gca,'fontname','century')
contourf(transpose(VOTESa),'edgecolor','none')
title("Household Dem's Voting Rule")
subtitle(' ')
yticks([50 100 150 200 250])
yticklabels(round(agrid([50 100 150 200 250]), 2))
ylabel('Household Wealth')
xticklabels(round(lgrid, 2))
xlabel('Idiosyncratic Productivity (ε)')
% set(gca, 'YScale', 'log')
colormap(t1, custom)

t2 = nexttile;
hold on
set(gca,'fontname','century')
contourf(transpose(VOTESb),'edgecolor','none')
title("Household Rep's Voting Rule")
subtitle(' ')
yticks([50 100 150 200 250])
yticklabels(round(agrid([50 100 150 200 250]), 2))
ylabel('Household Wealth')
xticklabels(round(lgrid, 2))
xlabel('Idiosyncratic Productivity (ε)')
% set(gca, 'YScale', 'log')
colormap(t2, custom)

load results_0.10_noG.mat

t3 = nexttile;
hold on
set(gca,'fontname','century')
contourf(transpose(VOTESa),'edgecolor','none')
title("Household Dem's Voting Rule")
subtitle(' ')
yticks([50 100 150 200 250])
yticklabels(round(agrid([50 100 150 200 250]), 2))
ylabel('Household Wealth')
xticklabels(round(lgrid, 2))
xlabel('Idiosyncratic Productivity (ε)')
% set(gca, 'YScale', 'log')
colormap(t3, custom)

t4 = nexttile;
hold on
set(gca,'fontname','century')
contourf(transpose(VOTESb),'edgecolor','none')
title("Household Rep's Voting Rule")
subtitle(' ')
yticks([50 100 150 200 250])
yticklabels(round(agrid([50 100 150 200 250]), 2))
ylabel('Household Wealth')
xticklabels(round(lgrid, 2))
xlabel('Idiosyncratic Productivity (ε)')
% set(gca, 'YScale', 'log')
colormap(t4, custom)

load results_0.13_noG.mat

t5 = nexttile;
hold on
set(gca,'fontname','century')
contourf(transpose(VOTESa),'edgecolor','none')
title("Household Dem's Voting Rule")
subtitle(' ')
yticks([50 100 150 200 250])
yticklabels(round(agrid([50 100 150 200 250]), 2))
ylabel('Household Wealth')
xticklabels(round(lgrid, 2))
xlabel('Idiosyncratic Productivity (ε)')
% set(gca, 'YScale', 'log')
colormap(t5, custom)

t6 = nexttile;
hold on
set(gca,'fontname','century')
contourf(transpose(VOTESb),'edgecolor','none')
title("Household Rep's Voting Rule")
subtitle(' ')
yticks([50 100 150 200 250])
yticklabels(round(agrid([50 100 150 200 250]), 2))
ylabel('Household Wealth')
xticklabels(round(lgrid, 2))
xlabel('Idiosyncratic Productivity (ε)')
% set(gca, 'YScale', 'log')
colormap(t6, custom)

cb = colorbar;
cb.Layout.Tile = 'east';

