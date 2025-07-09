cd C:\VAASAVI\Dropbox\Education\OSU\Ongoing_Research\Populism\political-polarization\p
restoredefaultpath;
clear all; clc;
addpath(genpath(pwd));cd ../d

%% transitions
load transition_P_to_L_jan2025

T = length(wt);
lag = 10;
cutoff = 50;

wt = aiyagari.getW(rt, alpha, delta);
impliedK = (rt + delta)/alpha;
impliedK = impliedK.^(1/(alpha-1));
impliedK = impliedK;
k0 = (r0 + delta)/alpha;
k0 = k0.^(1/(alpha-1));
k0 = k0.*(lt(1));

k1 = (r1 + delta)/alpha;
k1 = k1.^(1/(alpha-1));

fig = figure;
set(fig, 'Units', 'inches', 'Position', [1, 1, 8.5, 4.8]);
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [8.5, 4.8], ...
         'PaperPosition', [0, 0, 8.5, 4.8]);
tiledlayout(1,3)

nexttile
hold on
plot([-10:0 1:cutoff-1], [repelem(r0, 10) rt(1:cutoff)], "LineWidth", 1.0, ...
    "Color",[0.141 0.243 0.212])
title("$r_t$",'Interpreter', 'latex')
xlabel("Time Period",'Interpreter', 'latex')
xline(1, '--', 'Shock', 'Color',[0.5 0.5 0.5],'Interpreter', 'latex');
yline(r0,'-','Liberalism', 'Interpreter', 'latex')
yline(rt(cutoff),'-','Populism', 'Interpreter', 'latex')
hold off

nexttile
hold on
plot([-10:0 1:cutoff-1], [repelem(wt(1), 10) wt(1:cutoff)], "LineWidth", 1.0, ...
    "Color",[0.486 0.663 0.51])
title("$w_t$", 'Interpreter', 'latex')
xlabel("Time Period",'Interpreter', 'latex')
xline(1, '--', 'Shock', 'Color',[0.5 0.5 0.5],'Interpreter', 'latex');
yline(wt(1),'-','Liberalism', 'Interpreter', 'latex')
yline(wt(cutoff),'-','Populism', 'Interpreter', 'latex')
hold off


nexttile
hold on
plot([-10:0 1:cutoff-1], [repelem(k0, 12) Kguess(2:cutoff-1)], "LineWidth", 1.0, ...
    "Color",[0.761 0.659 0.243])
title("$K_t$",'Interpreter', 'latex')
xlabel("Time Period",'Interpreter', 'latex')
xline(1, '--', 'Shock', 'Color',[0.5 0.5 0.5],'Interpreter', 'latex');
yline(k0,'-','Liberalism', 'Interpreter', 'latex')
yline(impliedK(cutoff),'-','Populism', 'Interpreter', 'latex')  
hold off

load transition_P_to_L_jan2025

T = length(wt);
lag = 10;
cutoff = 50;
impliedK = (rt + delta)/alpha;
impliedK = impliedK.^(1/(alpha-1));
impliedK = impliedK.*lt;
k0 = (r0 + delta)/alpha;
k0 = k0.^(1/(alpha-1));
k0 = k0.*(lt(1));
wt = aiyagari.getW(rt, alpha, delta);

fig = figure;
set(fig, 'Units', 'inches', 'Position', [1, 1, 8.5, 4.8]);
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [8.5, 4.8], ...
         'PaperPosition', [0, 0, 8.5, 4.8]);
tiledlayout(1,3)

nexttile
hold on
plot([-10:0 1:cutoff-1], [repelem(r0, 10) rt(1:cutoff)], "LineWidth", 1.0, ...
    "Color",[0.141 0.243 0.212])
title("$r_t$",'Interpreter', 'latex')
xlabel("Time Period",'Interpreter', 'latex')
xline(1, '--', 'Shock', 'Color',[0.5 0.5 0.5],'Interpreter', 'latex');
yline(r0,'-','Populism', 'Interpreter', 'latex')
yline(rt(cutoff),'-','Liberalism', 'Interpreter', 'latex')
hold off

nexttile
hold on
plot([-10:0 1:cutoff-1], [repelem(wt(1), 10) wt(1:cutoff)], "LineWidth", 1.0, ...
    "Color",[0.486 0.663 0.51])
title("$w_t$", 'Interpreter', 'latex')
xlabel("Time Period",'Interpreter', 'latex')
xline(1, '--', 'Shock', 'Color',[0.5 0.5 0.5],'Interpreter', 'latex');
yline(wt(1),'-','Populism', 'Interpreter', 'latex')
yline(wt(cutoff),'-','Liberalism', 'Interpreter', 'latex')
hold off


nexttile
hold on
plot([-10:0 1:cutoff-1], [repelem(k0, 12) Kguess(2:cutoff-1)], "LineWidth", 1.0, ...
    "Color",[0.761 0.659 0.243])
title("$K_t$",'Interpreter', 'latex')
xlabel("Time Period",'Interpreter', 'latex')
xline(1, '--', 'Shock', 'Color',[0.5 0.5 0.5],'Interpreter', 'latex');
yline(k0,'-','Populism', 'Interpreter', 'latex')
yline(impliedK(cutoff),'-','Liberalism', 'Interpreter', 'latex')  
hold off
%% savings distributions

xticks = [0 amu(100:100:1000)];
xticks = round(xticks, 0);

tiledlayout(3,1)
% t = 0 (old world)

nexttile
distr = Farray{1};
distr = sum(distr,1);
distr = distr(1:1000);

hold on
plot(distr, "LineWidth", 1.0, ...
    "Color",[0.141 0.243 0.212])
title("Asset Distribution Before Shock (t=0)")
xlabel("Assets")
xticklabels(xticks)
hold off

nexttile
distr = Farray{2};
distr = sum(distr,1);
distr = distr(1:1000);

hold on
plot(distr, "LineWidth", 1.0, ...
    "Color",[0.141 0.243 0.212])
title("Asset Distribution at Shock (t=1)")
xlabel("Assets")
xticklabels(xticks)
hold off

nexttile
distr = Farray{T};
distr = sum(distr,1);
distr = distr(1:1000);

hold on
plot(distr, "LineWidth", 1.0, ...
    "Color",[0.141 0.243 0.212])
title("Asset Distribution In New Steady State (t=T)")
xlabel("Assets")
xticklabels(xticks)
hold off

%% savings rules

tiledlayout(1,2)

nexttile
rule = Garray{1};
nmu = length(amu);
rule2500 = zeros(nl, nmu);

for im = 1:nmu
    kval = amu(im);
    for il = 1:nl  
        [ix, we] = compute.weight(agrid, kval);
            
        %split between rep and dem capital choices
        kdval = rule(il,ix)*we + rule(il,ix+1)*(1.0 - we);

        rule2500(il, im) = kdval;
    end
end

rule = rule2500([1 4 7], :);
rule = rule(:,1:25:2500)';

hold on
plot(rule, "LineWidth", 1.0)
title("Savings Decision Rule Before Shock (t=0)")
xlabel("Assets")
plot(1:100, 1:100, '--', "Color", [0 0 0])
ylim([0 100])
hold off

nexttile
rule = Garray{T};
rule2500 = zeros(nl, nmu);

for im = 1:nmu
    kval = amu(im);
    for il = 1:nl  
        [ix, we] = compute.weight(agrid, kval);
            
        %split between rep and dem capital choices
        kdval = rule(il,ix)*we + rule(il,ix+1)*(1.0 - we);

        rule2500(il, im) = kdval;
    end
end

rule = rule2500([1 4 7], :);
rule = rule(:,1:25:2500)';

hold on
plot(rule, "LineWidth", 1.0)
title("Savings Decision Rule In New Steady State (t=T)")
xlabel("Assets")
plot(1:100, 1:100, '--', "Color", [0 0 0])
ylim([0 100])
hold off

hL = legend('\epsilon = 0.25', '\epsilon = 1.00', '\epsilon = 3.94', ...
    ['45' char(176) 'line']);
hL.Layout.Tile = 'East';

%% a' - a rules

tiledlayout(1,3)

nexttile
rule = Garray{1};
nmu = length(amu);
rule2500 = zeros(nl, nmu);

for im = 1:nmu
    kval = amu(im);
    for il = 1:nl  
        [ix, we] = compute.weight(agrid, kval);
            
        %split between rep and dem capital choices
        kdval = rule(il,ix)*we + rule(il,ix+1)*(1.0 - we) - kval;

        rule2500(il, im) = kdval;
    end
end

rule = rule2500([1 4 7], :);
rule = rule(:,1:25:2500)';

hold on
plot(rule, "LineWidth", 1.0)
title("Investment Decision Rule Before Shock (t=0)")
xlabel("a")
ylabel("a' - a")
plot(0:4, 0:4, '--', "Color", [0 0 0])
hold off

nexttile
rule = Garray{2};
nmu = length(amu);
rule2500 = zeros(nl, nmu);

for im = 1:nmu
    kval = amu(im);
    for il = 1:nl  
        [ix, we] = compute.weight(agrid, kval);
            
        %split between rep and dem capital choices
        kdval = rule(il,ix)*we + rule(il,ix+1)*(1.0 - we)-kval;

        rule2500(il, im) = kdval;
    end
end

rule = rule2500([1 4 7], :);
rule = rule(:,1:25:2500)';

hold on
plot(rule, "LineWidth", 1.0)
title("Investment Decision Rule at Shock (t=1)")
xlabel("a")
ylabel("a' - a")
plot(0:4, 0:4, '--', "Color", [0 0 0])
hold off

nexttile
rule = Garray{T};
nmu = length(amu);
rule2500 = zeros(nl, nmu);

for im = 1:nmu
    kval = amu(im);
    for il = 1:nl  
        [ix, we] = compute.weight(agrid, kval);
            
        %split between rep and dem capital choices
        kdval = rule(il,ix)*we + rule(il,ix+1)*(1.0 - we)-kval;

        rule2500(il, im) = kdval;
    end
end

rule = rule2500([1 4 7], :);
rule = rule(:,1:25:2500)';

hold on
plot(rule, "LineWidth", 1.0)
title("Investment Decision Rule in New Steady State (t=T)")
xlabel("a")
ylabel("a' - a")
plot(0:4, 0:4, '--', "Color", [0 0 0])
hold off


hL = legend('\epsilon = 0.25', '\epsilon = 1.00', '\epsilon = 3.94', ...
    ['45' char(176) 'line']);
hL.Layout.Tile = 'East';


%% 

hold on;
colormap(gray)
mesh((VPL<VP)')
ylabel('Assets $a$', 'FontSize', 16, 'Interpreter', 'latex');
xlabel('Idiosyncratic Productivity $\varepsilon$', 'FontSize', 16, 'Interpreter', 'latex');
title('Votes for Liberlism:', ...
      'FontSize', 18, 'Interpreter', 'latex');
subtitle('Vote for in Gray, Vote Against in White', ...
      'FontSize', 12, 'Interpreter', 'latex')
set(gca, 'FontSize', 14);
yticklabels(round(agrid(1:50:250),2))
xticklabels(round(lgrid,2))
hold off;

hold on;
colormap(gray)
mesh((VLP<VL)')
ylabel('Assets $a$', 'FontSize', 16, 'Interpreter', 'latex');
xlabel('Idiosyncratic Productivity $\varepsilon$', 'FontSize', 16, 'Interpreter', 'latex');
title('Votes for Populism:', ...
      'FontSize', 18, 'Interpreter', 'latex');
subtitle('Vote for in Gray, Vote Against in White', ...
      'FontSize', 12, 'Interpreter', 'latex')
set(gca, 'FontSize', 14);
yticklabels(round(agrid(1:50:250),2))
xticklabels(round(lgrid,2))
hold off;

hold on;
colormap(gray(3))
mesh((VPL>VP)' + (VLP>VL)')
ylabel('Assets $a$', 'FontSize', 16, 'Interpreter', 'latex');
xlabel('Idiosyncratic Productivity $\varepsilon$', 'FontSize', 16, 'Interpreter', 'latex');
title('A Middle Class', ...
      'FontSize', 18, 'Interpreter', 'latex');
set(gca, 'FontSize', 14);
yticklabels(round(agrid(1:50:250),2))
xticklabels(round(lgrid,2))
hold off;



