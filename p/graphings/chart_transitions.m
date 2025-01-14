cd ../d

load transition_results_jan2025

T = length(wt);
lag = 10;
cutoff = 50;

wt = aiyagari.getW(rt, alpha, delta);

%% transitions
tiledlayout(1,3)

nexttile
hold on
plot([-10:0 1:cutoff-1], [repelem(r0, 10) rt(1:cutoff)], "LineWidth", 1.0, ...
    "Color",[0.141 0.243 0.212])
title("r_t")
xlabel("Time Period")
xline(1, '--', 'Shock', 'Color',[0.5 0.5 0.5]);
yline(r0,'-','Steady State 0')
yline(r1,'-','Steady State 1')
hold off

nexttile
hold on
plot([-10:0 1:cutoff-1], [repelem(wt(1), 10) wt(1:cutoff)], "LineWidth", 1.0, ...
    "Color",[0.486 0.663 0.51])
title("w_t")
xlabel("Time Period")
xline(1, '--', 'Shock', 'Color',[0.5 0.5 0.5]);
yline(wt(1),'-','Steady State 0')
yline(wt(T),'-','Steady State 1')
hold off


nexttile
hold on
plot([-10:0 1:cutoff-1], [repelem(impliedK(1), 12) Kguess(2:cutoff-1)], "LineWidth", 1.0, ...
    "Color",[0.761 0.659 0.243])
title("K_t")
xlabel("Time Period")
xline(1, '--', 'Shock', 'Color',[0.5 0.5 0.5]);
yline(impliedK(1),'-','Steady State 0')
yline(impliedK(T),'-','Steady State 1')  
hold off

%% savings distributions

xticks = [0 amu([100:100:1000])];
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

tiledlayout(3,1)

nexttile
rule = Garray{1};
rule = rule([1 4 7], :);

hold on
plot(rule', "LineWidth", 1.0)
title("Savings Decision Rule Before Shock (t=0)")
xlabel("Assets")
xticklabels(xticks)
hold off

nexttile
rule = Garray{2};
rule = rule([1 4 7], :);

hold on
plot(rule, "LineWidth", 1.0, ...
    "Color",[0.761 0.659 0.243])
title("Savings Decision Rule at Shock (t=1)")
xlabel("Assets")
xticklabels(xticks)
hold off

nexttile
rule = Garray{T};
rule = rule([1 4 7], :);

hold on
plot(rule, "LineWidth", 1.0)
title("Savings Decision Rule In New Steady State (t=T)")
xlabel("Assets")
xticklabels(xticks)
hold off