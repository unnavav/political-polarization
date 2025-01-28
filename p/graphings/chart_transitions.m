addpath(genpath(pwd));
cd ../d

load transition_L_to_P_jan2025

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