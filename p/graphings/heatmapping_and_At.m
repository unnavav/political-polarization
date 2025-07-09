%% Charting Heat Maps
% 

cd C:\VAASAVI\Dropbox\Education\OSU\Ongoing_Research\Populism\political-polarization\p
restoredefaultpath;
clear all; clc;
addpath(genpath(pwd));

n = 256; % Number of colors
purpleMap = [linspace(1, 0.3, n)', linspace(1, 0.3, n)', linspace(1, 0.6, n)'];

% Shifted low color: pale blue-purple (slightly more red and blue than #BBE1FA)
startColor = [200, 200, 255] / 255;  % soft periwinkle

% Shifted high color: navy with slight purple tint
endColor   = [50, 30, 100] / 255;    % midnight indigo

% Linear interpolation
bluePurpleMap = [linspace(startColor(1), endColor(1), n)', ...
                 linspace(startColor(2), endColor(2), n)', ...
                 linspace(startColor(3), endColor(3), n)'];

rawTurbo = turbo(n);

% Apply a desaturation blend toward gray
gray = ones(n,3);  % white-ish gray anchor
blend = 0.5;       % 0 = full turbo, 1 = fully gray

softTurbo = (1 - blend) * rawTurbo + blend * gray;

c = parula(n);
white = ones(n, 3);
pastel_strength = 0.25;  % 0 = original color, 1 = full white
pastelTurbo = (1 - pastel_strength) * c + pastel_strength * white;

hexColors = ["#11302a", "#036264", "#8f5774", "#e0829d", "#dac4d0"];
rgbColors = sscanf(join(hexColors,''),'#%2x%2x%2x',[3,length(hexColors)]).'/255;

% Interpolate to 256 steps
darkAcademiaChic = interp1(linspace(0,1,size(rgbColors,1)), rgbColors, linspace(0,1,n));


hexColors = ["#FC7E01", "#D9A1A7", "#A3CADB", "#70AABD", "#212121"];
rgbColors = sscanf(join(hexColors,''),'#%2x%2x%2x',[3,length(hexColors)]).'/255;
sunswept = interp1(linspace(0,1,size(rgbColors,1)), rgbColors, linspace(0,1,n));


%% the actual graphing

cd ../d/steadystates/

load results_rho85sig2_t0.9000_eta0.5684.mat

eta_grid_labels = string(round(etagrid(2:13),2));
tau_grid_labels = string(round(taugrid,2));

P = cell2mat(parray);

P = P(2:end,:);

fig = figure;
set(fig, 'Units', 'inches', 'Position', [1, 1, 4.8, 4.8]);
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [4.8, 4.8], ...
         'PaperPosition', [0, 0, 4.8, 4.8]);

[X, Y] = meshgrid(etagrid(2:13), taugrid);
pcolor(X, Y, P');
shading interp;                % for smooth color blending
clim([.45 .7]);
colormap(sunswept);
colorbar;
hold on;
[C, h] = contour(X, Y, P', [0.5 0.5], 'k', 'LineWidth', 1.5);
clabel(C, h, 'Interpreter', 'latex', 'FontSize', 16, 'Color', [0.1 0.1 0.1], ...
    'FontWeight', 'bold');
% Add a text annotation
text(0.15, 0.5, 'Populist Majority', 'Interpreter', 'latex', ...
     'Color', [0.2 0.2 0.2], 'FontSize', 12, 'Rotation', 0);

% Labels and styling
xlabel('Migration Rate $(\eta)$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Tax Progressivity $(\tau)$', 'Interpreter', 'latex', 'FontSize', 22);
title('Percentage Voting for Populists, $G = 0$', 'Interpreter', 'latex', 'FontSize', 28);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
set(gcf, 'Color', 'w');
hold off

cd ../../v/Paper/
print(fig, 'm_vote_share_contour.pdf', '-dpdf', '-vector');
cd ../../d/steadystates/

% changing comparisons, now to baseline populist taxation
tvotes = parray;
for i = 1:13
    for j = 1:ntau
        acond = compute.condense(Warray{i,j}, amu, agrid);
            if j ~= ntau
                Votes_Pop = EVarray{i,j} < EVarray{i,ntau};
                p = Votes_Pop.*acond;
                ev_diff = EVarray{i,j} - EVarray{1,j};
                meanEV_diffs(i,j) = mean(ev_diff(:));
                stdEV_diffs(i,j)  = std(ev_diff(:));
                maxEV_diffs(i,j)  = max(ev_diff(:));
                minEV_diffs(i,j)  = min(ev_diff(:));
            else
                p = zeros(size(acond));
            end
    
            tvotes{i,j} = sum(sum(p));
    end
end

T = cell2mat(tvotes);
T = T(:,1:(ntau-1));

fig = figure;
set(fig, 'Units', 'inches', 'Position', [1, 1, 4.8, 4.8]);
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [4.8, 4.8], ...
         'PaperPosition', [0, 0, 4.8, 4.8]);


[X, Y] = meshgrid(etagrid(1:13), taugrid(1:ntau-1));
pcolor(X, Y, T');
shading interp;                % for smooth color blending
clim([min(T(:)) .75]);
colormap(sunswept);
colorbar;
hold on;
[C, h] = contour(X, Y, T', [0.5 0.5], 'k', 'LineWidth', 1.5);
clabel(C, h, 'Interpreter', 'latex', 'FontSize', 16, 'Color', [0.1 0.1 0.1], ...
    'FontWeight', 'bold');

% Labels and styling
xlabel('Migration Rate $(\eta)$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Tax Progressivity $(\tau)$', 'Interpreter', 'latex', 'FontSize', 22);
title('Percentage Voting for Populists', 'Interpreter', 'latex', 'FontSize', 28);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
set(gcf, 'Color', 'w');
hold off

cd ../../v/Paper/
print(fig, 't_vote_share_contour.pdf', '-dpdf', '-vector');
cd ../../d/steadystates/

% changing comparisons, now to baseline literally nothing
bvotes = parray;
for i = 1:13
    for j = 1:ntau
        acond = compute.condense(Warray{i,j}, amu, agrid);
            if i ~= 1 & j~=1
                Votes_Pop = EVarray{i,j} < EVarray{1,1};
                p = Votes_Pop.*acond;
                ev_diff = EVarray{i,j} - EVarray{1,j};
                meanEV_diffs(i,j) = mean(ev_diff(:));
                stdEV_diffs(i,j)  = std(ev_diff(:));
                maxEV_diffs(i,j)  = max(ev_diff(:));
                minEV_diffs(i,j)  = min(ev_diff(:));
            else
                p = zeros(size(acond));
            end
    
            bvotes{i,j} = sum(sum(p));
    end
end

B = cell2mat(bvotes);
B = B(2:13,2:ntau);

fig = figure;
set(fig, 'Units', 'inches', 'Position', [1, 1, 4.8, 4.8]);
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [4.8, 4.8], ...
         'PaperPosition', [0, 0, 4.8, 4.8]);


[X, Y] = meshgrid(etagrid(2:13), taugrid(2:ntau));
pcolor(X, Y, B');
shading interp;                % for smooth color blending
% clim([min(B(:)) .75]);
colormap(sunswept);
colorbar;
hold on;
[C, h] = contour(X, Y, B', [0.5 0.5], 'k', 'LineWidth', 1.5);
clabel(C, h, 'Interpreter', 'latex', 'FontSize', 16, 'Color', [0.1 0.1 0.1], ...
    'FontWeight', 'bold');

% Labels and styling
xlabel('Migration Rate $(\eta)$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Tax Progressivity $(\tau)$', 'Interpreter', 'latex', 'FontSize', 22);
title('Percentage Voting for Populists', 'Interpreter', 'latex', 'FontSize', 28);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
set(gcf, 'Color', 'w');
hold off

cd ../../v/Paper/
print(fig, 'b_vote_share_contour.pdf', '-dpdf', '-vector');
cd ../../d/steadystates/

%% now in less redistribution

load results_rho85sig2_lamscale_0.9_t0.9000_eta0.5684.mat

eta_grid_labels = string(round(etagrid(2:13),2));
tau_grid_labels = string(round(taugrid,2));
set(gca, 'YDir', 'normal');
set(gca, 'TickLabelInterpreter', 'latex');

P90 = cell2mat(parray);

P90 = P90(2:end,:);

fig = figure;
set(fig, 'Units', 'inches', 'Position', [1, 1, 4.8, 4.8]);
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [4.8, 4.8], ...
         'PaperPosition', [0, 0, 4.8, 4.8]);

[X, Y] = meshgrid(etagrid(2:13), taugrid);
pcolor(X, Y, P90');
shading interp;                % for smooth color blending
clim([.45 .7]);
colormap(sunswept);
colorbar;
hold on;
[C, h] = contour(X, Y, P90', [0.5 0.5], 'k', 'LineWidth', 1.5);
clabel(C, h, 'Interpreter', 'latex', 'FontSize', 16, 'Color', [0.1 0.1 0.1], ...
    'FontWeight', 'bold');
% Add a text annotation
text(0.1, 0.65, 'Populist Majority', 'Interpreter', 'latex', ...
     'Color', [0.2 0.2 0.2], 'FontSize', 12, 'Rotation', 0);

% Labels and styling
xlabel('Migration Rate $(\eta)$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Tax Progressivity $(\tau)$', 'Interpreter', 'latex', 'FontSize', 22);
title('Percentage Voting for Populists, $G = 10\%$', 'Interpreter', 'latex', 'FontSize', 28);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
set(gcf, 'Color', 'w');
hold off

cd ../../v/Paper/
print(fig, 'm_vote_share_l90.pdf', '-dpdf', '-vector');
cd ../../d/steadystates/

% changing comparisons, now to baseline populist taxation
tvotes = parray;
for i = 1:13
    for j = 1:ntau
        acond = compute.condense(Warray{i,j}, amu, agrid);
            if j ~= ntau
                Votes_Pop = EVarray{i,j} < EVarray{i,ntau};
                p = Votes_Pop.*acond;
                ev_diff = EVarray{i,j} - EVarray{1,j};
                meanEV_diffs(i,j) = mean(ev_diff(:));
                stdEV_diffs(i,j)  = std(ev_diff(:));
                maxEV_diffs(i,j)  = max(ev_diff(:));
                minEV_diffs(i,j)  = min(ev_diff(:));
            else
                p = zeros(size(acond));
            end
    
            tvotes{i,j} = sum(sum(p));
    end
end

T90 = cell2mat(tvotes);
T90 = T90(:,1:(ntau-1));

fig = figure;
set(fig, 'Units', 'inches', 'Position', [1, 1, 4.8, 4.8]);
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [4.8, 4.8], ...
         'PaperPosition', [0, 0, 4.8, 4.8]);

[X, Y] = meshgrid(etagrid(1:13), taugrid(1:ntau-1));
pcolor(X, Y, T90');
shading interp;                % for smooth color blending
clim([min(T90(:)) .75]);
colormap(sunswept);
colorbar;
hold on;
[C, h] = contour(X, Y, T90', [0.5 0.5], 'k', 'LineWidth', 1.5);
clabel(C, h, 'Interpreter', 'latex', 'FontSize', 16, 'Color', [0.1 0.1 0.1], ...
    'FontWeight', 'bold');

% Labels and styling
xlabel('Migration Rate $(\eta)$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Tax Progressivity $(\tau)$', 'Interpreter', 'latex', 'FontSize', 22);
title('Percentage Voting for Populists, $G = 10\%$', 'Interpreter', 'latex', 'FontSize', 28);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
set(gcf, 'Color', 'w');
hold off

cd ../../v/Paper/
print(fig, 't_vote_share_l90.pdf', '-dpdf', '-vector');
cd ../../d/steadystates/

% changing comparisons, now to baseline literally nothing
bvotes = parray;
for i = 1:13
    for j = 1:ntau
        acond = compute.condense(Warray{i,j}, amu, agrid);
            if i ~= 1 & j~=1
                Votes_Pop = EVarray{i,j} < EVarray{1,1};
                p = Votes_Pop.*acond;
                ev_diff = EVarray{i,j} - EVarray{1,j};
                meanEV_diffs(i,j) = mean(ev_diff(:));
                stdEV_diffs(i,j)  = std(ev_diff(:));
                maxEV_diffs(i,j)  = max(ev_diff(:));
                minEV_diffs(i,j)  = min(ev_diff(:));
            else
                p = zeros(size(acond));
            end
    
            bvotes{i,j} = sum(sum(p));
    end
end

B90 = cell2mat(bvotes);
B90 = B90(2:13,2:ntau);

fig = figure;
set(fig, 'Units', 'inches', 'Position', [1, 1, 4.8, 4.8]);
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [4.8, 4.8], ...
         'PaperPosition', [0, 0, 4.8, 4.8]);


[X, Y] = meshgrid(etagrid(2:13), taugrid(2:ntau));
pcolor(X, Y, B90');
shading interp;                % for smooth color blending
% clim([min(B(:)) .75]);
colormap(sunswept);
colorbar;
hold on;
[C, h] = contour(X, Y, B90', [0.5 0.5], 'k', 'LineWidth', 1.5);
clabel(C, h, 'Interpreter', 'latex', 'FontSize', 16, 'Color', [0.1 0.1 0.1], ...
    'FontWeight', 'bold');

% Labels and styling
xlabel('Migration Rate $(\eta)$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Tax Progressivity $(\tau)$', 'Interpreter', 'latex', 'FontSize', 22);
title('Percentage Voting for Populists', 'Interpreter', 'latex', 'FontSize', 28);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
set(gcf, 'Color', 'w');
hold off

cd ../../v/Paper/
print(fig, 'b_vote_share_l90.pdf', '-dpdf', '-vector');
cd ../../d/steadystates/

%% the diffs

eta_grid_labels = string(round(etagrid(2:13),2));
tau_grid_labels = string(round(taugrid,2));
set(gca, 'YDir', 'normal');
set(gca, 'TickLabelInterpreter', 'latex');

Pdiff = (P-P90);

fig = figure;
set(fig, 'Units', 'inches', 'Position', [1, 1, 4.8, 4.8]);
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [4.8, 4.8], ...
         'PaperPosition', [0, 0, 4.8, 4.8]);


[X, Y] = meshgrid(etagrid(2:13), taugrid);
pcolor(X, Y, Pdiff');
shading interp;                % for smooth color blending
% clim([.45 .7]);
colormap(sunswept);
colorbar;
hold on;
[C, h] = contour(X, Y, Pdiff', [0.5 0.5], 'k', 'LineWidth', 1.5);
clabel(C, h, 'Interpreter', 'latex', 'FontSize', 16, 'Color', [0.1 0.1 0.1], ...
    'FontWeight', 'bold');

% Labels and styling
xlabel('Migration Rate $(\eta)$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Tax Progressivity $(\tau)$', 'Interpreter', 'latex', 'FontSize', 22);
title('Change in Populist Votes', 'Interpreter', 'latex', 'FontSize', 28);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
set(gcf, 'Color', 'w');
hold off

cd ../../v/Paper/
print(fig, 'm_diff_vote_share.pdf', '-dpdf', '-vector');
cd ../../d/steadystates/

Tdiff = T-T90;

fig = figure;
set(fig, 'Units', 'inches', 'Position', [1, 1, 4.8, 4.8]);
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [4.8, 4.8], ...
         'PaperPosition', [0, 0, 4.8, 4.8]);


[X, Y] = meshgrid(etagrid(1:13), taugrid(1:ntau-1));
pcolor(X, Y, Tdiff');
shading interp;                % for smooth color blending
% clim([min(T(:)) .75]);
colormap(sunswept);
colorbar;
hold on;
[C, h] = contour(X, Y, Tdiff', [0.5 0.5], 'k', 'LineWidth', 1.5);
clabel(C, h, 'Interpreter', 'latex', 'FontSize', 16, 'Color', [0.1 0.1 0.1], ...
    'FontWeight', 'bold');

% Labels and styling
xlabel('Migration Rate $(\eta)$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Tax Progressivity $(\tau)$', 'Interpreter', 'latex', 'FontSize', 22);
title('Change in Populist Votes', 'Interpreter', 'latex', 'FontSize', 28);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
set(gcf, 'Color', 'w');
hold off

cd ../../v/Paper/
print(fig, 't_diff_vote_share.pdf', '-dpdf', '-vector');
cd ../../d/steadystates/

Bdiff = B - B90;

fig = figure;
set(fig, 'Units', 'inches', 'Position', [1, 1, 4.8, 4.8]);
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [4.8, 4.8], ...
         'PaperPosition', [0, 0, 4.8, 4.8]);


[X, Y] = meshgrid(etagrid(2:13), taugrid(2:ntau));
pcolor(X, Y, Bdiff');
shading interp;                % for smooth color blending
% clim([min(B(:)) .75]);
colormap(sunswept);
colorbar;
hold on;
[C, h] = contour(X, Y, Bdiff', [0.5 0.5], 'k', 'LineWidth', 1.5);
clabel(C, h, 'Interpreter', 'latex', 'FontSize', 16, 'Color', [0.1 0.1 0.1], ...
    'FontWeight', 'bold');

% Labels and styling
xlabel('Migration Rate $(\eta)$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Tax Progressivity $(\tau)$', 'Interpreter', 'latex', 'FontSize', 22);
title('Change in Populist Votes', 'Interpreter', 'latex', 'FontSize', 28);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
set(gcf, 'Color', 'w');
hold off

cd ../../v/Paper/
print(fig, 'b_diff_vote_share.pdf', '-dpdf', '-vector');
cd ../../d/steadystates/

%% A(\tau) chart

e_vals = lgrid;                      % Wage earnings axis (x)
K = Karray{1,1};                    % no distortions ss K
w = vaas.calcw(alpha,K,0);

% stealing this code from chatGPT because I am truly dissociating at this
% point
A = @(tau) sum(w .* lgrid .* ldist) ./ sum((w .* lgrid).^(1 - tau) .* ldist);
tau_vals = linspace(0, 1, 100);  
A_vals = arrayfun(A, tau_vals);

tax_collected = w*lgrid(4) - A_vals.*(w*lgrid(4)).^(1-tau_vals);
tax_collected = w*lgrid(1) - A_vals.*(w*lgrid(1)).^(1-tau_vals);

hold on;
tax_collected = w*lgrid(7) - A_vals.*(w*lgrid(7)).^(1-tau_vals);
plot(tau_vals, tax_collected, 'LineWidth', 2, 'Color', '#11302a'); 
tax_collected = w*lgrid(4) - A_vals.*(w*lgrid(4)).^(1-tau_vals);
plot(tau_vals, tax_collected, 'LineWidth', 2, 'Color', '#8f5774'); 
tax_collected = w*lgrid(1) - A_vals.*(w*lgrid(1)).^(1-tau_vals);
plot(tau_vals, tax_collected, 'LineWidth', 2, 'Color', '#dac4d0'); 
xlabel('Taxation Progressivity ($\tau$)', ...
      'FontSize', 18, 'Interpreter', 'latex');
ylabel('Balanced Budget Constant ($A(\tau)$)', ...
      'FontSize', 18, 'Interpreter', 'latex');
title('Change in Balanced Budget Constant $A(\tau)$', ...
      'FontSize', 22, 'Interpreter', 'latex');
hold off

%% taxation collected chart

e_vals = lgrid;                      % Wage earnings axis (x)
K = Karray{1,1};                    % no distortions ss K
w = vaas.calcw(alpha,K,0);

% stealing this code from chatGPT because I am truly dissociating at this
% point
TaxCollected = @(tau) sum(w .* lgrid .* ldist) ./ sum((w .* lgrid).^(1 - tau) .* ldist) ...
                     * sum((w .* lgrid).^(1 - tau) .* ldist);
Tvals = arrayfun(TaxCollected, tau_vals);

hold on;
plot(tau_vals, Tvals, 'LineWidth', 2, 'Color', '#034C53'); 
xlabel('Taxation Progressivity ($\tau$)', ...
      'FontSize', 18, 'Interpreter', 'latex');
ylabel('Balanced Budget Constant ($A(\tau)$)', ...
      'FontSize', 18, 'Interpreter', 'latex');
title('Taxation Collected', ...
      'FontSize', 22, 'Interpreter', 'latex');
hold off

% this is meaningless. starting over by quantiles

%%

y = w .* lgrid;                     % Pre-tax income vector
n_tau = length(tau_vals);            % Number of τ values
n_grid = length(lgrid);            % Number of types

% Initialize output matrices [n_grid × n_tau]
net_tax_paid_mat = zeros(n_grid, n_tau);
mtr_mat          = zeros(n_grid, n_tau);
redistribution_mat = zeros(n_grid, n_tau);
net_income_mat = zeros(n_grid, n_tau);

% Loop over tau values
for j = 1:n_tau
    tau = tau_vals(j);

    % Balanced budget multiplier A(tau)
    A = sum(y .* ldist) / sum((y).^(1 - tau) .* ldist);

    % Post-tax income
    net_income = A .* y.^(1 - tau);

    % Net tax paid
    net_tax_paid = y - net_income;

    % Marginal tax rate
    mtr = A .* (1 - tau) .* y.^(-tau);

    % Redistribution share
    share_received = net_income .* ldist / sum(net_income .* ldist);
    share_earned   = y .* ldist / sum(y .* ldist);
    redistribution = share_received - share_earned;

    % Store results
    net_income_mat(:,j) = net_income;
    net_tax_paid_mat(:, j) = net_tax_paid;
    mtr_mat(:, j) = mtr;
    redistribution_mat(:, j) = redistribution;
end

%changing the colormap
n_types = length(lgrid);
idx = round(linspace(1, size(darkAcademiaChic, 1), n_types));
seven_colors = darkAcademiaChic(idx, :);  % Now [7 × 3]

figure;
hold on;

% Use a colormap with as many colors as ε types
n_grid = length(lgrid);
legendEntries = cell(length(1:2:n_types),1);

% Plot every other ε-type (1, 3, 5, 7)
legix = 1;
for i = 1:2:n_types
    plot(tau_vals, net_tax_paid_mat(i, :), 'Color', seven_colors(i, :), 'LineWidth', 1.5);
    legendEntries{legix} = sprintf('$\\varepsilon = %.2f$', lgrid(i));
    legix = legix + 1;
end


% Label the colorbar with ε values
cb.Ticks = linspace(0, 1, 5);  % 5 tick marks
cb.TickLabels = arrayfun(@(x) sprintf('%.2f', x), ...
                         linspace(min(lgrid), max(lgrid), 5), 'UniformOutput', false);

xlabel('Tax Progressivity ($\tau$)', ...
      'FontSize', 18, 'Interpreter', 'latex');
ylabel('Net Tax Paid', ...
      'FontSize', 18, 'Interpreter', 'latex');
title('Net Tax Paid by Household Type ($\varepsilon$)', ...
      'FontSize', 22, 'Interpreter', 'latex');

legend(legendEntries, 'Location', 'bestoutside', 'Interpreter', 'latex', ...
    'FontSize', 18);
hold off

%
figure;
hold on;

% Plot every other ε-type (1, 3, 5, 7)
legix = 1;
for i = 1:2:n_types
    plot(tau_vals, mtr_mat(i, :), 'Color', seven_colors(i, :), 'LineWidth', 1.5);
    legendEntries{legix} = sprintf('$\\varepsilon = %.2f$', lgrid(i));
    legix = legix + 1;
end


% Label the colorbar with ε values
cb.Ticks = linspace(0, 1, 5);  % 5 tick marks
cb.TickLabels = arrayfun(@(x) sprintf('%.2f', x), ...
                         linspace(min(lgrid), max(lgrid), 5), 'UniformOutput', false);

xlabel('Tax Progressivity ($\tau$)', ...
      'FontSize', 18, 'Interpreter', 'latex');
ylabel('Net Tax Paid', ...
      'FontSize', 18, 'Interpreter', 'latex');
title('Marginal Tax Rate by HH Type ($\varepsilon$)', ...
      'FontSize', 22, 'Interpreter', 'latex');

legend(legendEntries, 'Location', 'bestoutside', 'Interpreter', 'latex', ...
    'FontSize', 18);
hold off

%
figure;
hold on;

% Plot every other ε-type (1, 3, 5, 7)
legix = 1;
for i = 1:2:n_types
    plot(tau_vals, redistribution_mat(i, :), 'Color', seven_colors(i, :), 'LineWidth', 1.5);
    legendEntries{legix} = sprintf('$\\varepsilon = %.2f$', lgrid(i));
    legix = legix + 1;
end


% Label the colorbar with ε values
cb.Ticks = linspace(0, 1, 5);  % 5 tick marks
cb.TickLabels = arrayfun(@(x) sprintf('%.2f', x), ...
                         linspace(min(lgrid), max(lgrid), 5), 'UniformOutput', false);

xlabel('Tax Progressivity ($\tau$)', ...
      'FontSize', 18, 'Interpreter', 'latex');
ylabel('Net Tax Paid', ...
      'FontSize', 18, 'Interpreter', 'latex');
title('Redistribution by HH Type ($\varepsilon$)', ...
      'FontSize', 22, 'Interpreter', 'latex');

legend(legendEntries, 'Location', 'bestoutside', 'Interpreter', 'latex', ...
    'FontSize', 18);
hold off


fig = figure;
set(fig, 'Units', 'inches', 'Position', [1, 1, 12.8, 7.2]);
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [12.8, 7.2], ...
         'PaperPosition', [0, 0, 12.8, 7.2]);
hold on;

% Plot every other ε-type (1, 3, 5, 7)
legix = 1;
for i = 1:2:n_types
    plot(tau_vals, net_income_mat(i, :), 'Color', seven_colors(i, :), 'LineWidth', 1.5);
    legendEntries{legix} = sprintf('$\\varepsilon = %.2f$', lgrid(i));
    legix = legix + 1;
end

xlabel('Tax Progressivity ($\tau$)', ...
      'FontSize', 18, 'Interpreter', 'latex');
ylabel('Net Income ($\lambda(R)(w \varepsilon)^{1-\tau}$)', ...
      'FontSize', 18, 'Interpreter', 'latex');
title('Net Income by Productivity Type ($\varepsilon$)', ...
      'FontSize', 22, 'Interpreter', 'latex');

legend(legendEntries, 'Location', 'bestoutside', 'Interpreter', 'latex', ...
    'FontSize', 18);
hold off
print(fig, 'net_income_tau.pdf', '-dpdf', '-vector');


% Not sure if I want to do anything with the code below, but it's nice to
% have it. 
% 
% eps_thresh = zeros(size(tau_vals));
% for j = 1:length(tau_vals)
%     tau = tau_vals(j);
%     if tau == 0
%         eps_thresh(j) = NaN;  % Avoid log(0) or div-by-zero
%     else
%         A = sum(w .* lgrid .* ldist) / sum((w .* lgrid).^(1 - tau) .* ldist);
%         eps_thresh(j) = A^(1 / tau) / w;
%     end
% end
% plot(tau_vals, eps_thresh, 'k--', 'LineWidth', 2);
% text(0.5, eps_thresh(find(tau_vals > 0.5, 1)), ...
%      '$\varepsilon^*(\tau)$ threshold', 'Interpreter', 'latex', 'FontSize', 14);