%% Charting Heat Maps
% 
cd C:\VAASAVI\Dropbox\Education\OSU\Ongoing_Research\Populism\political-polarization\p
restoredefaultpath;
clear all; clc;
addpath(genpath(pwd));

n = 3; % Number of colors
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

cd ../../v/
mkdir thresholds
cd thresholds

smallylabels = round(agrid(100:5:150),2);

liberalism = EVarray{7,1};
populist = EVarray{1,1};
A = populist>liberalism;

fig = figure;
set(fig, 'Units', 'inches', 'Position', [1, 1, 8.5, 6.4]);
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [8.5, 6.4], ...
         'PaperPosition', [0, 0, 8.5, 6.4]);

hold on
imagesc(lgrid, agrid(100:150), A')
% Labels and styling
colormap(darkAcademiaChic)
h1 = plot(nan, nan, 's', 'MarkerFaceColor', darkAcademiaChic(1,:), ...
    'MarkerEdgeColor', darkAcademiaChic(1,:), 'DisplayName', 'Liberalist');
h3 = plot(nan, nan, 's', 'MarkerFaceColor', darkAcademiaChic(3,:), ...
    'MarkerEdgeColor', darkAcademiaChic(3,:), 'DisplayName', 'Populist');

legend([h1 h3], 'Location', 'eastoutside', 'FontSize', 12, 'Interpreter', 'latex')
% legend('Position', [0.2 0.75 0.2 0.1]);

xlabel('Idiosyncratic Productivity $(\varepsilon)$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Asset Position $(a)$', 'Interpreter', 'latex', 'FontSize', 22);
% xticklabels(round(lgrid,2));
yticklabels(smallylabels);
title('Households Voting for Populists', 'Interpreter', 'latex', 'FontSize', 28);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
set(gcf, 'Color', 'w');
hold off
print(fig, 'baseline_migration.pdf', '-dpdf', '-vector');

liberalism = EVarray{1,1};
populist = EVarray{1,7};
A = populist>liberalism;

fig = figure;
set(fig, 'Units', 'inches', 'Position', [1, 1, 8.5, 6.4]);
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [8.5, 6.4], ...
         'PaperPosition', [0, 0, 8.5, 6.4]);

hold on
imagesc(lgrid, agrid(100:150), A')
% Labels and styling
colormap(darkAcademiaChic)
h1 = plot(nan, nan, 's', 'MarkerFaceColor', darkAcademiaChic(1,:), ...
    'MarkerEdgeColor', darkAcademiaChic(1,:), 'DisplayName', 'Liberalist');
h3 = plot(nan, nan, 's', 'MarkerFaceColor', darkAcademiaChic(3,:), ...
    'MarkerEdgeColor', darkAcademiaChic(3,:), 'DisplayName', 'Populist');

legend([h1 h3], 'Location', 'eastoutside', 'FontSize', 12, 'Interpreter', 'latex')
% legend('Position', [0.2 0.75 0.2 0.1]);

xlabel('Idiosyncratic Productivity $(\varepsilon)$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Asset Position $(a)$', 'Interpreter', 'latex', 'FontSize', 22);
% xticklabels(round(lgrid,2));
yticklabels(smallylabels);
title('Households Voting for Populists', 'Interpreter', 'latex', 'FontSize', 28);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
set(gcf, 'Color', 'w');
hold off
print(fig, 'baseline_taxation.pdf', '-dpdf', '-vector');


xlabels = round(lgrid,2);
ylabels = round(agrid(50:50:250),2);

heatmap(cell2mat(parray(2:13,:)))

liberalism = EVarray{13,1};
populist = EVarray{1,1};
A = populist>liberalism;

liberalism = EVarray{3,1};
populist = EVarray{1,1};
B = populist>liberalism;

liberalism = EVarray{13,1};
populist = EVarray{3,1};
D = populist > liberalism;

C = A+B;
E = C + D;
C = C(:,100:150);
E = E(:,100:150);

fig = figure;
set(fig, 'Units', 'inches', 'Position', [1, 1, 6.375, 4.8]);
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [6.375, 4.8], ...
         'PaperPosition', [0, 0, 6.375, 4.8]);
% set(gca, 'Position', [0 0 1 1])  % fills the figure canvas

hold on
imagesc(lgrid, agrid(100:150), C')
% Labels and styling
colormap(darkAcademiaChic)
h1 = plot(nan, nan, 's', 'MarkerFaceColor', darkAcademiaChic(1,:), ...
    'MarkerEdgeColor', darkAcademiaChic(1,:), 'DisplayName', 'A');
h2 = plot(nan, nan, 's', 'MarkerFaceColor', darkAcademiaChic(2,:), ...
    'MarkerEdgeColor', darkAcademiaChic(2,:), 'DisplayName', 'B');
h3 = plot(nan, nan, 's', 'MarkerFaceColor', darkAcademiaChic(3,:), ...
    'MarkerEdgeColor', darkAcademiaChic(3,:), 'DisplayName', 'C');

legend([h1 h2 h3], 'Location', 'eastoutside', 'FontSize', 12, 'Interpreter', 'latex')

xlabel('Idiosyncratic Productivity $(\varepsilon)$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Asset Position $(a)$', 'Interpreter', 'latex', 'FontSize', 22);
% xticklabels(round(lgrid,2));
yticklabels(smallylabels);
title('Results of Different Pairwise Votes', 'Interpreter', 'latex', 'FontSize', 28);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
set(gcf, 'Color', 'w');
hold off
print(fig, 'migration_shifts_56_09.pdf', '-dpdf', '-vector');

liberalism = EVarray{1,1};
populist = EVarray{1,13};
mesh(populist>liberalism)
A = populist>liberalism;

liberalism = EVarray{1,1};
populist = EVarray{1,7};
B = populist>liberalism;

C = A+B;
mesh(C)

fig = figure;
set(fig, 'Units', 'inches', 'Position', [1, 1, 6.375, 4.8]);
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [6.375, 4.8], ...
         'PaperPosition', [0, 0, 6.375, 4.8]);
hold on
imagesc(lgrid, agrid(100:150), C')
% Labels and styling
colormap(darkAcademiaChic)
h1 = plot(nan, nan, 's', 'MarkerFaceColor', darkAcademiaChic(1,:), ...
    'MarkerEdgeColor', darkAcademiaChic(1,:), 'DisplayName', 'A');
h2 = plot(nan, nan, 's', 'MarkerFaceColor', darkAcademiaChic(2,:), ...
    'MarkerEdgeColor', darkAcademiaChic(2,:), 'DisplayName', 'B');
h3 = plot(nan, nan, 's', 'MarkerFaceColor', darkAcademiaChic(3,:), ...
    'MarkerEdgeColor', darkAcademiaChic(3,:), 'DisplayName', 'C');

legend([h1 h2 h3], 'Location', 'eastoutside', 'FontSize', 12, 'Interpreter', 'latex')

xlabel('Idiosyncratic Productivity $(\varepsilon)$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Asset Position $(a)$', 'Interpreter', 'latex', 'FontSize', 22);
% xticklabels(round(lgrid,2));
yticklabels(smallylabels);
title('Households Voting for Populists', 'Interpreter', 'latex', 'FontSize', 28);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
set(gcf, 'Color', 'w');
hold off
print(fig, 'taxation_shifts_56_09_.pdf', '-dpdf', '-vector');

liberalism = EVarray{13,7};
populist = EVarray{7,13};
A = populist>liberalism;

liberalism = EVarray{3,1};
populist = EVarray{1,3};
B = populist>liberalism;

C = A+B;

fig = figure;
set(fig, 'Units', 'inches', 'Position', [1, 1, 6.375, 4.8]);
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [6.375, 4.8], ...
         'PaperPosition', [0, 0, 6.375, 4.8]);

hold on
imagesc(lgrid, agrid(100:150), C')
% Labels and styling
colormap(darkAcademiaChic)
h1 = plot(nan, nan, 's', 'MarkerFaceColor', darkAcademiaChic(1,:), ...
    'MarkerEdgeColor', darkAcademiaChic(1,:), 'DisplayName', 'Never Votes');
h2 = plot(nan, nan, 's', 'MarkerFaceColor', darkAcademiaChic(2,:), ...
    'MarkerEdgeColor', darkAcademiaChic(2,:), 'DisplayName', 'Votes Sometimes');
h3 = plot(nan, nan, 's', 'MarkerFaceColor', darkAcademiaChic(3,:), ...
    'MarkerEdgeColor', darkAcademiaChic(3,:), 'DisplayName', 'Always Votes');

legend([h1 h2 h3], 'Location', 'eastoutside', 'FontSize', 12, 'Interpreter', 'latex')

xlabel('Idiosyncratic Productivity $(\varepsilon)$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Asset Position $(a)$', 'Interpreter', 'latex', 'FontSize', 22);
% xticklabels(round(lgrid,2));
yticklabels(smallylabels);
title('Households Voting for Populists', 'Interpreter', 'latex', 'FontSize', 28);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
set(gcf, 'Color', 'w');
hold off
print(fig, 'both_shifts_56_09_.pdf', '-dpdf', '-vector');

%%

% Loop through EVarray dimensions (rows: τ, cols: η)
[nEta, nTau] = size(EVarray(1:13,:));
agrid_span = 100:150;

baseline_EV = EVarray{5, 3}; %kinda close to US labor force but idk
customColors = [1 1 1; darkAcademiaChic(1,:); darkAcademiaChic(2,:); darkAcademiaChic(3,:)];

% Directory setup
mkdir results_tikz_inner;
cd results_tikz_inner;

% Loop through valid populist-liberal pairs
for low_i = 1:floor(nTau/2)
    for high_j = ceil(nEta/2):nEta

        high_i = nTau - low_i + 1;  % liberal τ
        low_j  = nEta - high_j + 1; % populist η

        % Skip if out of bounds
        if high_i < 1 || low_j < 1 || high_i > nTau || low_j > nEta
            continue;
        end

        % Voting map under this regime
        populist = EVarray{low_j, high_i};  % low τ, high η
        liberal  = EVarray{high_j, low_i};  % high τ, low η

        voteP = populist > liberal;
        voteL = liberal  > populist;

        % Region code map: 0 = never, 1 = populist only, 2 = liberal only, 3 = both
        regionMap = uint8(voteP) + 2 * uint8(voteL);
  
        % === PLOT ===
        fig = figure();
        set(fig, 'Units', 'inches', 'Position', [1, 1, 6.375, 4.8]);
        set(fig, 'PaperUnits', 'inches', ...
                 'PaperSize', [6.375, 4.8], ...
                 'PaperPosition', [0, 0, 6.375, 4.8]);

        hold on
        imagesc(lgrid, agrid(agrid_span), regionMap(:,agrid_span)');
        colormap(customColors);

        h1 = plot(nan, nan, 's', 'MarkerFaceColor', darkAcademiaChic(1,:), ...
            'MarkerEdgeColor', darkAcademiaChic(1,:), 'DisplayName', 'Populist Only');
        h2 = plot(nan, nan, 's', 'MarkerFaceColor', darkAcademiaChic(2,:), ...
            'MarkerEdgeColor', darkAcademiaChic(2,:), 'DisplayName', 'Liberal Only');
        h3 = plot(nan, nan, 's', 'MarkerFaceColor', darkAcademiaChic(3,:), ...
            'MarkerEdgeColor', darkAcademiaChic(3,:), 'DisplayName', 'Both');

        legend([h1 h2 h3], 'Location', 'eastoutside', 'FontSize', 12, 'Interpreter', 'latex');

        xlabel('Idiosyncratic Productivity $(\varepsilon)$', 'Interpreter', 'latex', 'FontSize', 22);
        ylabel('Asset Position $(a)$', 'Interpreter', 'latex', 'FontSize', 22);
        yticklabels(smallylabels);
        title(sprintf('$\\tau_L = %d$, $\\eta_P = %d$', high_i, high_j), ...
              'Interpreter', 'latex', 'FontSize', 24);

        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
        set(gcf, 'Color', 'w');
        hold off

        % Save
        fname = sprintf('innerspan_tauL%d_etaP%d.pdf', high_i, high_j);
        print(fig, fname, '-dpdf', '-vector');
        close(fig);

    end
end

cd ..