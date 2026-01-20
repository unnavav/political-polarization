%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Charting Exogenous R Forecasting Rules
% vaasavi | nov 2025
%
% this code charts the forecast rules for exogenous regime switching
% it's basically debugging by staring at how the model behaves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

restoredefaultpath;
clear all; clc;
cd ../../ % get back to parent code directory
addpath(genpath(pwd));
cd ../v
mkdir Forecast_Rules

cd ../d/
load ks_rmseed_in_pr_exo_logit_R_tau086600_eta0000.mat
cd ../v/Forecast_Rules
mkdir tau_086600
cd tau_086600

for i = 1:length(Kforearray)

    Kfore = Kforearray{i};

    Kfore1 = Kfore(1,:);
    Kfore2 = Kfore(2,:);
    Kfore1 = exp(Kfore1(1) + Kfore1(2).*log(Kgrid));
    Kfore2 = exp(Kfore2(1) + Kfore2(2).*log(Kgrid));

    graph_title = sprintf('Capital Forecast Round %i', i);
    
    f = figure;

    set(gcf, 'Position', [100, 100, 800, 500]); 

    hold on;
    plot(Kgrid,Kfore1 , 'LineWidth', 1, 'Color', '#034C53');  % overlay smooth cutoff curve
    plot(Kgrid,Kfore2, 'LineWidth', 1, 'Color', '#F38C79');
    plot(Kgrid,Kgrid, '--k',   'LineWidth', 1)
    ylim([min(Kgrid)*.85 max(Kgrid)*1.15])

    xlabel('Today $K$', ...
          'FontSize', 18, 'Interpreter', 'latex');
    ylabel('Future $K$', ...
          'FontSize', 18, 'Interpreter', 'latex');
    title(graph_title, ...
          'FontSize', 22, 'Interpreter', 'latex');
    subtitle('$\tau_p = .6; \tau_l = .09; \eta_p = \eta_l = 0$', ...
          'FontSize', 18, 'Interpreter', 'latex')
    
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
    legendEntries{1} = "Populism";
    legendEntries{2} = "Liberalism";
    legendEntries{3} = "45 Degree Line";
    legend(legendEntries, 'Location', 'bestoutside', 'Interpreter', 'latex');
    hold off

    exportgraphics(f,sprintf("%i.png",i))
end

cd ../../../d/
load ks_rmseed_in_pr_exo_logit_R_tau086018_eta0000.mat
cd ../v/Forecast_Rules
mkdir tau_086018
cd tau_086018

for i = 1:length(Kforearray)

    Kfore = Kforearray{i};

    Kfore1 = Kfore(1,:);
    Kfore2 = Kfore(2,:);
    Kfore1 = exp(Kfore1(1) + Kfore1(2).*log(Kgrid));
    Kfore2 = exp(Kfore2(1) + Kfore2(2).*log(Kgrid));

    graph_title = sprintf('Capital Forecast Round %i', i);
    f = figure;

    set(gcf, 'Position', [100, 100, 800, 500]); 

    hold on;
    plot(Kgrid,Kfore1 , 'LineWidth', 1, 'Color', '#034C53');  % overlay smooth cutoff curve
    plot(Kgrid,Kfore2, 'LineWidth', 1, 'Color', '#F38C79');
    plot(Kgrid,Kgrid, '--k',   'LineWidth', 1)
    ylim([min(Kgrid)*.85 max(Kgrid)*1.15])
    xlabel('Today $K$', ...
          'FontSize', 18, 'Interpreter', 'latex');
    ylabel('Future $K$', ...
          'FontSize', 18, 'Interpreter', 'latex');
    title(graph_title, ...
          'FontSize', 22, 'Interpreter', 'latex');
    subtitle('$ \tau_p = .181; \tau_l = .086; \eta_p = \eta_l = 0 $', ...
          'FontSize', 16, 'Interpreter', 'latex')
    
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
    legendEntries{1} = "Populism";
    legendEntries{2} = "Liberalism";
    legendEntries{3} = "45 Degree Line";
    legend(legendEntries, 'Location', 'bestoutside', 'Interpreter', 'latex');
    hold off

    exportgraphics(f,sprintf("%i.png",i))
end

cd ../../../d/
load ks_rmseed_in_pr_exo_logit_R_tau000000_eta5000.mat
cd ../v/Forecast_Rules
mkdir eta_500000
cd eta_500000

for i = 1:length(Kforearray)

    Kfore = Kforearray{i};

    Kfore1 = Kfore(1,:);
    Kfore2 = Kfore(2,:);
    Kfore1 = exp(Kfore1(1) + Kfore1(2).*log(Kgrid));
    Kfore2 = exp(Kfore2(1) + Kfore2(2).*log(Kgrid));

    graph_title = sprintf('Capital Forecast Round %i', i);
    f = figure;

    set(gcf, 'Position', [100, 100, 800, 500]); 

    hold on;
    plot(Kgrid,Kfore1 , 'LineWidth', 1, 'Color', '#034C53');  % overlay smooth cutoff curve
    plot(Kgrid,Kfore2, 'LineWidth', 1, 'Color', '#F38C79');
    plot(Kgrid,Kgrid, '--k',   'LineWidth', 1)
    ylim([min(Kgrid)*.85 max(Kgrid)*1.15])
    xlabel('Today $K$', ...
          'FontSize', 18, 'Interpreter', 'latex');
    ylabel('Future $K$', ...
          'FontSize', 18, 'Interpreter', 'latex');
    title(graph_title, ...
          'FontSize', 22, 'Interpreter', 'latex');
    subtitle('$ \tau_p = \tau_l = 0; \eta_p = 0;  \eta_l = .5 $', ...
          'FontSize', 16, 'Interpreter', 'latex')
    
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
    legendEntries{1} = "Populism";
    legendEntries{2} = "Liberalism";
    legendEntries{3} = "45 Degree Line";
    legend(legendEntries, 'Location', 'bestoutside', 'Interpreter', 'latex');
    hold off

    exportgraphics(f,sprintf("%i.png",i))
end


cd ../../../d/
load ks_rmseed_in_pr_exo_logit_R_tau000000_eta2500.mat
cd ../v/Forecast_Rules
mkdir eta_2500000
cd eta_2500000

for i = 1:length(Kforearray)

    Kfore = Kforearray{i};

    Kfore1 = Kfore(1,:);
    Kfore2 = Kfore(2,:);
    Kfore1 = exp(Kfore1(1) + Kfore1(2).*log(Kgrid));
    Kfore2 = exp(Kfore2(1) + Kfore2(2).*log(Kgrid));

    graph_title = sprintf('Capital Forecast Round %i', i);
    f = figure;

    set(gcf, 'Position', [100, 100, 800, 500]); 

    hold on;
    plot(Kgrid,Kfore1 , 'LineWidth', 1, 'Color', '#034C53');  % overlay smooth cutoff curve
    plot(Kgrid,Kfore2, 'LineWidth', 1, 'Color', '#F38C79');
    plot(Kgrid,Kgrid, '--k', 'LineWidth', 1)
    ylim([min(Kgrid)*.85 max(Kgrid)*1.15])
    xlabel('Today $K$', ...
          'FontSize', 18, 'Interpreter', 'latex');
    ylabel('Future $K$', ...
          'FontSize', 18, 'Interpreter', 'latex');
    title(graph_title, ...
          'FontSize', 22, 'Interpreter', 'latex');
    subtitle('$ \tau_p = \tau_l = 0; \eta_p = 0;  \eta_l = .25 $', ...
          'FontSize', 16, 'Interpreter', 'latex')
    
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
    legendEntries{1} = "Populism";
    legendEntries{2} = "Liberalism";
    legendEntries{3} = "45 Degree Line";
    legend(legendEntries, 'Location', 'bestoutside', 'Interpreter', 'latex');
    hold off

    exportgraphics(f,sprintf("%i.png",i))
end