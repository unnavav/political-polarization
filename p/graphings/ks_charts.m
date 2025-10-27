load ks_rmseed_in_pr_exo_R

Kfore1 = Kfore(1,:);
Kfore2 = Kfore(2,:);

x = 5:.2:20;
y1 = exp(Kfore1(1) + Kfore1(2)*log(x));
y2 = exp(Kfore2(1) + Kfore2(2)*log(x));
hold on
plot(x, y1, 'LineWidth', 2, 'Color', '#034C53')
plot(x, y2', 'LineWidth', 2, 'Color', '#F38C79')

xlabel('Current Capital $K$', ...
      'FontSize', 18, 'Interpreter', 'latex');
ylabel('Future Capital $K$', ...
      'FontSize', 18, 'Interpreter', 'latex');
title('Capital Law of Motion (50/50 Switching)', ...
      'FontSize', 22, 'Interpreter', 'latex');

legendEntries = cell(2,1); % Preallocate legend text
legendEntries{1} = "Populist";
legendEntries{2} = "Liberalist";
legend(legendEntries, 'Location', 'bestoutside', 'FontSize', 18,'Interpreter', 'latex');
hold off


load ks_rmseed_in_pr_exo_R
y3 = exp(Kfore1(1) + Kfore1(2)*log(x));
y4 = exp(Kfore2(1) + Kfore2(2)*log(x));

hold on;
plot(x, y1, 'LineWidth', 2, 'Color', '#034C53')
plot(x, y2', 'LineWidth', 2, 'Color', '#F38C79')
plot(x, y3, 'LineWidth', 2,  'Color', 'red')
plot(x, y4', 'LineWidth', 2,  'Color', 'blue')
xlabel('Current Capital $K$', ...
      'FontSize', 18, 'Interpreter', 'latex');
ylabel('Future Capital $K$', ...
      'FontSize', 18, 'Interpreter', 'latex');
title('Capital Law of Motion (50/50 Switching)', ...
      'FontSize', 22, 'Interpreter', 'latex');

legendEntries = cell(4,1); % Preallocate legend text
legendEntries{1} = "Populist 50-50";
legendEntries{2} = "Liberalist 50-50";
legendEntries{3} = "Populist (High Persistence)";
legendEntries{4} = "Liberalist (High Persistence)";
legend(legendEntries, 'Location', 'bestoutside', 'FontSize', 18,'Interpreter', 'latex');
hold off