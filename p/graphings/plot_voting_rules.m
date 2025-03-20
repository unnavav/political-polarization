% Initialize figure
figure;
hold on;
colors = lines(3); % Generate distinct colors
legendEntries = cell(3,1); % Preallocate legend text

Votes = VL>VP;
Votes = Votes([1 4 7],:);
egrid = lgrid([1 4 7]);
% Loop over each productivity level
for i = 1:3
    e = egrid(i); % Get productivity level
     
    % Plot each line
    ve = Votes(i,1:170);
    plot(agrid(1:170), ve, 'Color', colors(i,:), 'LineWidth', 2);
    legendEntries{i} = sprintf('$\\varepsilon = %.2f$', e);

end

% Format plot
yticks([0 1]); % Set tick positions at 0 and 1
yticklabels({'Vote against', 'Vote for'}); % Custom labels

xlabel('Assets $a$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel(' ', 'FontSize', 16, 'Interpreter', 'latex');
title('Votes for High Migration ($\eta$)', ...
      'FontSize', 18, 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', 14);

legend(legendEntries, 'Location', 'bestoutside', 'Interpreter', 'latex');
