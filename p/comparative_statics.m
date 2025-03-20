eta_grid = linspace(0, 1, 100);
ks = linspace(5,20,100);

R = zeros(100,100);
W = R;
dR = R;
dW = W;

for ie = 1:100
    for ik = 1:100
        R(ie, ik) = vaas.calcr(alpha,delta,ks(ik),eta_grid(ie));
        W(ie, ik) = vaas.calcw(alpha,ks(ik),eta_grid(ie));
        dR(ie, ik) = vaas.drde(alpha,ks(ik),eta_grid(ie));
        dW(ie, ik) = vaas.dwde(alpha,ks(ik),eta_grid(ie));        
    end
end

K_values = [5, 10, 15]; % Three different capital levels
% Initialize figures
figure;
tiledlayout(1,2); % Two plots side by side
% Set Beamer-friendly aesthetics
set(gcf, 'Color', 'w'); % White background for Beamer
set(gcf, 'Units', 'inches', 'Position', [1, 1, 8.4, 4.8]);
set(gca, 'FontSize', 14, 'FontName', 'CMU Serif');

% Interest Rate Plot
nexttile;
hold on;
colors = lines(length(K_values)); % Generate distinct colors
for i = 1:length(K_values)
    K = K_values(i);
    r = alpha .* (K).^(alpha - 1) .* (1 + eta_grid).^(1 - alpha) - delta;
    plot(eta_grid, r, 'Color', colors(i,:), 'LineWidth', 2);
    % Label each line directly at a fixed position
    text(eta_grid(end-10), r(end-20), sprintf('K = %.1f', K), 'FontSize', 14, ...
         'FontName', 'CMU Serif', 'Interpreter', 'latex', 'Color', colors(i,:));
end
hold off;
xlabel('$\eta$ (Migration in Percent)', 'FontSize', 16, 'FontName', 'CMU Serif', 'Interpreter', 'latex');
ylabel('Interest Rate $r$', 'FontSize', 16, 'FontName', 'CMU Serif', 'Interpreter', 'latex');
title('Effect of Migration on Interest Rate', 'FontSize', 18, 'FontName', 'CMU Serif', 'Interpreter', 'latex');
grid on;

% Wage Plot
nexttile;
hold on;
for i = 1:length(K_values)
    K = K_values(i);
    w = (1 - alpha) .* (K).^alpha .* (1 + eta_grid).^(-alpha);
    plot(eta_grid, w, 'Color', colors(i,:), 'LineWidth', 2);
    % Label each line directly at a fixed position
    text(eta_grid(end-10), w(end-20), sprintf('K = %.1f', K), 'FontSize', 14, ...
         'FontName', 'CMU Serif', 'Interpreter', 'latex', 'Color', colors(i,:));
end
hold off;
xlabel('$\eta$ (Migration in Percent)', 'FontSize', 16, 'FontName', 'CMU Serif', 'Interpreter', 'latex');
ylabel('Wage $w$', 'FontSize', 16, 'FontName', 'CMU Serif', 'Interpreter', 'latex');
title('Effect of Migration on Wages', 'FontSize', 18, 'FontName', 'CMU Serif', 'Interpreter', 'latex');
grid on;

cd ../v/
print(gcf, 'comparative_statics.png', '-dpng', '-r300');

%% Initialize figures for partial derivatives
figure;
tiledlayout(1,2); % Two plots side by side
% Set Beamer-friendly aesthetics
set(gcf, 'Color', 'w'); % White background for Beamer
set(gcf, 'Units', 'inches', 'Position', [1, 1, 8.4, 4.8]);
set(gca, 'FontSize', 14, 'FontName', 'CMU Serif');

% Interest Rate Plot
nexttile;
hold on;
colors = lines(length(K_values)); % Generate distinct colors
for i = 1:length(K_values)
    K = K_values(i);
    dr = vaas.drde(alpha, K, eta_grid);
    plot(eta_grid, dr, 'Color', colors(i,:), 'LineWidth', 2);
    % Label each line directly at a fixed position
    text(eta_grid(end-10), dr(end-20), sprintf('K = %.1f', K), 'FontSize', 14, ...
         'FontName', 'CMU Serif', 'Interpreter', 'latex', 'Color', colors(i,:));
end
hold off;
xlabel('$\eta$ (Migration in Percent)', 'FontSize', 16, 'FontName', 'CMU Serif', 'Interpreter', 'latex');
ylabel('Interest Rate $r$', 'FontSize', 16, 'FontName', 'CMU Serif', 'Interpreter', 'latex');
title('$\frac{\partial r}{\partial \eta}$', 'FontSize', 24, 'FontName', 'CMU Serif', 'Interpreter', 'latex');
grid on;

% Wage Plot
nexttile;
hold on;
for i = 1:length(K_values)
    K = K_values(i);
    dw = vaas.dwde(alpha, K, eta_grid);
    plot(eta_grid, dw, 'Color', colors(i,:), 'LineWidth', 2);
    % Label each line directly at a fixed position
    text(eta_grid(end-10), dw(end-20), sprintf('K = %.1f', K), 'FontSize', 14, ...
         'FontName', 'CMU Serif', 'Interpreter', 'latex', 'Color', colors(i,:));
end
hold off;
xlabel('$\eta$ (Migration in Percent)', 'FontSize', 16, 'FontName', 'CMU Serif', 'Interpreter', 'latex');
ylabel('Wage $w$', 'FontSize', 16, 'FontName', 'CMU Serif', 'Interpreter', 'latex');
title('$\frac{\partial w}{\partial \eta}$', 'FontSize', 24, 'FontName', 'CMU Serif', 'Interpreter', 'latex');
grid on;

cd ../v/
print(gcf, 'comparative_statics_partials.png', '-dpng', '-r300');
