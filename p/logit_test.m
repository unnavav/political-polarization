x1 = linspace(0,.03,10000);
% y1 = randi(2,[1000, 1])-1;
y1 = rand([10000, 1])>=.95;


nlmtest = fitnlm(x1, y1, mymodelfun, beta0);
nlmtest.Coefficients.Estimate

scatter(x1 , y1)
hold on
fplot(@(kk) mymodelfun(nlmtest.Coefficients.Estimate, kk), [0,1])
hold off


%% kernel distributins

ksdensity(K_curr(R_curr == 1))
hold on;
ksdensity(K_curr(R_curr == 2))


%% coeffs

x1 = -log(.95^-1);
x2 = -log(.05^-1);

scatter(X1 , Y1)
hold on
fplot(@(kk) mymodelfun(nlm1.Coefficients.Estimate, kk), [min(X1), max(X1)])
hold off

%% maximum likelihood

betagrid = linspace(-10, 10, 100);

% Evaluate the likelihood over the beta grid
likelihood1 = zeros(100);
likelihood2 = likelihood1;
likelihood3 = likelihood1;
likelihood4 = likelihood1;

for i = 1:100
    for j = 1:100
        beta = [betagrid(i) betagrid(j)];

        p1 = mymodelfun(beta, exp(X1));
        ptest = mymodelfun(beta, x1);
      
        % Raw likelihood (product form)
        likelihood3(i, j) = prod(p1.^Y1 .* (1-p1).^(1-Y1));
        likelihood4(i, j) = prod(ptest'.^y1 .* (1-ptest').^(1-y1));

        % calculating likelihood
        likelihood1(i, j) = sum(log(mymodelfun(beta, X1)) .* Y1 + ...
            log(1 - mymodelfun(beta, X1)) .* (1 - Y1));
        likelihood2(i, j) = sum(log(mymodelfun(beta, X2)) .* Y2 + ...
            log(1 - mymodelfun(beta, X2)) .*(1- Y2));
        likelihoodexp(i, j) = sum(log(mymodelfun(beta, exp(X1))) .* Y1 + ...
            log(1 - mymodelfun(beta, X1)) .*(1- Y1));
    end
end



scatter(X1 , Y1)
hold on
fplot(@(kk) mymodelfun([betagrid(81) betagrid(81)], kk), [min(X1), max(X1)])
hold off

[x,y] = meshgrid(betagrid, betagrid);
mesh(x, y, log(likelihood4))

%% Contour plot (shows the ridge clearly)
figure('Position', [100 100 700 600]);
contourf(x, y, log(likelihood3), 30, 'LineColor', 'none');
hold on;
contour(x, y, log(likelihood3), 20, 'LineWidth', 1.5, 'LineColor', 'k');
colorbar;
colormap('turbo');

% Mark the theoretical maximum
p_hat = mean(Y1);
beta0_star = log(p_hat/(1-p_hat));
plot(beta0_star, 0, 'rx', 'MarkerSize', 15, 'LineWidth', 3);
text(beta0_star + 0.5, 0.5, sprintf('Theoretical max\n(\\beta_0=%.2f, \\beta_1=0)', beta0_star), ...
    'FontSize', 11, 'BackgroundColor', 'white', 'EdgeColor', 'k');

xlabel('\beta_0 (Intercept)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('\beta_1 (Slope on log K)', 'FontSize', 14, 'FontWeight', 'bold');
title({'Log-Likelihood Surface: Regime 1 Transition Model', ...
       sprintf('Weak identification: \\beta_0 and \\beta_1 are nearly unidentified')}, ...
       'FontSize', 13);
grid on;
set(gca, 'FontSize', 12);
hold off;
