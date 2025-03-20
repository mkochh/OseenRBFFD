clear
close all
format compact
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

poisson_mixed_2D_irregular_data;
x = positions(1, :);
y = positions(2, :);
scatter(x, y, 10, solution, 'filled');
title('$u(x,y)$');
xlabel('$x$')
ylabel('$y$')
colorbar
colormap jet