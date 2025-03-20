clear
close all
format compact
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

poisson_neumann_2D_data;
x = positions(1, :);
y = positions(2, :);
scatter(x, y, 10, solution, 'filled');
title('$u(x,y)$');
xlabel('$x$')
ylabel('$y$')
axis('equal')
xlim([0 1])
ylim([0 1])
colorbar
colormap jet