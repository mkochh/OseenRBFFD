clear
close all
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

run parametric_domain_2D_data.m
x = positions(1, :);
y = positions(2, :);
scatter(x, y, 10, solution, 'filled');
title('$u(x, y)$');
axis('equal')
xlabel('$x$')
ylabel('$y$')
colorbar
colormap jet