clear
close all
format compact
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

custom_operators_biharmonic_data;
x = positions(1, :);
y = positions(2, :);
% z = positions(3, :);

analytical = 1/2*(1+sum(positions.^2))';

err = norm(solution - analytical) / norm(analytical)

scatter(x, y, 10, solution, 'filled');
% scatter3(x, y, z, 10, solution, 'filled');
daspect([1 1 1])
title('$u(x,y,z)$');
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
xlim([-1 1])
ylim([-1 1])
zlim([-1 1])
colorbar
colormap jet
