clear
close all
format compact
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% h5disp('poisson_dirichlet_3D_irregular.h5')
solution = h5read('poisson_dirichlet_3D_irregular.h5','/solution');
positions = h5read('poisson_dirichlet_3D_irregular.h5','/positions');

x = positions(:, 1);
y = positions(:, 2);
z = positions(:, 3);
scatter3(x, y, z, 10, solution, 'filled');
title('$u(x,y,z)$');
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
% xlim([0 1])
% ylim([0 1])
% ylim([0 1])
colorbar
colormap jet