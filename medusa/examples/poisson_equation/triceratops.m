clear
close all
format compact
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

solution = h5read('triceratops.h5','/solution');
positions = h5read('triceratops_domain.h5', '/domain/pos');

x = positions(:, 1);
y = positions(:, 2);
z = positions(:, 3);
scatter3(x, y, z, 10, solution, 'filled');
title('$u(x,y,z)$');
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
daspect([1 1 1])
colorbar
colormap jet
caxis([0, inf])