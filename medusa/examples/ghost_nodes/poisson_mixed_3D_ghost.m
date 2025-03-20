clear
close all
format compact
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

poisson_mixed_3D_ghost_data;
I = types ~= 0;  % ignore ghost nodes
x = positions(1, I);
y = positions(2, I);
z = positions(3, I);

scatter3(x, y, z, 10, solution(I), 'filled');
view([42.2 26.4])
title('$u(x, y, z)$');
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
o = 1.1;
xlim([-o o])
ylim([-o o])
zlim([-o o])
daspect([1 1 1])
colorbar
colormap jet