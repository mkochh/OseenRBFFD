clear
close all
format compact
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

point_contact2d_data;
x = positions(1, :)';
y = positions(2, :)';
u = displacement(:, 1);
v = displacement(:, 2);

sxx = stress(:, 1);
syy = stress(:, 2);
sxy = stress(:, 3);

mises = sqrt(sxx.^2 - sxx.*syy + syy.^2 + 3*sxy.^2);

f = 5e6;
scatter(x+f*u, y+f*v, 15, mises/1000, 'filled');
title('Von Mises stress and displacements of point contact.');
daspect([1 1 1])
xlabel('$x$')
ylabel('$y$')
xlim([-1.2, 1.2])
ylim([-1.2, 0.1])
c = colorbar;
title(c, 'kPa');
caxis([0, inf])
colormap jet
grid on
box on

