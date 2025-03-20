clear
close all
format compact
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

point_contact3d_data;
x = positions(1, :)';
y = positions(2, :)';
z = positions(3, :)';
u = displacement(:, 1);
v = displacement(:, 2);
w = displacement(:, 3);

sxx = stress(:, 1);
syy = stress(:, 2);
szz = stress(:, 3);
sxy = stress(:, 4);
sxz = stress(:, 5);
syz = stress(:, 6);

mises = sqrt(0.5*((sxx - syy).^2 + (syy - szz).^2 + (szz - sxx).^2 ...
                + 6*(sxy.^2 + sxz.^2 + syz.^2)));

f = 0;
scatter3(x+f*u, y+f*v, z+f*w, 10, mises, 'filled');
view([81.3000   14.0000])
title('Von Mises stress and displacements of point contact.');
daspect([1 1 1])
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
% xlim([-1.2, 1.2])
% ylim([-1.2, 0.1])
c = colorbar;
caxis([0, inf])
colormap jet
grid on
box on

