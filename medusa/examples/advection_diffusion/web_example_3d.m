clear
close all
format compact
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

web_example_3d_data;
x = positions(1, :);
y = positions(2, :);
z = positions(3, :);

if 0

[X, Y, Z] = meshgrid(deal(linspace(-1, 1, 200)));
V = griddata(x, y, z, solution, X, Y, Z);

V((X-0.1).^2 + (Y-0.1).^2 + (Z-0.1).^2 < 0.3^2) = nan;

sx = 0;
sy = 0;
sz = 0;
h = slice(X, Y, Z, V, sx, sy, sz);
set(h, 'edgecolor', 'none')
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
title('$u(x, y, z)$')
xlim([-1 1])
ylim([-1 1])
zlim([-1 1])
daspect([1 1 1])
grid on
box on
xticks(-1:0.5:1)
yticks(-1:0.5:1)
colormap jet
colorbar
caxis([0, inf])

% exportf(gcf, 'adv-dif-3d-slice.png', '-m3')

end

if 0

sx = 0;
sy = 0;
sz = 0;
contourslice(X, Y, Z, V, sx, sy, sz, 20);
view(3)
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
title('$u(x, y, z)$')
xlim([-1 1])
ylim([-1 1])
zlim([-1 1])
daspect([1 1 1])
grid on
box on
xticks(-1:0.5:1)
yticks(-1:0.5:1)
colormap jet
colorbar

end

if 1

n = 100;
[X, Y, Z] = meshgrid(linspace(-1, 0, n), linspace(-1, 1, 2*n), linspace(-1, 1, 2*n));
V = griddata(x, y, z, solution, X, Y, Z);

V((X-0.1).^2 + (Y-0.1).^2 + (Z-0.1).^2 < 0.3^2) = nan;

vals = [0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03];
nval = length(vals);
p = zeros(size(vals));
for i = 1:nval
    isosurface(X, Y, Z, V, vals(i))
end

colormap jet
box on; grid on; daspect([1 1 1])
view([53.7000   14.8000]); lighting gouraud; % camlight; % alpha(0.3);
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
title('$u(x, y, z)$')
xlim([-1 0])
ylim([-1 1])
zlim([-1 1])
colorbar

% exportf(gcf, 'adv-dif-3d-iso.png', '-m3')

end