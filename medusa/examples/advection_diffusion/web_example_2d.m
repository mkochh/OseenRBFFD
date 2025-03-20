clear
close all
format compact
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

web_example_2d_data;
x = positions(1, :);
y = positions(2, :);

if 1

scatter(x, y, 5, solution, 'filled');
title('$u(x, y)$');
xlabel('$x$')
ylabel('$y$')
zlabel('$u(x, y)$')
xlim([-1 1])
ylim([-1 1])
axis square
grid on
box on
xticks(-1:0.5:1)
yticks(-1:0.5:1)
colormap jet
colorbar

% exportf(gcf, 'adv-dif-2d-scat.png', '-m3')

end

if 0

tri = delaunay(x,y);
trisurf(tri, x, y, solution, 'EdgeColor', 'None');
lighting phong
shading interp
view([145.3000   40.4000])
xlabel('$x$')
ylabel('$y$')
zlabel('$u(x, y)$')
title('$u(x, y)$');
xlim([-1 1])
ylim([-1 1])
axis square
grid on
box on
xticks(-1:0.5:1)
yticks(-1:0.5:1)
colormap jet

% exportf(gcf, 'adv-dif-2d-surf.png', '-m3')

end