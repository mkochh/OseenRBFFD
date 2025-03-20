clear
close all
format compact
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

custom_rbf_data;
x = positions(1, :)';
y = positions(2, :)';

I = abs(x.^2+y.^2 - 1) < 1e-5;

analytical = sin(x).*sin(y);

err2 = norm(solution(I) - analytical(I))

err = norm(solution - analytical) / norm(analytical)

scatter(x, y, 10, solution, 'filled');
daspect([1 1 1])
title('$u(x,y,z)$');
xlabel('$x$')
ylabel('$y$')
xlim([-1 1])
ylim([-1 1])
colorbar
colormap jet
box on
grid on
