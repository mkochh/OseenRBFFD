clear
close all
format compact
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

custom_stencils_visibility_data;
x = positions(1, :);
y = positions(2, :);
% z = positions(3, :);
N = length(x);
s = length(stencils_closest);

sc = reshape(stencils_closest, [s/N N])';
sv = reshape(stencils_visibility, [s/N N])';

target = [0.941; 0.317];
[~, I] = min(sqrt(sum((positions - target).^2, 1)));

figure
hold on; box on; grid on; daspect([1 1 1]);
plot(x, y, 'x');
scatter(x(sc(I, :)+1), y(sc(I, :)+1), 15, 'k', 'filled');
scatter(x(I), y(I), 15, 'r', 'filled');
daspect([1 1 1])
title('closest stencil');
xlabel('$x$')
ylabel('$y$')

figure
hold on; box on; grid on; daspect([1 1 1]);
plot(x, y, 'x');
vis = sv(I, :)+1; vis = vis(vis > 0);
scatter(x(vis), y(vis), 15, 'k', 'filled');
scatter(x(I), y(I), 15, 'r', 'filled');
daspect([1 1 1])
title('visibility stencil');
xlabel('$x$')
ylabel('$y$')