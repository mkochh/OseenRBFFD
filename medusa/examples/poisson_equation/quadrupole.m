clear
close all
format compact
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

domain = h5read('quadrupole.h5', '/domain/pos');
sol = h5read('quadrupole.h5', '/sol');
scatter(domain(:,1), domain(:,2), 12, sol, 'filled');

title('Four charged circles');
xlabel('$x$')
ylabel('$y$')
axis equal;
xlim([0 1])
ylim([0 1])
colorbar
colormap(b2r(-1,1))