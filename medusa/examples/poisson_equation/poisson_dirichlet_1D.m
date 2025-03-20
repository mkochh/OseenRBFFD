clear
close all
format compact
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

poisson_dirichlet_1D_data;
plot(positions(1,:), solution,'.')
title('$u(x)$');
xlabel('$x$')
ylabel('$u$')
xlim([0 1])
ylim([0 1])