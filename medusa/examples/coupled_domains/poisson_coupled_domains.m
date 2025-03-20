clear
close all
format compact
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

in = hdf5read('poisson_coupled_domains.h5', '/inner/pos');
out = hdf5read('poisson_coupled_domains.h5', '/outer/pos');
sol = hdf5read('poisson_coupled_domains.h5', '/sol');
lam1 = hdf5read('poisson_coupled_domains.h5', '/', 'lam1');
lam2 = hdf5read('poisson_coupled_domains.h5', '/', 'lam2');

M = hdf5read('poisson_coupled_domains.h5', '/M');
M = sparse(M(:, 1), M(:, 2), M(:, 3));

r = [in; out];
scatter3(r(:,1), r(:,2), sol, 10, sol, 'filled')
title("Coupled domain poisson equations, $\lambda_{in} = $ "  + num2str(lam1) + ", $\lambda_{out} = $ "  + num2str(lam2))
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')

