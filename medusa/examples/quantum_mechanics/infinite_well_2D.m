clear
close all
format compact
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

pos = hdf5read('infinite_well_2D.h5', '/pos');
rsol = hdf5read('infinite_well_2D.h5', '/rsol');
csol = hdf5read('infinite_well_2D.h5', '/csol');
x = pos(:, 1);
y = pos(:, 2);

f1 = figure(1);
scatter(x, y, 10, rsol, 'filled');
title("$\Re{\Psi(x,y)}$");
axis('equal')
xlabel('$x$')
ylabel('$y$')
xlim([0 1])
ylim([0 1])
colorbar
colormap jet
saveas(f1, 'infinite_well_2D_rsol.pdf')

f2 = figure(2);
scatter(x, y, 10, csol, 'filled');
title("$\Im{\Psi(x,y)}$");
axis('equal')
xlabel('$x$')
ylabel('$y$')
xlim([0 1])
ylim([0 1])
colorbar
colormap jet
saveas(f2, 'infinite_well_2D_csol.pdf')