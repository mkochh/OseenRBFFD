clear
close all
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

run inside_by_oversampling_data.m

scatter(positions_leaking(1, :), positions_leaking(2, :), 10, 'filled');
title('Leaking domain');
axis('equal')
xlabel('$x$')
ylabel('$y$')

figure()
scatter(positions(1, :), positions(2, :), 10, 'filled');
title('Domain filled via oversampling');
axis('equal')
xlabel('$x$')
ylabel('$y$')