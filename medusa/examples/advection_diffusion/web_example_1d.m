clear
close all
format compact
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

web_example_1d_data;
x = positions(1, :);
scatter(x, solution, 5, solution, 'filled');
title('$u(x)$');
grid on
box on
colormap jet
xlabel('$x$')
ylabel('$u(x)$')
xlim([-1.1 1.1])
ylim([-0.005 0.04])
yticks(0:0.01:0.04)


exportf(gcf, 'adv-dif-1d.png', '-m3')