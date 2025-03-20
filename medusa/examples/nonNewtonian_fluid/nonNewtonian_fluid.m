clear 
close all
format compact
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');  

nonNewtonian_fluid_data;
x = positions(1, :)';
y = positions(2, :)';
hold on 
quiver(x, y, u(:,1), u(:,2), 5);
scatter(x, y, 10, T, 'filled');
title('$u(x,y)$');
axis('equal')
xlabel('$x$')
ylabel('$y$')
xlim([0 1])
ylim([0 1])
colorbar
colormap jet

Ra
Pr
