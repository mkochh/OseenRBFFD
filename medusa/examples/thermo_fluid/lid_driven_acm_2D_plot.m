clear 
close all
format compact
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');  

lid_driven_acm_2D_data;
x = positions(1, :)';
y = positions(2, :)';
quiver(x, y, u(:,1), u(:,2));
title('$u(x,y)$');
axis('equal')
xlabel('$x$')
ylabel('$y$')
xlim([0 1])
ylim([0 1])
colorbar
colormap jet