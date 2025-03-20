clear
close all
format compact
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

cahnHilliard_visualization_data;
x = positions(1, :);
y = positions(2, :);

xBounds = [min(x), max(x)];
yBounds = [min(y), max(y)];

corners = [-border, -border];
for i = 1:4
    corners = cat(1, corners,  [-corners(end, 2), corners(end, 1)]);
end
    
figure(1)
for i = 1:4
    subplot(2,2,i)
    hold on
    scatter(x, y, 10, variables(:, i), 'filled');
    plot(corners(:, 1), corners(:, 2), 'k', 'LineWidth', 1)
    hold off
    caxis([-1.1, 1.1]);
    daspect([1 1 1])
    title(sprintf('$c(x, y), t=%.2f$', times(i)));
    xlabel('$x$')
    ylabel('$y$')
    xlim(xBounds)
    ylim(yBounds)
end
sgtitle("Time evolution of Cahn-Hilliard equation")

%colorbar
h = colorbar('southoutside');
set(h, 'Position', [.1 .88 .8 .02]);
colormap jet

set(gcf,'Position',[0 0 700 750])
set(gcf,'color','w');
ha=get(gcf,'children');
set(ha(3),'position',[.55 .1 .35 .35])
set(ha(4),'position',[.1 .1 .35 .35])
set(ha(5),'position',[.55 .5 .35 .35])
set(ha(6),'position',[.1 .5 .35 .35])
