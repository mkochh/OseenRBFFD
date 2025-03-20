clear
close all

set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');


fn='wave_equation_2D.h5'
pause on
pos = h5read(fn,'/pos');

for i=1:999
   num = int2str(i);
   eval(['E',num,'=','h5read(''',fn,''',','''/step',num,'/E''',');']);
end

for i=1 : 999

    num = int2str(i);
    var=eval(['E',num]);

    plot3(pos(:,1),pos(:,2),var,'.-');
    tri = delaunay(pos(:,1),pos(:,2));
    plot(pos(:,1),pos(:,2),'.');
    [r,c] = size(tri);

    h = trisurf(tri,pos(:,1),pos(:,2),var);
    axis vis3d;
    caxis manual;
    caxis([-0.3 0.3]);
    axis equal
    xlim([-1 1]);
    ylim([-1 1]);
    zlim([-0.4 0.4]);
    axis on;
    l = light('Position',[-20 -15 29]);
    lighting phong
    shading faceted
    xlabel('$x$');
    ylabel('$y$');
    zlabel('$z$')
    pause(0.001);
end