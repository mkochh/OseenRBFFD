close all
clear

model = '../../test/testdata/bunny.off';
file = 'bunny_poisson.h5';

[f, v] = read_off(model);
T = triangulation(f, v);

pos = h5read(file, '/pos');
sol = h5read(file, '/sol');
x = pos(:, 1);
y = pos(:, 2);
z = pos(:, 3);

I = y > -20;
x = x(I);
y = y(I);
z = z(I);
sol = sol(I);

figure;
hold on
box on
grid on
daspect([1 1 1])
view(3);
gray = 0.5*[1 1 1];
trisurf(T, 'EdgeColor', gray, 'FaceColor', 'none');
scatter3(x, y, z, 3, sol,'filled');
colorbar
colormap jet
title('Cross-section of the solution')

