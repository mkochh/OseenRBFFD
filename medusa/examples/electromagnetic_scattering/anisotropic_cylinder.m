clear
close all
format compact
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

r1 = hdf5read('anisotropic_cylinder.h5', '/inner/pos');
r2 = hdf5read('anisotropic_cylinder.h5', '/outer/pos');
rsol = hdf5read('anisotropic_cylinder.h5', '/rsol');
csol = hdf5read('anisotropic_cylinder.h5', '/csol');

r = [r2; r1];
insol = rsol(length(r2)+1:end) + 1.0j*csol(length(r2)+1:end); %inside the scatterer
outsol = rsol(1:length(r2)) + 1.0j*csol(1:length(r2)); %outiside the scatterer

f1 = setfig('b1', [840 840]);
scatter(r1(:, 1), r1(:, 2), 50, abs(insol),  'filled')
xlabel('$x$')
ylabel('$y$')
xlim([-0.3 0.3])
ylim([-0.3 0.3])
title("Field $v$  inside the cylinder")
axis equal
colorbar
exportf(f1, 'InsideField.png')

f2 = setfig('b2', [840 840]);
scatter(r2(:, 1), r2(:, 2), 10, abs(outsol),  'filled')
title("Scattered field $u^s$  outside the cylinder")
xlabel('$x$')
ylabel('$y$')
xlim([-0.8 0.8])
ylim([-0.8 0.8])
axis equal
colorbar
exportf(f2, 'OutsideField.png')

f3 = setfig('b3', [840 840]);
scatter(r1(:, 1), r1(:, 2), 50, real(insol),  'filled')
xlabel('$x$')
ylabel('$y$')
xlim([-0.3 0.3])
ylim([-0.3 0.3])
title("Field $\Re{v}$  inside the cylinder")
axis equal
colorbar
exportf(f3, 'RealInsideField.png')

f4 = setfig('b4', [840 840]);
scatter(r2(:, 1), r2(:, 2), 10, real(outsol),  'filled')
title("Scattered field $\Re{u^s}$  outside the cylinder")
xlabel('$x$')
ylabel('$y$')
xlim([-0.8 0.8])
ylim([-0.8 0.8])
axis equal
colorbar
exportf(f4, 'RealOutsideField.png')

f5 = setfig('b1', [840 840]);
scatter(r1(:, 1), r1(:, 2), 50, imag(insol),  'filled')
xlabel('$x$')
ylabel('$y$')
xlim([-0.3 0.3])
ylim([-0.3 0.3])
title("Field $\Im{v}$  inside the cylinder")
axis equal
colorbar
exportf(f5, 'ImagInsideField.png')

f6 = setfig('b2', [840 840]);
scatter(r2(:, 1), r2(:, 2), 10, imag(outsol),  'filled')
title("Scattered field $\Im{u^s}$  outside the cylinder")
xlabel('$x$')
ylabel('$y$')
xlim([-0.8 0.8])
ylim([-0.8 0.8])
axis equal
colorbar
exportf(f6, 'ImagOutsideField.png')
