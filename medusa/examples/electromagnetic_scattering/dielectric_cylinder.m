clear all;
close all;

filename = "dielectric_cylinder.h5";

r1 = h5read(filename, "/outer_dom/pos");
r2 = h5read(filename, "/cyl_dom/pos");
t1 = h5read(filename, "/outer_dom/types");
t2 = h5read(filename, "/cyl_dom/types");
inter = h5read(filename, "/outer_inter_nodes");

lambda = h5readatt(filename, "/", "wavelength");
rad = h5readatt(filename, '/config/', "geometry.r");
a = h5readatt(filename, '/config/', "geometry.a");
d = h5readatt(filename, '/config/', "PML.d");
epsr = h5readatt(filename, '/config/', "case.epsr");

phi = linspace(0, 2*pi, 100);
cilx = a*cos(phi);
cily = a*sin(phi);

k = (2*pi/lambda);

r = [r1; r2];
t = [t1; t2];
EIK = exp(1.0i*k*r(:, 1));

rsol = h5read(filename, "/rsol");
csol = h5read(filename, "/csol");

tsol = rsol+1.0i*csol;
tsol(inter+1) = tsol(inter+1)+ EIK(inter+1);
rsol = real(tsol);
csol = imag(tsol);

f1 = setfig('a1', [900, 900]);
title("Imaginary part of $E_z$");
scatter(r(:, 1), r(:, 2), 40, csol, 'filled', 'o')
text(0-a/2, 0, "$\varepsilon = $ "+epsr,'FontSize', 25)
annotation('arrow', [0.7 0.6], [0.51 0.51])
hold on
plot(cilx, cily, 'k', 'LineWidth', 2)
colorbar

f2 = setfig('a2', [900, 900]);
title("Real part of $E_z$");
scatter(r(:, 1), r(:, 2), 40, rsol, 'filled', 'o')
hold on
plot(cilx, cily, 'k', 'LineWidth', 2)
text(0-a/2, 0, "$\varepsilon = $ "+epsr,'FontSize', 25)
annotation('arrow', [0.7 0.6], [0.51 0.51])
colorbar
