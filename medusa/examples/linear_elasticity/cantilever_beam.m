clear
close all
format compact
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

cantilever_beam_data;
x = positions(1, :)';
y = positions(2, :)';
u = displacement(:, 1);
v = displacement(:, 2);

sxx = stress(:, 1);
syy = stress(:, 2);
sxy = stress(:, 3);

mises = sqrt(sxx.^2 - sxx.*syy + syy.^2 + 3*sxy.^2);

f = 1e5;
scatter(x+f*u, y+f*v, 15, mises/1000, 'filled');
title('Von Mises stress and deflection of catilever beam.');
daspect([1 1 1])
xlabel('$x$')
ylabel('$y$')
xlim([-1, 31])
ylim([-5, 5])
c = colorbar;
title(c, 'kPa');
caxis([0, inf])
colormap jet
grid on
box on

% Analytical solution
I = D^3 / 12;
asxx = I.^(-1).*P.*x.*y;
asyy = zeros(size(x));
asxy = (1/2).*I.^(-1).*P.*((1/4).*D.^2+(-1).*y.^2);

au = (1/24).*E.^(-1).*I.^(-1).*P.*y.*(3.*D.^2.*(1+nu)-...
    4.*(3.*L.^2+(-3).*x.^2+(2+nu).*y.^2));
av = (-1/24).*E.^(-1).*I.^(-1).*P.*(3.*D.^2.*(1+nu).*(L+(-1).*x)+...
    4.*(2.*L.^3+(-3).*L.^2.*x+x.^3+3.*nu.*x.*y.^2));

errInfDispl = norm([au-u; av-v], inf) / norm([au av], inf)
errInfStress = norm([asxx - sxx, asyy - syy, asxy - sxy], inf) / norm([sxx syy sxy], inf)
