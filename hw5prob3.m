clear all
close all
home

g = 9.81; 
T = 12; 
w = 2 * pi / T;
d = 4; 
D = 0.3; 
x = D/2;
H = 0.78 * d; 
L = 73.8;
rho = 1040;
mu = 0.001002;
nu = 0.000001004;
Cd = 0.65;
Cm = 1.6;
k = 2 * pi / L;

u = (H * g * T / 2 / L) * cosh(2*pi*(H/2+d)/L) / cosh(2*pi*d/L);
Re = rho * u * D / mu;
KC = u * T / D;
beta = D^2 / nu / T;

n = 0.5 * (1 + (2*k*d) / sinh(2*k*d)); 
E = 1/8 * rho * g * H^2;

t = linspace(0,T,1000);

Fd = Cd * D * n * E .* cos(k*x-w*t) .* abs(cos(k*x-w*t));
Fi = Cm * pi * D * E * (D/H) .* tanh(k*d) .* sin(k*x-w*t);
Ft = Fd + Fi;

maxDrag = max(abs(Fd))
maxInertia = max(abs(Fi))
maxTotal = max(abs(Ft))

figureHandle = figure('Position',[25,55,900,900]);
plot(t,Fd/1000,'-r', 'LineWidth', 2); hold on;
plot(t,Fi/1000,'-b', 'LineWidth', 2); hold on;
plot(t,Ft/1000,'-m', 'LineWidth', 2); hold on;
axis([0 12 -5 5]); set(gca,'YTickLabel',sprintf('%.1f|',get(gca,'YTick')));
title('Drag, Inertia, and Total Forces over one Wave Period')
ylabel('Wave Force, kN'); xlabel('Time, s');
legend('Drag Forces', 'Inertial Forces', 'Total Forces', 'location', 'best');

Ftmax = 1070.74 * cos(-.566*t) .* abs(cos(-.566*t)) + 2003.64 * sinh(-.566*t);
maxTotal2 = max(abs(Ftmax))