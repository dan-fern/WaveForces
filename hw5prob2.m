clear all
close all
home

h = 6;
gamma = 64; 
g = 32.2;
Cd = 2.0;
L = 16;
W = 6; 

Fr = linspace(0.5,1,100);
v = Fr * sqrt(g * h);

Fdrag = 1/2 * Cd * gamma / g * v.^2 * L * W;

figureHandle = figure('Position',[25,55,1080,720]);
plot(Fr, Fdrag, 'LineWidth', 2); 
axis([0.5 1 0.5 40000]);
set(gca,'YTickLabel', sprintf('%.0f|',get(gca,'YTick')));
title('2b.) Drag Force vs Froude Number using Quadratic Drag Law')
ylabel('Drag Force, lbf'); xlabel('Froude Number');
%legend('Drag Forces', 'location', 'best');
