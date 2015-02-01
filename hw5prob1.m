clear all
close all
home

gamma = 64; 
Cs = 2;
w = 50;
h = 6;
H = linspace(0.1,6,60);
%H = 0.78 * h;
eta = H / 2;

FB = h + eta - 8;
index = find(FB < 0);
FB(index) = 0;

Fslam = 1/2 * gamma * Cs .* (H.^2/h) .* (2 - (1.22*H.^2/H)).^(-2) .* FB * w;

figureHandle = figure('Position',[25,55,1080,720]);
plot(H, Fslam, 'LineWidth', 2); 
%axis([0 6 0 600])
title('1g.) Wave Slam Equation vs Wave Height (Tamiczek et al. 2014)')
ylabel('Wave Slam Force, lbf'); xlabel('Wave Height, ft');
%legend('Drag Forces', 'location', 'best');
