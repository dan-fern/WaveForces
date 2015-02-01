clear all
close all
home

g = 9.81; 
rho = 1040;
T = 12; 
w = 2 * pi / T;
d = 4; 
H = 0.78 * d; 
L = 73.8;
k = 2 * pi / L;
Cp = 1.6;

Ffema = 1.1 * Cp * rho * g * d^2 + 1.9 * rho * g * d^2;

lambda1 = 1;
lambda2 = 1;
eta = 0.75 * (1 + lambda1) * H;
hb = H / 0.78;
alpha1 = 0.6 + 1/2*((4*pi*d/L)/sinh(4*pi*d/L))^2;                    
alpha2 = (hb-d)/(3*hb) .* (H./d)^2;
%alpha2b = 2*d./H;
alpha3 = 1 - (hb/d) * (1 - 1/cosh(2*pi*d/L));

z = 0;          
p1 = 1/2 * (1 + 1) * (alpha1*lambda1 + alpha2*lambda2) * rho * g * H;
p2 = p1 / cosh(k*d);          
p3 = alpha3 * p1;
p4 = 0; 
P = 0.5 .* ((p1 + p3) .* d + (p1 + p4) .* eta); 

%%
Bm = 2;

del11 = 0.93 * (Bm/L - 0.12) + 0.36 * (0.4 - d/hb);
del22 = -.36 * (Bm/L - 0.12) + 0.93 * (0.4 - d/hb);

if del11 > 0 del1 = 15 * del11; 
else del1 = 20 * del11; end
if del22 > 0 del2 = 3 * del22; 
else del2 = 4.9 * del22; end

if del2 > 0 alphaIB = 1/(cosh(del1) * sqrt(cosh(del2))); 
else alphaIB = cos(del2) / cosh(del1); end
if H/d > 2 alphaIH = 2; 
else alphaIH = H/d; end

alphaI = alphaIH * alphaIB;

if alphaI > alpha2 alphastar = slphaI; 
else alphastar = alpha2; end

p1new = 1/2 * (1 + 1) * (alpha1*lambda1 + alphastar) * rho * g * H;
Pnew = 0.5 .* ((p1new + p3) .* d + (p1new + p4) .* eta); 

