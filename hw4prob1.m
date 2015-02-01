clear all
close all
home

g = 9.81;  
rho = 1040;                 
kd = [.5 1 2];
d = 5; 

k = kd / d;          
L = 2 * pi ./ k;                
T = sqrt((2.*pi.*L)./(g.*tanh(2.*pi.*d./L)));

Hmax = [2.34 2.34 1.31]; 
hb = Hmax ./ 0.78;
eta = [3.51 3.51 1.98];
alpha1 = [0.96 0.75 0.61];                    
alpha2 = (hb-d) ./ (3*hb) .* (Hmax./d).^2;
%alpha2b = 2*h./Hmax;

% d) Dynamic Pressure at SWL
z = 0;          
p1_GODA = 1/2 .* (1 + 1) .* (alpha1 + alpha2) .* rho .* g .* Hmax;
p1_LWT = (cosh(kd+k.*z)./cosh(kd)) .* rho .* g .* Hmax; 
compare_d = p1_GODA - p1_LWT; 

% e) Dynamic Pressure at the seafloor
z = -5;    
p2_GODA = p1_GODA ./ cosh(kd);          
p2_LWT = (cosh(kd+k.*z)./cosh(kd)) .* rho .* g .* Hmax; 
compare_e = abs(p2_GODA - p2_LWT);
         
% f) Total Force for Goda/LWT minus hydrostatic

alpha3 = 1 - (d ./ d) * (1 - 1./cosh(kd));
p3_GODA = alpha3 .* p1_GODA;
p4_GODA = 0;
P = 0.5 .* ((p1_GODA + p2_GODA) .* d + (p1_GODA + p4_GODA) .* eta); 
Fmax = rho .* g .* (((4.*d).^2 + (Hmax).^2) ./ 8 ...
     + d .* (Hmax)./2.*(tanh(kd)./kd) - d.^2); 
Fdiff = abs(P - Fmax);
hydrostatic = rho .* g .* d;
compare_f = Fdiff ./ hydrostatic .* 100;
