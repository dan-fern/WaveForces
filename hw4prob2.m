%% PROBLEM 2
clear all
Hoprime = 6.3;
T = 11.4;
B = 15;
g = 9.81;
tidelevel = 0.6;
seaslope = 1/100;
L0 = g.*T.^2./(2*pi);
k = 2*pi/L0;

d = 10.1;
kd = k*d;
hprime = 7.1;
d = 5.6;
hcrit = 3.4;

Lo = 1.56.*T.^2;

Bostar = 0.052*(Hoprime/Lo).^-0.38*exp(20*(seaslope).^1.5);
B1star = 0.63*exp(3.8*seaslope);
Bstarmax = max([1.65, 0.53*(Hoprime./Lo).^-0.29*exp(2.4*seaslope)]);

Ksi = (tanh(kd) + kd.*(1-tanh(kd).^2)).^-.5;
Ks = Ksi;

Hmax = ([Bostar*Hoprime+B1star*d, Bstarmax*Hoprime, 1.8*Ks*Hoprime])

%% method 2 
syms  h
h_30 = solve((((((2*pi)/30)*(Hoprime/Lo)*(tanh(wavenumber(T,d,g)*d) + wavenumber(T,d,g)*d.*(1-tanh(wavenumber(T,d,g)*d).^2)).^-.5)^.5)*Lo)-d == 0,d,'Real', true);

syms h
h_50 = solve((((((2*pi)/50)*(Hoprime/Lo)*(tanh(wavenumber(T,d,g)*d) + wavenumber(T,d,g)*d.*(1-tanh(wavenumber(T,d,g)*d).^2)).^-.5)^.5)*Lo)-d == 0,d,'Real', true);

Ks_50 = (tanh(wavenumber(T,h_50,g)*h_50) + wavenumber(T,h_50,g)*h_50.*(1-tanh(wavenumber(T,h_50,g)*h_50).^2)).^-.5;

d = 10.1;
B = ((2*sqrt(3))/sqrt(2*pi*Hoprime/Lo))*(d/Lo);
C50 = (Ks_50*(h_50/Lo)^(3/2))*(sqrt(2*pi*(Hoprime/Lo)*Ks_50)-2*sqrt(3)*(h_50/Lo));
C = (C50/(sqrt(2*pi*(Hoprime/Lo))))*((Lo/d)^(3/2));

syms Ks_temp
Ks = solve(Ks_temp*(sqrt(Ks_temp) - B) - C == 0,Ks_temp,'Real', true );

H = Ks * Hoprime

%% Method 3

Ksi = (tanh(kd) + kd.*(1-tanh(kd).^2)).^-.5;
Ks = Ksi +0.0015*(d/Lo)^(-2.87)*(Hoprime/Lo)^1.27;
H = Ks*Hoprime



