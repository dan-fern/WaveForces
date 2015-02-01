% CE 645 Homework 2 Problem 5

clear all
close all
home

L = 95.6; % m
T = 9; % s
H = 4; % m
h = 15; % m
D = 0.3; % m
l = 60; % m
U = 2.3; % m/s
k = (2*pi)/L; % 1/m
omega = (2*pi)/T; % 1/s
rho = 1040; % kg/m^3
g = 9.81; % m/s^2

% use coefficients given in class, can't justify coefficients for inertial
% and lift given in D & D
C_d = 0.6;
C_m = 1.8;

% from D & D
% C_m = 3.29;
% C_l = 4.493;

E = (1/2)*rho*g*H^2;
n = (1/2)*(1+(2*k*h)/sinh(2*k*h));

% assume point zero is over the pipe
x = 0.15; % m
t = linspace(0,T,18); % s

for j = 1:length(t)
    F_d(:,j) = C_d*D*n*E*cos(k*x-omega*t(j))*abs(cos(k*x-omega*t(j)))/l; % N/m
    F_i(:,j) = C_m*(pi/4)*D^2*E*tanh(k*h)*sin(k*x-omega*t(j))/l; % N/m
    
    % don't include lift because it's an upward force, not a horizontal
    % force
    
    % u(:,j) = (H/2)*(g*T/L)*cosh(2*pi*((15-0.3)+h)/L)/cosh(2*pi*h/L)*sin(k*x-omega*t(j));
    % F_l(:,j) = C_l*rho*D*(1/2)*u(j)^2/l;
    
    F_t = F_d + F_i;
end

figure(1)
subplot(3,1,1)
plot(t,F_d)
title('Wave Forces on Exposed Pipeline')
ylabel('Drag Force per unit Length (N/m)')
subplot(3,1,2)
plot(t,F_i)
ylabel('Inertial Force per unit Length (N/m)')
subplot(3,1,3)
plot(t,F_t)
ylabel('Total Forces per unit Length (N/m)')
xlabel('Time (s)')

F_t_max = max(F_t)
F_t_min = min(F_t)
