clc;
clear;
% clf;

L = 220; % the length of cable is 220m;
n = 5; % number of modes used to represent cable;
x_d = 4; % damper location is 4m;
m = 85; % Unit: kg/m;
T = 7115*1e3; % Unit: N;
Diameter = 0.118; % Unit: m;
angle = 22.5; % Unit: degrees;
c = 0.05; % visous damping per unit;

numberofModes = (1:n);
w0 = sqrt(T/m).*(numberofModes)./(2*L);

[K,M,C,modes,w] = getMatrix(L, x_d, m, n, T, c);

