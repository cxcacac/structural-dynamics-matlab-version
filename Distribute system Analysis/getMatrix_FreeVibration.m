clc;
clear;
% clf;
% parpool('local',4);
L = 220; % the length of cable is 220m;
n = 5; % number of modes used to represent cable;
dLocation = 4; % damper location is 4m;
m = 85; % Unit: kg/m;
T = 7115*1e3; % Unit: KN;
Diameter = 0.118; % Unit: m;
angle = 22.5; % Unit: degrees;

% from Pacheco, phi = sin(pi*i*x/L)
M = zeros(n);
K = zeros(n);

for i = 1:n
    for j = 1:n
        M(i,j) = m*L/2*(i==j);
    end
end

for i = 1:n
    for j = 1:n
        K(i,j) = T*pi^2*i^2/(2*L)*(i==j);
    end
end

[modes, w2] = eig(K,M);

w = diag(sqrt(w2))./(2*pi);