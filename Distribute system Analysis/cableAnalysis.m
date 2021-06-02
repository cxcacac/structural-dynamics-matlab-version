% the state-space relations based on Li(2007)
% Vibration control of stay cables of the shandong binzhou rivers;

clc;
clear;
clf;
%% basic informations
L = 220; % the length of cable is 220m;
n = 5; % number of modes used to represent cable;
x_d = 4; % damper location is 4m;
m = 85; % Unit: kg/m;
T = 7115*1e3; % Unit: N;
Diameter = 0.118; % Unit: m;
angle = 22.5; % Unit: degrees;
c = 0.1; % visous damping per unit;

%% check if freq is right;
numberofModes = (1:n);
w0 = sqrt(T/m).*(numberofModes)./(2*L);

%% get K,M,C matrix;
[K,M,C,modes,freq] = getMatrix(L, x_d, m, n, T, c);

%% get A,B,D matrix;
A = [zeros(n), eye(n); -inv(M)*K,-inv(M)*C];
phiD = zeros(n,1);

for i = 1:n
    phiD(i) = shapeFunction(x_d, i, x_d, L);
end

B = [zeros(n,1); inv(M)*phiD];
D0 = [zeros(n); inv(M)];

t = (0:0.02:50);
dt = 0.02;
slots = length(t);
f = 30*sin(2*pi*freq(1).*t); % load

F = zeros(n, slots);
syms x;
for i = 1:n
    F(i,:) = int(shapeFunction(x,i,x_d, L),x,0, L).*f;
end
u = 0;

X0 = zeros(3*n,1);
R = zeros(n,1);

%% use wilson-theta method to solve equation;
F0 = [F, zeros(n,1)];

[q,dq,dq] = wilsonSeta(M,K,C,dt,F0,X0,K,R);
%% get cable response
 
d = phiD'*q; % at damper location;

figure(1);
plot(t,d);
