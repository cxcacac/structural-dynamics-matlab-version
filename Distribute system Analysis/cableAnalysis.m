% the state-space relations based on Li(2007)
% Vibration control of stay cables of the shandong binzhou rivers;
clc;
clear;
clf;

%% basic informations
L = 220; % the length of cable is 220m;
n = 5; % number of modes used to represent cable;
x_d = 4; % damper location is 4m, about 2% away from cable;
m = 85; % Unit: kg/m;
T = 7115*1e3; % Unit: N;
Diameter = 0.118; % Unit: m;
angle = 22.5; % Unit: degrees;
cRatio = 0.01; % damping ratio;

%% use mathematical method check if freq is right;
numberofModes = (1:n);
w0 = sqrt(T/m).*(numberofModes)./(2*L);

%% get K,M,C matrix;
[K,M,C,modes,freq] = getMatrix(L, x_d, m, n, T, cRatio);

%% get A,B,D matrix;
A = [zeros(n), eye(n); -inv(M)*K,-inv(M)*C];
% location 
ALocation = L/3;
BLocation = L/2; % mid point location
CLocation = 2*L/3;
DLocation = x_d; % damper location
ELocation = L/6; % excitation location

phiA = zeros(n,1); % phi(x_{L/2});
phiB = zeros(n,1);
phiC = zeros(n,1);
phiD = zeros(n,1);
phiE = zeros(n,1);

for i = 1:n
    phiA(i) = shapeFunction(ALocation, i, x_d, L);
    phiB(i) = shapeFunction(BLocation, i, x_d, L);
    phiC(i) = shapeFunction(CLocation, i, x_d, L);
    phiD(i) = shapeFunction(DLocation, i, x_d, L);
    phiE(i) = shapeFunction(ELocation, i, x_d, L);
end

B = [zeros(n,1); M\phiD];
D0 = [zeros(n); inv(M)];

%% get load vector;
t = (0:0.02:50);
dt = 0.02;
slots = length(t);
f1 = 30*sin(2*pi*freq(2).*t); % second freq;
f2 = 50*sin(2*pi*freq(1).*t); % first freq;
F1 = zeros(n, slots);
F2 = zeros(n, slots);
syms x;
for i = 1:n
    F1(i,:) = int(shapeFunction(x,i,x_d, L),x,0, L).*f1;
    F2(i,:) = int(shapeFunction(x,i,x_d, L),x,0, L).*f2;
end
u = 0; % control force
X0 = zeros(3*n,1); % initial state;
R = zeros(n,1); % wilson-seta parameters;

%% use NewMark-beta method to solve equation;
[x0,v0,a0] = NewMarkBeta(M,K,C,F1,dt,X0,phiD, 'free-vibration');
[x1,v1,a1] = NewMarkBeta(M,K,C,F1,dt,X0,phiD, 'passive');
[x2,v2,a2] = NewMarkBeta(M,K,C,F1,dt,X0,phiD, 'semi-active');

%% use Phi matrix to get cable response(Free Vibration)
A_d0 = phiA'*x0; % displacement at point A(L/3)
B_d0 = phiB'*x0; % displacement at point B(L/2)
C_d0 = phiC'*x0; % displacement at point C(2*L/3)
D_d0 = phiD'*x0; % displacement at point Damper(x_d)
E_d0= phiE'*x0; % displacement at point Excitation(L/6)

A_v0 = phiA'*v0; % veloctiy at point A(L/3)
B_v0 = phiB'*v0; % velocity at point B(L/2)
C_v0 = phiC'*v0; % velocity at point C(2*L/3)
D_v0 = phiD'*v0; % velocity at point Damper(x_d)
E_v0 = phiE'*v0; % velocity at point Excitation(L/6)

A_a0 = phiA'*a0; % veloctiy at point A(L/3)
B_a0 = phiB'*a0; % velocity at point B(L/2)
C_a0 = phiC'*a0; % velocity at point C(2*L/3)
D_a0 = phiD'*a0; % velocity at point Damper(x_d)
E_a0 = phiE'*a0; % velocity at point Excitation(L/6)

%% use Phi matrix to get cable response(Controlled Vibration)
A_d1 = phiA'*x1; % displacement at point A(L/3)
B_d1 = phiB'*x1; % displacement at point B(L/2)
C_d1 = phiC'*x1; % displacement at point C(2*L/3)
D_d1 = phiD'*x1; % displacement at point Damper(x_d)
E_d1= phiE'*x1; % displacement at point Excitation(L/6)

A_v1 = phiA'*v1; % veloctiy at point A(L/3)
B_v1 = phiB'*v1; % velocity at point B(L/2)
C_v1 = phiC'*v1; % velocity at point C(2*L/3)
D_v1 = phiD'*v1; % velocity at point Damper(x_d)
E_v1 = phiE'*v1; % velocity at point Excitation(L/6)

A_a1 = phiA'*a1; % acceleration at point A(L/3)
B_a1 = phiB'*a1; % acceleration at point B(L/2)
C_a1 = phiC'*a1; % acceleration at point C(2*L/3)
D_a1 = phiD'*a1; % acceleration at point Damper(x_d)
E_a1 = phiE'*a1; % acceleration at point Excitation(L/6)

%% use Phi matrix to get cable response(Passive Controlled Vibration)
A_d2 = phiA'*x2; % displacement at point A(L/3)
B_d2 = phiB'*x2; % displacement at point B(L/2)
C_d2 = phiC'*x2; % displacement at point C(2*L/3)
D_d2 = phiD'*x2; % displacement at point Damper(x_d)
E_d2= phiE'*x2; % displacement at point Excitation(L/6)

A_v2 = phiA'*v2; % veloctiy at point A(L/3)
B_v2 = phiB'*v2; % velocity at point B(L/2)
C_v2 = phiC'*v2; % velocity at point C(2*L/3)
D_v2 = phiD'*v2; % velocity at point Damper(x_d)
E_v2 = phiE'*v2; % velocity at point Excitation(L/6)

A_a2 = phiA'*a2; % acceleration at point A(L/3)
B_a2 = phiB'*a2; % acceleration at point B(L/2)
C_a2 = phiC'*a2; % acceleration at point C(2*L/3)
D_a2 = phiD'*a2; % acceleration at point Damper(x_d)
E_a2 = phiE'*a2; % acceleration at point Excitation(L/6)

%% make graph
% B_d1 = zeros(size(B_d0));
figure(1);
plot(t,B_d0,'b-', t,B_d1, 'k:', t, B_d2, 'r-');
