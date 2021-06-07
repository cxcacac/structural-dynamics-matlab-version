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
cRatio = 0.02; % damping ratio;

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

B = [zeros(n,1); inv(M)*phiE];
D0 = [zeros(n); inv(M)];

t = (0:0.02:50);
dt = 0.02;
slots = length(t);
f = 30*sin(2*pi*freq(2).*t); % load
% f = 50*sin(2*pi*freq(1).*t); % load first freq;
F = zeros(n, slots);
syms x;
for i = 1:n
    F(i,:) = int(shapeFunction(x,i,x_d, L),x,0, L).*f;
end
u = 0; % control force
X0 = zeros(3*n,1); % initial state;
R = zeros(n,1); % wilson-seta parameters;

%% use NewMark-beta method to solve equation;
cn = length(M); % dof;
n = length(F); % number of dt;

% when beta = 1/6, it is linear acceleration;
% when beta >= 1/4, it is unconditional stable;
beta = 1/4;
gamma = 1/2; 
a0=1/beta/dt^2; a1=gamma/beta/dt; a2=1/beta/dt; a3=1/2/beta-1;
a4=gamma/beta-1; a5=dt*(gamma/beta-2)/2; a6=dt*(1-gamma); a7=gamma*dt;

% initial state
x = zeros(cn, n); 
v = zeros(cn, n); 
a = zeros(cn, n); 
x(:,1) = X0(1:cn); 
v(:,1) = X0(cn+1:2*cn);
a(:,1) = X0(2*cn+1:3*cn);

K_ = K + a0*M + a1*C;
iK = inv(K_);

% x - modal displacement
% v - modal velocitie
% a - modal acceleration
for i = 1:length(F)-1
    F_(:,i+1)=F(:,i+1)+M*(a0*x(:,i)+a2*v(:,i)+a3*a(:,i))+C*(a1*x(:,i)+a4*v(:,i)+a5*a(:,i));
    x(:,i+1)=iK*F_(:,i+1);
    a(:,i+1)=a0*(x(:,i+1)-x(:,i))-a2*v(:,i)-a3*a(:,i);
    v(:,i+1)=v(:,i)+a6*a(:,i)+a7*a(:,i+1);
end

%% use Phi matrix to get cable response
A_d = phiA'*x; % displacement at point A(L/3)
B_d = phiB'*x; % displacement at point B(L/2)
C_d = phiC'*x; % displacement at point C(2*L/3)
D_d = phiD'*x; % displacement at point Damper(x_d)
E_d= phiE'*x; % displacement at point Excitation(L/6)

A_v = phiA'*v; % veloctiy at point A(L/3)
B_v = phiB'*v; % velocity at point B(L/2)
C_v = phiC'*v; % velocity at point C(2*L/3)
D_v = phiD'*v; % velocity at point Damper(x_d)
E_v= phiE'*v; % velocity at point Excitation(L/6)

A_a = phiA'*a; % veloctiy at point A(L/3)
B_a = phiB'*a; % velocity at point B(L/2)
C_a = phiC'*a; % velocity at point C(2*L/3)
D_a = phiD'*a; % velocity at point Damper(x_d)
E_a= phiE'*a; % velocity at point Excitation(L/6)

%% make graph
figure(1);
plot(t,B_d);
