clc;
clear;
% clf;

length = 220; % the length of cable is 220m;
n = 20; % number of modes used to represent cable;
dLocation = 4;
mass = 85; % Unit: kg/m;
T = 7115*1e3; % Unit: KN;
Diameter = 0.118; % Unit: m;
angle = 22.5; % Unit: degrees;

M = zeros(n);
K = zeros(n);

syms x;

for i = 1:n
    for j = 1:n
        % using the dot (.) operator for element-wise multiplication (.*).
        a = shapeFunction(x,i,dLocation,length);
        b = shapeFunction(x,j,dLocation,length);
        aa= matlabFunction(diff(diff(b)).*a);
        h = @(x) shapeFunction(x,i,dLocation,length).*shapeFunction(x,j,dLocation,length);
        g = @(x) aa(x);
        M(i,j) = mass.*integral(h, 0, length);
        K(i,j) = (-T).*integral(g, 0, length);
    end
end

[mode, w2] = eig(K,M);
w = diag(sqrt(w2));

