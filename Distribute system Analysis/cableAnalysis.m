clc;
clear;
% clf;

length = 220; % the length of cable is 220m;
n = 5; % number of modes used to represent cable;
dLocation = 4; % damper location is 4m;
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
        bb = diff(diff(b)).*a;
        M(i,j) = mass.*int(a.*b,0, length);
        K(i,j) = (-T).*int(bb, 0, length);
    end
end

[mode, w2] = eig(K,M);
w = diag(sqrt(w2));

