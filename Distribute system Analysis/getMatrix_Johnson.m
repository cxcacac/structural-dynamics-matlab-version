% The procedure is based on Johnson 's paper(2007)
% Semiactive damping of stay cables;
clc;
clear;
% clf;

L = 220; % the length of cable is 220m;
n = 5; % number of modes used to represent cable;
x_d = 0.02; % damper location is 4m;
m = 85; % Unit: kg/m;
T = 7115*1e3; % Unit: N;
Diameter = 0.118; % Unit: m;
angle = 22.5; % Unit: degrees;
c = 0.05; % visous damping per unit;

numberofModes = (1:n);
w0 = sqrt(T/m).*(numberofModes)./(2*L);

M = zeros(n);
K = zeros(n);

for i = 1:n
    for j = 1:n
       if(i==1 && j==1)
           M(i,j) = 1/3;
       elseif(i > 1 && j > 1)
           M(i,j) = 1/2*(i==j);
       else
           k = max(i,j)-1;
           M(i,j) = sin(k*x_d*pi)/((pi*k)^2*x_d*(1-x_d));
       end
    end
end

for i = 1:n
    for j = 1:n
       if(i==1 && j==1)
           K(i,j) = 1/(x_d*(1-x_d)*pi^2);
       elseif(i > 1 && j > 1)
           K(i,j) = (i-1)^2*(i==j)/2;
       else
           k = max(i,j)-1;
           K(i,j) = sin(k*x_d*pi)/(x_d*(1-x_d)*pi^2);
       end
    end
end

[modes, w2] = eig(K, M);
w = diag(sqrt(w2))./(2*pi);

