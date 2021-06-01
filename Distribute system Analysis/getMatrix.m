function [K,M,C,modes,w] = getMatrix(L, x_d, m, n, T,c)
% calculate the stiffness matrix, mass matrix and damping matrix;
% the script is based on Li's paper (2007,2008)
% (1)Vibration control of stay cables of the shandong binzhou rivers;
% (2)Negative stiffness characteristic of active and semi-active control
% considering the static deformation as shape function (Johnson, 2007)s

M = zeros(n);
K = zeros(n);

for i = 1:n
    for j = 1:n
       if(i==1 && j==1)
           M(i,j) = m*L/3;
       elseif(i > 1 && j > 1)
           M(i,j) = m*L/2*(i==j);
       else
           k = max(i,j)-1;
           M(i,j) = (m*L^3*sin(k*x_d*pi/L))/((pi*k)^2*x_d*(L-x_d));
       end
    end
end

for i = 1:n
    for j = 1:n
       if(i==1 && j==1)
           K(i,j) = L*T/(x_d*(L-x_d));
       elseif(i > 1 && j > 1)
           K(i,j) = T*pi^2*(i-1)^2*(i==j)/(2*L);
       else
           k = max(i,j)-1;
           K(i,j) = (T*L*sin(k*x_d*pi/L))/(x_d*(L-x_d));
       end
    end
end

[modes, w2] = eig(K, M);
w = diag(sqrt(w2))./(2*pi); % check w == w0;

C = (M*c)./m;