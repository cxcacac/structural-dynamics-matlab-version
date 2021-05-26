function [K] = stiffnessShear(k)
% input: stiffness at each layer;
% output: stiffness matrix;
cn = length(k);
% zeros(n) means create (n,n) matrix;
% zeros(m,n) means create (m,n) matrix;
K = zeros(cn);
for i=1:cn-1
    K(i,i) = k(i)+k(i+1);
    K(i,i+1) = -k(i+1);
    K(i+1, i) = -k(i+1);
end
K(cn,cn) = k(cn);