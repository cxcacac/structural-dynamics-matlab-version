function [K] = relativeK(k0)
cn = length(k0);
K = zeros(cn);
for i = 1:cn-1
    K(i,i) = k0(i);
    K(i,i+1) = -k0(i+1);
end
K(cn, cn) = k0(cn);