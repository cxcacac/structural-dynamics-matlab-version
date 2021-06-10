function [F] = getLoadVector(load, n, slots, x_d, L)
% n: number of dof;
F = zeros(n, slots);
syms x;
for i = 1:n
    F(i,:) = int(shapeFunction(x,i,x_d, L),x,0, L).*load;
end