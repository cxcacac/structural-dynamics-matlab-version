function [d,v] = stateSpaceMatrixTransfer(A,B,F,X0,dt)
cn = length(X0)/2;
F0 = B*F;

As = expm(A.*dt);
Y = zeros(length(X0), length(F));

Y(:, 1) = X0;
for i = 1:length(F)-1
    Y(:,i+1) = As*Y(:,i)+(As*F0(:,i)).*dt;
end
d = Y(1:cn,:);
v = Y(cn+1:2*cn,:);
end