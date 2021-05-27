function [dMOF, vMOF, aMOF] = modalSuperposition(M,K,C,F,dt)
[z,d] = eig(K,M);

m_ = diag(z'*M*z);
k_ = diag(z'*K*z);
c_ = diag(z'*C*z);
F_ = z'*F;

X0 = zeros(3,1);

for i = 1:length(m_)
    [d,v,a] = NewMarkBeta(m_(i), k_(i), c_(i), dt, F_(i,:), X0, k_(i),0);
    % directly adding lines on D, V, A;
    D(i,:) = d;
    V(i,:) = v;
    A(i,:) = a;
end
dMOF = z*D;
vMOF = z*V;
aMOF = z*A;
end