function [C,T,z] = dampC(K,M,es,n)
% input: 
% es: damping ratio; n: num of frequencies needed;
% ouput:
% C: damping matrix; T: periods; z:vibration type
[z, d] = eig(K,M);
w = diag(sqrt(d));

cn = length(w);
C = zeros(cn);

esn = 2.*es*ones(n,1);
for i = 1:n
    wn(:,i) = w(1:n).^(2i-3);
end
a = inv(wn)*esn;

for i = 1:n
    C = C + a(i).*M*(inv(M)*K)^(i-1);
end

T = (2*pi)./w;