function [C,T,z] = dampR(K,M,E,es,flag)
% input:
% K: stiffness matrix, M: mass matrix;
% E: loading location matrix;
% es: damping ratio;
% flag = 1: seismic loads, flag = 0: normal loads;
% output:
% C: damping matrix, T: periods, z: vibration modes;
[z, d] = eig(K, M);
w = diag(sqrt(d));

cn = length(w);
C = zeros(cn);

if(flag==0)
    ix1 = 1;
    ix2 = 2;
else
    for i=1:cn
        M_star(i) = z(:,i)'*M*z(:,i);
        % declear specific index
        r(i) = sum(z(:,i)'*M*E)./M_star(i);
        meq(i) = r(i)*r(i)*M_star(i);
    end
    [Meq, ix] = sort(meq);
    ix1 = ix(length(ix));
    ix2 = ix(length(ix)-1);
end
a = 2*es*w(ix1)*w(ix2)/(w(ix1)+w(ix2));
b = 2*es/(w(ix1)+w(ix2));
C = a.*M + b.*K;
T = (2*pi)./w;
