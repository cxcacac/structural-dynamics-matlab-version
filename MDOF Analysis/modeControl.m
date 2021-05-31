function [G] = modeControl(M, C, K, D, F,s)
% s: controlled modes;
cn = length(M);
sn = length(s);
[mode, w2] = eig(K,M);
TT = mode(:,s);

m = TT'*M*TT;
k = TT'*K*TT;
c = TT'*C*TT;
d = TT'*D;
f = TT'*F;

for i=1:sn
    A = [0,1; -inv(m(s(i),s(i)))*k(i,i), -inv(k(s(i),s(i)))*c(s(i),s(i))];
    B = [0; inv(m(s(i), s(i)))];
    Q = [k(s(i), s(i)), 0; 0, m(s(i),s(i))];
    R = 0.161;
    Gc(i,:) = lqr(A,B,Q,R);
end

if(sn==cn)
    Lc = inv(d);
elseif(sn>cn)
    Lc = d'*pinv(d'*d);
elseif(sn<cn)
    Lc = pinv(d'*d)*d';
end

Tc = [pinv(TT), zeros(size(pinv(TT))); zeros(size(pinv(TT))), pinv(TT)];
Gc_1 = [diag([Gc(:,1)]), diag([Gc(:,2)])];
G = Lc*Gc_1*Tc;