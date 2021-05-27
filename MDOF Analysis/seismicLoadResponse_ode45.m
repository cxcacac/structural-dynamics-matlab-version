clc;
clf;
m = [2762,2760,2300];
k = [2.485,1.921,1.522].*1e4;
es = 0.05;

wavefile = char('elcentro_NS.dat');
ugmax = 0.7;
[ug, t, tf,dt] = wave(wavefile, ugmax);

[M] = lumpMass(m);
[K] = stiffnessShear(k);
[E,F] = waveForce(ug, M);
flag = 1;
[C,T,z] = dampR(K,M,E,es,flag);

[A, B, D, L] = ssLinear(M, K, C);
cn = length(m);
X0 = zeros(2*cn, 1);
X = X0;
[d1,v1] = rk4('stateSpaceEq', [0, tf], dt, X0, A, B, F);
% in rk4, use B, F as coefficient to calculate F0;
F0 = B*F;
for i = 1:length(t)-1
    [tout, Y] = ode45(@stateSpaceEq_ode45, [t(i), t(i)+dt], X0, [], F0(:,i), A);
    X0 = Y(length(Y),:)';
    X = [X,X0];
end
d = X(1:cn,:);
v = X(cn:2*cn,:);

figure(1);
plot(t,d(1,:),'r-', t, d1(1,:), 'k:');

xlabel('time/s');
ylabel('displacement/m');
