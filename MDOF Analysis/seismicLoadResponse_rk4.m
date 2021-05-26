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

[d,v] = rk4('stateSpaceEq', [0, tf], dt, X0, A, B, F);

figure(1);
plot(t,d);

xlabel('time/s');
ylabel('displacement/m');
