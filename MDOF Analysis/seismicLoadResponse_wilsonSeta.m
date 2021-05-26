clear;
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

cn = length(m);
X0 = zeros(3*cn, 1);
R = zeros(cn,1);
F = [F, zeros(cn,1)];
[d,v,a] = wilsonSeta(M,K,C,dt,F,X0,K,R);

figure(1);
plot(t,d);

xlabel('time/s');
ylabel('displacement/m');